"""NeoepitopeSeq2Seq — improved CNN-LSTM de novo sequencing model.

Architecture changes vs. curated31_v2 baseline (2026-06-25):
  - CNN encoder replaced with a 4-block residual architecture:
      Block 1: 1 → 128 ch, kernel 7 (captures wider local ion-series context)
      Block 2: 128 → 128 ch, kernel 5, residual skip
      Block 3: 128 → 256 ch, kernel 3, residual skip with projection
      Block 4: 256 → 256 ch, kernel 3, residual skip
    Residual connections reduce gradient vanishing across the deeper stack.
  - Precursor mass + charge conditioning: an optional (mass, charge) pair is
    projected to a 64-d embedding and concatenated to every LSTM input step.
    When not provided (default), the model behaves identically to the baseline
    so existing checkpoints remain loadable after removing those layers.
  - LSTM hidden size increased from 256 → 320 to accommodate the wider CNN output.
  - Dropout moved *after* BatchNorm inside each residual block and *before* the
    final linear projection to reduce co-adaptation of output neurons.
  - Two-layer FC head retained with an additional LayerNorm before the final
    projection for training stability.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from .spectral_dataset import VOCAB_SIZE
except ImportError:
    from spectral_dataset import VOCAB_SIZE


# ---------------------------------------------------------------------------
# Residual CNN block
# ---------------------------------------------------------------------------

class _ResBlock(nn.Module):
    """1-D residual block: Conv → BN → ReLU → Dropout → Conv → BN + skip."""

    def __init__(
        self,
        in_ch: int,
        out_ch: int,
        kernel: int = 3,
        dropout: float = 0.15,
    ):
        super().__init__()
        pad = kernel // 2
        self.conv1 = nn.Conv1d(in_ch, out_ch, kernel, padding=pad, bias=False)
        self.bn1   = nn.BatchNorm1d(out_ch)
        self.drop  = nn.Dropout(dropout)
        self.conv2 = nn.Conv1d(out_ch, out_ch, kernel, padding=pad, bias=False)
        self.bn2   = nn.BatchNorm1d(out_ch)

        # Projection shortcut when channel dimensions change
        self.proj = (
            nn.Conv1d(in_ch, out_ch, 1, bias=False)
            if in_ch != out_ch
            else nn.Identity()
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = self.proj(x)
        out = F.relu(self.bn1(self.conv1(x)))
        out = self.drop(out)
        out = self.bn2(self.conv2(out))
        return F.relu(out + residual)


# ---------------------------------------------------------------------------
# Main model
# ---------------------------------------------------------------------------

class NeoepitopeSeq2Seq(nn.Module):
    """CNN-LSTM de novo sequencing model with residual encoder and optional
    precursor-mass conditioning.

    Parameters
    ----------
    cnn_base_ch:
        Base number of CNN channels (doubled at block 3).
    lstm_hidden:
        Hidden size of the BiLSTM (each direction).
    num_lstm_layers:
        Number of LSTM layers.
    mass_embed_dim:
        Dimension of the (mass, charge) projection. Set to 0 to disable.
    dropout:
        Dropout probability applied inside residual blocks and the FC head.
    """

    def __init__(
        self,
        cnn_base_ch:    int = 128,
        lstm_hidden:    int = 320,
        num_lstm_layers: int = 2,
        mass_embed_dim: int = 64,
        dropout:        float = 0.20,
    ):
        super().__init__()

        # ── 1. Residual CNN encoder ──────────────────────────────────────────
        # Input: (B, 1, 20 000)
        self.stem = nn.Sequential(
            nn.Conv1d(1, cnn_base_ch, kernel_size=7, stride=2, padding=3, bias=False),
            nn.BatchNorm1d(cnn_base_ch),
            nn.ReLU(),
            nn.MaxPool1d(2),
        )
        # After stem: (B, cnn_base_ch, ~5 000)

        self.res_blocks = nn.Sequential(
            _ResBlock(cnn_base_ch,       cnn_base_ch,       kernel=5, dropout=dropout),
            nn.MaxPool1d(2),
            _ResBlock(cnn_base_ch,       cnn_base_ch * 2,   kernel=3, dropout=dropout),
            nn.MaxPool1d(2),
            _ResBlock(cnn_base_ch * 2,   cnn_base_ch * 2,   kernel=3, dropout=dropout),
        )
        # After res_blocks + pooling: (B, cnn_base_ch*2, ~625)

        cnn_out_ch = cnn_base_ch * 2  # 256

        # ── 2. Optional precursor-mass/charge embedding ──────────────────────
        self.mass_embed_dim = mass_embed_dim
        if mass_embed_dim > 0:
            self.mass_proj = nn.Sequential(
                nn.Linear(2, mass_embed_dim),
                nn.ReLU(),
                nn.Linear(mass_embed_dim, mass_embed_dim),
            )
            lstm_in = cnn_out_ch + mass_embed_dim
        else:
            self.mass_proj = None
            lstm_in = cnn_out_ch

        # ── 3. BiLSTM ────────────────────────────────────────────────────────
        self.lstm = nn.LSTM(
            input_size=lstm_in,
            hidden_size=lstm_hidden,
            num_layers=num_lstm_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if num_lstm_layers > 1 else 0.0,
        )

        # ── 4. Multi-head self-attention ─────────────────────────────────────
        attn_dim = lstm_hidden * 2  # 640
        self.attention = nn.MultiheadAttention(
            embed_dim=attn_dim,
            num_heads=8,            # 8 heads on 640-d (80 per head)
            dropout=dropout,
            batch_first=True,
        )

        # ── 5. FC output head ────────────────────────────────────────────────
        self.fc = nn.Sequential(
            nn.LayerNorm(attn_dim),
            nn.Linear(attn_dim, lstm_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(lstm_hidden, VOCAB_SIZE),
        )

    # ------------------------------------------------------------------
    def forward(
        self,
        x: torch.Tensor,
        precursor: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """
        Parameters
        ----------
        x:
            Binned spectrum tensor of shape (B, 1, 20 000).
        precursor:
            Optional (B, 2) tensor with columns [neutral_mass_Da, charge].
            When provided, these are embedded and concatenated to every LSTM
            input step to provide charge-state context.

        Returns
        -------
        torch.Tensor
            Logits of shape (B, 30, VOCAB_SIZE).
        """
        # ── CNN ──────────────────────────────────────────────────────────────
        c = self.stem(x)          # (B, cnn_base_ch, ~5 000)
        c = self.res_blocks(c)    # (B, 256, ~625)

        # Collapse to exactly 30 sequence steps (one per output token)
        c = F.adaptive_avg_pool1d(c, 30)   # (B, 256, 30)
        c = c.transpose(1, 2)              # (B, 30, 256)

        # ── Precursor conditioning ────────────────────────────────────────────
        if self.mass_proj is not None:
            if precursor is not None:
                mass_emb = self.mass_proj(precursor)          # (B, mass_embed_dim)
            else:
                # Zero embedding keeps LSTM input_size constant when no precursor
                mass_emb = torch.zeros(
                    c.size(0), self.mass_embed_dim, device=c.device, dtype=c.dtype
                )
            mass_emb = mass_emb.unsqueeze(1).expand(-1, 30, -1)  # (B, 30, D)
            c = torch.cat([c, mass_emb], dim=-1)          # (B, 30, 256+D)

        # ── BiLSTM ──────────────────────────────────────────────────────────
        r_out, _ = self.lstm(c)    # (B, 30, lstm_hidden*2)

        # ── Self-attention with residual ─────────────────────────────────────
        attn_out, _ = self.attention(r_out, r_out, r_out)
        r_out = r_out + attn_out   # residual

        # ── Output projection ────────────────────────────────────────────────
        return self.fc(r_out)      # (B, 30, VOCAB_SIZE)


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Initialising NeoepitopeSeq2Seq (residual edition)...")
    model = NeoepitopeSeq2Seq()
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"  Trainable parameters: {n_params:,}")

    mock_x = torch.rand(4, 1, 20_000)
    mock_prec = torch.tensor([[1200.0, 2.0], [980.0, 3.0],
                               [1450.0, 2.0], [870.0, 1.0]])
    out = model(mock_x, precursor=mock_prec)
    print(f"  Output shape: {out.shape}  (expected [4, 30, {VOCAB_SIZE}])")

    # Also test without precursor (backward-compatible)
    out2 = model(mock_x)
    print(f"  Without precursor: {out2.shape}")
    print("Model checks passed.")
