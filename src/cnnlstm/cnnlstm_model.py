import torch
import torch.nn as nn

try:
    from .spectral_dataset import VOCAB_SIZE
except ImportError:
    from spectral_dataset import VOCAB_SIZE

class NeoepitopeSeq2Seq(nn.Module):
    """
    CNN-LSTM Deep Learning Architecture defined by Phase 3 (Workplan 3).
    Translates raw MGF mass-spectrometry signal into a natural Amino Acid sequence mapping 
    to identify and prioritize high-affinity biological neoantigens de novo.
    """
    def __init__(self, cnn_channels=128, lstm_hidden=256, num_lstm_layers=2):
        super().__init__()
        
        # 1. Feature Extractor (1D Convolutional Neural Network)
        # Because we binned M/Z into 20,000 discrete spaces, local spectral proximity matters.
        # CNN finds local correlation sequences (b and y ion ladders).
        self.cnn = nn.Sequential(
            nn.Conv1d(in_channels=1, out_channels=cnn_channels, kernel_size=5, stride=2, padding=2),
            nn.BatchNorm1d(cnn_channels),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),
            
            nn.Conv1d(in_channels=cnn_channels, out_channels=cnn_channels*2, kernel_size=5, stride=2, padding=2),
            nn.BatchNorm1d(cnn_channels*2),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2)
        )
        
        # 2. Sequence Interpreter (Bi-Directional Long Short-Term Memory)
        # Translates abstracted CNN feature steps into temporal amino acid predictions
        self.lstm = nn.LSTM(
            input_size=cnn_channels*2,
            hidden_size=lstm_hidden,
            num_layers=num_lstm_layers,
            batch_first=True,
            bidirectional=True
        )
        
        # 3. Dense Classification Output
        # Maps the combined Bi-directional sequence matrices to the 20 biological amino acids + tokens
        # Two-layer head: (hidden*2 -> hidden) -> VOCAB_SIZE for better generalisation
        self.fc = nn.Sequential(
            nn.Linear(lstm_hidden * 2, lstm_hidden),
            nn.ReLU(),
            nn.Dropout(p=0.3),
            nn.Linear(lstm_hidden, VOCAB_SIZE)
        )

    def forward(self, x):
        # 1. CNN extraction: (Batch, 1, 20000) -> (Batch, 128, 1250)
        c_features = self.cnn(x) 
        
        # Add pooling to fix sequence length to exactly 30 steps
        # This maps the 1250 spectral features to the 30 sequence slots
        c_features = nn.functional.adaptive_avg_pool1d(c_features, 30)
        
        # Reshape for LSTM (Batch, 30, 128)
        c_features = c_features.transpose(1, 2)
        
        # 2. LSTM Sequencing
        r_out, (h_n, c_n) = self.lstm(c_features)
        
        # 3. Output prediction: (Batch, 30, VOCAB_SIZE)
        predictions = self.fc(r_out)
        
        return predictions

if __name__ == "__main__":
    # Test instantiation and forward pass
    print("Initializing NeoepitopeSeq2Seq Model...")
    model = NeoepitopeSeq2Seq()
    print(model)
    print("Testing forward pass logic with mock tensor (Batch=4, M/Z_Bins=20000)...")
    mock_input = torch.rand(4, 1, 20000)
    output = model(mock_input)
    print(f"Output Matrix Shape: {output.shape} (Batch, SeqLocs, ProbabilityDist)")
    print("Model tests successfully passed parameters.")
