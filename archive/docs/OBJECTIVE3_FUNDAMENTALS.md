# Objective 3: The Ground-Level Guide to Neoepitope Discovery

This guide explains the "what, why, and how" of our hybrid immunopeptidomics pipeline, intended for researchers and engineers setting up the system manually.

---

## 1. The Big Picture (The Bio Context)
**The Goal**: To find "Neoepitopes"—tiny protein fragments (peptides) that exist on the surface of cancer cells but NOT on healthy cells. These are the "targets" used to create personalized cancer vaccines.

**The Challenge**: Standard databases (like UniProt) only contain "normal" human proteins. Cancer proteins are **mutated**, meaning they have typos in their sequences. A standard search will miss them.

---

## 2. Our "Hybrid" Solution
We use two different methods and combine them to find the "hidden" mutations.

### Method A: Database Search (MaxQuant)
*   **Role**: Identifies the "known" and "normal" peptides.
*   **Math**: Matches your experimental data against a known library (UniProt).
*   **Analogy**: Looking for words in a standard dictionary.

### Method B: De Novo Sequencing (Deep Learning)
*   **Role**: Predicts sequences directly from the signal, even if they aren't in any dictionary.
*   **Math**: Uses a **Seq2Seq Neural Network** (similar to Google Translate) to translate Mass Spectrometry "waves" into Amino Acid letters (A, C, D, E...).
*   **Analogy**: Sounding out a word you've never heard before based on phonetics.

---

## 3. Step-by-Step Logic (The Pipeline)

### Step 1: Data Retrieval
*   **Input**: A Project ID (e.g., `PXD005231`).
*   **Action**: Downloads massive `.raw` files containing millions of "spectral signatures."

### Step 2: Spectral Conversion
*   **Input**: `.raw` (Proprietary Thermo format).
*   **Output**: `.mgf` (Plain text list of peaks).
*   **Why?**: Neural networks cannot read binary proprietary files; they need the raw mass-to-charge ($m/z$) numbers.

### Step 3: The Subtraction Math
This is the core of our discovery logic:
1.  **Run MaxQuant**: Identify everything that is "Normal."
2.  **Run De Novo**: Predict everything that is "Possible."
3.  **Subtract**: `Potential Neoepitopes = (De Novo Results) - (MaxQuant Results)`
*   **Example**: If De Novo predicts `TYR-MUT-SEQ` and MaxQuant finds nothing, it's a high-priority candidate.

---

## 4. Sample Input & Output Examples

### A. Input Data (The Spectrum)
A single "scan" in an MGF file looks like this:
```text
BEGIN IONS
TITLE=Scan_1234_RT_45.2
PEPMASS=450.231  # The mass of the peptide
110.071 1205.4   # [Mass] [Intensity]
129.102 4500.1   # This is the "signal"
END IONS
```

### B. Intermediate Output (The Prediction)
The Deep Learning model takes that scan and outputs:
*   **Predicted Sequence**: `LLGRNSFEV`
*   **Confidence Score**: `0.92` (92% sure)

### C. Final Output (The Ranked List)
The final `ranked_neoantigens.tsv` combines biology and math:

| Peptide | Source | Mutation | Binding Score (nM) | Final Rank |
| :--- | :--- | :--- | :--- | :--- |
| **SLYNTVATL** | De Novo | V -> A | 12.5 (Strong) | **1** |
| **LLGRNSFEV** | De Novo | Novel | 450.0 (Weak) | **14** |

---

## 5. Manual Execution Check-list
When running this manually, always check:
1.  **Paths**: Is `dotnet` in the right `bin/` folder?
2.  **Fasta**: Is the `uniprot_human.fasta` absolute path correct in the XML?
3.  **Symlinks**: Are the large `.raw` files "pointed to" correctly so you don't run out of disk space?

---
**Summary**: We find the "knowns" to reveal the "unknowns." By subtracting standard biology from deep-learning predictions, we isolate the mutations that matter for cancer therapy.
