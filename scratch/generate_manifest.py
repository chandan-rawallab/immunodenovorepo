import os
import pandas as pd
import re

psms_dir = "/home/amity/Documents/experiments/data/psms"
files = [f for f in os.listdir(psms_dir) if f.startswith("msms_") and f.endswith(".txt")]

manifest_data = []

for f in sorted(files):
    # Pattern: msms_YYYYMMDD_PATIENT_RUN.raw.txt
    # Example: msms_20160513_TIL1_R2.raw.txt
    match = re.search(r"msms_(.*)\.raw\.txt", f)
    if match:
        run_id = match.group(1)
        # Try to extract patient_id. In these samples, it's often the part before _R1/_R2 or similar.
        # Bassani-Sternberg 2016 samples: TIL1, TIL2, DC1W6, GD149-1, etc.
        parts = run_id.split('_')
        if len(parts) >= 2:
            # The date is usually the first part. Let's look for the patient name.
            # 20160513_TIL1_R2 -> TIL1
            # 20160823_QEh1_LC2_HuPa_SA_HLApI_CM647_2_MG_1 -> CM647 (probably)
            
            # Simple heuristic for now:
            if "TIL" in run_id:
                patient_id = [p for p in parts if "TIL" in p][0]
            elif "GD149" in run_id:
                patient_id = "GD149"
            elif "DC" in run_id and "W6" in run_id:
                patient_id = [p for p in parts if "DC" in p][0]
            elif "CM647" in run_id:
                patient_id = "CM647"
            elif "RA957" in run_id:
                patient_id = "RA957"
            elif "MD155" in run_id:
                patient_id = "MD155"
            elif "Apher" in run_id:
                 patient_id = [p for p in parts if "Apher" in p][0]
            else:
                patient_id = "Unknown"
        else:
            patient_id = "Unknown"
            
        manifest_data.append({
            "run_id": run_id,
            "patient_id": patient_id,
            "filename": f,
            "hla_alleles": "TBD",
            "cohort": "PXD005231"
        })

df = pd.DataFrame(manifest_data)
# Group by patient_id to see run distributions
# patient_runs = df.groupby("patient_id")["run_id"].apply(list).reset_index()

df.to_csv("/home/amity/Documents/experiments/configs/sample_manifest.tsv", sep="\t", index=False)
print(f"Created manifest with {len(df)} runs.")
