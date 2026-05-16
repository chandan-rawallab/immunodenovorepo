# RAW to MGF Conversion Guide

This guide explains how to set up the environment and run the spectral conversion pipeline to transform Thermo `.raw` files into `.mgf` format for deep learning training.

## 1. Prerequisites

### .NET Runtime
The conversion tool requires the **.NET 8.0 Runtime**.
- **Linux**: Install via your package manager (e.g., `sudo apt install dotnet-runtime-8.0`).
- **Local/Portable**: If you don't have root access, you can use a local `dotnet` binary (as currently configured in this project in `bin/dotnet8`).

### ThermoRawFileParser
We use the [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) by Compomics.
- The binaries should be located in `bin/ThermoRawFileParser/`.

---

## 2. Installation & Troubleshooting

### Shared Library Issues
On some Linux systems, the tool may fail to find `.NET` shared libraries. If you see errors like `libhostpolicy.so not found`, ensure the following files are in the same folder as the `ThermoRawFileParser.dll`:
- `libhostpolicy.so`
- `libcoreclr.so`
- `Mono.Options.dll`

### Runtime Configuration
Ensure `ThermoRawFileParser.runtimeconfig.json` targets the version of .NET installed on your system. 
Example for .NET 8:
```json
{
  "runtimeOptions": {
    "tfm": "net8.0",
    "framework": {
      "name": "Microsoft.NETCore.App",
      "version": "8.0.0"
    }
  }
}
```

---

## 3. Running the Pipeline

We provide a wrapper script that automates the conversion for multiple files and ensures they align with your MaxQuant results.

### Automated Mode (Recommended)
This script scans your raw data folder and converts files that have matching MaxQuant results.

```bash
# Run from the project root
python3 src/utils/conversion_pipeline.py
```

### Script Logic
- **Input Folder**: `data/raw/` (scanned recursively).
- **Output Folder**: `/home/amity/hla_data_mgf` (centralized storage).
- **Synchronization**: By default, it only converts files that have an existing `msms.txt` in `production_results/`.

### Manual/Standalone Usage
To convert a single file manually:
```bash
# Using the dotnet multiplexer
dotnet exec bin/ThermoRawFileParser/ThermoRawFileParser.dll -i /path/to/data.raw -o /path/to/output/ -f 0
```
*(Flag `-f 0` specifies MGF format)*

---

## 4. Configuration
You can modify the paths in `src/utils/conversion_pipeline.py` to point to new data locations:
- `RAW_DIR`: Where your `.raw` files live.
- `MGF_OUTPUT_DIR`: Where you want the `.mgf` files saved.
- `RESULTS_DIR`: Where the script looks for MaxQuant results to decide what to convert.
