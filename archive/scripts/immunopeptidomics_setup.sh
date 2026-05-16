#!/bin/bash

# =========================================================================
# PROJECT: IMMUNOPEPTIDOMICS DISCOVERY PIPELINE
# ARCHITECTURE: HIGH-PERFORMANCE SEQUENTIAL/PARALLEL PROCESSING
# =========================================================================

echo "[INIT] Initializing Modern Academic Research Environment..."

# Smart Path Detection: Prevent "immunopeptidomics_pipeline/immunopeptidomics_pipeline" nesting
if [[ "$PWD" == *"/immunopeptidomics_pipeline" ]]; then
    echo "[INFO] Already inside project root. Deploying to current directory."
    PROJECT_ROOT="."
else
    PROJECT_ROOT="immunopeptidomics_pipeline"
fi

# 1. Standardized Directory Structure
mkdir -p $PROJECT_ROOT/bin
mkdir -p $PROJECT_ROOT/data/raw
mkdir -p $PROJECT_ROOT/data/reference
mkdir -p $PROJECT_ROOT/src/configs
mkdir -p $PROJECT_ROOT/src/utils
mkdir -p $PROJECT_ROOT/results/identifications
mkdir -p $PROJECT_ROOT/logs

cd $PROJECT_ROOT

# 2. Local Computational Runtimes
# 2.1 .NET 8 Runtime (Required for MaxQuant and v2.0 Conversion)
if [ ! -d "bin/dotnet_runtime" ]; then
    echo "[DEPLOY] Installing .NET 8 Runtime Core..."
    mkdir -p bin/dotnet_runtime
    wget -q https://dot.net/v1/dotnet-install.sh -O bin/dotnet-install.sh
    chmod +x bin/dotnet-install.sh
    ./bin/dotnet-install.sh --channel 8.0 --install-dir ./bin/dotnet_runtime
fi

# 2.2 Conversion Tool (v2.0.0-dev Native Linux)
if [ ! -d "bin/ThermoRawFileParser" ]; then
    echo "[DEPLOY] Downloading ThermoRawFileParser v2.0.0-dev (Native Linux)..."
    wget -c "https://github.com/CompOmics/ThermoRawFileParser/releases/download/v.2.0.0-dev/ThermoRawFileParser-v.2.0.0-dev-linux.zip" -O bin/ThermoRawFileParser.zip
    unzip -q bin/ThermoRawFileParser.zip -d bin/ThermoRawFileParser
    rm bin/ThermoRawFileParser.zip
    chmod +x bin/ThermoRawFileParser/ThermoRawFileParser
fi

# 2.3 Reference Foundation (UniProt Human Reviewed)
if [ ! -f "data/reference/uniprot_human_reviewed.fasta" ]; then
    echo "[DEPLOY] Downloading UniProt Human Reviewed Proteome..."
    wget -c "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:9606)" -O data/reference/uniprot_human_reviewed.fasta
fi

# 3. Virtual Environment Provisions
echo "[ENV] Conda Deployment Protocol:"
echo "------------------------------------------------"
echo "conda create -n hla_deep_learning python=3.10 -y"
echo "conda activate hla_deep_learning"
echo "pip install torch numpy pandas pyteomics glob2"
echo "------------------------------------------------"

# 4. Orchestration Scripts
cat <<EOF > src/execute_production_search.sh
#!/bin/bash
# High-Performance Execution Entrypoint
export PATH=\$PWD/bin/dotnet_runtime:\$PATH
./bin/dotnet_runtime/dotnet ./bin/MaxQuant_v2.8.0.0/bin/MaxQuantCmd.dll src/configs/mqpar_production.xml
EOF
chmod +x src/execute_production_search.sh

echo "[SUCCESS] Modern Foundation Validated."
echo "!! MANUAL ACTION REQUIRED: Copy MaxQuant into bin/MaxQuant_v2.8.0.0 !!"
