#!/usr/bin/env bash
# =============================================================================
# masc-utils: GNPS CMN analysis environment setup
# Run from the root of your extracted GNPS zip (GNPS_CMN_ALL/)
# Usage: bash setup_env.sh
# =============================================================================

set -e

ENV_NAME="masc-utils"
PYTHON_VERSION="3.11"

echo ">>> Removing Zone.Identifier junk files (WSL artifacts)..."
find . -name "*.Zone.Identifier" -type f -delete
echo "    Done."

echo ">>> Creating conda environment: $ENV_NAME (Python $PYTHON_VERSION)..."
conda create -y -n "$ENV_NAME" python="$PYTHON_VERSION"

echo ">>> Activating environment..."
# Note: 'conda activate' doesn't work in subshells — after this script completes,
# run: conda activate masc-utils
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo ">>> Installing core packages..."
conda install -y -c conda-forge \
    pandas \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    networkx \
    scikit-learn \
    jupyter \
    jupyterlab \
    ipykernel \
    openpyxl \
    tqdm \
    pyyaml

echo ">>> Installing pip-only packages..."
pip install \
    pyteomics \
    matchms \
    spec2vec

echo ">>> Registering kernel for Jupyter..."
python -m ipykernel install --user --name "$ENV_NAME" --display-name "masc-utils"

echo ""
echo "============================================================"
echo " Setup complete!"
echo " To activate:  conda activate $ENV_NAME"
echo " To launch:    jupyter lab"
echo "============================================================"
