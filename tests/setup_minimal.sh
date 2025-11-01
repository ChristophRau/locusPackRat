#!/bin/bash

# Setup script for minimal R environment to test dependency installation
# This creates a barebones conda environment with ONLY R

echo "============================================"
echo "PackRat Minimal Environment Setup"
echo "============================================"
echo ""
echo "This creates a minimal environment with ONLY R"
echo "to test how PackRat handles dependency installation."
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH"
    exit 1
fi

ENV_NAME="packrat-minimal"

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Removing existing '${ENV_NAME}' environment..."
    conda env remove -n ${ENV_NAME} -y
fi

# Create minimal environment
echo "Creating minimal conda environment '${ENV_NAME}'..."
echo "This will have ONLY R, no packages pre-installed."
conda env create -f tests/environment_minimal.yml

echo ""
echo "============================================"
echo "Minimal environment created!"
echo "============================================"
echo ""
echo "To test dependency installation:"
echo "  conda activate ${ENV_NAME}"
echo "  Rscript tests/test_dependency_install.R"
echo ""
echo "This will show which dependencies are:"
echo "  - Missing and need installation"
echo "  - Successfully installed automatically"
echo "  - Failed to install"