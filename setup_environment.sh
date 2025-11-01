#!/bin/bash

# Setup script for PackRat development environment
# This creates a minimal conda environment with R and necessary packages

echo "=================================="
echo "PackRat Environment Setup"
echo "=================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH"
    echo "Please install conda or load the appropriate module"
    exit 1
fi

# Set environment name
ENV_NAME="packrat-dev"

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '${ENV_NAME}' already exists."
    echo -n "Do you want to remove and recreate it? (y/n): "
    read answer
    if [ "$answer" = "y" ]; then
        echo "Removing existing environment..."
        conda env remove -n ${ENV_NAME} -y
    else
        echo "Using existing environment. To activate:"
        echo "  conda activate ${ENV_NAME}"
        exit 0
    fi
fi

# Create environment
echo "Creating conda environment '${ENV_NAME}'..."
conda env create -f environment.yml

echo ""
echo "=================================="
echo "Environment created successfully!"
echo "=================================="
echo ""
echo "To use this environment:"
echo "  conda activate ${ENV_NAME}"
echo ""
echo "To test the installation:"
echo "  conda activate ${ENV_NAME}"
echo "  Rscript tests/test_installation.R"
echo ""
echo "To develop interactively:"
echo "  conda activate ${ENV_NAME}"
echo "  R"
echo "  > devtools::load_all()"