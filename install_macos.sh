#!/bin/bash
# =============================================================================
# Cell-Specific Epigenetic dPCR Assay Design Pipeline — macOS Installation
# =============================================================================
# Installs all dependencies needed to run the pipeline on macOS.
# Does NOT require biomni — everything runs locally.
#
# Usage:
#   chmod +x install_macos.sh
#   ./install_macos.sh
#
# Prerequisites:
#   - macOS (Intel or Apple Silicon)
#   - Xcode Command Line Tools (xcode-select --install)
# =============================================================================

set -e

echo "=========================================="
echo "Cell-Specific Epigenetic dPCR Assay Design"
echo "Pipeline — macOS Installation"
echo "=========================================="
echo ""

# Check for Homebrew
if ! command -v brew &> /dev/null; then
    echo "Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
    eval "$(/opt/homebrew/bin/brew shellenv)"
fi

echo "✓ Homebrew installed"

# Install system dependencies
echo ""
echo "Installing system dependencies..."
brew install python@3.11 samtools htslib bedtools

echo "✓ System dependencies installed (samtools, htslib/tabix, bedtools)"

# Create virtual environment
echo ""
echo "Creating Python virtual environment..."
python3.11 -m venv .venv
source .venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install Python packages
echo ""
echo "Installing Python packages..."
pip install primer3-py==2.3.0
pip install openpyxl pandas numpy
pip install reportlab
pip install pysam
pip install scipy

echo "✓ Python packages installed (including scipy for find_markers t-tests)"

# Install wgbstools (from GitHub — not on PyPI, not pip-installable)
echo ""
echo "Installing wgbstools (from GitHub source)..."
if [ ! -d "wgbs_tools" ]; then
    git clone https://github.com/nloyfer/wgbs_tools.git
fi

# Step 1: Compile C++ tools (patter, segmentor, etc.)
cd wgbs_tools
python3 setup.py
cd ..

# Step 2: Make wgbstools available on PATH.
# wgbstools is a symlink to src/python/wgbs_tools.py. When Python runs it,
# it resolves the symlink and adds src/python/ to sys.path[0], which lets
# it find find_markers.py and other modules via importlib.
# pip install -e . does NOT work — pyproject.toml has the packages line
# commented out and there is no __init__.py package directory.
ln -sf "$(pwd)/wgbs_tools/wgbstools" .venv/bin/wgbstools

echo "✓ wgbstools installed"

# Verify wgbstools is on PATH
if ! command -v wgbstools &> /dev/null; then
    echo "  WARNING: wgbstools not on PATH. Add .venv/bin to your PATH or run:"
    echo "    ln -sf \"$(pwd)/wgbs_tools/wgbstools\" .venv/bin/wgbstools"
fi

echo "✓ wgbstools available: $(which wgbstools 2>/dev/null || echo 'NOT FOUND')"

# Re-verify scipy after wgbstools install — the failed pip install -e .
# in older versions of this script could roll back scipy.
if ! python3 -c "import scipy" 2>/dev/null; then
    echo "  scipy was lost during wgbstools install, reinstalling..."
    pip install scipy
fi
echo "✓ scipy verified"

# NOTE: init_genome is NOT needed for the pipeline.
# find_markers works directly with the blocks file and beta files.
# The pipeline finds CpG positions by scanning the genomic sequence
# in step 2, so no CpG dictionary is required.

# Verify installation
echo ""
echo "=========================================="
echo "Verifying installation..."
echo "=========================================="

echo -n "Python: "; python3 --version
echo -n "samtools: "; samtools --version | head -1
echo -n "tabix: "; tabix --version 2>/dev/null | head -1 || echo "included with htslib/samtools"
echo -n "primer3-py: "; python3 -c "import primer3; print(primer3.__version__)"
echo -n "openpyxl: "; python3 -c "import openpyxl; print(openpyxl.__version__)"
echo -n "pandas: "; python3 -c "import pandas; print(pandas.__version__)"
echo -n "reportlab: "; python3 -c "import reportlab; print(reportlab.Version)"
echo -n "scipy: "; python3 -c "import scipy; print(scipy.__version__)"
echo -n "wgbstools: "; wgbstools --version 2>/dev/null || echo "WARNING: wgbstools not working"

echo ""
echo "=========================================="
echo "Installation complete!"
echo "=========================================="
echo ""
echo "To activate the environment:"
echo "  source .venv/bin/activate"
echo ""
echo "To run the full pipeline (DMR discovery + primer design):"
echo "  python -m methyl_panel.pipeline --steps all --discover-dmrs \\"
echo "    --cell-type MONO \\"
echo "    --genome data/hg19/hg19.fa.gz \\"
echo "    --min-tm 58 --opt-tm 60 --max-tm 62 \\"
echo "    --output-dir results/"
echo ""
echo "To run for all 7 cell types, repeat for each:"
echo "  for CT in MONO BCELL NK GRAN CD3T CD4T CD8T; do"
echo "    python -m methyl_panel.pipeline --steps all --discover-dmrs \\"
echo "      --cell-type \$CT \\"
echo "      --genome data/hg19/hg19.fa.gz \\"
echo "      --min-tm 58 --opt-tm 60 --max-tm 62 \\"
echo "      --output-dir results/\$CT"
echo "  done"
echo ""
echo "To load DMR blocks from a pre-computed Excel file (step 1):"
echo "  python -m methyl_panel.pipeline --steps all \\"
echo "    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \\"
echo "    --genome data/hg19/hg19.fa.gz \\"
echo "    --min-tm 58 --opt-tm 60 --max-tm 62 \\"
echo "    --output-dir results/"
