#!/bin/bash
# ==============================================================================
# Cell-Specific Epigenetic dPCR Assay Design Pipeline — All-in-One Local Setup
#
# Handles everything: Homebrew deps, Python venv, wgbstools, data download,
# genome decompression.
#
# Idempotent: safe to re-run. Skips completed steps.
#
# USAGE:
#   chmod +x setup_local.sh
#   ./setup_local.sh
#
# One-liner (from existing repo):
#   cd ~/cell-specific-epigenetic-dpcr-assay-design-pipeline && git pull origin main && chmod +x setup_local.sh && ./setup_local.sh
#
# ESTIMATED TIME: 2-3 hours (dominated by data download)
# DISK SPACE: ~10 GB (without dbSNP) or ~20 GB (with dbSNP)
# ==============================================================================

set -euo pipefail
cd "$(dirname "$0")"
REPO_ROOT="$(pwd)"

# --- Colors ---
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'
ok()   { echo -e "${GREEN}✓ $1${NC}"; }
fail() { echo -e "${RED}✗ $1${NC}"; exit 1; }
warn() { echo -e "${YELLOW}⚠ $1${NC}"; }
step() { echo -e "\n${GREEN}=== $1 ===${NC}"; }

# --- Detect platform ---
if [ -d /opt/homebrew ]; then
    HOMEBREW_BIN="/opt/homebrew/bin"
elif [ -d /usr/local/Homebrew ]; then
    HOMEBREW_BIN="/usr/local/bin"
else
    HOMEBREW_BIN="/usr/local/bin"
fi
export PATH="$HOMEBREW_BIN:$PATH"

THREADS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
if [ "$THREADS" -gt 8 ]; then THREADS=8; fi

echo "============================================"
echo "Cell-Specific Epigenetic dPCR Assay Design"
echo "Pipeline — Local Setup"
echo "============================================"
echo "  Repo:     $REPO_ROOT"
echo "  Homebrew: $HOMEBREW_BIN"
echo "  Threads:  $THREADS"
echo ""

# ==============================================================================
# STEP 1: Homebrew dependencies
# ==============================================================================
step "Step 1/6: Homebrew dependencies"

if ! command -v brew &> /dev/null; then
    echo "  Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi

BREW_PKGS="python@3.11 samtools htslib bedtools"
NEED_INSTALL=""
for pkg in $BREW_PKGS; do
    if ! brew list "$pkg" &>/dev/null; then
        NEED_INSTALL="$NEED_INSTALL $pkg"
    fi
done
if [ -n "$NEED_INSTALL" ]; then
    echo "  Installing:$NEED_INSTALL"
    brew install $NEED_INSTALL
else
    echo "  All Homebrew packages already installed"
fi
ok "Homebrew dependencies ready"

# ==============================================================================
# STEP 2: Python virtual environment + packages
# ==============================================================================
step "Step 2/6: Python environment"

if [ ! -f ".venv/bin/activate" ]; then
    echo "  Creating .venv..."
    "$HOMEBREW_BIN/python3.11" -m venv .venv
fi
source .venv/bin/activate

echo "  Checking Python packages..."
NEED_PIP=""
for pkg in primer3 openpyxl pandas numpy reportlab pysam scipy; do
    python -c "import $pkg" 2>/dev/null || NEED_PIP="$NEED_PIP $pkg"
done
if [ -n "$NEED_PIP" ]; then
    echo "  Installing missing packages:$NEED_PIP"
    pip install --upgrade pip -q
    pip install primer3-py openpyxl pandas numpy reportlab pysam scipy -q
else
    echo "  All Python packages already installed"
fi
ok "Python environment ready"

# ==============================================================================
# STEP 3: wgbstools
# ==============================================================================
step "Step 3/6: wgbstools"

if [ -f ".venv/bin/wgbstools" ] && .venv/bin/wgbstools --version 2>/dev/null; then
    echo "  wgbstools already installed"
else
    echo "  Installing wgbstools from GitHub..."
    if [ ! -d "wgbs_tools" ]; then
        git clone https://github.com/nloyfer/wgbs_tools.git
    fi
    (cd wgbs_tools && python3 setup.py 2>/dev/null)
    ln -sf "$(pwd)/wgbs_tools/wgbstools" .venv/bin/wgbstools
fi

if .venv/bin/wgbstools --version 2>/dev/null; then
    ok "wgbstools ready"
else
    warn "wgbstools may not be fully functional (needed for DMR discovery)"
fi

# ==============================================================================
# STEP 4: Data download
# ==============================================================================
step "Step 4/6: Data download"

chmod +x download_data.sh 2>/dev/null || true
if [ -f "download_data.sh" ]; then
    echo "  Running download_data.sh (skips existing files)..."
    bash download_data.sh || warn "download_data.sh had errors (some optional files may be missing)"
else
    warn "download_data.sh not found — skipping data download"
fi

# Verify critical data files
for f in data/GSE186458_blocks.s205.bed.gz data/sample_groups.csv data/sample_manifest.csv; do
    [ -f "$f" ] || fail "Missing required file: $f"
done

BETA_COUNT=$(ls data/beta_files/*.beta 2>/dev/null | wc -l | tr -d ' ')
if [ "$BETA_COUNT" -lt 100 ]; then
    warn "Only $BETA_COUNT beta files found (expected 207). Some cell types may fail."
else
    ok "Beta files: $BETA_COUNT"
fi

# ==============================================================================
# STEP 5: Genome — verify integrity, decompress, index
# ==============================================================================
step "Step 5/6: Genome (decompress + index)"

GENOME_FA="data/hg19/hg19.fa"
GENOME_GZ="data/hg19/hg19.fa.gz"

# Check if we already have a valid decompressed genome
if [ -f "$GENOME_FA" ] && [ -f "$GENOME_FA.fai" ]; then
    CHROM_COUNT=$(grep -c "^>" "$GENOME_FA" 2>/dev/null || echo 0)
    if [ "$CHROM_COUNT" -ge 25 ]; then
        echo "  Genome already decompressed and indexed ($CHROM_COUNT sequences)"
        ok "Genome ready"
    else
        warn "Genome file exists but has only $CHROM_COUNT sequences — may be incomplete"
        echo "  Re-decompressing..."
        rm -f "$GENOME_FA" "$GENOME_FA.fai"
        gunzip -c "$GENOME_GZ" > "$GENOME_FA" 2>/dev/null || {
            echo "  gunzip failed — re-downloading genome..."
            rm -f "$GENOME_GZ"
            wget -O "$GENOME_GZ" "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
            gunzip -c "$GENOME_GZ" > "$GENOME_FA"
        }
        samtools faidx "$GENOME_FA"
        ok "Genome re-decompressed and indexed"
    fi
else
    # Need to decompress
    if [ ! -f "$GENOME_GZ" ]; then
        echo "  Downloading hg19 genome (~900 MB)..."
        wget -O "$GENOME_GZ" "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    fi

    # Verify gzip integrity
    if ! gzip -t "$GENOME_GZ" 2>/dev/null; then
        warn "Genome file is corrupted — re-downloading..."
        rm -f "$GENOME_GZ"
        wget -O "$GENOME_GZ" "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
        gzip -t "$GENOME_GZ" || fail "Genome download still corrupted. Check your network."
    fi

    echo "  Decompressing genome..."
    gunzip -c "$GENOME_GZ" > "$GENOME_FA"
    echo "  Indexing with samtools faidx..."
    samtools faidx "$GENOME_FA"
    ok "Genome decompressed and indexed ($(grep -c '^>' "$GENOME_FA") sequences)"
fi

# ==============================================================================
# STEP 6: Create env_setup.sh + verify everything
# ==============================================================================
step "Step 6/6: Environment setup file + verification"

# Create env_setup.sh that run_cell_type.sh will source
cat > env_setup.sh << ENVEOF
#!/bin/bash
# Auto-generated by setup_local.sh — do not edit manually.
# Sources: .venv (Python + wgbstools), Homebrew (samtools etc.)
REPO_ROOT="$REPO_ROOT"
HOMEBREW_BIN="$HOMEBREW_BIN"

export PATH="$HOMEBREW_BIN:\$PATH"
source "$REPO_ROOT/.venv/bin/activate"
ENVEOF
chmod +x env_setup.sh
ok "Created env_setup.sh"

# --- Final verification ---
echo ""
echo "============================================"
echo "Verification"
echo "============================================"

PASS=0
FAIL=0
check() {
    if eval "$2" >/dev/null 2>&1; then
        ok "$1"
        PASS=$((PASS + 1))
    else
        fail "$1"
        FAIL=$((FAIL + 1))
    fi
}

check "samtools"            'samtools --version 2>/dev/null | head -1 | grep -q samtools'
check "tabix"               'tabix --version 2>/dev/null | head -1 | grep -q tabix'
check "primer3-py"          'python -c "import primer3"'
check "openpyxl"            'python -c "import openpyxl"'
check "pandas"              'python -c "import pandas"'
check "scipy"               'python -c "import scipy"'
check "reportlab"           'python -c "import reportlab"'
check "pysam"               'python -c "import pysam"'
check "wgbstools"           '.venv/bin/wgbstools --version 2>/dev/null'
check "genome FASTA"        '[ -f "$GENOME_FA" ]'
check "genome index"        '[ -f "$GENOME_FA.fai" ]'
check "beta files"          '[ $(ls data/beta_files/*.beta 2>/dev/null | wc -l) -ge 100 ]'
check "blocks file"         '[ -f data/GSE186458_blocks.s205.bed.gz ]'
check "groups CSV"          '[ -f data/sample_groups.csv ]'
check "manifest CSV"        '[ -f data/sample_manifest.csv ]'
check "blood Excel"         '[ -f dmr_regions_v2.2.8_260714.xlsx ]'
check "env_setup.sh"        '[ -f env_setup.sh ]'

echo ""
echo "============================================"
echo "Setup complete: $PASS passed, $FAIL failed"
echo "============================================"

if [ "$FAIL" -gt 0 ]; then
    warn "Some checks failed. Review the output above."
    exit 1
fi

echo ""
echo "To run a cell type:"
echo "  ./run_cell_type.sh Blood-Monocytes"
echo ""
echo "To list all available cell types:"
echo "  ./run_cell_type.sh --list"
echo ""
echo "To run a blood cell type (uses pre-computed DMRs):"
echo "  ./run_cell_type.sh MONO"
echo ""
