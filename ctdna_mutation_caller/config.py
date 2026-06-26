"""
config.py
---------
Configuration file for the LLM-Assisted ctDNA Somatic Mutation Caller.
Defines global directories, thresholds for variant filtering, LLM API
parameters, and cancer driver gene lists for annotations.
"""

import os

# ── Directories ────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
DOWNLOADS_DIR = os.path.join(BASE_DIR, "downloads")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs")
PLOTS_DIR = os.path.join(OUTPUTS_DIR, "plots")
REPORTS_DIR = os.path.join(OUTPUTS_DIR, "reports")
LLM_CACHE_DIR = os.path.join(OUTPUTS_DIR, "llm_cache")

# Ensure required directories exist
for directory in [DATA_DIR, DOWNLOADS_DIR, OUTPUTS_DIR, PLOTS_DIR, REPORTS_DIR, LLM_CACHE_DIR]:
    os.makedirs(directory, exist_ok=True)

# ── GDC API Configurations ─────────────────────────────────────────────────────
# Filters for TCGA-LUAD somatic MAF downloads
GDC_API_BASE = "https://api.gdc.cancer.gov"
TCGA_PROJECT = "TCGA-LUAD"

# ── Statistical and Biological Variant Filtering Thresholds ────────────────────
# ctDNA requires deep sequencing, so we filter out low-coverage sites.
MIN_DEPTH = 20           # Minimum total read coverage (depth) at variant site
MIN_ALT_COUNT = 3        # Minimum number of alt-supporting reads
MIN_VAF = 0.005          # Minimum Variant Allele Frequency (0.5% threshold)
MAX_STRAND_BIAS = 0.8    # Maximum strand bias ratio (0.0 perfectly balanced, 1.0 extreme bias)
MAX_NORMAL_VAF = 0.01    # Maximum Normal VAF (to exclude germline variants, i.e. normal sample must be ~0 VAF)
MIN_AF_DIFF = 0.01       # Minimum VAF(tumor) - VAF(normal) somatic signal strength

# ── LLM API Settings ───────────────────────────────────────────────────────────
# Standard OpenAI-compatible API configurations
LLM_API_KEY_ENV = "OPENAI_API_KEY"
LLM_BASE_URL = os.environ.get("OPENAI_BASE_URL", "https://api.openai.com/v1")
LLM_MODEL = os.environ.get("OPENAI_MODEL_NAME", "gpt-4o-mini")
LLM_TEMPERATURE = 0.0
LLM_MAX_TOKENS = 600

# Optimization & Cost Controls
LLM_BATCH_SIZE = 10      # Number of variants to process per batch (prompt consolidation)
LLM_MAX_WORKERS = 4      # Multi-threading workers for concurrent LLM requests
LLM_CACHE_ENABLED = True # Cache LLM calls locally to avoid repeated API bills

# ── Biological Annotation Fallbacks ────────────────────────────────────────────
# Used if the official COSMIC Cancer Gene Census (CGC) file is not provided.
KNOWN_ONCOGENES = {
    "EGFR", "KRAS", "BRAF", "PIK3CA", "ALK", "RET", "MET", "ERBB2", 
    "MAP2K1", "NRAS", "HRAS", "MYC", "CCND1", "FGFR1", "FGFR3", "JAK2"
}

KNOWN_TUMOR_SUPPRESSORS = {
    "TP53", "RB1", "STK11", "KEAP1", "PTEN", "APC", "BRCA1", "BRCA2", 
    "NF1", "SMAD4", "CDKN2A", "ARID1A", "VHL", "WT1", "ATM", "CREBBP",
    "KMT2D", "FAT1"
}

# Large, highly mutated passenger genes that often represent sequencing artifacts
# or passenger noise rather than drivers in ultra-low VAF settings.
ARTIFACT_PRONE_GENES = {
    "TTN", "MUC16", "RYR2", "OBSCN", "USH2A", "SYNE1", "NEB", "DST", 
    "LRP1B", "PCLO", "CSMD3", "ZFHX4"
}
