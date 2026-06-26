"""
variant_filter.py
-----------------
Implements statistical and biological filters to identify candidate somatic mutations
in ctDNA. Includes a rescue mechanism for ultra-low VAF variants in known driver genes.
"""

import logging
import pandas as pd
import numpy as np

import config

logger = logging.getLogger(__name__)


class VariantFilter:
    """Applies rule-based statistical and biological filters to variants."""

    def __init__(self, 
                 min_depth: int = config.MIN_DEPTH, 
                 min_alt_count: int = config.MIN_ALT_COUNT, 
                 min_vaf: float = config.MIN_VAF, 
                 max_strand_bias: float = config.MAX_STRAND_BIAS, 
                 max_normal_vaf: float = config.MAX_NORMAL_VAF, 
                 min_af_diff: float = config.MIN_AF_DIFF):
        
        self.min_depth = min_depth
        self.min_alt_count = min_alt_count
        self.min_vaf = min_vaf
        self.max_strand_bias = max_strand_bias
        self.max_normal_vaf = max_normal_vaf
        self.min_af_diff = min_af_diff

    def filter_variants(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Applies sequential filters and tags variants with detailed status logs.
        Supports a "Driver Rescue Heuristic" for ultra-low VAF mutations (0.1% to 0.5%)
        found in known driver genes.
        
        Args:
            df (pd.DataFrame): Preprocessed DataFrame.
            
        Returns:
            pd.DataFrame: DataFrame with filtering columns appended ('filter_status', 'is_candidate').
        """
        logger.info("Applying statistical and biological filters to variants...")
        
        out = df.copy()
        
        filter_status = []
        is_candidate = []
        rescue_count = 0
        
        for idx, row in out.iterrows():
            failures = []
            
            # Extract features
            depth = row.get("t_depth", 0)
            alt_count = row.get("t_alt_count", 0)
            vaf = row.get("VAF", 0.0)
            normal_vaf = row.get("normal_VAF", 0.0)
            af_diff = row.get("af_diff", 0.0)
            sb = row.get("strand_bias", 0.0)
            gene = str(row.get("Hugo_Symbol", "")).upper()
            cosmic_hit = row.get("cosmic_hit", 0)
            
            # Rule 1: Read Depth Coverage
            if depth < self.min_depth:
                failures.append("low_depth")
                
            # Rule 2: Minimum Alt Supporting Reads
            if alt_count < self.min_alt_count:
                failures.append("low_alt_count")
                
            # Rule 3: Strand Bias Heuristic (reject artifacts)
            if sb > self.max_strand_bias:
                failures.append("strand_bias_artifact")
                
            # Rule 4: Germline Leakage / Matched Normal checks
            if normal_vaf > self.max_normal_vaf or af_diff < self.min_af_diff:
                failures.append("germline_leakage")
            
            # Check VAF threshold (Rule 5)
            vaf_failed = vaf < self.min_vaf
            
            # Check if this gene is a known driver for the Rescue Heuristic
            is_driver = (
                gene in config.KNOWN_ONCOGENES or 
                gene in config.KNOWN_TUMOR_SUPPRESSORS or 
                cosmic_hit == 1
            )
            
            # Rescue Heuristic:
            # If the variant passed all other quality rules (depth, alt, strand bias, germline)
            # but failed the VAF threshold (e.g., VAF is ultra-low, say 0.1% to 0.5%),
            # AND it is in a known cancer driver gene, we rescue it for LLM consideration.
            is_rescued = False
            if len(failures) == 0 and vaf_failed:
                if vaf >= 0.001 and is_driver: # Ultra-low VAF threshold
                    is_rescued = True
                    rescue_count += 1
                else:
                    failures.append("low_vaf")
            elif vaf_failed:
                failures.append("low_vaf")
                
            # Compile results
            if is_rescued:
                status = "RESCUED_DRIVER"
                candidate = True
            elif len(failures) == 0:
                status = "PASS"
                candidate = True
            else:
                status = ";".join(failures)
                candidate = False
                
            filter_status.append(status)
            is_candidate.append(candidate)
            
        out["filter_status"] = filter_status
        out["is_candidate"] = is_candidate
        
        n_total = len(out)
        n_pass = sum(out["filter_status"] == "PASS")
        n_rescued = sum(out["filter_status"] == "RESCUED_DRIVER")
        n_filtered = n_total - n_pass - n_rescued
        
        logger.info(f"Filtering complete summary:")
        logger.info(f"  Total input variants   : {n_total:,}")
        logger.info(f"  Passed filters (PASS)  : {n_pass:,} ({n_pass/n_total*100:.1f}%)")
        logger.info(f"  Rescued drivers        : {n_rescued:,} ({n_rescued/n_total*100:.1f}%)")
        logger.info(f"  Filtered out (Noise)   : {n_filtered:,} ({n_filtered/n_total*100:.1f}%)")
        
        return out
