"""
preprocessing.py
----------------
Preprocesses input raw datasets (MAF, VCF, or CSV) and engineers features
specifically relevant for low-VAF ctDNA somatic mutation detection.
"""

import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class Preprocessor:
    """Preprocesses mutation data and engineers statistical features for caller/LLM."""

    def __init__(self):
        pass

    def preprocess(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Runs the full preprocessing pipeline on the input DataFrame.
        
        Args:
            df (pd.DataFrame): Raw parsed mutations.
            
        Returns:
            pd.DataFrame: Feature-engineered and standardized DataFrame.
        """
        logger.info("Starting preprocessing and feature engineering...")
        
        # 1. Handle missing values for critical read count columns
        df = df.copy()
        
        # Ensure depth columns are numeric and filled
        df["t_depth"] = pd.to_numeric(df.get("t_depth", pd.Series(np.nan, index=df.index)), errors="coerce").fillna(0)
        df["t_alt_count"] = pd.to_numeric(df.get("t_alt_count", pd.Series(np.nan, index=df.index)), errors="coerce").fillna(0)
        df["t_ref_count"] = pd.to_numeric(df.get("t_ref_count", pd.Series(np.nan, index=df.index)), errors="coerce")
        
        # Reconstruct ref count if missing
        df["t_ref_count"] = df["t_ref_count"].fillna(df["t_depth"] - df["t_alt_count"])
        
        # Handle matched normal columns
        df["n_depth"] = pd.to_numeric(df.get("n_depth", pd.Series(np.nan, index=df.index)), errors="coerce").fillna(0)
        df["n_alt_count"] = pd.to_numeric(df.get("n_alt_count", pd.Series(np.nan, index=df.index)), errors="coerce").fillna(0)
        df["n_ref_count"] = pd.to_numeric(df.get("n_ref_count", pd.Series(np.nan, index=df.index)), errors="coerce").fillna(df["n_depth"] - df["n_alt_count"])

        # 2. Calculate VAF and somatic allele frequency difference
        logger.debug("Engineering Variant Allele Frequency (VAF) features...")
        
        # Tumor VAF
        df["VAF"] = (df["t_alt_count"] / df["t_depth"].replace(0, np.nan)).fillna(0.0).clip(0.0, 1.0)
        
        # Normal VAF
        n_vaf = (df["n_alt_count"] / df["n_depth"].replace(0, np.nan)).fillna(0.0).clip(0.0, 1.0)
        df["normal_VAF"] = n_vaf
        
        # Somatic signal strength (af_diff)
        df["af_diff"] = (df["VAF"] - n_vaf).clip(-1.0, 1.0)

        # 3. Calculate Strand Bias
        # strand_bias = |fwd_alt_fraction - 0.5| * 2  -> 0.0 is balanced, 1.0 is extreme bias (all alt on one strand)
        logger.debug("Engineering strand bias features...")
        has_strand_counts = "t_alt_count_forward" in df.columns and "t_alt_count_reverse" in df.columns
        
        if has_strand_counts:
            fwd = pd.to_numeric(df["t_alt_count_forward"], errors="coerce").fillna(0.0)
            rev = pd.to_numeric(df["t_alt_count_reverse"], errors="coerce").fillna(0.0)
            total_alt = fwd + rev
            
            # Use 0.5 fraction (neutral) if no alt reads are present to avoid Division by Zero
            fwd_frac = (fwd / total_alt.replace(0, np.nan)).fillna(0.5)
            df["strand_bias"] = (np.abs(fwd_frac - 0.5) * 2.0).clip(0.0, 1.0)
        else:
            # Fallback if strand counts are not present
            df["strand_bias"] = df.get("strand_bias", 0.0)
            df["strand_bias"] = pd.to_numeric(df["strand_bias"], errors="coerce").fillna(0.0)

        # 4. Categorical Encodings (Variant Type & Classification)
        logger.debug("Encoding categorical mutation properties...")
        
        # Variant Type (SNP, DNP, TNP, INS, DEL)
        vtype_map = {"SNP": 0, "DNP": 1, "TNP": 2, "INS": 3, "DEL": 4}
        if "Variant_Type" in df.columns:
            df["variant_type_enc"] = df["Variant_Type"].map(vtype_map).fillna(0).astype(int)
        else:
            df["variant_type_enc"] = df.get("variant_type_enc", 0).astype(int)

        # Classification mapping
        if "Variant_Classification" in df.columns:
            missense_classes = {"Missense_Mutation"}
            nonsense_classes = {
                "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Del", 
                "Frame_Shift_Ins", "Nonstop_Mutation", "Translation_Start_Site"
            }
            silent_classes = {"Silent"}

            df["is_missense"] = df["Variant_Classification"].isin(missense_classes).astype(int)
            df["is_nonsense"] = df["Variant_Classification"].isin(nonsense_classes).astype(int)
            df["is_silent"] = df["Variant_Classification"].isin(silent_classes).astype(int)
        else:
            df["is_missense"] = df.get("is_missense", 0).astype(int)
            df["is_nonsense"] = df.get("is_nonsense", 0).astype(int)
            df["is_silent"] = df.get("is_silent", 0).astype(int)

        # 5. VEP Functional Impact
        impact_map = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
        raw_impact_col = "IMPACT" if "IMPACT" in df.columns else ("impact" if "impact" in df.columns else None)
        
        if raw_impact_col:
            df["impact_enc"] = df[raw_impact_col].astype(str).str.upper().map(impact_map).fillna(0).astype(int)
        else:
            df["impact_enc"] = df.get("impact_enc", 0).astype(int)

        # 6. PolyPhen/SIFT proxies
        if "PolyPhen" in df.columns:
            def _encode_polyphen(val):
                v = str(val).lower()
                if "probably_damaging" in v or "damaging" in v:
                    return 2
                if "possibly_damaging" in v:
                    return 1
                return 0
            df["polyphen_enc"] = df["PolyPhen"].apply(_encode_polyphen).astype(int)
        else:
            df["polyphen_enc"] = df.get("polyphen_enc", 0).astype(int)

        if "SIFT" in df.columns:
            df["sift_enc"] = df["SIFT"].apply(
                lambda v: 1 if "deleterious" in str(v).lower() else 0
            ).astype(int)
        else:
            df["sift_enc"] = df.get("sift_enc", 0).astype(int)

        # 7. dbSNP Novelty Flag
        if "dbSNP_RS" in df.columns:
            df["dbsnp_novel"] = (
                (df["dbSNP_RS"].isna()) | 
                (df["dbSNP_RS"].astype(str).str.lower() == "novel") |
                (df["dbSNP_RS"].astype(str) == ".")
            ).astype(int)
        else:
            df["dbsnp_novel"] = df.get("dbsnp_novel", 1).astype(int)

        # 8. Ensure essential metadata fields exist
        meta_fields = {
            "Hugo_Symbol": "GENE_UNK",
            "Variant_Classification": "Unknown",
            "Chromosome": "chrUn",
            "Start_Position": 0,
            "Reference_Allele": "N",
            "Tumor_Seq_Allele2": "N",
            "Tumor_Sample_Barcode": "SAMPLE_UNK",
            "FILTER_raw": "Unknown"
        }
        for field, default in meta_fields.items():
            if field not in df.columns:
                df[field] = default
                
        logger.info(f"Preprocessing completed. Engineered {df.shape[1]} columns for {len(df):,} mutations.")
        return df
