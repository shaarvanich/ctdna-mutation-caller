"""
annotation.py
-------------
Annotates variants with biological metadata, matching them against the COSMIC Cancer
Gene Census (CGC) and lists of oncogenes, tumor suppressors, and artifact-prone genes.
"""

import os
import logging
import pandas as pd

import config

logger = logging.getLogger(__name__)


class VariantAnnotator:
    """Annotates variants with biological and clinical cancer markers."""

    def __init__(self, cosmic_path: str = None):
        self.cosmic_path = cosmic_path if cosmic_path else os.path.join(config.DATA_DIR, "cancer_gene_census.csv")
        self.cosmic_genes = set()
        self._load_cosmic()

    def _load_cosmic(self):
        """Loads COSMIC Cancer Gene Census genes from local CSV file if it exists."""
        if os.path.exists(self.cosmic_path):
            try:
                cgc = pd.read_csv(self.cosmic_path, low_memory=False)
                # Find the gene symbol column (usually 'Gene Symbol' or contains 'gene')
                gene_col = next((c for c in cgc.columns 
                                 if "gene" in c.lower() and "symbol" in c.lower()), None)
                if gene_col is None:
                    gene_col = cgc.columns[0] # Fallback
                    
                self.cosmic_genes = set(cgc[gene_col].dropna().str.strip().str.upper())
                logger.info(f"Loaded {len(self.cosmic_genes):,} genes from COSMIC CGC CSV.")
            except Exception as e:
                logger.warning(f"Could not load COSMIC CGC file: {e}. Fallback to built-in drivers.")
        else:
            logger.info("COSMIC CGC CSV not found. Using built-in database fallback.")

    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotates variants with gene classification (driver, tumor suppressor, artifact-prone).
        
        Args:
            df (pd.DataFrame): Dataframe containing variants (with 'Hugo_Symbol').
            
        Returns:
            pd.DataFrame: Dataframe with annotation columns added.
        """
        logger.info("Annotating mutations with cancer genetics metadata...")
        out = df.copy()
        
        cosmic_hits = []
        gene_roles = []
        is_driver = []
        is_artifact_prone = []
        
        for idx, row in out.iterrows():
            gene = str(row.get("Hugo_Symbol", "")).strip().upper()
            
            # 1. Check COSMIC Census
            in_cosmic = gene in self.cosmic_genes
            
            # Fallback if no COSMIC census loaded: check if in our driver lists
            in_oncogene = gene in config.KNOWN_ONCOGENES
            in_tsg = gene in config.KNOWN_TUMOR_SUPPRESSORS
            in_artifact = gene in config.ARTIFACT_PRONE_GENES
            
            is_cgc_hit = 1 if (in_cosmic or in_oncogene or in_tsg) else 0
            cosmic_hits.append(is_cgc_hit)
            
            # 2. Determine biological role
            if in_oncogene:
                role = "oncogene"
                driver_flag = 1
                art_flag = 0
            elif in_tsg:
                role = "tumor_suppressor"
                driver_flag = 1
                art_flag = 0
            elif in_artifact:
                role = "artifact_prone"
                driver_flag = 0
                art_flag = 1
            else:
                role = "other"
                driver_flag = 1 if in_cosmic else 0
                art_flag = 0
                
            gene_roles.append(role)
            is_driver.append(driver_flag)
            is_artifact_prone.append(art_flag)
            
        out["cosmic_hit"] = cosmic_hits
        out["gene_role"] = gene_roles
        out["is_driver"] = is_driver
        out["is_artifact_prone"] = is_artifact_prone
        
        n_drivers = sum(out["is_driver"] == 1)
        n_artifacts = sum(out["is_artifact_prone"] == 1)
        
        logger.info(f"Annotations complete: {n_drivers} mutations labeled driver; {n_artifacts} in artifact-prone genes.")
        return out
