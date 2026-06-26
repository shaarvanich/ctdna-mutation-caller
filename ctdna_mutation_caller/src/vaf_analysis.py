"""
vaf_analysis.py
---------------
Analyzes variant allele frequencies (VAF), stratifies variants into clinical tiers
(ultra-low <1%, low 1-5%, medium/high >5%), and evaluates filter performance.
"""

import logging
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class VAFAnalyzer:
    """Analyzes and stratifies variants by VAF, evaluating filtering efficiency."""

    def __init__(self):
        pass

    def stratify_by_vaf(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Adds a 'vaf_tier' column to the DataFrame based on Variant Allele Frequency.
        
        Tiers:
            - '<1%'  : Ultra-low frequency ctDNA mutations
            - '1-5%' : Low-frequency ctDNA mutations
            - '>5%'  : Standard somatic mutations
        """
        df = df.copy()
        
        def _get_tier(vaf):
            if vaf < 0.01:
                return "<1%"
            elif vaf <= 0.05:
                return "1-5%"
            else:
                return ">5%"
                
        df["vaf_tier"] = df["VAF"].apply(_get_tier)
        return df

    def compute_tier_metrics(self, df: pd.DataFrame) -> dict:
        """
        Computes detailed filtering sensitivity and precision metrics stratified by VAF tiers.
        Requires a 'label' column in the dataset (representing GDC validation status) to compute metrics.
        """
        df = self.stratify_by_vaf(df)
        
        tiers = ["<1%", "1-5%", ">5%", "Overall"]
        metrics = {}
        
        logger.info("Analyzing filtering performance across VAF tiers...")
        
        # Check if 'label' column exists for performance evaluation
        has_label = "label" in df.columns
        
        for tier in tiers:
            if tier == "Overall":
                tier_df = df
            else:
                tier_df = df[df["vaf_tier"] == tier]
            total_count = len(tier_df)
            
            if total_count == 0:
                metrics[tier] = {
                    "total_count": 0,
                    "candidates_called": 0,
                    "calling_rate": 0.0
                }
                continue
                
            candidates_called = int(tier_df["is_candidate"].sum())
            calling_rate = candidates_called / total_count
            
            tier_metrics = {
                "total_count": total_count,
                "candidates_called": candidates_called,
                "calling_rate": calling_rate
            }
            
            if has_label:
                # True Positives: is_candidate=True & label=1
                tp = int(((tier_df["is_candidate"] == True) & (tier_df["label"] == 1)).sum())
                # False Positives: is_candidate=True & label=0
                fp = int(((tier_df["is_candidate"] == True) & (tier_df["label"] == 0)).sum())
                # True Negatives: is_candidate=False & label=0
                tn = int(((tier_df["is_candidate"] == False) & (tier_df["label"] == 0)).sum())
                # False Negatives: is_candidate=False & label=1
                fn = int(((tier_df["is_candidate"] == False) & (tier_df["label"] == 1)).sum())
                
                # Calculations
                sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
                precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
                f1_score = (2 * precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0.0
                
                # Advanced Evaluation Metrics
                balanced_accuracy = (sensitivity + specificity) / 2.0
                npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0
                
                mcc_denom = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
                mcc = (tp * tn - fp * fn) / mcc_denom if mcc_denom > 0 else 0.0
                
                tier_metrics.update({
                    "true_positives": tp,
                    "false_positives": fp,
                    "true_negatives": tn,
                    "false_negatives": fn,
                    "sensitivity": sensitivity,
                    "specificity": specificity,
                    "precision": precision,
                    "f1_score": f1_score,
                    "balanced_accuracy": balanced_accuracy,
                    "npv": npv,
                    "mcc": mcc
                })
                
                logger.info(
                    f"VAF Tier {tier:<6} | N={total_count:<4} | Sens: {sensitivity*100:5.1f}% | "
                    f"Prec: {precision*100:5.1f}% | F1: {f1_score:.3f} | BalAcc: {balanced_accuracy*100:5.1f}% | MCC: {mcc:.3f}"
                )
            else:
                logger.info(f"VAF Tier {tier:<6} | N={total_count:<4} | Calling Rate: {calling_rate*100:5.1f}%")
                
            metrics[tier] = tier_metrics
            
        return metrics
