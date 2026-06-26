"""
visualization.py
----------------
Generates publication-quality charts and plots for ctDNA somatic mutations,
including VAF distributions, TMB, gene frequencies, and VAF tier comparisons.
"""

import os
import logging
import gc
import matplotlib
matplotlib.use("Agg")  # Run headless for server/PyCharm terminal executions
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import config

logger = logging.getLogger(__name__)


class PipelineVisualizer:
    """Generates analytical visualizations for the ctDNA somatic caller pipeline."""

    def __init__(self, out_dir: str = config.PLOTS_DIR):
        self.out_dir = out_dir
        os.makedirs(self.out_dir, exist_ok=True)
        # Apply clean visual styles
        sns.set_theme(style="whitegrid")
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

        # Theme Colors: HSL-inspired premium palette
        self.colors = {
            "somatic_high": "#1D9E75",     # Emerald green
            "somatic_low": "#7F77DD",      # Soft purple
            "germline": "#4A90E2",         # Steel blue
            "artifact": "#D85A30",         # Warm coral
            "filtered": "#9B9B9B",         # Slate gray
            "oncogene": "#1D9E75",
            "tumor_suppressor": "#E04F5F", # Crimson
            "artifact_prone": "#D85A30",
            "other": "#7F77DD"
        }

    def plot_vaf_distribution(self, df: pd.DataFrame):
        """Plots the VAF distribution of variants, categorized by LLM classification."""
        logger.info("Generating Plot: VAF distribution...")
        fig, ax = plt.subplots(figsize=(10, 5))
        
        classes_to_plot = ["somatic_high_conf", "somatic_low_conf", "germline", "sequencing_artifact"]
        
        # Plot kernel density estimates or histograms for each class
        has_data = False
        for cls in classes_to_plot:
            subset = df[df["llm_classification"] == cls]
            if len(subset) > 0:
                has_data = True
                sns.histplot(
                    data=subset,
                    x="VAF",
                    bins=50,
                    element="step",
                    stat="density",
                    alpha=0.4,
                    color=self.colors.get(cls.split("_")[0] if "somatic" in cls else cls, "#333333"),
                    label=cls.replace("_", " ").title(),
                    ax=ax
                )
                
        if not has_data:
            logger.warning("No classified variants to plot for VAF distribution.")
            plt.close(fig)
            return

        ax.set_xlabel("Variant Allele Frequency (VAF)", fontsize=12)
        ax.set_ylabel("Density", fontsize=12)
        ax.set_title("VAF Distribution Stratified by LLM-Assisted Call Category", fontsize=14, fontweight="bold", pad=15)
        ax.set_xlim(0, max(0.25, df["VAF"].max() if len(df) > 0 else 0.25)) # Focus on low range
        ax.legend(title="LLM Call Class", frameon=True)
        fig.tight_layout()
        
        out_path = os.path.join(self.out_dir, "01_vaf_distribution.png")
        fig.savefig(out_path, dpi=100)
        plt.clf()
        plt.close(fig)
        gc.collect()
        logger.info(f"Saved: {out_path}")

    def plot_tumor_mutation_burden(self, df: pd.DataFrame):
        """Plots the somatic mutation count (TMB proxy) for each patient sample."""
        logger.info("Generating Plot: Tumor mutation burden...")
        
        # Filter for called somatic mutations (high and low confidence)
        somatic = df[df["llm_classification"].str.startswith("somatic", na=False)]
        
        if len(somatic) == 0:
            logger.warning("No somatic mutations detected. Skipping TMB plot.")
            return
            
        tmb_counts = somatic["Tumor_Sample_Barcode"].value_counts().head(25) # Top 25 samples
        
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.barplot(
            x=tmb_counts.index,
            y=tmb_counts.values,
            color="#5C6BC0",
            edgecolor="#3F51B5",
            ax=ax
        )
        
        ax.set_xticks(range(len(tmb_counts)))
        ax.set_xticklabels(tmb_counts.index, rotation=45, ha="right", fontsize=9)
        ax.set_xlabel("Tumor Sample Barcode", fontsize=12)
        ax.set_ylabel("Somatic Mutation Count (TMB Proxy)", fontsize=12)
        ax.set_title("Top 25 Samples by Somatic Mutation Burden", fontsize=14, fontweight="bold", pad=15)
        fig.tight_layout()
        
        out_path = os.path.join(self.out_dir, "02_tumor_mutation_burden.png")
        fig.savefig(out_path, dpi=100)
        plt.clf()
        plt.close(fig)
        gc.collect()
        logger.info(f"Saved: {out_path}")

    def plot_gene_frequency(self, df: pd.DataFrame):
        """Plots mutation frequency for the top mutated genes, colored by biological driver role."""
        logger.info("Generating Plot: Gene mutation frequency...")
        
        # Filter for somatic mutations
        somatic = df[df["llm_classification"].str.startswith("somatic", na=False)]
        if len(somatic) == 0:
            logger.warning("No somatic mutations detected. Skipping gene frequency plot.")
            return
            
        top_genes = somatic["Hugo_Symbol"].value_counts().head(15)
        
        # Create a mapping of genes to roles in top_genes
        gene_roles = []
        for g in top_genes.index:
            role = df[df["Hugo_Symbol"] == g]["gene_role"].values[0]
            gene_roles.append(role)
            
        temp_df = pd.DataFrame({
            "Gene": top_genes.index,
            "Count": top_genes.values,
            "Role": gene_roles
        })
        
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.barplot(
            data=temp_df,
            y="Gene",
            x="Count",
            hue="Role",
            palette=self.colors,
            dodge=False,
            ax=ax
        )
        
        ax.set_xlabel("Number of Somatic Mutations", fontsize=12)
        ax.set_ylabel("Gene Symbol", fontsize=12)
        ax.set_title("Top 15 Most Frequently Mutated Genes in Cohort", fontsize=14, fontweight="bold", pad=15)
        ax.legend(title="Biological Role", frameon=True)
        fig.tight_layout()
        
        out_path = os.path.join(self.out_dir, "03_gene_frequency.png")
        fig.savefig(out_path, dpi=100)
        plt.clf()
        plt.close(fig)
        gc.collect()
        logger.info(f"Saved: {out_path}")

    def plot_low_vaf_histogram(self, df: pd.DataFrame):
        """Plots a detailed histogram of called somatic mutations below 5% VAF (clinical ctDNA focus)."""
        logger.info("Generating Plot: Low VAF histogram...")
        
        # Filter for low-VAF somatic calls (VAF < 5%)
        low_vaf_somatic = df[
            (df["llm_classification"].str.startswith("somatic", na=False)) & 
            (df["VAF"] < 0.05)
        ]
        
        if len(low_vaf_somatic) == 0:
            logger.warning("No low-VAF (<5%) somatic mutations found to plot.")
            return
            
        fig, ax = plt.subplots(figsize=(10, 5))
        
        sns.histplot(
            data=low_vaf_somatic,
            x="VAF",
            hue="llm_classification",
            bins=30,
            multiple="stack",
            palette={
                "somatic_high_conf": self.colors["somatic_high"],
                "somatic_low_conf": self.colors["somatic_low"]
            },
            alpha=0.8,
            ax=ax
        )
        
        ax.set_xlabel("Variant Allele Frequency (VAF)", fontsize=12)
        ax.set_ylabel("Mutation Count", fontsize=12)
        ax.set_title("Low-VAF (VAF < 5%) Called Somatic Mutations", fontsize=14, fontweight="bold", pad=15)
        ax.set_xlim(0, 0.05)
        
        # Adjust legend labels
        legend = ax.get_legend()
        if legend:
            legend.set_title("Call Category")
            for t in legend.get_texts():
                t.set_text(t.get_text().replace("_", " ").title())
                
        fig.tight_layout()
        out_path = os.path.join(self.out_dir, "04_low_vaf_histogram.png")
        fig.savefig(out_path, dpi=100)
        plt.clf()
        plt.close(fig)
        gc.collect()
        logger.info(f"Saved: {out_path}")

    def plot_sensitivity_tiers(self, tier_metrics: dict):
        """Generates a sensitivity and calling rate chart across the three VAF tiers."""
        logger.info("Generating Plot: Sensitivity tiers chart...")
        
        tiers = ["<1%", "1-5%", ">5%"]
        calling_rates = [tier_metrics[t]["calling_rate"] * 100 for t in tiers]
        
        has_sens = "sensitivity" in tier_metrics["<1%"]
        
        fig, ax = plt.subplots(figsize=(8, 5))
        x_indexes = np.arange(len(tiers))
        width = 0.35
        
        if has_sens:
            sensitivities = [tier_metrics[t]["sensitivity"] * 100 for t in tiers]
            rects1 = ax.bar(x_indexes - width/2, sensitivities, width, label="Sensitivity (Recall)", color="#1D9E75")
            rects2 = ax.bar(x_indexes + width/2, calling_rates, width, label="Pipeline Calling Rate", color="#7F77DD")
            
            # Label bars
            ax.bar_label(rects1, fmt='%.1f%%', padding=3, fontsize=9)
            ax.bar_label(rects2, fmt='%.1f%%', padding=3, fontsize=9)
        else:
            rects = ax.bar(x_indexes, calling_rates, width * 1.5, label="Pipeline Calling Rate", color="#7F77DD")
            ax.bar_label(rects, fmt='%.1f%%', padding=3, fontsize=10)
            
        ax.set_xticks(x_indexes)
        ax.set_xticklabels(tiers)
        ax.set_ylim(0, 110)
        ax.set_xlabel("VAF Tier", fontsize=12)
        ax.set_ylabel("Percentage (%)", fontsize=12)
        ax.set_title("Detection Performance and Calling Rates by VAF Tier", fontsize=14, fontweight="bold", pad=15)
        ax.legend(frameon=True)
        fig.tight_layout()
        
        out_path = os.path.join(self.out_dir, "05_vaf_tier_metrics.png")
        fig.savefig(out_path, dpi=100)
        plt.clf()
        plt.close(fig)
        gc.collect()
        logger.info(f"Saved: {out_path}")

    def plot_confusion_matrix(self, df: pd.DataFrame):
        """Generates a 2x2 confusion matrix heatmap displaying called vs filtered counts."""
        if "label" not in df.columns:
            logger.warning("No ground-truth 'label' column found. Skipping confusion matrix plot.")
            return
            
        logger.info("Generating Plot: Confusion matrix...")
        
        # True Positives: is_candidate=True & label=1
        tp = int(((df["is_candidate"] == True) & (df["label"] == 1)).sum())
        # False Positives: is_candidate=True & label=0
        fp = int(((df["is_candidate"] == True) & (df["label"] == 0)).sum())
        # True Negatives: is_candidate=False & label=0
        tn = int(((df["is_candidate"] == False) & (df["label"] == 0)).sum())
        # False Negatives: is_candidate=False & label=1
        fn = int(((df["is_candidate"] == False) & (df["label"] == 1)).sum())
        
        cm = np.array([[tn, fp],
                       [fn, tp]])
                       
        total = len(df)
        labels = np.array([
            [f"True Neg (Noise)\n{tn}\n({tn/total*100:.1f}%)", f"False Pos (Leakage)\n{fp}\n({fp/total*100:.1f}%)"],
            [f"False Neg (Missed)\n{fn}\n({fn/total*100:.1f}%)", f"True Pos (Somatic)\n{tp}\n({tp/total*100:.1f}%)"]
        ])
        
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.heatmap(
            cm, 
            annot=labels, 
            fmt="", 
            cmap="Blues", 
            cbar=False, 
            linewidths=1.5, 
            linecolor="#DDDDDD", 
            xticklabels=["Filtered (Neg)", "Called (Pos)"], 
            yticklabels=["Actual Noise (0)", "Actual Somatic (1)"],
            annot_kws={"size": 11, "weight": "bold"},
            ax=ax
        )
        
        ax.set_ylabel("Ground Truth (GDC Label)", fontsize=12, fontweight="bold")
        ax.set_xlabel("Pipeline Call Decision", fontsize=12, fontweight="bold")
        ax.set_title("Pipeline Somatic Mutation Confusion Matrix", fontsize=13, fontweight="bold", pad=15)
        fig.tight_layout()
        
        out_path = os.path.join(self.out_dir, "06_confusion_matrix.png")
        fig.savefig(out_path, dpi=100)
        plt.clf()
        plt.close(fig)
        gc.collect()
        logger.info(f"Saved: {out_path}")

    def plot_roc_pr_curves(self, df: pd.DataFrame):
        """Generates ROC and Precision-Recall curves using numpy from scratch."""
        if "label" not in df.columns or "llm_confidence" not in df.columns:
            logger.warning("Missing 'label' or 'llm_confidence' columns. Skipping curves plot.")
            return
            
        logger.info("Generating Plot: ROC and PR curves...")
        
        # Prediction score = candidate flag * llm confidence
        y_scores = (df["is_candidate"].astype(float) * df["llm_confidence"]).values
        y_true = df["label"].values
        
        # Sort descending
        desc_score_indices = np.argsort(y_scores)[::-1]
        y_scores = y_scores[desc_score_indices]
        y_true = y_true[desc_score_indices]
        
        # Find unique threshold locations
        distinct_value_indices = np.where(np.diff(y_scores))[0]
        threshold_idxs = np.r_[distinct_value_indices, y_true.size - 1]
        
        # Calculate TPs and FPs at each threshold
        tps = np.cumsum(y_true)[threshold_idxs]
        fps = 1 + threshold_idxs - tps
        
        # Add point for threshold = max_score + epsilon (where everything is negative)
        tps = np.r_[0, tps]
        fps = np.r_[0, fps]
        
        n_pos = y_true.sum()
        n_neg = y_true.size - n_pos
        
        tpr = tps / n_pos if n_pos > 0 else np.zeros_like(tps)
        fpr = fps / n_neg if n_neg > 0 else np.zeros_like(fps)
        
        # Precision-Recall
        with np.errstate(divide="ignore", invalid="ignore"):
            precision = tps / (tps + fps)
        precision[0] = 1.0 # 0/0 edge case
        recall = tpr
        
        # Custom AUC calculation helper using trapezoidal rule (safe for NumPy 2.0)
        def _calculate_auc(x_vals, y_vals):
            sort_idx = np.argsort(x_vals)
            x_s = x_vals[sort_idx]
            y_s = y_vals[sort_idx]
            return np.sum((x_s[1:] - x_s[:-1]) * (y_s[1:] + y_s[:-1]) / 2.0)
            
        roc_auc = _calculate_auc(fpr, tpr)
        pr_auc = _calculate_auc(recall, precision)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Panel 1: ROC Curve
        ax1.plot(fpr, tpr, color="#1D9E75", lw=2, label=f"ROC Curve (AUC = {roc_auc:.3f})")
        ax1.plot([0, 1], [0, 1], color="#9B9B9B", lw=1, linestyle="--")
        ax1.set_xlim([0.0, 1.0])
        ax1.set_ylim([0.0, 1.05])
        ax1.set_xlabel("False Positive Rate (FPR)", fontsize=11)
        ax1.set_ylabel("True Positive Rate (TPR)", fontsize=11)
        ax1.set_title("Receiver Operating Characteristic (ROC)", fontsize=12, fontweight="bold")
        ax1.legend(loc="lower right", frameon=True)
        
        # Panel 2: PR Curve
        ax2.plot(recall, precision, color="#7F77DD", lw=2, label=f"PR Curve (AUC = {pr_auc:.3f})")
        baseline = n_pos / len(df) if len(df) > 0 else 0.5
        ax2.plot([0, 1], [baseline, baseline], color="#9B9B9B", lw=1, linestyle="--", label=f"Baseline ({baseline:.2f})")
        ax2.set_xlim([0.0, 1.0])
        ax2.set_ylim([0.0, 1.05])
        ax2.set_xlabel("Recall (Sensitivity)", fontsize=11)
        ax2.set_ylabel("Precision (PPV)", fontsize=11)
        ax2.set_title("Precision-Recall (PR) Curve", fontsize=12, fontweight="bold")
        ax2.legend(loc="lower left", frameon=True)
        
        fig.tight_layout()
        out_path = os.path.join(self.out_dir, "07_performance_curves.png")
        fig.savefig(out_path, dpi=100)
        plt.clf()
        plt.close(fig)
        gc.collect()
        logger.info(f"Saved: {out_path}")

    def generate_all_plots(self, df: pd.DataFrame, tier_metrics: dict):
        """Helper method to run the entire plotting suite."""
        logger.info("Generating all pipeline plots...")
        self.plot_vaf_distribution(df)
        self.plot_tumor_mutation_burden(df)
        self.plot_gene_frequency(df)
        self.plot_low_vaf_histogram(df)
        self.plot_sensitivity_tiers(tier_metrics)
        self.plot_confusion_matrix(df)
        self.plot_roc_pr_curves(df)
        logger.info(f"All visualizations successfully written to {self.out_dir}/")
