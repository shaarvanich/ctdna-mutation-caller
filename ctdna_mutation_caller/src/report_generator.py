"""
report_generator.py
--------------------
Generates CSV spreadsheets of called mutations, structured JSON stats, 
and premium clinical interpretation Markdown reports.
"""

import os
import json
import logging
from datetime import datetime
import pandas as pd

import config

logger = logging.getLogger(__name__)


class ReportGenerator:
    """Compiles and exports pipeline calling outputs and clinical reports."""

    def __init__(self, out_dir: str = config.REPORTS_DIR):
        self.out_dir = out_dir
        os.makedirs(self.out_dir, exist_ok=True)

    def export_csv_report(self, df: pd.DataFrame) -> str:
        """
        Exports a filtered CSV of called somatic mutations (high & low confidence).
        """
        logger.info("Exporting somatic mutation list to CSV...")
        
        # Filter for called somatic mutations
        somatic_df = df[df["llm_classification"].str.startswith("somatic", na=False)].copy()
        
        # Select and order key columns
        export_cols = [
            "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
            "t_depth", "t_alt_count", "VAF", "normal_VAF", "af_diff", "strand_bias", 
            "cosmic_hit", "gene_role", "llm_classification", "llm_confidence", "llm_relevance", 
            "llm_explanation"
        ]
        
        # Keep columns that are actually present
        present_cols = [c for c in export_cols if c in somatic_df.columns]
        somatic_df = somatic_df[present_cols]
        
        # Sort by gene and VAF descending
        if "VAF" in somatic_df.columns:
            somatic_df = somatic_df.sort_values(by=["Hugo_Symbol", "VAF"], ascending=[True, False])
            
        out_path = os.path.join(self.out_dir, "somatic_mutations_called.csv")
        somatic_df.to_csv(out_path, index=False)
        logger.info(f"Saved: {out_path} (Called somatic count: {len(somatic_df):,})")
        return out_path

    def export_json_summary(self, df: pd.DataFrame, tier_metrics: dict, source_file: str) -> str:
        """
        Exports structured summary metrics of the pipeline run to a JSON file.
        """
        logger.info("Exporting run summary stats to JSON...")
        
        n_total = len(df)
        n_candidates = int(df["is_candidate"].sum())
        n_pass = int((df["filter_status"] == "PASS").sum())
        n_rescued = int((df["filter_status"] == "RESCUED_DRIVER").sum())
        
        # Class distributions
        classes = df["llm_classification"].value_counts().to_dict()
        
        # Top mutated genes in called somatic mutations
        somatic = df[df["llm_classification"].str.startswith("somatic", na=False)]
        top_genes = somatic["Hugo_Symbol"].value_counts().head(10).to_dict() if len(somatic) > 0 else {}
        
        summary = {
            "metadata": {
                "timestamp": datetime.now().isoformat(),
                "input_source": os.path.basename(source_file),
                "pipeline_version": "1.0.0-LLM-Assisted",
                "rules": {
                    "min_depth": config.MIN_DEPTH,
                    "min_alt_count": config.MIN_ALT_COUNT,
                    "min_vaf": config.MIN_VAF,
                    "max_strand_bias": config.MAX_STRAND_BIAS,
                    "max_normal_vaf": config.MAX_NORMAL_VAF
                }
            },
            "cohort_statistics": {
                "total_variants_analyzed": n_total,
                "variants_passing_rules_as_candidates": n_candidates,
                "rule_pass_count": n_pass,
                "driver_rescue_count": n_rescued,
                "candidate_vaf_distribution": {
                    "min_vaf": float(somatic["VAF"].min()) if len(somatic) > 0 else 0.0,
                    "max_vaf": float(somatic["VAF"].max()) if len(somatic) > 0 else 0.0,
                    "mean_vaf": float(somatic["VAF"].mean()) if len(somatic) > 0 else 0.0
                }
            },
            "llm_classifications": {
                k: int(v) for k, v in classes.items()
            },
            "vaf_tier_breakdown": tier_metrics,
            "top_mutated_somatic_genes": top_genes
        }
        
        out_path = os.path.join(self.out_dir, "pipeline_run_summary.json")
        with open(out_path, "w") as f:
            json.dump(summary, f, indent=4)
        logger.info(f"Saved: {out_path}")
        return out_path

    def export_clinical_report(self, df: pd.DataFrame, source_file: str) -> str:
        """
        Generates a premium Markdown Clinical Interpretation Report.
        """
        logger.info("Generating Clinical Interpretation Report...")
        
        somatic = df[df["llm_classification"].str.startswith("somatic", na=False)].copy()
        somatic_high = df[df["llm_classification"] == "somatic_high_conf"]
        somatic_low = df[df["llm_classification"] == "somatic_low_conf"]
        germline = df[df["llm_classification"] == "germline"]
        artifacts = df[df["llm_classification"] == "sequencing_artifact"]
        
        # Sort somatic mutations by clinical relevance and VAF
        relevance_order = {"pathogenic": 0, "uncertain_significance": 1, "benign": 2}
        somatic["rel_order"] = somatic["llm_relevance"].map(relevance_order).fillna(3)
        somatic = somatic.sort_values(by=["rel_order", "VAF"], ascending=[True, False])
        
        # Exec stats
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        report_md = []
        report_md.append(f"# Clinical Interpretation Report: ctDNA Somatic Mutation Call")
        report_md.append(f"**Pipeline**: LLM-Assisted ctDNA Mutation Caller for Low-VAF Detection")
        report_md.append(f"**Date Generated**: {timestamp}")
        report_md.append(f"**Dataset Analysed**: `{os.path.basename(source_file)}` (N={len(df):,})")
        report_md.append(f"\n---\n")
        
        # 1. Executive Summary
        report_md.append(f"## 1. Executive Summary")
        report_md.append(
            f"This analysis uses a hybrid rule-based statistical filter and large language model (LLM) reasoning "
            f"to identify somatic mutations in circulating tumor DNA (ctDNA). A total of **{len(df):,}** variant "
            f"calls were processed. Rule-based filters flagged **{len(df) - int(df['is_candidate'].sum()):,}** variants as "
            f"obvious noise or germline. The remaining **{int(df['is_candidate'].sum()):,}** candidates were processed using "
            f"LLM chain-of-thought analysis.\n"
        )
        
        report_md.append(f"### Mutation Classification Summary")
        report_md.append(f"| Call Category | Count | Description |")
        report_md.append(f"| :--- | :---: | :--- |")
        report_md.append(f"| **Somatic (High Confidence)** | {len(somatic_high)} | High-confidence somatic calls, potential therapeutic targets. |")
        report_md.append(f"| **Somatic (Low Confidence)** | {len(somatic_low)} | Passenger mutations or variants in long passenger genes. |")
        report_md.append(f"| **Germline Leakage** | {len(germline)} | Inherited variants detected in matched normal sample. |")
        report_md.append(f"| **Sequencing Artifacts** | {len(artifacts)} | PCR errors, chemistry noise, or strand-biased anomalies. |")
        report_md.append(f"| *Filtered out by Rules* | {len(df) - len(somatic_high) - len(somatic_low) - len(germline) - len(artifacts)} | Failed depth, alt count, or other hard filters. |")
        
        # 2. Pipeline Performance Evaluation
        from vaf_analysis import VAFAnalyzer
        analyzer = VAFAnalyzer()
        tier_metrics = analyzer.compute_tier_metrics(df)
        has_label = "label" in df.columns

        if has_label:
            report_md.append(f"\n## 2. Pipeline Performance Evaluation")
            report_md.append(
                "Because a ground-truth dataset label was available, we evaluated the pipeline's detection accuracy, "
                "stratified by Variant Allele Frequency (VAF) tiers. Below is a comprehensive diagnostic metrics table:"
            )
            report_md.append(
                "\n| VAF Tier | Total Count | Called Somatic | Calling Rate | Sensitivity (Recall) | Specificity | Precision (PPV) | NPV | F1 Score | Balanced Acc | MCC |"
            )
            report_md.append(
                "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |"
            )
            for tier in ["<1%", "1-5%", ">5%", "Overall"]:
                m = tier_metrics[tier]
                t_count = m["total_count"]
                called = m["candidates_called"]
                crate = f"{m['calling_rate']*100:.1f}%"
                
                sens = f"{m.get('sensitivity', 0.0)*100:.1f}%" if "sensitivity" in m else "N/A"
                spec = f"{m.get('specificity', 0.0)*100:.1f}%" if "specificity" in m else "N/A"
                prec = f"{m.get('precision', 0.0)*100:.1f}%" if "precision" in m else "N/A"
                npv = f"{m.get('npv', 0.0)*100:.1f}%" if "npv" in m else "N/A"
                f1 = f"{m.get('f1_score', 0.0):.3f}" if "f1_score" in m else "N/A"
                bal_acc = f"{m.get('balanced_accuracy', 0.0)*100:.1f}%" if "balanced_accuracy" in m else "N/A"
                mcc = f"{m.get('mcc', 0.0):.3f}" if "mcc" in m else "N/A"
                
                report_md.append(
                    f"| **{tier}** | {t_count:,} | {called:,} | {crate} | {sens} | {spec} | {prec} | {npv} | {f1} | {bal_acc} | {mcc} |"
                )
        else:
            report_md.append(f"\n## 2. Pipeline Calling Statistics")
            report_md.append(
                "Below are the pipeline's somatic calling rates stratified by Variant Allele Frequency (VAF) tiers:"
            )
            report_md.append(
                "\n| VAF Tier | Total Count | Called Somatic | Calling Rate |"
            )
            report_md.append(
                "| :--- | :---: | :---: | :---: |"
            )
            for tier in ["<1%", "1-5%", ">5%", "Overall"]:
                m = tier_metrics[tier]
                t_count = m["total_count"]
                called = m["candidates_called"]
                crate = f"{m['calling_rate']*100:.1f}%"
                report_md.append(f"| **{tier}** | {t_count:,} | {called:,} | {crate} |")
                
        # 3. Key Findings
        report_md.append(f"\n## 3. Key Somatic Driver Findings")
        report_md.append(
            "Below are the identified somatic mutations, prioritized by clinical pathogenicity and allele fraction. "
            "Special emphasis is placed on low-VAF mutations (<5%) salvaged by driver gene heuristics."
        )
        
        if len(somatic) == 0:
            report_md.append("\n> [!NOTE]\n> No somatic mutations were called in this cohort.")
        else:
            report_md.append("\n| Gene | Coord | Change | VAF | Normal VAF | Class | Relevance | Action |")
            report_md.append("| :--- | :--- | :---: | :---: | :---: | :--- | :--- | :--- |")
            for idx, row in somatic.iterrows():
                coord = f"Chr{row['Chromosome']}:{int(row['Start_Position'])}"
                ref = row.get("Reference_Allele", "N")
                alt = row.get("Tumor_Seq_Allele2", "N")
                change = f"{ref}→{alt}"
                vaf_str = f"{row['VAF']*100:.2f}%"
                nvaf_str = f"{row['normal_VAF']*100:.2f}%" if row['normal_VAF'] > 0 else "0.0%"
                cls_clean = row['llm_classification'].replace("_", " ").title()
                rel_clean = row['llm_relevance'].replace("_", " ").upper()
                action_clean = row['llm_action'].upper()
                
                # Highlight pathogenic variants
                gene_str = f"**{row['Hugo_Symbol']}**" if rel_clean == "PATHOGENIC" else row['Hugo_Symbol']
                
                report_md.append(f"| {gene_str} | {coord} | {change} | {vaf_str} | {nvaf_str} | {cls_clean} | {rel_clean} | {action_clean} |")
                
        # 4. Pathologist Detailed Explanations
        report_md.append(f"\n## 4. Bioinformatic & Clinical Explanations")
        report_md.append("Detailed LLM explanations for selected actionable somatic mutations:")
        
        # Focus on pathogenic driver genes
        drivers_reported = somatic[somatic["is_driver"] == 1].head(5)
        
        if len(drivers_reported) == 0:
            # Fallback to top somatic variants if no driver genes
            drivers_reported = somatic.head(3)
            
        if len(drivers_reported) == 0:
            report_md.append("\n*No somatic variants available for report details.*")
        else:
            for idx, row in drivers_reported.iterrows():
                report_md.append(f"\n### {row['Hugo_Symbol']} (Chr{row['Chromosome']}:{int(row['Start_Position'])})")
                
                # Add alert color based on relevance
                if row['llm_relevance'] == "pathogenic":
                    alert_type = "> [!IMPORTANT]"
                elif row['llm_classification'] == "somatic_high_conf":
                    alert_type = "> [!NOTE]"
                else:
                    alert_type = "> [!WARNING]"
                    
                report_md.append(f"{alert_type}")
                report_md.append(f"> **Classification**: {row['llm_classification'].replace('_', ' ').title()} (Confidence: {row['llm_confidence']:.2f})")
                report_md.append(f"> **Relevance**: {row['llm_relevance'].replace('_', ' ').upper()}")
                report_md.append(f"> **VAF**: {row['VAF']*100:.3f}% (Tumor depth: {int(row['t_depth'])}x)")
                report_md.append(f"> **Explanation**: {row['llm_explanation']}")
                report_md.append(f"\n**Bioinformatic Reasoning**:\n{row['llm_reasoning']}\n")
                
        # 5. Pipeline Parameters
        report_md.append(f"\n## 5. Pipeline Parameters and Architecture")
        report_md.append(
            f"This analysis was performed using a two-stage filter. First-stage statistical filters "
            f"enforced a minimum coverage depth of **{config.MIN_DEPTH}x**, minimum alt counts of **{config.MIN_ALT_COUNT}**, "
            f"and rejected variants with forward/reverse strand bias exceeding **{config.MAX_STRAND_BIAS}**. "
            f"Matched-normal control filtering rejected alleles with normal VAF > **{config.MAX_NORMAL_VAF*100:.1f}%**.\n\n"
            f"Second-stage classification utilized the LLM model `{config.LLM_MODEL}` on an OpenAI-compatible interface, "
            f"configured with temperature={config.LLM_TEMPERATURE}."
        )
        
        out_path = os.path.join(self.out_dir, "clinical_interpretation_report.md")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n".join(report_md))
            
        logger.info(f"Saved Markdown Clinical Report: {out_path}")
        return out_path
