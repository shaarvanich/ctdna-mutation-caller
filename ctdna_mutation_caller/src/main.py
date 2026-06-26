"""
main.py
-------
Main entry point and orchestrator for the LLM-Assisted ctDNA Somatic Mutation Caller.
Integrates loaders, preprocessors, filters, LLM interpreters, annotators,
visualizers, and report builders into a single CLI tool.
"""

import os
import sys
import argparse
import logging
from datetime import datetime
import pandas as pd

# Adjust path to import config and local modules if run from root or src
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import config
from data_loader import DataLoader
from preprocessing import Preprocessor
from variant_filter import VariantFilter
from vaf_analysis import VAFAnalyzer
from annotation import VariantAnnotator
from llm_interpreter import LLMInterpreter
from visualization import PipelineVisualizer
from report_generator import ReportGenerator

# Setup logging
log_file = os.path.join(config.OUTPUTS_DIR, "pipeline.log")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file, encoding="utf-8")
    ]
)
logger = logging.getLogger("main_pipeline")


def run_pipeline(args):
    logger.info("=" * 70)
    logger.info("  LLM-ASSISTED ctDNA SOMATIC MUTATION CALLER PIPELINE")
    logger.info("=" * 70)
    logger.info(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    loader = DataLoader()
    preprocessor = Preprocessor()
    analyzer = VAFAnalyzer()
    
    input_file = args.input
    is_simulated = False
    
    # ── STEP 1: Data Acquisition / Loader ─────────────────────────────────────
    # If mode is explicitly simulate
    if args.mode == "simulate":
        logger.info(f"Mode: Simulate. Generating {args.simulate_size} synthetic variants...")
        sim_prefix = os.path.join(config.DATA_DIR, "simulated_ctdna")
        maf_path, vcf_path = loader.generate_simulated_dataset(sim_prefix, size=args.simulate_size)
        input_file = maf_path
        is_simulated = True
        
    elif args.mode == "download":
        logger.info("Mode: Download. Fetching TCGA-LUAD dataset from GDC API...")
        try:
            downloaded_maf = loader.download_tcga_luad_maf()
            logger.info(f"Successfully downloaded GDC MAF to: {downloaded_maf}")
            logger.info("You can now run the full pipeline using this file as --input.")
            return
        except Exception as e:
            logger.error(f"Download failed: {e}")
            sys.exit(1)
            
    # Default full pipeline or other steps
    else:
        # If no input is specified, search for files
        if not input_file:
            logger.info("No input file specified. Searching for data files...")
            
            # Check 1: Provided features.csv in extracted directory or data directory
            extracted_features = os.path.join(config.BASE_DIR, "extracted_files", "features.csv")
            data_features = os.path.join(config.DATA_DIR, "features.csv")
            local_maf = os.path.join(config.DATA_DIR, "TCGA-LUAD.maf")
            
            if os.path.exists(extracted_features):
                logger.info(f"Found features.csv in extracted_files. Copying to data directory...")
                shutil_copy = True
                try:
                    import shutil
                    shutil.copy2(extracted_features, data_features)
                    input_file = data_features
                except Exception as ex:
                    logger.warning(f"Could not copy file: {ex}. Using extracted directory path directly.")
                    input_file = extracted_features
            elif os.path.exists(data_features):
                logger.info(f"Found existing features.csv at: {data_features}")
                input_file = data_features
            elif os.path.exists(local_maf):
                logger.info(f"Found local TCGA-LUAD.maf at: {local_maf}")
                input_file = local_maf
                
            # If no local files found, try GDC download
            else:
                logger.info("No local datasets found. Attempting to download TCGA-LUAD somatic MAF via GDC API...")
                try:
                    input_file = loader.download_tcga_luad_maf()
                except Exception as e:
                    logger.warning("GDC API download failed or was interrupted. Fallback to generating simulated dataset...")
                    sim_prefix = os.path.join(config.DATA_DIR, "simulated_ctdna")
                    maf_path, _ = loader.generate_simulated_dataset(sim_prefix, size=args.simulate_size)
                    input_file = maf_path
                    is_simulated = True

    # Load data based on extension
    logger.info(f"Loading input file: {input_file}")
    if input_file.endswith(".maf") or ".maf" in input_file:
        raw_df = loader.parse_maf(input_file)
    elif input_file.endswith(".vcf") or ".vcf" in input_file:
        raw_df = loader.parse_vcf(input_file)
    elif input_file.endswith(".csv"):
        raw_df = loader.parse_csv(input_file)
    else:
        logger.error("Unsupported file format. Must be MAF, VCF, or CSV.")
        sys.exit(1)

    if len(raw_df) == 0:
        logger.error("Input file contains 0 records. Exiting.")
        sys.exit(1)

    # ── STEP 2 & 3: Preprocessing & Biological Annotation ─────────────────────
    # Preprocess features
    preprocessed_df = preprocessor.preprocess(raw_df)
    
    # Annotate driver role metadata
    annotator = VariantAnnotator(cosmic_path=args.cosmic)
    annotated_df = annotator.annotate(preprocessed_df)

    # ── STEP 4: Rule-based Filtering ──────────────────────────────────────────
    # Apply hard filters & Driver Rescue Heuristic
    v_filter = VariantFilter(
        min_depth=args.min_depth,
        min_alt_count=args.min_alt_count,
        min_vaf=args.min_vaf,
        max_strand_bias=args.max_strand_bias,
        max_normal_vaf=args.max_normal_vaf
    )
    filtered_df = v_filter.filter_variants(annotated_df)

    # Exit early if only running filtering
    if args.mode == "filter":
        out_csv = os.path.join(config.OUTPUTS_DIR, "filtered_variants.csv")
        filtered_df.to_csv(out_csv, index=False)
        logger.info(f"Mode filter finished. Results written to: {out_csv}")
        return

    # ── STEP 5: LLM Interpretation ────────────────────────────────────────────
    # Classify candidate mutations using LLM (OpenAI/Local/Mock)
    if args.run_llm:
        interpreter = LLMInterpreter(
            api_key=args.api_key,
            model=args.model,
            temperature=args.temperature,
            cache_enabled=not args.no_cache,
            max_workers=args.llm_workers
        )
        # Process candidates (allow limit for testing/cost control)
        final_df = interpreter.interpret_candidates(filtered_df, limit=args.llm_limit)
    else:
        logger.info("LLM interpretation disabled via --no-llm. Tagging candidates as somatic_unverified.")
        final_df = filtered_df.copy()
        final_df["llm_classification"] = final_df["is_candidate"].apply(
            lambda c: "somatic_unverified" if c else "filtered"
        )
        final_df["llm_confidence"] = final_df["is_candidate"].apply(lambda c: 0.5 if c else 0.0)
        final_df["llm_relevance"] = "uncertain_significance"
        final_df["llm_reasoning"] = "LLM interpretation disabled."
        final_df["llm_explanation"] = "LLM interpretation disabled."
        final_df["llm_action"] = final_df["is_candidate"].apply(lambda c: "report" if c else "filter")

    # ── STEP 6: VAF Tier Analysis ──────────────────────────────────────────────
    # Stratify calls by VAF tiers (<1%, 1-5%, >5%) and compute detection stats
    tier_metrics = analyzer.compute_tier_metrics(final_df)

    # ── STEP 7: Visualization ─────────────────────────────────────────────────
    visualizer = PipelineVisualizer()
    visualizer.generate_all_plots(final_df, tier_metrics)

    # ── STEP 8: Report Generation ─────────────────────────────────────────────
    reporter = ReportGenerator()
    csv_path = reporter.export_csv_report(final_df)
    json_path = reporter.export_json_summary(final_df, tier_metrics, input_file)
    md_path = reporter.export_clinical_report(final_df, input_file)

    logger.info("=" * 70)
    logger.info("  PIPELINE EXECUTION COMPLETE")
    logger.info("=" * 70)
    logger.info(f"CSV Mutation Report  : [somatic_mutations_called.csv](file:///{csv_path.replace(os.sep, '/')})")
    logger.info(f"JSON Statistics      : [pipeline_run_summary.json](file:///{json_path.replace(os.sep, '/')})")
    logger.info(f"Clinical Markdown    : [clinical_interpretation_report.md](file:///{md_path.replace(os.sep, '/')})")
    logger.info(f"Visualizations plots : [plots/](file:///{config.PLOTS_DIR.replace(os.sep, '/')})")
    logger.info("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="LLM-Assisted ctDNA Somatic Mutation Caller for Low Variant Allele Frequency Detection"
    )
    
    # Execution modes
    parser.add_argument(
        "--mode", 
        choices=["full", "download", "simulate", "filter"], 
        default="full",
        help="Pipeline execution mode. 'full' runs all steps; 'download' fetches MAF from GDC; "
             "'simulate' generates synthetic data; 'filter' runs statistical rules only."
    )
    
    # Input/Output paths
    parser.add_argument("--input", help="Path to input mutation dataset (.maf, .vcf, or .csv)")
    parser.add_argument("--cosmic", help="Path to COSMIC Cancer Gene Census CSV")
    
    # Simulator settings
    parser.add_argument("--simulate-size", type=int, default=200, help="Number of variants to generate in simulation mode")
    
    # Filter thresholds override
    parser.add_argument("--min-depth", type=int, default=config.MIN_DEPTH, help="Minimum total coverage depth")
    parser.add_argument("--min-alt-count", type=int, default=config.MIN_ALT_COUNT, help="Minimum alt supporting reads")
    parser.add_argument("--min-vaf", type=float, default=config.MIN_VAF, help="Minimum Variant Allele Frequency")
    parser.add_argument("--max-strand-bias", type=float, default=config.MAX_STRAND_BIAS, help="Maximum strand bias (0-1)")
    parser.add_argument("--max-normal-vaf", type=float, default=config.MAX_NORMAL_VAF, help="Maximum normal control VAF")
    
    # LLM Settings
    parser.add_argument("--no-llm", dest="run_llm", action="store_false", help="Disable LLM interpretation step")
    parser.add_argument("--api-key", help="OpenAI-compatible API key")
    parser.add_argument("--model", default=config.LLM_MODEL, help="Model name (e.g. gpt-4o-mini)")
    parser.add_argument("--temperature", type=float, default=config.LLM_TEMPERATURE, help="LLM temperature")
    parser.add_argument("--llm-limit", type=int, help="Max number of candidate mutations to query (saves API cost)")
    parser.add_argument("--llm-workers", type=int, default=config.LLM_MAX_WORKERS, help="Number of concurrent LLM threads")
    parser.add_argument("--no-cache", action="store_true", help="Disable local LLM response caching")
    
    parser.set_defaults(run_llm=True)
    args = parser.parse_args()
    
    # Check if directories exist
    os.makedirs(config.DATA_DIR, exist_ok=True)
    os.makedirs(config.DOWNLOADS_DIR, exist_ok=True)
    os.makedirs(config.OUTPUTS_DIR, exist_ok=True)
    
    run_pipeline(args)


if __name__ == "__main__":
    main()
