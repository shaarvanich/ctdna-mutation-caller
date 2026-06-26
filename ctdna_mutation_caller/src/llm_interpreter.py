"""
llm_interpreter.py
------------------
Uses an OpenAI-compatible API to perform expert bioinformatics classification
of ctDNA mutations. Implements local caching, multi-threading, and a mock fallback.
"""

import os
import json
import hashlib
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

from openai import OpenAI
import config

logger = logging.getLogger(__name__)


class LLMInterpreter:
    """Interprets candidate mutations using LLM reasoning and clinical rules."""

    def __init__(self, 
                 api_key: str = None, 
                 base_url: str = config.LLM_BASE_URL, 
                 model: str = config.LLM_MODEL, 
                 temperature: float = config.LLM_TEMPERATURE, 
                 max_tokens: int = config.LLM_MAX_TOKENS,
                 cache_enabled: bool = config.LLM_CACHE_ENABLED,
                 max_workers: int = config.LLM_MAX_WORKERS):
        
        self.api_key = api_key if api_key else os.environ.get(config.LLM_API_KEY_ENV)
        self.base_url = base_url
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.cache_enabled = cache_enabled
        self.max_workers = max_workers
        
        self.mock_mode = False
        if not self.api_key:
            logger.warning(
                f"No API key found in environment variable '{config.LLM_API_KEY_ENV}'. "
                "Switching to local rule-based Mock LLM Mode. "
                "Set the key to run real LLM queries."
            )
            self.mock_mode = True
            self.client = None
        else:
            try:
                self.client = OpenAI(api_key=self.api_key, base_url=self.base_url)
                logger.info(f"Initialized OpenAI-compatible client connecting to {self.base_url} (Model: {self.model})")
            except Exception as e:
                logger.error(f"Failed to initialize LLM client: {e}. Falling back to Mock mode.")
                self.mock_mode = True
                self.client = None

    def _get_cache_path(self, variant_features: dict) -> str:
        """Generates a unique cache file path using MD5 hash of variant features."""
        # Use key genomic attributes for hashing to ensure stable cache hits
        key_features = {
            "gene": variant_features.get("Hugo_Symbol"),
            "chrom": variant_features.get("Chromosome"),
            "pos": variant_features.get("Start_Position"),
            "ref": variant_features.get("Reference_Allele"),
            "alt": variant_features.get("Tumor_Seq_Allele2"),
            "t_depth": variant_features.get("t_depth"),
            "t_alt": variant_features.get("t_alt_count"),
            "n_depth": variant_features.get("n_depth"),
            "n_alt": variant_features.get("n_alt_count"),
            "sb": variant_features.get("strand_bias")
        }
        feature_string = json.dumps(key_features, sort_keys=True)
        feature_hash = hashlib.md5(feature_string.encode('utf-8')).hexdigest()
        return os.path.join(config.LLM_CACHE_DIR, f"{feature_hash}.json")

    def _generate_mock_response(self, features: dict) -> dict:
        """
        Generates a high-fidelity mock response based on biological and statistical rules.
        Used as a fallback when no API key is provided.
        """
        gene = str(features.get("Hugo_Symbol", "")).upper()
        vaf = features.get("VAF", 0.0)
        depth = features.get("t_depth", 0)
        sb = features.get("strand_bias", 0.0)
        normal_vaf = features.get("normal_VAF", 0.0)
        cosmic_hit = features.get("cosmic_hit", 0)
        impact = features.get("impact_enc", 0)
        dbsnp_novel = features.get("dbsnp_novel", 1)
        
        is_driver = gene in config.KNOWN_ONCOGENES or gene in config.KNOWN_TUMOR_SUPPRESSORS or cosmic_hit == 1
        is_artifact_gene = gene in config.ARTIFACT_PRONE_GENES
        
        reasoning = []
        
        # 1. Evaluate somatic vs germline status
        if normal_vaf > 0.30:
            classification = "germline"
            confidence = 0.95
            relevance = "benign" if dbsnp_novel == 0 else "uncertain_significance"
            reasoning.append(f"Variant VAF in matched-normal is {normal_vaf*100:.1f}%, indicating a germline allele.")
        elif normal_vaf > 0.01:
            classification = "germline"
            confidence = 0.70
            relevance = "uncertain_significance"
            reasoning.append(f"Variant detected in normal sample (VAF={normal_vaf*100:.1f}%), likely germline leakage.")
        
        # 2. Evaluate sequencing artifact
        elif sb > 0.85:
            classification = "sequencing_artifact"
            confidence = 0.88
            relevance = "benign"
            reasoning.append(f"Variant shows extreme forward-reverse strand bias ({sb:.2f}). This is characteristic of sequencing chemistry artifacts.")
        elif depth < 30 and vaf < 0.01:
            classification = "sequencing_artifact"
            confidence = 0.75
            relevance = "benign"
            reasoning.append(f"Variant coverage is low ({depth}x) and variant allele frequency is ultra-low ({vaf*100:.2f}%), likely PCR/sequencing background noise.")
            
        # 3. Classify somatic variants
        else:
            if is_driver:
                classification = "somatic_high_conf"
                confidence = 0.92 if vaf >= 0.01 else 0.82
                relevance = "pathogenic"
                reasoning.append(f"Variant occurs in key cancer driver gene {gene} and lacks matched normal allele or strand bias.")
                reasoning.append(f"High probability of being an active driver somatic mutation (VAF={vaf*100:.2f}%).")
            elif is_artifact_gene:
                classification = "somatic_low_conf"
                confidence = 0.65
                relevance = "uncertain_significance"
                reasoning.append(f"Variant resides in {gene}, a large passenger gene known for high background mutation rates. Although quality metrics look clean, it could be a passenger somatic event or localized artifact.")
            else:
                classification = "somatic_high_conf" if vaf >= 0.02 else "somatic_low_conf"
                confidence = 0.80 if vaf >= 0.02 else 0.70
                relevance = "uncertain_significance" if impact < 2 else "pathogenic"
                reasoning.append(f"Clean somatic signal detected in {gene} (VAF={vaf*100:.2f}%) with no normal sample leakage.")
                if impact >= 2:
                    reasoning.append("VEP predicts high/moderate functional impact, suggesting pathogenic biological significance.")
                else:
                    reasoning.append("VEP predicts low or silent impact, suggesting a passenger mutation.")
                    
        return {
            "classification": classification,
            "confidence_score": confidence,
            "clinical_relevance": relevance,
            "reasoning_steps": " ".join(reasoning),
            "explanation": f"Mock classified variant as {classification} in {gene} based on coverage ({depth}x), VAF ({vaf*100:.2f}%), and driver status.",
            "recommended_action": "report" if "somatic" in classification else "filter"
        }

    def _query_llm(self, features: dict) -> dict:
        """Sends a single mutation details to the LLM and retrieves classification."""
        if self.mock_mode:
            # Simulate network latency for mock calls
            time.sleep(0.01)
            return self._generate_mock_response(features)
            
        # Compile prompts
        system_prompt = (
            "You are a clinical cancer genomics expert and bioinformatician. Your task is to evaluate a candidate "
            "circulating tumor DNA (ctDNA) somatic mutation call. You will review statistical sequencing metrics "
            "and biological markers, and output a structured JSON classification.\n\n"
            "Classes:\n"
            "- 'somatic_high_conf': Clear somatic signal with solid coverage, balanced strand reads, and absence in normal control.\n"
            "- 'somatic_low_conf': Likely real somatic variant, but in a passenger/artifact gene, or has borderline quality counts.\n"
            "- 'germline': Variant represents inherited germline allele (usually normal VAF ~50% or ~100%).\n"
            "- 'sequencing_artifact': Sequencing error or PCR noise. Indicated by extreme strand bias, poor coverage, or presence in panel of normals.\n\n"
            "Return EXACTLY a JSON object with these keys (no extra text outside the JSON):\n"
            "{\n"
            "  \"classification\": \"somatic_high_conf\" | \"somatic_low_conf\" | \"germline\" | \"sequencing_artifact\",\n"
            "  \"confidence_score\": float (0.0 to 1.0),\n"
            "  \"clinical_relevance\": \"pathogenic\" | \"benign\" | \"uncertain_significance\",\n"
            "  \"reasoning_steps\": \"Detailed explanation of your biological and statistical reasoning\",\n"
            "  \"explanation\": \"A 1-2 sentence summary for clinicians\",\n"
            "  \"recommended_action\": \"report\" | \"follow_up\" | \"filter\"\n"
            "}"
        )
        
        user_prompt = (
            f"Please evaluate this variant:\n"
            f"- Gene Symbol: {features.get('Hugo_Symbol')}\n"
            f"- Coordinates: Chromosome {features.get('Chromosome')}, Position {features.get('Start_Position')}\n"
            f"- Alleles: Ref={features.get('Reference_Allele')}, Alt={features.get('Tumor_Seq_Allele2')}\n"
            f"- Tumor Coverage: Depth={features.get('t_depth')}x, Alt Reads={features.get('t_alt_count')}, VAF={features.get('VAF')*100:.3f}%\n"
            f"- Matched Normal Coverage: Depth={features.get('n_depth')}x, Alt Reads={features.get('n_alt_count')}, VAF={features.get('normal_VAF')*100:.3f}%\n"
            f"- Strand Bias Index (0=balanced, 1=skewed): {features.get('strand_bias'):.2f}\n"
            f"- VEP Impact Score: {features.get('impact_enc')} / 3 (Missense={features.get('is_missense')}, Nonsense={features.get('is_nonsense')}, Silent={features.get('is_silent')})\n"
            f"- Functional Prediction: SIFT={features.get('sift_enc')} (1=del, 0=tol), PolyPhen={features.get('polyphen_enc')} (2=damaging, 0=benign)\n"
            f"- dbSNP Database Status: {'Novel (Not in germline databases)' if features.get('dbsnp_novel') == 1 else 'Known Common Polymorphism'}\n"
            f"- Cancer driver databases: {'COSMIC Cancer Gene Census member' if features.get('cosmic_hit') == 1 else 'No driver record'}"
        )
        
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                response_format={"type": "json_object"},
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                temperature=self.temperature,
                max_tokens=self.max_tokens,
                timeout=15
            )
            
            content = response.choices[0].message.content
            
            # Clean possible markdown wrap-around in case model ignored JSON instruction
            if content.startswith("```"):
                content = content.replace("```json", "").replace("```", "").strip()
                
            return json.loads(content)
            
        except Exception as e:
            logger.error(f"Error querying LLM API for {features.get('Hugo_Symbol')} chr{features.get('Chromosome')}:{features.get('Start_Position')}: {e}")
            # Return fallback mock response if API fails mid-run
            return self._generate_mock_response(features)

    def process_variant(self, variant_row: dict) -> dict:
        """Processes a single variant, checking the local cache first before calling the LLM API."""
        cache_path = self._get_cache_path(variant_row)
        
        # Load from cache if enabled and file exists
        if self.cache_enabled and os.path.exists(cache_path):
            try:
                with open(cache_path, "r") as f:
                    cached_res = json.load(f)
                logger.debug(f"Cache HIT for {variant_row.get('Hugo_Symbol')} {variant_row.get('Chromosome')}:{variant_row.get('Start_Position')}")
                return cached_res
            except Exception as e:
                logger.warning(f"Failed to read cache file {cache_path}: {e}")
                
        # Cache miss, call query
        logger.debug(f"Cache MISS for {variant_row.get('Hugo_Symbol')} {variant_row.get('Chromosome')}:{variant_row.get('Start_Position')}. Querying LLM...")
        result = self._query_llm(variant_row)
        
        # Write to cache
        if self.cache_enabled and result:
            try:
                with open(cache_path, "w") as f:
                    json.dump(result, f, indent=2)
            except Exception as e:
                logger.warning(f"Failed to write cache file {cache_path}: {e}")
                
        return result

    def interpret_candidates(self, df: pd.DataFrame, limit: int = None) -> pd.DataFrame:
        """
        Runs LLM interpretation on all variants where 'is_candidate' is True.
        Utilizes ThreadPoolExecutor for fast concurrent API queries.
        
        Args:
            df (pd.DataFrame): Dataframe containing variants.
            limit (int, optional): Max number of variants to interpret (for testing/cost savings).
            
        Returns:
            pd.DataFrame: Dataframe with LLM classification columns appended.
        """
        # Filter candidate variants
        candidates = df[df["is_candidate"] == True].copy()
        
        if len(candidates) == 0:
            logger.warning("No candidate variants found for LLM interpretation. Skip.")
            # Create empty columns
            df["llm_classification"] = "filtered"
            df["llm_confidence"] = 0.0
            df["llm_relevance"] = "benign"
            df["llm_reasoning"] = "Filtered out by statistical rules."
            df["llm_explanation"] = "Filtered out by statistical rules."
            df["llm_action"] = "filter"
            return df
            
        if limit:
            logger.info(f"Limiting LLM interpretation to first {limit} candidate variants (out of {len(candidates)}).")
            candidates = candidates.head(limit)
            
        logger.info(f"Running LLM interpretation on {len(candidates)} candidates using {self.max_workers} threads...")
        
        # List to collect records
        variant_records = candidates.to_dict(orient="records")
        results_map = {}
        
        start_time = time.time()
        
        # Threaded execution
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_idx = {
                executor.submit(self.process_variant, record): idx 
                for idx, record in enumerate(variant_records)
            }
            
            # Retrieve completed futures
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                record = variant_records[idx]
                key = (
                    record.get("Chromosome", "chrUn"),
                    record.get("Start_Position", 0),
                    record.get("Reference_Allele", "N"),
                    record.get("Tumor_Seq_Allele2", "N")
                )
                try:
                    llm_res = future.result()
                    results_map[key] = llm_res
                except Exception as exc:
                    logger.error(f"Variant index {idx} generated an exception during LLM execution: {exc}")
                    results_map[key] = self._generate_mock_response(record)
                    
        elapsed_time = time.time() - start_time
        logger.info(f"LLM Interpretation finished in {elapsed_time:.2f} seconds.")
        
        # Integrate LLM results back into original DataFrame
        llm_classifications = []
        llm_confidences = []
        llm_relevances = []
        llm_reasonings = []
        llm_explanations = []
        llm_actions = []
        
        for idx, row in df.iterrows():
            key = (
                row.get("Chromosome", "chrUn"),
                row.get("Start_Position", 0),
                row.get("Reference_Allele", "N"),
                row.get("Tumor_Seq_Allele2", "N")
            )
            
            # If variant was selected and processed by LLM
            if key in results_map:
                res = results_map[key]
                llm_classifications.append(res.get("classification", "unknown"))
                llm_confidences.append(float(res.get("confidence_score", 0.0)))
                llm_relevances.append(res.get("clinical_relevance", "uncertain_significance"))
                llm_reasonings.append(res.get("reasoning_steps", ""))
                llm_explanations.append(res.get("explanation", ""))
                llm_actions.append(res.get("recommended_action", "filter"))
            else:
                # If variant was filtered out by rule-based filters
                llm_classifications.append("filtered")
                llm_confidences.append(0.0)
                llm_relevances.append("benign")
                llm_reasonings.append("Filtered by statistical heuristics.")
                llm_explanations.append("Filtered by statistical heuristics.")
                llm_actions.append("filter")
                
        df["llm_classification"] = llm_classifications
        df["llm_confidence"] = llm_confidences
        df["llm_relevance"] = llm_relevances
        df["llm_reasoning"] = llm_reasonings
        df["llm_explanation"] = llm_explanations
        df["llm_action"] = llm_actions
        
        # Summary log
        n_high = sum(df["llm_classification"] == "somatic_high_conf")
        n_low = sum(df["llm_classification"] == "somatic_low_conf")
        n_germline = sum(df["llm_classification"] == "germline")
        n_artifact = sum(df["llm_classification"] == "sequencing_artifact")
        
        logger.info(f"LLM Classification Summary:")
        logger.info(f"  Somatic High Confidence : {n_high}")
        logger.info(f"  Somatic Low Confidence  : {n_low}")
        logger.info(f"  Germline Leakage        : {n_germline}")
        logger.info(f"  Sequencing Artifact     : {n_artifact}")
        
        return df
