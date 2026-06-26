"""
data_loader.py
--------------
Loads public cancer mutation datasets (MAF, VCF, CSV), queries the GDC API
to download real TCGA-LUAD somatic MAFs, and contains a simulation engine
for generating low-VAF somatic mutations, germline variants, and artifacts.
"""

import os
import gzip
import shutil
import logging
import random
import pandas as pd
import numpy as np
import requests

from config import DATA_DIR, DOWNLOADS_DIR, GDC_API_BASE, TCGA_PROJECT

logger = logging.getLogger(__name__)


class DataLoader:
    """Handles data downloading, parsing, and simulation for ctDNA mutation analysis."""

    def __init__(self):
        pass

    # ── 1. GDC API Downloader ──────────────────────────────────────────────────
    def download_tcga_luad_maf(self) -> str:
        """
        Queries the GDC API to search for the public open-access Masked Somatic Mutation MAF
        file for TCGA-LUAD, downloads it, and extracts it to the downloads directory.
        
        Returns:
            str: Path to the extracted MAF file.
        """
        logger.info(f"Querying GDC API for {TCGA_PROJECT} Masked Somatic Mutation MAF...")
        
        # Search filter for GDC files
        filters = {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": [TCGA_PROJECT]}},
                {"op": "in", "content": {"field": "files.data_type", "value": ["Masked Somatic Mutation"]}},
                {"op": "in", "content": {"field": "files.data_format", "value": ["MAF"]}}
            ]
        }
        
        params = {
            "filters": pd.json_normalize(filters).to_json(orient='records')[0] if hasattr(pd, "json_normalize") else str(filters).replace("'", '"'),
            "format": "JSON",
            "size": "5"
        }
        
        # In python, we can query GDC files endpoint
        files_url = f"{GDC_API_BASE}/files"
        try:
            # Let's construct a direct filter payload
            import json
            response = requests.get(files_url, params={"filters": json.dumps(filters), "format": "JSON", "size": "5"}, timeout=30)
            response.raise_for_status()
            res_data = response.json()
            
            hits = res_data.get("data", {}).get("hits", [])
            if not hits:
                raise ValueError(f"No MAF files found for project {TCGA_PROJECT} on GDC.")
            
            # Select the first file (usually the latest release)
            file_info = hits[0]
            file_id = file_info["id"]
            file_name = file_info["file_name"]
            file_size = file_info["file_size"]
            
            logger.info(f"Found GDC file: {file_name} (ID: {file_id}, Size: {file_size / (1024*1024):.2f} MB)")
            
            # Paths
            download_path = os.path.join(DOWNLOADS_DIR, file_name)
            extracted_path = os.path.join(DOWNLOADS_DIR, file_name.replace(".gz", ""))
            
            # If already extracted, return it
            if os.path.exists(extracted_path):
                logger.info(f"MAF file already exists and is decompressed: {extracted_path}")
                return extracted_path
                
            # If gzip exists but not decompressed, decompress it
            if os.path.exists(download_path):
                logger.info("Compressed MAF file exists. Decompressing...")
                self._decompress_gz(download_path, extracted_path)
                return extracted_path
            
            # Download file
            data_url = f"{GDC_API_BASE}/data/{file_id}"
            logger.info(f"Downloading MAF from GDC: {data_url}")
            
            with requests.get(data_url, stream=True, timeout=120) as stream:
                stream.raise_for_status()
                with open(download_path, "wb") as f_out:
                    shutil.copyfileobj(stream.raw, f_out)
            
            logger.info(f"Successfully downloaded to {download_path}")
            logger.info("Decompressing MAF...")
            self._decompress_gz(download_path, extracted_path)
            return extracted_path
            
        except Exception as e:
            logger.error(f"Failed to download real TCGA-LUAD dataset: {e}")
            logger.info("Falling back to local data checking or simulation...")
            raise e

    def _decompress_gz(self, source: str, target: str):
        """Decompresses a .gz file."""
        with gzip.open(source, 'rb') as f_in:
            with open(target, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        logger.info(f"Extracted to {target}")

    # ── 2. MAF Parser ──────────────────────────────────────────────────────────
    def parse_maf(self, path: str) -> pd.DataFrame:
        """
        Parses a Mutation Annotation Format (MAF) file.
        Skips comments starting with '#' and standardizes columns.
        """
        logger.info(f"Parsing MAF file: {path}")
        if not os.path.exists(path):
            raise FileNotFoundError(f"MAF file not found: {path}")
            
        # Read the file skipping comment lines
        df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
        logger.info(f"Loaded {len(df):,} rows from MAF.")
        
        # Normalize and map columns to common naming
        column_mapping = {
            "Hugo_Symbol": "Hugo_Symbol",
            "Chromosome": "Chromosome",
            "Start_Position": "Start_Position",
            "Reference_Allele": "Reference_Allele",
            "Tumor_Seq_Allele2": "Tumor_Seq_Allele2",
            "Variant_Classification": "Variant_Classification",
            "Variant_Type": "Variant_Type",
            "Tumor_Sample_Barcode": "Tumor_Sample_Barcode",
            "FILTER": "FILTER_raw",
            "dbSNP_RS": "dbSNP_RS",
            "IMPACT": "IMPACT",
            "PolyPhen": "PolyPhen",
            "SIFT": "SIFT",
            # Read counts
            "t_depth": "t_depth",
            "t_alt_count": "t_alt_count",
            "t_ref_count": "t_ref_count",
            "n_depth": "n_depth",
            "n_alt_count": "n_alt_count",
            "n_ref_count": "n_ref_count",
            # Forward/reverse counts if available for strand bias
            "t_alt_count_forward": "t_alt_count_forward",
            "t_alt_count_reverse": "t_alt_count_reverse",
            "t_ref_count_forward": "t_ref_count_forward",
            "t_ref_count_reverse": "t_ref_count_reverse"
        }
        
        # Find matching columns
        rename_dict = {}
        for raw_col, std_col in column_mapping.items():
            if raw_col in df.columns:
                rename_dict[raw_col] = std_col
                
        df = df.rename(columns=rename_dict)
        
        # Ensure minimum key columns exist, fill if missing
        required_cols = ["Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", 
                         "Tumor_Seq_Allele2", "t_depth", "t_alt_count"]
        for col in required_cols:
            if col not in df.columns:
                if col == "t_depth" and "t_ref_count" in df.columns and "t_alt_count" in df.columns:
                    df["t_depth"] = df["t_ref_count"] + df["t_alt_count"]
                else:
                    df[col] = np.nan
                    
        return df

    # ── 3. Pure-Python VCF Parser ──────────────────────────────────────────────
    def parse_vcf(self, path: str) -> pd.DataFrame:
        """
        Parses a Variant Call Format (VCF) file natively in Python.
        Handles standard fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
        and parses FORMAT/sample columns (e.g., AD, DP, AF, SB).
        """
        logger.info(f"Parsing VCF file: {path}")
        if not os.path.exists(path):
            raise FileNotFoundError(f"VCF file not found: {path}")
            
        variants = []
        samples = []
        
        with open(path, "r") as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    # Header line
                    headers = line.strip().split("\t")
                    if len(headers) > 9:
                        samples = headers[9:]
                    continue
                
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                    
                chrom = parts[0]
                pos = int(parts[1])
                db_id = parts[2]
                ref = parts[3]
                alt = parts[4]
                qual = parts[5]
                filt = parts[6]
                info_raw = parts[7]
                
                # Parse INFO column
                info_dict = {}
                for item in info_raw.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info_dict[k] = v
                    else:
                        info_dict[item] = True
                
                # Default values
                t_depth, t_alt, t_ref = np.nan, np.nan, np.nan
                n_depth, n_alt, n_ref = 0, 0, 0
                t_fwd_alt, t_rev_alt = np.nan, np.nan
                
                # Parse sample format columns if sample data is present
                if len(parts) > 8 and samples:
                    format_keys = parts[8].split(":")
                    
                    # We assume the first sample is Tumor and the second is Normal, 
                    # or we try to identify them based on header/INFO.
                    # Commonly: Tumor is first or named TUMOR.
                    tumor_idx = 0
                    normal_idx = 1 if len(samples) > 1 else None
                    
                    # Read tumor sample
                    tumor_vals = parts[9].split(":")
                    tumor_data = dict(zip(format_keys, tumor_vals))
                    
                    # Extract Depth (DP)
                    if "DP" in tumor_data:
                        t_depth = float(tumor_data["DP"])
                    elif "DP" in info_dict:
                        t_depth = float(info_dict["DP"])
                        
                    # Extract Allele Depth (AD) containing [Ref, Alt] counts
                    if "AD" in tumor_data:
                        ad_parts = tumor_data["AD"].split(",")
                        if len(ad_parts) >= 2:
                            t_ref = float(ad_parts[0])
                            t_alt = float(ad_parts[1])
                            if np.isnan(t_depth) or t_depth == 0:
                                t_depth = t_ref + t_alt
                    
                    # Extract VAF (AF)
                    t_vaf = np.nan
                    if "AF" in tumor_data:
                        t_vaf = float(tumor_data["AF"].split(",")[0])
                    elif t_depth > 0:
                        t_vaf = t_alt / t_depth
                        
                    # Extract strand counts if available in format (e.g. SB, DP4)
                    # SB is usually 4 integers: ref_fwd, ref_rev, alt_fwd, alt_rev
                    if "SB" in tumor_data:
                        sb_parts = tumor_data["SB"].split(",")
                        if len(sb_parts) == 4:
                            t_fwd_alt = float(sb_parts[2])
                            t_rev_alt = float(sb_parts[3])
                    elif "DP4" in info_dict:
                        dp4_parts = info_dict["DP4"].split(",")
                        if len(dp4_parts) == 4:
                            t_fwd_alt = float(dp4_parts[2])
                            t_rev_alt = float(dp4_parts[3])
                            
                    # Read normal sample if present
                    if normal_idx is not None and len(parts) > 10:
                        normal_vals = parts[10].split(":")
                        normal_data = dict(zip(format_keys, normal_vals))
                        if "DP" in normal_data:
                            n_depth = float(normal_data["DP"])
                        if "AD" in normal_data:
                            nad_parts = normal_data["AD"].split(",")
                            if len(nad_parts) >= 2:
                                n_ref = float(nad_parts[0])
                                n_alt = float(nad_parts[1])
                                if n_depth == 0:
                                    n_depth = n_ref + n_alt
                                    
                # Map INFO annotations like Gene Name if present (e.g., from VEP, SnpEff)
                gene = info_dict.get("GENE", info_dict.get("Gene", info_dict.get("HUGO", ".")))
                if gene == "." and "CSQ" in info_dict:
                    # CSQ is VEP format: Allele|Consequence|IMPACT|Symbol|Gene|...
                    csq_fields = info_dict["CSQ"].split(",")[0].split("|")
                    if len(csq_fields) > 3:
                        gene = csq_fields[3] # Hugo Symbol is usually 4th field
                        
                v_class = info_dict.get("VC", info_dict.get("Variant_Classification", "."))
                impact = info_dict.get("IMPACT", ".")
                
                variants.append({
                    "Hugo_Symbol": gene if gene != "." else f"GENE_{chrom}_{pos}",
                    "Chromosome": chrom,
                    "Start_Position": pos,
                    "Reference_Allele": ref,
                    "Tumor_Seq_Allele2": alt,
                    "t_depth": t_depth,
                    "t_alt_count": t_alt,
                    "t_ref_count": t_ref,
                    "n_depth": n_depth,
                    "n_alt_count": n_alt,
                    "n_ref_count": n_ref,
                    "t_alt_count_forward": t_fwd_alt,
                    "t_alt_count_reverse": t_rev_alt,
                    "FILTER_raw": filt,
                    "dbSNP_RS": db_id if db_id != "." else "novel",
                    "Variant_Classification": v_class,
                    "IMPACT": impact
                })
                
        df = pd.DataFrame(variants)
        logger.info(f"Loaded {len(df):,} variants from VCF.")
        return df

    # ── 4. CSV Loader ──────────────────────────────────────────────────────────
    def parse_csv(self, path: str) -> pd.DataFrame:
        """Loads a preprocessed feature CSV (e.g. features.csv) and ensures correct column types."""
        logger.info(f"Loading CSV file: {path}")
        if not os.path.exists(path):
            raise FileNotFoundError(f"CSV file not found: {path}")
        df = pd.read_csv(path)
        logger.info(f"Loaded {len(df):,} rows from CSV.")
        return df

    # ── 5. Simulation Engine ───────────────────────────────────────────────────
    def generate_simulated_dataset(self, out_prefix: str, size: int = 200) -> tuple:
        """
        Generates simulated high-fidelity somatic mutations, sequencing artifacts, 
        and germline variants, outputting both a MAF and a VCF file.
        
        Args:
            out_prefix: File path prefix for output (creates out_prefix.maf and out_prefix.vcf)
            size: Total number of variants to simulate
            
        Returns:
            tuple: (maf_path, vcf_path)
        """
        logger.info(f"Simulating {size} variants with ctDNA low-VAF characteristics...")
        
        driver_genes = ["EGFR", "TP53", "KRAS", "STK11", "KEAP1", "RB1", "SMAD4", "CDKN2A", "BRAF", "PIK3CA", "PTEN", "ALK"]
        large_passenger_genes = ["TTN", "MUC16", "RYR2", "OBSCN", "USH2A", "SYNE1", "NEB", "DST"]
        other_genes = ["APOB", "LRP2", "PCLO", "CSMD1", "ADGRV1", "ZFHX4", "FAT4", "DMD", "SPTA1", "XIRP2"]
        
        variant_classifications = ["Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Silent"]
        impacts = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
        
        variants = []
        
        for i in range(size):
            # Select variant category:
            # 35% somatic, 45% sequencing artifact/noise, 20% germline
            choice = random.random()
            
            chrom = f"chr{random.choice(list(range(1, 23)) + ['X'])}"
            pos = random.randint(100000, 200000000)
            ref_allele = random.choice(["A", "C", "G", "T"])
            alt_allele = random.choice([a for a in ["A", "C", "G", "T"] if a != ref_allele])
            
            # Matched sample barcode
            sample_barcode = f"SIM-PATIENT-{random.randint(101, 150)}"
            
            if choice < 0.35:
                # ── Category 1: Somatic ctDNA Variant (True Somatic) ──
                # Low VAF (0.1% to 5.0%), high depth (characteristic of ctDNA sequencing, e.g. 500x - 5000x)
                t_depth = random.randint(500, 3000)
                # Ensure VAF is low: 0.1% to 5.0%
                vaf = random.uniform(0.001, 0.05)
                t_alt = max(3, int(round(vaf * t_depth)))
                t_ref = t_depth - t_alt
                vaf = t_alt / t_depth # recalculate exact
                
                # Normal sample (matched control) has zero or near-zero VAF (somatic signal)
                n_depth = random.randint(100, 500)
                n_alt = 0 if random.random() < 0.95 else 1 # occasional sequencing noise in normal
                n_ref = n_depth - n_alt
                
                # Strand bias: balanced (true mutations appear on both strands)
                t_alt_fwd = int(round(t_alt * random.uniform(0.35, 0.65)))
                t_alt_rev = t_alt - t_alt_fwd
                t_ref_fwd = int(round(t_ref * 0.5))
                t_ref_rev = t_ref - t_ref_fwd
                
                # Biology
                gene = random.choice(driver_genes) if random.random() < 0.8 else random.choice(other_genes)
                v_class = random.choice(variant_classifications[:-1]) # rarely silent
                impact = "HIGH" if v_class in ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"] else "MODERATE"
                db_id = "novel" if random.random() < 0.95 else f"rs{random.randint(1000000, 99999999)}"
                polyphen = "probably_damaging" if impact == "HIGH" else ("possibly_damaging" if random.random() < 0.7 else "benign")
                sift = "deleterious" if impact in ["HIGH", "MODERATE"] else "tolerated"
                filt_raw = "PASS"
                label = 1
                
            elif choice < 0.80:
                # ── Category 2: Sequencing Artifact / Noise ──
                # Very low alt count, low depth, extreme strand bias
                t_depth = random.randint(50, 400) # lower coverage often yields artifacts
                vaf = random.uniform(0.005, 0.04)
                t_alt = random.choice([1, 2, 3, 4]) # very few reads
                t_ref = t_depth - t_alt
                vaf = t_alt / t_depth
                
                # Normal sample
                n_depth = random.randint(30, 200)
                n_alt = 0
                n_ref = n_depth
                
                # Strand bias: extremely skewed (e.g. all alt reads on one strand due to dye/sequencing chemistry error)
                t_alt_fwd = t_alt if random.random() < 0.85 else 0
                t_alt_rev = t_alt - t_alt_fwd
                t_ref_fwd = int(round(t_ref * 0.5))
                t_ref_rev = t_ref - t_ref_fwd
                
                # Biology
                # Often occurs in massive passenger genes or random other genes
                gene = random.choice(large_passenger_genes) if random.random() < 0.6 else random.choice(other_genes)
                v_class = random.choice(variant_classifications)
                impact = "LOW" if v_class == "Silent" else random.choice(impacts)
                db_id = "novel" if random.random() < 0.8 else f"rs{random.randint(1000000, 99999999)}"
                polyphen = "benign" if random.random() < 0.8 else "possibly_damaging"
                sift = "tolerated" if random.random() < 0.8 else "deleterious"
                
                # Often filtered out by raw caller pipelines
                filt_raw = random.choice(["panel_of_normals", "weak_evidence", "strand_bias", "germline"])
                label = 0
                
            else:
                # ── Category 3: Germline Variant (Leakage) ──
                # High VAF in BOTH tumor and normal (heterozygous ~50%, homozygous ~100%)
                is_homo = random.random() < 0.15
                target_vaf = 1.0 if is_homo else 0.50
                
                t_depth = random.randint(100, 500)
                t_alt = int(round(t_depth * random.uniform(target_vaf - 0.1, min(1.0, target_vaf + 0.1))))
                t_ref = t_depth - t_alt
                
                n_depth = random.randint(50, 300)
                n_alt = int(round(n_depth * random.uniform(target_vaf - 0.1, min(1.0, target_vaf + 0.1))))
                n_ref = n_depth - n_alt
                
                # Strand bias: balanced
                t_alt_fwd = int(round(t_alt * 0.5))
                t_alt_rev = t_alt - t_alt_fwd
                t_ref_fwd = int(round(t_ref * 0.5))
                t_ref_rev = t_ref - t_ref_fwd
                
                # Biology
                gene = random.choice(other_genes)
                v_class = random.choice(variant_classifications)
                impact = "LOW" if v_class == "Silent" else "MODERATE"
                db_id = f"rs{random.randint(1000000, 99999999)}" # germline variants are mostly known in dbSNP
                polyphen = "benign" if random.random() < 0.9 else "possibly_damaging"
                sift = "tolerated" if random.random() < 0.9 else "deleterious"
                
                filt_raw = "germline" if random.random() < 0.9 else "PASS"
                label = 0
            
            variants.append({
                "Hugo_Symbol": gene,
                "Chromosome": chrom,
                "Start_Position": pos,
                "Reference_Allele": ref_allele,
                "Tumor_Seq_Allele2": alt_allele,
                "Variant_Classification": v_class,
                "Variant_Type": "SNP",
                "Tumor_Sample_Barcode": sample_barcode,
                "FILTER_raw": filt_raw,
                "dbSNP_RS": db_id,
                "IMPACT": impact,
                "PolyPhen": polyphen,
                "SIFT": sift,
                "t_depth": t_depth,
                "t_alt_count": t_alt,
                "t_ref_count": t_ref,
                "n_depth": n_depth,
                "n_alt_count": n_alt,
                "n_ref_count": n_ref,
                "t_alt_count_forward": t_alt_fwd,
                "t_alt_count_reverse": t_alt_rev,
                "t_ref_count_forward": t_ref_fwd,
                "t_ref_count_reverse": t_ref_rev,
                "label": label
            })
            
        df = pd.DataFrame(variants)
        
        # Write MAF format file
        maf_path = f"{out_prefix}.maf"
        with open(maf_path, "w") as f:
            f.write("#version 2.4\n")
            f.write("#GDC Masked Somatic Mutation simulated dataset\n")
            df.to_csv(f, sep="\t", index=False)
            
        # Write VCF format file
        vcf_path = f"{out_prefix}.vcf"
        with open(vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=ctDNASomaticMutationSimulator\n")
            f.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene Symbol\">\n")
            f.write("##INFO=<ID=VC,Number=1,Type=String,Description=\"Variant Classification\">\n")
            f.write("##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"VEP Functional Impact\">\n")
            f.write("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Depth counts: RefFwd, RefRev, AltFwd, AltRev\">\n")
            f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
            f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth: Ref, Alt\">\n")
            f.write("##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Strand bias counts\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n")
            
            for index, row in df.iterrows():
                info = f"GENE={row['Hugo_Symbol']};VC={row['Variant_Classification']};IMPACT={row['IMPACT']};DP4={int(row['t_ref_count_forward'])},{int(row['t_ref_count_reverse'])},{int(row['t_alt_count_forward'])},{int(row['t_alt_count_reverse'])}"
                
                # Tumor Genotype
                t_gt = "0/1"
                t_format = f"{t_gt}:{int(row['t_depth'])}:{int(row['t_ref_count'])},{int(row['t_alt_count'])}:{int(row['t_ref_count_forward'])},{int(row['t_ref_count_reverse'])},{int(row['t_alt_count_forward'])},{int(row['t_alt_count_reverse'])}"
                
                # Normal Genotype
                n_gt = "0/0" if row['n_alt_count'] == 0 else "0/1"
                n_format = f"{n_gt}:{int(row['n_depth'])}:{int(row['n_ref_count'])},{int(row['n_alt_count'])}:0,0,0,0"
                
                f.write(f"{row['Chromosome']}\t{row['Start_Position']}\t{row['dbSNP_RS']}\t{row['Reference_Allele']}\t{row['Tumor_Seq_Allele2']}\t100\t{row['FILTER_raw']}\t{info}\tGT:DP:AD:SB\t{t_format}\t{n_format}\n")
                
        logger.info(f"Generated simulated datasets:\n  MAF: {maf_path}\n  VCF: {vcf_path}")
        return maf_path, vcf_path
