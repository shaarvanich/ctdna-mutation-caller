

# LLM Assisted Somatic Mutation Calling for CtDNA Detection in Computational Oncology

A hybrid bioinformatics + LLM pipeline for classifying low-frequency somatic variants in circulating tumor DNA (ctDNA) liquid biopsy data, using large language model reasoning to integrate biological context that conventional threshold-based callers miss.

## Overview

Detecting true cancer mutations in ctDNA is hard because the signal is faint. In early-stage cancer or minimal residual disease (MRD) monitoring, tumor-derived DNA can make up less than 1%, sometimes less than 0.1%, of the total cell-free DNA circulating in a patient's blood. At these frequencies, real mutations are statistically almost indistinguishable from PCR errors, oxidative DNA damage, and other sequencing artifacts.

Conventional variant callers (MuTect2, VarScan2, Strelka2, etc.) handle this by applying fixed VAF (variant allele frequency) thresholds and quality filters. That approach forces an impossible trade-off: thresholds loose enough to catch real low-frequency mutations also let through a flood of false positives, while thresholds tight enough to control noise end up discarding genuine driver mutations.

This project explores a different approach: instead of treating a variant's VAF as a single scalar cutoff, it treats each candidate variant as a small body of evidence (VAF, sequencing depth, strand bias, gene identity, hotspot status, mutational context) and asks an LLM to reason through that evidence the way a molecular pathologist would, producing both a classification and a written justification for it.

## Why This Project Exists

Statistical filters are good at one thing: rejecting variants based on numbers. They are not good at incorporating biological knowledge, such as "this mutation sits in a well-known cancer hotspot, so it deserves the benefit of the doubt even at a borderline VAF." A pathologist reviewing the same variant brings exactly that kind of contextual judgment, weighing gene identity, known driver status, and clinical history alongside the raw numbers.

This project was built to test whether an LLM can approximate that judgment, acting as a reasoning layer on top of standard variant-calling statistics rather than replacing them, and to produce classifications that come with a human-readable explanation instead of an opaque confidence score. The goal isn't to outperform a deep learning classifier on raw accuracy; it's to build something interpretable, auditable, and biologically grounded, suited to the way ctDNA review actually happens in practice.

## How It Works

The pipeline runs candidate variants through three stages, each handled by the tool best suited to it rather than throwing everything at the LLM:

**1. Statistical pre-filtering**
Cheap, rule-based filters remove the unambiguous noise first: variants with too few supporting reads, extreme strand bias, or evidence of the same variant in the matched normal sample. This step needs no LLM calls and clears out the bulk of artifacts before anything expensive happens.

**2. Driver rescue**
A biologically motivated heuristic re-examines variants that the statistical filter just rejected. If a rejected variant falls in a well-characterized cancer driver gene (cross-referenced against the COSMIC Cancer Gene Census), it gets pulled back in for review instead of being discarded outright. This is the step that lets the pipeline recover real low-VAF driver mutations that a pure threshold filter would always miss.

**3. LLM chain-of-thought classification**
Every variant that survives stages 1–2 is passed to an LLM with a structured prompt encoding its quality metrics (VAF, depth, read counts, strand bias, matched normal VAF) and biological annotations (gene role, driver status, hotspot prevalence, trinucleotide context). The model is prompted to reason step by step: sequencing quality, biological plausibility, mutational signature consistency, clinical context, before assigning the variant to one of four tiers:

- Somatic, high confidence
- Somatic, low confidence
- Germline
- Sequencing artifact

Each classification is returned with a natural-language justification, not just a score, so the result can be read and sanity-checked by a human reviewer.

## Modular Design & Architecture

The pipeline is deliberately split so that LLM calls are reserved for the cases that actually need reasoning, keeping cost and latency manageable:

```
Raw variant calls (VCF / tabular input)
        │
        ▼
┌───────────────────────┐
│ Statistical pre-filter │  → discards clear artifacts (low depth, strand bias, normal contamination)
└───────────────────────┘
        │
        ▼
┌───────────────────────┐
│   Driver rescue        │  → re-admits low-VAF variants in known driver genes (COSMIC CGC)
└───────────────────────┘
        │
        ▼
┌───────────────────────┐
│  LLM CoT classifier    │  → structured prompt → 4-tier classification + written justification
└───────────────────────┘
        │
        ▼
Classified output (tier, confidence, justification) + benchmark/visualization reports
```

Design principles behind the split:
- **Cheap filters first.** Arithmetic-based rejection rules run before any API call, so the LLM only ever sees a pre-screened candidate set.
- **Biological priors are explicit, not learned.** The driver-rescue step encodes domain knowledge (driver gene membership) directly, rather than relying on the LLM to infer it implicitly every time.
- **Reasoning is reserved for ambiguity.** The LLM stage is the only part of the pipeline that requires nuanced, multi-factor judgment, so it's the only part that costs an API call per variant.
- **Caching for reproducibility.** LLM responses are cached against a hash of the variant's input features, so re-running the pipeline on the same data doesn't reissue identical queries or introduce unnecessary run-to-run variability.

## Tech Stack

- **Python**: core pipeline logic and orchestration
- **Pandas / NumPy**: variant table processing, statistical filtering, feature engineering
- **LLM API (chain-of-thought prompting)**: structured reasoning and classification of pre-filtered candidate variants
- **COSMIC Cancer Gene Census**: reference annotation for driver gene status and hotspot prevalence
- Benchmarking utilities for confusion matrix, sensitivity/precision/F1, ROC/PR-AUC, and VAF-tier-stratified performance analysis

## Current Results

On a 2,000-variant synthetic benchmark dataset:

| Metric | Value |
|---|---|
| Sensitivity | 81.3% |
| Precision | 56.3% |
| F1 score | 0.665 |
| 1–5% VAF tier sensitivity | 79.4% |
| 1–5% VAF tier precision | 61.7% (F1 = 0.694) |
| ROC-AUC | 0.500 |
| PR-AUC | 0.544 (vs. 0.550 baseline) |

These numbers reflect a deliberately sensitivity-biased design: in a clinical liquid biopsy context, missing a real driver mutation is far more costly than flagging a false positive that a clinician can filter out on review. The pipeline is explicitly conservative rather than tuned for raw discriminative accuracy.

## Known Limitations

- **Hallucination risk.** The LLM can occasionally fabricate or misstate biological claims (e.g., incorrect hotspot attribution), so outputs should be treated as decision support, not a final diagnostic call.
- **Reproducibility.** Even at low temperature, LLM outputs aren't perfectly deterministic across model versions or hardware; response caching gives run-level reproducibility but not cross-version guarantees.
- **Benchmark gaps.** There is currently no community-standard ground-truth dataset for ultra-low-VAF ctDNA variants, so evaluation relies on synthetic data that can't fully capture real-world sequencing artifacts and tumor heterogeneity.
- **Not clinically validated.** This is a research/prototype pipeline, not a regulatory-cleared diagnostic tool, and is intended to augment expert review rather than replace it.

## What's Next

- **Retrieval-augmented generation (RAG):** dynamically pull current ClinVar, COSMIC, and OncoKB annotations into the reasoning prompt instead of relying solely on the model's parametric knowledge.
- **Protein language model integration:** incorporate evolutionary substitution scores (e.g., ESM-2) as an additional functional-impact signal for variants lacking strong database annotation.
- **Mutational signature context:** feed whole-sample signature decomposition results into the per-variant reasoning step to better separate signature-consistent mutations from artifacts.
- **Agentic reasoning:** move from single-pass prompting toward an agent that can autonomously query databases and refine its classification across multiple steps for ambiguous variants.
- **Improved benchmarking:** build out evaluation against orthogonally validated low-VAF variants (e.g., ddPCR-confirmed) rather than synthetic data alone.
- **Output validation:** add automated post-hoc fact-checking of LLM justifications against canonical databases to catch hallucinated claims before they reach a report.

## Disclaimer

This project is a research prototype for exploring LLM-assisted variant interpretation. It is not validated for clinical use and should not be used to inform patient care decisions.
