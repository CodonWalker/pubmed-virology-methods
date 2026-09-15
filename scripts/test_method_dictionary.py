"""Regression tests for the frozen virology-method dictionary."""
from __future__ import annotations

import render_virology_methods_updated as methods


OVERVIEW_CASES = {
    "pcr": "A reverse-transcription PCR assay with LAMP-based detection was evaluated.",
    "serology": "We measured antibodies using an enzyme-linked immunosorbent assay.",
    "transcriptomics": "Host responses were profiled by RNA-seq and transcriptomics.",
    "ai_all": "A machine learning classifier predicted viral host range.",
    "phylogenetics": "Phylogenetic and phylodynamic analyses reconstructed transmission.",
    "animal_models": "The candidate was tested in a hamster model.",
    "genome_seq": "Whole-genome sequencing was performed using an NGS platform.",
    "cell_culture": "Virus was propagated in primary cell culture.",
    "metagenomics": "Metagenomic sequencing characterized the virome.",
    "viral_vectors": "An AAV-based vector delivered the transgene.",
    "molecular_modeling": "Molecular docking and molecular dynamics identified ligands.",
    "cytometry": "Flow cytometry was used for immunophenotyping.",
    "proteomics": "Quantitative proteomics mapped the interactome.",
    "single_cell": "Single-cell RNA sequencing resolved infected populations.",
    "rna_vaccines": "A nucleoside-modified mRNA–based vaccination strategy was tested.",
    "math_epi": "A compartmental SEIR model estimated transmission.",
    "systems_biology": "A systems biology multi-omics analysis was conducted.",
    "metabolomics": "Untargeted metabolomics and lipidomics were performed.",
    "biosensors": "A microfluidic electrochemical biosensor detected antigen.",
    "vlp_nanovaccines": "A virus-like particle vaccine induced neutralizing antibodies.",
    "crispr_diagnostics": "The SHERLOCK platform enabled CRISPR-based detection.",
    "cryo_em": "Cryo-electron microscopy resolved the capsid structure.",
    "wastewater": "Wastewater-based epidemiology supported sewage surveillance.",
    "organoids": "Airway organoids were cultured in an organ-on-a-chip device.",
    "cell_microscopy": "Live-cell imaging and confocal microscopy tracked entry.",
    "epigenomics": "Epitranscriptomic m6A modification was profiled.",
    "long_reads": "Oxford Nanopore direct RNA sequencing produced long reads.",
    "crystallography_nmr": "X-ray crystallography and NMR spectroscopy determined structure.",
    "annotation_homology": "Genome annotation used multiple sequence alignment and BLAST searches.",
    "spatial_omics": "Spatial transcriptomics mapped infected tissue.",
    "compound_screening": "High-throughput compound screening identified inhibitors.",
    "immune_repertoire": "T-cell receptor repertoire sequencing measured clonal expansion.",
    "deep_mutational": "Multiplexed assays of variant effects enabled deep mutational scanning.",
    "functional_screens": "A genome-wide CRISPR screen identified host factors.",
}

AI_CASES = {
    "classic_ml": "A random forest machine learning model was trained.",
    "deep_learning": "A transformer-based model and CNN classifier were evaluated.",
    "sequence_lm": "A protein foundation model and ESM-2 embeddings were used.",
    "text_llm": "GPT-4o and a large language model summarized the reports.",
    "ai_structure": "AlphaFold predicted the protein structure.",
    "nlp": "Natural language processing and named entity recognition mined abstracts.",
    "generative": "A denoising diffusion probabilistic model generated candidate sequences.",
}

NEGATIVE_CASES = {
    "LAMP-1 lysosomal membrane protein expression was measured.": "pcr",
    "Replication protein A (RPA) bound the viral genome.": "pcr",
    "Messenger RNA expression increased after infection.": "rna_vaccines",
    "A random diffusion model described particle motion.": "generative",
    "The CNN news report discussed the outbreak.": "deep_learning",
}


def main() -> int:
    for key, text in OVERVIEW_CASES.items():
        overview, _ = methods.classify(text)
        assert key in overview, (key, text, overview)
    for key, text in AI_CASES.items():
        overview, ai = methods.classify(text)
        assert key in ai and "ai_all" in overview, (key, text, overview, ai)
    for text, forbidden in NEGATIVE_CASES.items():
        overview, ai = methods.classify(text)
        assert forbidden not in overview | ai, (forbidden, text, overview, ai)
    print(f"OK: {len(OVERVIEW_CASES)} overview, {len(AI_CASES)} AI, {len(NEGATIVE_CASES)} negative cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
