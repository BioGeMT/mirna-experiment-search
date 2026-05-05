"""
Built-in controlled vocabulary for miRNA experiment search.

Users choose canonical terms such as:
  RNA-seq
  CLIP
  CLASH
  qPCR
  overexpression
  knockdown

The pipeline expands these internally into common synonyms.
"""


EXPERIMENT_TYPES = {
    "RNA-seq": [
        "RNA-seq",
        "RNA seq",
        "RNA sequencing",
        "mRNA-seq",
        "transcriptome",
        "transcriptomic",
        "gene expression profiling by high throughput sequencing",
    ],
    "small_RNA-seq": [
        "small RNA-seq",
        "small RNA seq",
        "small RNA sequencing",
        "miRNA-seq",
        "microRNA-seq",
        "miRNA sequencing",
        "microRNA sequencing",
    ],
    "CLIP": [
        "CLIP-seq",
        "CLIP seq",
        "eCLIP",
        "HITS-CLIP",
        "PAR-CLIP",
        "iCLIP",
        "AGO CLIP",
        "AGO2 CLIP",
        "Argonaute CLIP",
    ],
    "CLASH": [
        "CLASH",
        "AGO CLASH",
        "AGO2 CLASH",
        "Argonaute CLASH",
        "CLEAR-CLIP",
        "CLEAR CLIP",
        "chimeric reads",
        "miRNA target chimera",
        "microRNA target chimera",
    ],
    "qPCR": [
        "qPCR",
        "RT-qPCR",
        "RT qPCR",
        "real-time PCR",
        "quantitative PCR",
        "PCR",
    ],
    "microarray": [
        "microarray",
        "expression array",
        "gene expression array",
        "miRNA array",
        "microRNA array",
    ],
    "proteomics": [
        "proteomics",
        "mass spectrometry",
        "LC-MS",
        "LC MS",
    ],
    "reporter_assay": [
        "luciferase",
        "reporter assay",
        "dual luciferase",
        "3' UTR reporter",
        "3 UTR reporter",
    ],
}


PERTURBATION_TYPES = {
    "overexpression": [
        "overexpression",
        "over-expression",
        "overexpressed",
        "overexpress",
        "mimic",
        "miRNA mimic",
        "microRNA mimic",
        "transfection",
    ],
    "knockdown": [
        "knockdown",
        "knock-down",
        "silencing",
        "silenced",
        "inhibitor",
        "miRNA inhibitor",
        "microRNA inhibitor",
        "antagomir",
        "siRNA",
        "shRNA",
    ],
    "knockout": [
        "knockout",
        "knock-out",
        "KO",
        "CRISPR",
        "deletion",
        "deleted",
    ],
    "perturbation_any": [
        "overexpression",
        "over-expression",
        "mimic",
        "knockdown",
        "knock-down",
        "inhibitor",
        "antagomir",
        "knockout",
        "knock-out",
        "CRISPR",
    ],
}


CONTROL_TERMS = [
    "control",
    "ctrl",
    "mock",
    "vehicle",
    "untreated",
    "wild type",
    "wild-type",
    "WT",
    "negative control",
    "scramble",
    "scrambled",
    "non-targeting",
    "non targeting",
]


SPECIES_TERMS = {
    "Homo sapiens": [
        "Homo sapiens",
        "human",
    ],
    "Mus musculus": [
        "Mus musculus",
        "mouse",
    ],
    "Rattus norvegicus": [
        "Rattus norvegicus",
        "rat",
    ],
}


DEFAULT_EXCLUDE_TERMS = [
    "pre-miRNA",
    "pre miRNA",
    "pre-miR",
    "pre miR",
    "pri-miRNA",
    "pri miRNA",
    "pri-miR",
    "pri miR",
    "precursor miRNA",
    "precursor microRNA",
    "primary miRNA",
    "primary microRNA",
]


COMMON_CELL_LINES = {
    "HEK293T": ["HEK293T", "293T"],
    "HEK293": ["HEK293", "HEK-293", "293 cells"],
    "HEK293FT": ["HEK293FT", "293FT"],
    "HeLa": ["HeLa"],
    "HepG2": ["HepG2", "Hep G2"],
    "A549": ["A549"],
    "MCF7": ["MCF7", "MCF-7"],
    "K562": ["K562"],
    "HCT116": ["HCT116"],
    "U2OS": ["U2OS"],
}