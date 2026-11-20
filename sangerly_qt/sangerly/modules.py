"""The steps of SangeRly in workflow order. `key` is the id of the tab in the Shiny app (app.R), `logo` the picture in www/img/.
Every step runs in the Shiny workspace inside the Qt window, so the tools are exactly the ones of the Shiny program."""
from dataclasses import dataclass, field


@dataclass
class Module:
    key: str
    title: str
    short: str                      # text under the ribbon icon
    step: int
    color: str
    glyph: str                      # fallback icon name in theme.icon (used when the logo picture is missing)
    summary: str
    logo: str = ""
    details: list = field(default_factory=list)
    ready: bool = True


MODULES = [
    Module("assembly", "Assembly", "ASSEMBLY", 1, "#2f6fed", "assembly",
           "Assemble Forward + Reverse .ab1 reads, check the chromatograms, edit the consensus.", "logo_assembly.png"),
    Module("blast", "BLAST", "BLAST", 2, "#0ea5a4", "blast",
           "Identify the species: search the consensus against the NCBI databases.", "logo_blast.png"),
    Module("msa", "MSA and trimming", "MSA", 3, "#8b5cf6", "align",
           "Align many sequences and trim the poorly aligned columns.", "logo_msa.png"),
    Module("tree", "Phylogeny", "TREE", 4, "#16a34a", "tree",
           "NJ and maximum-likelihood (IQ-TREE) trees, interactive viewer, tree comparison.", "logo_phylo.png"),
    Module("vis", "Visualisation", "VISUAL", 5, "#0891b2", "tree",
           "Prepare Newick files for iTOL and ChiPlot and draw publication figures.", "logo_visual.png"),
    Module("popgen", "Population structure", "POPULATION", 6, "#ef4444", "popgen",
           "Haplotypes, diversity, Tajima's D, Fu's Fs, networks, AMOVA / PhiST, export for DnaSP, Arlequin and PopART.", "logo_popup.png"),
    Module("barcode", "Genetic distances", "DISTANCES", 7, "#475569", "barcode",
           "Distance matrices, barcode gap, PCoA, statistical tests and species identification.", "logo_barcode.png"),
    Module("primer", "Primer mining and design", "PRIMER", 8, "#d99a06", "primer",
           "Find target sequences in BOLD and NCBI, design primers, check their quality and BLAST them.", "logo_primer.png"),
]
BY_KEY = {m.key: m for m in MODULES}

# the "More" tools (tab ids as in app.R)
EXTRAS = [
    ("open", "Open and view files", "Just look at files: AB1 chromatograms, FASTA, FASTQ, GenBank, Newick trees and tables."),
    ("project", "Project and report", "Save or open the whole project, marker presets, Methods paragraph and HTML report."),
    ("genbank", "GenBank export", "FASTA with source modifiers and a feature table for BankIt / table2asn."),
    ("tutorial", "Tutorial", "Step-by-step lessons for every tool."),
    ("help", "Help", "Installation, use, antivirus, slow start and other problems."),
    ("about", "About and citation", "Versions of the programs, how to cite, and the references of the methods."),
]
TAB_TITLES = {m.key: m.title for m in MODULES} | {k: t for k, t, _ in EXTRAS} | {"home": "Home"}


# Entries of the drop-down menu of every tool (menu text, id or title of the tab inside the Shiny tool) -- "-" is a separator.
# A list opens nested tabs one after the other. The MSA alignment, trimming and tree dialogs are in dialogs.py.
SECTIONS = {
    "assembly": [("Input files", "Input files"), ("Samples", "Samples"), ("Overview", "Overview"), "-",
                 ("Alignment", "Alignment"), ("Read and QC", "Read and QC"), ("Needs review", "Needs review"),
                 ("SNPs and variants", "SNPs and variants"), ("Translation", "Translation"), "-",
                 ("Sequence analysis", "Sequence analysis"), ("Distances and dendrogram", "Distances and dendrogram"),
                 ("Consensus and export", "Consensus and export")],
    "blast": [("Descriptions", "Descriptions"), ("Graphic Summary", "Graphic Summary"), ("Alignments", "Alignments"), ("Taxonomy", "Taxonomy")],
    "msa": [("Editor", "Editor"), ("Sequences", "Sequences"), ("Trimmed result", "Trimmed result"), "-",
            ("Mutations", "Mutations"), ("Site statistics", "Site statistics"), ("Columns", "Columns")],
    "tree": [("NJ tree", "NJ tree"), ("IQ-TREE tree", "IQ-TREE tree"), ("ModelFinder", "ModelFinder"), "-",
             ("Interactive view", "Interactive view"), ("Compare trees", "Compare trees"), ("Distances", "Distances"), "-",
             ("IQ-TREE log and report", "IQ-TREE log and report"), ("Save results", "Save results")],
    "vis": [("Tree", "Tree"), ("Edit names", "Edit names"), ("External viewers (iTOL, ChiPlot)", "External viewers")],
    "popgen": [("Haplotype", "Haplotype"), ("Network", "Network"), ("Map", "Map"), "-",
               ("AMOVA", "AMOVA"), ("Structure", "Structure"), ("Sliding window", "Sliding window"), "-",
               ("SFS and linkage (LD)", ["More tests", "SFS and linkage (LD)"]), ("Mantel test", ["More tests", "Mantel test"]),
               ("Hardy-Weinberg", ["More tests", "Hardy-Weinberg"]), ("Phasing (EM)", ["More tests", "Phasing (EM)"]),
               ("Fst outliers (FDIST)", ["More tests", "Fst outliers (FDIST)"]), "-", ("Result", "Result")],
    "barcode": [("Result (distances and group means)", "Result"), "-", ("Distribution", "Distribution"), ("Box plot", "Box plot"), ("Statistics", "Statistics"), ("PCoA", "PCoA"), "-",
                ("Distance matrix", "Distance matrix"), ("Per species", "Per species"), ("Max intra vs nearest inter", "Max intra vs nearest inter"),
                ("Threshold", "Threshold"), "-", ("Identification", "Identification"), ("Tree", "Tree"), ("Pairs", "Pairs")],
    "primer": [("Find sequences and primers", "mine"), ("Design", "design"), ("Quality", "ana"), ("BLAST", "blast")],
}


# Tools > Genetic distances > Compute ... (menu text, analysis code of the R page; the codes are in R/barcode_module.R, .bc_analyses)
DIST_ANALYSES = [
    ("Compute pairwise distances", "pairwise"),
    ("Compute overall mean distance", "overall"),
    ("Compute within group mean distance", "within"),
    ("Compute between group mean distance", "between"),
    ("Compute net between group mean distances", "net"),
    ("Barcode gap and species identification", "barcode"),
]


def sub_id(sub):
    """The id of a section: the last title of a nested list (popgen: [\"More tests\", \"Mantel test\"])."""
    return sub[-1] if isinstance(sub, (list, tuple)) else sub


# (tool, section id) -> menu text; used for the titles of the document tabs
SUB_LABELS = {(k, sub_id(e[1])): e[0] for k, items in SECTIONS.items() for e in items if e != "-"}
