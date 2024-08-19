from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter

from wf.pathway_enrichment import Reference_Type, pathway_enrichment

metadata = LatchMetadata(
    display_name="Pathway Enrichment",
    author=LatchAuthor(
        name="LatchBio",
    ),
    parameters={
        "run_name": LatchParameter(
            display_name="Name of Run",
            batch_table_column=True,
        ),
        "contrast_file": LatchParameter(
            display_name="Input Contrast File",
            batch_table_column=True,
        ),
        "genome": LatchParameter(
            display_name="Genome",
            batch_table_column=True,
        ),
        "output_directory": LatchParameter(
            display_name="Output Directory",
            batch_table_column=True,
        ),
    },
)


@workflow(metadata)
def pathway_enrichment_wf(
    run_name: str,
    contrast_file: LatchFile,
    genome: Reference_Type,
    output_directory: LatchOutputDir = LatchOutputDir("latch:///Pathway_Enrichment"),
) -> LatchOutputDir:
    """Pathway Enrichment

    ## Pathway Enrichment

    Under the hood, the workflow uses clusterprofiler for Gene Set Enrichment Analysis (GSEA) and Over-representation analysis (ORA) with mSigDB, KEGG, WikiPathway, Disease Ontology (DO) and Gene Ontology (GO) pathway databases.
    For more information, see our
    [wiki]("https://wiki.latch.bio/workflows/bulk-rna-seq#bulk-rnaseq-quantification").
    """

    return pathway_enrichment(
        run_name=run_name,
        contrast_file=contrast_file,
        genome=genome,
        output_directory=output_directory,
    )


LaunchPlan(
    pathway_enrichment_wf,
    "Human Test Data",
    {
        "run_name": "human_pe_demo",
        "contrast_file": LatchFile(
            "s3://latch-public/test-data/18440/human_pe_demo/deseq2_contrast_results.csv"
        ),
        "genome": Reference_Type.homo_sapiens,
        "output_directory": LatchOutputDir("latch:///Pathway_Enrichment"),
    },
)
