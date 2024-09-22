from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter

from wf.pathway_enrichment import Reference_Type, pathway_enrichment

metadata = LatchMetadata(
    display_name="Pathway Enrichment",
    documentation="https://wiki.latch.bio/workflows/bulk-rna-seq#bulk-rnaseq-quantification",
    wiki_url="https://wiki.latch.bio/workflows/bulk-rna-seq#bulk-rnaseq-quantification",
    author=LatchAuthor(
        name="LatchBio",
        github="https://github.com/latchbio-workflows/wf-verified-pathway-enrichment",
    ),
    parameters={
        "run_name": LatchParameter(
            display_name="Name of Run",
            description="Name of run.",
            batch_table_column=True,
        ),
        "contrast_file": LatchParameter(
            display_name="Input Contrast File",
            description="CSV file containing log-fold change, p-values, etc. from a comparison between conditions from DESeq2.",
            batch_table_column=True,
        ),
        "genome": LatchParameter(
            display_name="Reference Genome",
            description="Reference genome for run.",
            batch_table_column=True,
        ),
        "output_directory": LatchParameter(
            display_name="Output Directory",
            description="The output directory where the results will be saved.",
            batch_table_column=True,
        ),
    },
)


@workflow(metadata)
def gene_ontology_pathway_analysis(
    run_name: str,
    contrast_file: LatchFile,
    genome: Reference_Type,
    output_directory: LatchOutputDir = LatchOutputDir("latch:///Pathway_Enrichment"),
) -> LatchOutputDir:
    """Implicate pathways and gene ontologies from a sample contrast comparison of differential gene results

    <html>
    <p align="center">
      <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    <p align="center">
      Implicate pathways and gene ontologies from a set of desired genes.
    </p>

    </html>

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
    gene_ontology_pathway_analysis,
    "Human Test Data",
    {
        "run_name": "human_pe_demo",
        "contrast_file": LatchFile(
            "s3://latch-public/nf-core/rnaseq/launchplan_data/human_pe_demo/deseq2_contrast_results.csv"
        ),
        "genome": Reference_Type.homo_sapiens,
        "output_directory": LatchOutputDir("latch:///Pathway_Enrichment"),
    },
)
