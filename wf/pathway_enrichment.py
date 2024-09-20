import subprocess
import sys
from enum import Enum
from pathlib import Path

from latch.executions import rename_current_execution
from latch.resources.tasks import medium_task
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile

sys.stdout.reconfigure(line_buffering=True)


class Reference_Type(Enum):
    homo_sapiens = "Homo sapiens (RefSeq GRCh38.p14)"
    mus_musculus = "Mus musculus (RefSeq GRCm39)"
    rattus_norvegicus = "Rattus norvegicus (RefSeq GRCr8)"
    # rhesus_macaque = "Macaca mulatta (RefSeq rheMac10/Mmul_10)"
    saccharomyces_cerevisiae = "Saccharomyces cerevisiae (RefSeq R64)"


@medium_task
def pathway_enrichment(
    run_name: str,
    contrast_file: LatchFile,
    genome: Reference_Type,
    output_directory: LatchOutputDir = LatchOutputDir("latch:///Pathway_Enrichment"),
) -> LatchOutputDir:
    rename_current_execution(str(run_name))

    print("Setting up local files and directories")
    local_output_directory = Path(f"/root/output/{run_name}")
    local_output_directory.mkdir(exist_ok=True, parents=True)

    contrast_file_p = Path(contrast_file)

    destination = local_output_directory / contrast_file_p.name
    contrast_file_p.rename(destination)

    contrast_file_p = destination

    pathways = ["GO", "KEGG", "DO", "WP", "MSIG"]

    for pathway in pathways:
        (local_output_directory / f"{pathway}/graphs").mkdir(
            parents=True, exist_ok=True
        )

    print("Running pathway enrichment")
    subprocess.run(
        [
            "Rscript",
            "/root/r_scripts/pathway_enrichment.R",
            str(contrast_file_p),
            str(local_output_directory),
            str(genome.name),
        ],
        check=True,
    )

    print("Uploading results")
    return LatchOutputDir(str("/root/output"), output_directory.remote_path)
