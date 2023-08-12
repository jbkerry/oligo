import functools
from importlib import metadata

import click

from .wrappers import run_capture, run_tiled, run_off_target

@click.group(invoke_without_command=True)
@click.option("--version", is_flag=True)
@click.option("--config", "-cfg", type=click.Path(exists=True), required=True)
@click.pass_context
def cli(ctx, version, config):
    if version:
        click.echo(f"oligo v{metadata.version('oligo-capture')}")
    elif ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
    else:
        ctx.obj["CONFIG"] = config

def common_options(f):
    @click.option(
        "--fasta",
        "-f",
        type=click.Path(exists=True),
        help="Path to reference genome fasta.",
        required=True,
    )
    @click.option(
        "--genome",
        "-g",
        type=click.Choice(["hg18", "hg19", "hg38", "mm9", "mm10"]),
        help="Genome build",
        required=True,
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options

def final_common_options(f):
    @click.option(
        "--oligo",
        "-o",
        type=click.INT,
        default=70,
        help="The size (in bp) of the oligo to design, default=70"
    )
    @click.option(
        "--blat",
        is_flag=True,
        help="Detect off-targets using BLAT instead of STAR.",
    )
    @click.option(
        "--star-index",
        "-s",
        type=click.Path(exists=True),
        help="Path to STAR index directory. Omit this option if running with BLAT (--blat)"
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options

@cli.command()
@common_options
@click.option(
    "--bed",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to bed file with capture viewpoint coordinates"
)
@click.option(
    "--enzyme",
    "-e",
    type=click.Choice(['DpnII', 'NlaIII', 'HindIII']),
    default="DpnII",
    help="Name of restriction enzyme, default=DpnII"
)
@final_common_options
@click.pass_context
def capture(ctx, **kwargs):
    click.echo('capture')
    run_capture(ctx.obj["CONFIG"], **kwargs)

@cli.command()
@common_options
@click.option(
    "--chr",
    "-c",
    "chrom",
    required=True,
    type=click.STRING,
    help="Chromosome number/letter on which to design the oligos.",
)
@click.option(
    "--region",
    "-r",
    type=click.STRING,
    help=(
        "The region in which to design the oligos; must be in the format 'start-stop' e.g. 10000-20000. Omit this "
        "option to design oligos across the entire chromosome."
    ),
)
@click.option(
    "--contig",
    is_flag=True,
    help="Run in contiguous mode (restriciton enzyme independent)."
)
@click.option(
    "--enzyme",
    "-e",
    type=click.Choice(['DpnII', 'NlaIII', 'HindIII']),
    default="DpnII",
    help="Name of restriction enzyme, default=DpnII. Omit this option if running in contiguous mode (--contig)"
)
@click.option(
    "--step-size",
    "-t",
    type=click.INT,
    default=70,
    help=(
        "Step size (in bp) between adjacent oligos when running in contiguous mode (--contig), default=70. Omit this "
        "option if you are not using the --contig flag"
    ),
)
@final_common_options
@click.pass_context
def tiled(ctx, **kwargs):
    click.echo('tiled')
    run_tiled(ctx.obj["CONFIG"], **kwargs)

@cli.command("off-target")
@common_options
@click.option(
    "--bed",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to bed file with capture viewpoint coordinates"
)
@click.option(
    "--step-size",
    "-t",
    type=click.INT,
    default=10,
    help="Step size (in bp) between adjacent oligos, default=10",
)
@click.option(
    "--max-dist",
    "-m",
    type=click.INT,
    default=200,
    help="The maximum distance away from the off-target site to design oligos to, default=200"
)
@final_common_options
@click.pass_context
def off_target(ctx, **kwargs):
    click.echo('off-target')
    run_off_target(ctx.obj["CONFIG"], **kwargs)