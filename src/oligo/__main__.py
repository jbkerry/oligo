import click

@click.group()
@click.option("--config", "-cfg", type=click.Path(exists=True))
@click.option("--genome", "-g", type=click.Choice(["hg18", "hg19", "hg38", "mm9", "mm10"]))
@click.option("--fasta", "-f", type=click.Path(exists=True))
@click.pass_context
def cli(ctx, config, genome, fasta):
    pass

@cli.command()
@click.option("--bed", "-b", required=True, type=click.STRING)
def capture(bed):
    click.echo('capture')

@cli.command()
def tiled():
    click.echo('tiled')

@cli.command("off-target")
def off_target():
    click.echo('tiled')

# print("Hello, World!")
if __name__ == "__main__":
    cli(obj={})