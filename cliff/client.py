import click

from .parser import SeqParser, MutParser, SeqArgs, MutArgs
from .epistasis import Epistasis
from .ruggness import Ruggness
from .metadata import MetaData


@click.group()
def cli():
    pass


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('-s', help='mutation label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-w', help='wild type sequence of dataset', type=str)
@click.option('-v', help='index offset of dataset', type=int, default=0)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
def rug_mut(filename: str, symbol: str, fitness: str, wild_type: str, vt_offset: int, chars: str):
    click.echo('[Mutation] Dataset -> [Ruggness] cauculation')
    click.echo(f'file: {filename}')
    click.echo('[Mutation] Args:')
    click.echo(
        f'mutation label: [{symbol}], fitness label:[{fitness}], offset:[{vt_offset}]')
    click.echo(f'variables: [{chars}]')
    click.echo(f'wile type: [{wild_type}]')

    args = MutArgs()
    args.mutation_label = symbol
    args.fitness_label = fitness
    args.wile_type = wild_type
    args.vt_offset = vt_offset
    scenery = MutParser.parse(filename, args)
    meta = MetaData(scenery, chars)

    calculator = Ruggness(meta)
    rug = calculator.calculate()
    click.echo("Ruggness: {}".format(rug))


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('-s', help='sequence label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
def rug_seq(filename: str, symbol: str, fitness: str, chars: str):
    click.echo('[Sequence] Dataset -> [Epistasis] cauculation')
    click.echo(f'file: {filename}')
    click.echo('[Sequence] Args:')
    click.echo(f'sequence label: [{symbol}], fitness label:[{fitness}]')
    click.echo(f'variables: [{chars}]')

    args = SeqArgs()
    args.sequence_label = s
    args.fitness_label = f
    scenery = SeqParser.parse(filename, args)
    meta = MetaData(scenery, c)

    calculator = Ruggness(meta)
    rug = calculator.calculate()
    click.echo("Ruggness: {}".format(rug))


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('-s', help='mutation label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-w', help='wild type sequence of dataset', type=str)
@click.option('-v', help='index offset of dataset', type=int, default=0)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
@click.option('-o', help='max order of epistasis calculation', type=int)
def epi_mut(filename: str, symbol: str, fitness: str, wild_type: str, vt_offset: int, chars: str, max_order: int):
    click.echo('[Mutation] Dataset -> [Epistasis] cauculation')
    click.echo(f'file: {filename}')
    click.echo('[Mutation] Args:')
    click.echo(
        f'mutation label: [{symbol}], fitness label:[{fitness}], offset:[{vt_offset}]')
    click.echo(f'variables: [{chars}]')
    click.echo(f'wile type: [{wild_type}]')
    click.echo('[Epistasis] Args:')
    click.echo(f'order range: [1-{max_order}]')

    args = MutArgs()
    args.mutation_label = symbol
    args.fitness_label = fitness
    args.wile_type = wild_type
    args.vt_offset = vt_offset
    scenery = MutParser.parse(filename, args)

    calculator = Epistasis(scenery, max_order, chars)
    epi = calculator.calculate()

    show_model = calculator.to_draw(epi)
    fig = show_model.plot()
    fig.savefig('output.png')
    click.echo('Epistasis probability: saved to output.png')


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('-s', help='sequence label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
@click.option('-o', help='max order of epistasis calculation', type=int)
def epi_seq(filename: str, symbol: str, fitness: str, chars: str, max_order: int):
    click.echo('[Sequence] Dataset -> [Epistasis] cauculation')
    click.echo(f'file: {filename}')
    click.echo('[Sequence] Args:')
    click.echo(f'sequence label: [{symbol}], fitness label:[{fitness}]')
    click.echo(f'variables: [{chars}]')
    click.echo(f'order range: [1-{max_order}]')

    args = SeqArgs()
    args.sequence_label = symbol
    args.fitness_label = fitness
    scenery = SeqParser.parse(filename, args)

    calculator = Epistasis(scenery, max_order, chars)
    epi = calculator.calculate()

    show_model = calculator.to_draw(epi)
    fig = show_model.plot()
    fig.savefig('output.png')
    click.echo('Epistasis probability: saved to output.png')


if __name__ == '__main__':
    cli()
