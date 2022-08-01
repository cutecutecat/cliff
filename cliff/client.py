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
def rug_mut(filename: str, s: str, f: str, w: str, v: int, c: str):
    click.echo('[Mutation] Dataset -> [Ruggness] cauculation')
    click.echo('file: {}'.format(filename))
    click.echo('[Mutation] Args:')
    click.echo(
        'mutation label: [{}], fitness label:[{}], offset:[{}]'.format(s, f, v))
    click.echo('variables: [{}]'.format(c))
    click.echo('wile type: [{}]'.format(w))

    args = MutArgs()
    args.mutation_label = s
    args.fitness_label = f
    args.wile_type = w
    args.vt_offset = v
    scenery = MutParser.parse(filename, args)
    meta = MetaData(scenery, c)

    calculator = Ruggness(meta)
    rug = calculator.calculate()
    click.echo("Ruggness: {}".format(rug))


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('-s', help='sequence label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
def rug_seq(filename: str, s: str, f: str, c: str):
    click.echo('[Sequence] Dataset -> [Epistasis] cauculation')
    click.echo('file: {}'.format(filename))
    click.echo('[Sequence] Args:')
    click.echo('sequence label: [{}], fitness label:[{}]'.format(s, f))
    click.echo('variables: [{}]'.format(c))

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
def epi_mut(filename: str, s: str, f: str, w: str, v: int, c: str, o: int):
    click.echo('[Mutation] Dataset -> [Epistasis] cauculation')
    click.echo('file: {}'.format(filename))
    click.echo('[Mutation] Args:')
    click.echo(
        'mutation label: [{}], fitness label:[{}], offset:[{}]'.format(s, f, v))
    click.echo('variables: [{}]'.format(c))
    click.echo('wile type: [{}]'.format(w))
    click.echo('[Epistasis] Args:')
    click.echo('order range: [1-{}]'.format(o))

    args = MutArgs()
    args.mutation_label = s
    args.fitness_label = f
    args.wile_type = w
    args.vt_offset = v
    scenery = MutParser.parse(filename, args)

    calculator = Epistasis(scenery, o, c)
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
def epi_seq(filename: str, s: str, f: str, c: str, o: int):
    click.echo('[Sequence] Dataset -> [Epistasis] cauculation')
    click.echo('file: {}'.format(filename))
    click.echo('[Sequence] Args:')
    click.echo('sequence label: [{}], fitness label:[{}]'.format(s, f))
    click.echo('variables: [{}]'.format(c))
    click.echo('order range: [1-{}]'.format(o))

    args = SeqArgs()
    args.sequence_label = s
    args.fitness_label = f
    scenery = SeqParser.parse(filename, args)

    calculator = Epistasis(scenery, o, c)
    epi = calculator.calculate()

    show_model = calculator.to_draw(epi)
    fig = show_model.plot()
    fig.savefig('output.png')
    click.echo('Epistasis probability: saved to output.png')


if __name__ == '__main__':
    cli()
