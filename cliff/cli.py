import click

from .parser import SeqParser, MutParser, SeqArgs, MutArgs
from .epistasis import Epistasis
from .ruggness import Ruggness
from .metadata import MetaData

@click.group()
def cli():
    pass

@cli.command()
@click.option('-s', help='mutation label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-w', help='wild type sequence of dataset', type=str)
@click.option('-v', help='index offset of dataset', type=int)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
def rug_mut(s:str, f:str, w:str, v:int, c:str):
    click.echo('[Mutation] Dataset -> [Ruggness] cauculation')
    click.echo('[Mutation] Args:')
    click.echo('mutation label: [{}], fitness label:[{}], offset:[{}]'.format(s, f, v))
    click.echo('variables: [{}]'.format(c))
    click.echo('wile type: [{}]'.format(w))

    args = MutArgs()
    args.mutation_label = s
    args.fitness_label = f
    args.wile_type = w
    args.vt_offset = v
    scenery = MutParser.parse(s, args)
    meta = MetaData(scenery, c)

    calculator = Ruggness(meta)
    rug = calculator.calculate()
    click.echo("Ruggness: {}".format(rug))
    

@cli.command()
@click.option('-s', help='sequence label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
def rug_seq(s:str, f:str, c:str):
    click.echo('[Sequence] Dataset -> [Epistasis] cauculation')
    click.echo('[Sequence] Args:')
    click.echo('sequence label: [{}], fitness label:[{}]'.format(s, f))
    click.echo('variables: [{}]'.format(c))

    args = SeqArgs()
    args.sequence_label = s
    args.fitness_label = f
    scenery = SeqParser.parse(s, args)
    meta = MetaData(scenery, c)

    calculator = Ruggness(meta)
    rug = calculator.calculate()
    click.echo("Ruggness: {}".format(rug))
    

@cli.command()
@click.option('-s', help='mutation label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-w', help='wild type sequence of dataset', type=str)
@click.option('-v', help='index offset of dataset', type=int)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
@click.option('-o', help='max order of epistasis calculation', type=int)
def epi_mut(s:str, f:str, w:str, v:int, c:str, o:int):
    click.echo('[Mutation] Dataset -> [Epistasis] cauculation')
    click.echo('[Mutation] Args:')
    click.echo('mutation label: [{}], fitness label:[{}], offset:[{}]'.format(s, f, v))
    click.echo('variables: [{}]'.format(c))
    click.echo('wile type: [{}]'.format(w))
    click.echo('[Epistasis] Args:')
    click.echo('order range: [1-{}]'.format(o))

    args = MutArgs()
    args.mutation_label = s
    args.fitness_label = f
    args.wile_type = w
    args.vt_offset = v
    scenery = MutParser.parse(s, args)

    calculator = Epistasis(scenery, o, c)
    epi = calculator.calculate()
    prob_model = calculator.to_prob(epi)
    table = prob_model.table_format()
    click.echo('Epistasis probability:')
    click.echo(table)

@cli.command()
@click.option('-s', help='sequence label of csv file', type=str)
@click.option('-f', help='fitness label of csv file', type=str)
@click.option('-c', help='input variables for sequence', default='ACDEFGHIKLMNPQRSTVWY', type=str)
@click.option('-o', help='max order of epistasis calculation', type=int)
def epi_seq(s:str, f:str, c:str, o: int):
    click.echo('[Sequence] Dataset -> [Epistasis] cauculation')
    click.echo('[Sequence] Args:')
    click.echo('sequence label: [{}], fitness label:[{}]'.format(s, f))
    click.echo('variables: [{}]'.format(c))
    click.echo('order range: [1-{}]'.format(o))
    
    args = SeqArgs()
    args.sequence_label = s
    args.fitness_label = f
    scenery = SeqParser.parse(s, args)

    calculator = Epistasis(scenery, o, c)
    epi = calculator.calculate()
    prob_model = calculator.to_prob(epi)
    table = prob_model.table_format()
    click.echo('Epistasis probability:')
    click.echo(table)

if __name__ == '__main__':
    cli()