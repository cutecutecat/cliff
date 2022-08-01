# cliff -- Epistasis and Ruggness calculation

[中文](/docs/README_CN.md) | [English](/README.md)

A Python API for calculating epistasis and ruggness of mutation dataset. We provide two ways of API calling or argument client.

Ruggedness is a measure of how far a dataset deviates from linear relationship. From a macro view, it is a measurement of space ruggedness of the fitness in the sequence space. The higher the ruggedness, the more complex patterns the dataset contains, and it is difficult to learn the relationship between its fitness relative to the sequence.

By our improved algorithm, the ruggedness can still be accurate enough from the original ruggedness when there are a large number of deletions in the mutation database.

Epistasis are used to measure the direction and magnitude of the relative phenotype effects of multiple inputs. The traditional Epistasis algorithm use input of mutations, and the mixed effect of multiple mutations is calculated.

Our algorithm provides a new perspective, with input in residue, it can calculate the mixed effect of multiple residue. We believe this will provide more information than traditional algorithms.

for 
```python
sequence = ["AAA", "AAT", "ATA", "TAA", 
            "ATT", "TAT", "TTA", "TTT"]
fitness = [0.1, 0.2, 0.4, 0.3, 0.3, 0.6, 0.8, 1.0]
```
The Epistasis:

![img](/docs/epistasis.png)

# install

Run
```shell
pip install git+https //github.com/cutecutecat/cliff
```

# usage

## use as a library

when calculating Ruggness:

```python
from cliff import MetaData, Ruggness
from cliff.parser import SeqArgs, SeqParser

# generate args for calculation
args = SeqArgs()
args.sequence_label = 'sequence'
args.fitness_label = 'fitness'
scenery = SeqParser.parse('input.csv', args)
meta = MetaData(scenery, 'ABCDEFGHI')

# begin calculation
calculator = Ruggness(meta)
rug = calculator.calculate()
```

when calculating Epistasis:

```python
from cliff import Epistasis
from cliff.parser import SeqArgs, SeqParser

# generate args for calculation
args = SeqArgs()
args.sequence_label = 'sequence'
args.fitness_label = 'fitness'
scenery = SeqParser.parse('input.csv', args)

# begin calculation
calculator = Epistasis(scenery, 3, 'ABCDEFGHI')
epi = calculator.calculate()

# plot and save result
show_model = calculator.to_draw(epi)
fig = show_model.plot()
fig.save_fig('output.png')
```

## use as a command line program

refer to help of `cliff --help`

# input file format

Input file should be a `csv` format file, which should at least contain two columns for different parser. We recommend using `sequence` parser for sake of convenience.

## sequence parser

* `sequence` for `str` input sequence like `ACGT`
* `fitness` for `float` sequence fitness like `1.25`

The name of these two columns need not to be `sequence` or `fitness`, tell parser the right name in `SeqArgs` is enough.

file example:
```csv
,Sequence,Fitness
0,AAAAA,0.49338512366434883
1,AAAAC,0.5096570842221034
2,AAAAD,0.5546400893159292
3,AAAAE,0.5265283438717137
4,AAAAF,0.4392058266974471
```

## mutation parser

* `mutation` for `str` input mutation from wild type like `A2T:A3T`
* `fitness` for `float` sequence fitness like `1.25`

The name of these two columns need not to be `mutation` or `fitness`, tell parser the right name in `MutArgs` is enough.

file example:
```csv
variant,score
,0.1
A3T,0.2
A2T,0.4
A1T,0.3
A2T:A3T,0.3
```