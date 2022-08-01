# cliff -- 粗糙度和上位效应计算工具

[中文](/docs/README_CN.md) | [English](/README.md)

用于计算蛋白/基因突变数据库粗糙度和上位效应的Python接口程序，提供调用库和命令行两种使用方式。

粗糙度用于衡量数据集偏离线性关联的程度。在宏观上的表现即为序列空间上适应度崎岖不平的程度。粗糙度越高，数据集越含有复杂的模式，难以去学习其适应度相对序列的关系。

使用我们改进的算法，在计算粗糙度时能够在突变数据库存在大量缺失时，仍然能够接近原始粗糙度。

上位效应用来衡量多个输入相对表现型的影响方向和幅度。传统的上位效应算法是以突变为单位输入，计算多个突变的混合作用。

我们的算法提供了一个新视角，以位点为单位输入，可以计算多个位点的混合作用。能够提供比传统算法更多的信息。

对于 
```python
sequence = ["AAA", "AAT", "ATA", "TAA", 
            "ATT", "TAT", "TTA", "TTT"]
fitness = [0.1, 0.2, 0.4, 0.3, 0.3, 0.6, 0.8, 1.0]
```
上位效应

![img](/docs/epistasis.png)

# 安装方法

执行
```shell
pip install git+https //github.com/cutecutecat/cliff
```

# 使用方法

## 用作库调用

计算粗糙度:

```python
from cliff import MetaData, Ruggness
from cliff.parser import SeqArgs, SeqParser

# 生成计算参数
args = SeqArgs()
args.sequence_label = 'sequence'
args.fitness_label = 'fitness'
scenery = SeqParser.parse('input.csv', args)
meta = MetaData(scenery, 'ABCDEFGHI')

# 开始计算
calculator = Ruggness(meta)
rug = calculator.calculate()
```

计算上位效应:

```python
from cliff import Epistasis
from cliff.parser import SeqArgs, SeqParser

# 生成计算参数
args = SeqArgs()
args.sequence_label = 'sequence'
args.fitness_label = 'fitness'
scenery = SeqParser.parse('input.csv', args)

# 开始计算
calculator = Epistasis(scenery, 3, 'ABCDEFGHI')
epi = calculator.calculate()

# 作图并保存结果
show_model = calculator.to_draw(epi)
fig = show_model.plot()
fig.save_fig('output.png')
```

## 用作命令行程序

参考 `cliff --help`

# input file format

输入文件需要符合`csv`格式，对于不同的解析器，至少需要包含两列。我们推荐使用`sequence parser`因为其更加简便。

## sequence 解析器

* `sequence` 格式 `str` 代表输入序列，例如 `ACGT`
* `fitness` 格式 `float` 代表序列适应度，例如 `1.25`

这两列的列名可以任意指定，不局限于 `sequence` 和 `fitness` , 只要在创建 `SeqArgs` 时正确指定参数即可。

示例文件:
```csv
,Sequence,Fitness
0,AAAAA,0.49338512366434883
1,AAAAC,0.5096570842221034
2,AAAAD,0.5546400893159292
3,AAAAE,0.5265283438717137
4,AAAAF,0.4392058266974471
```

## mutation 解析器

* `mutation` 格式 `str` 代表突变对象，例如 `A2T:A3T`
* `fitness` 格式 `float` 代表序列适应度，例如 `1.25`

这两列的列名可以任意指定，不局限于 `mutation` 和 `fitness` , 只要在创建 `MutArgs` 时正确指定参数即可。

示例文件:
```csv
variant,score
,0.1
A3T,0.2
A2T,0.4
A1T,0.3
A2T:A3T,0.3
```