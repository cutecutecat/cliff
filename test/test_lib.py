"""do the unit test of the API calling."""

import unittest
from os.path import join, dirname

from cliff import Ruggness, MetaData, Epistasis
from cliff.parser import SeqArgs, SeqParser, MutArgs, MutParser, Scenery


class TestLibCall(unittest.TestCase):
    """do the unit test of the API calling."""

    def test_load_sequence(self):
        """test load a sequence format csv"""
        args = SeqArgs()
        args.fitness_label = "Fitness"
        args.sequence_label = "Sequence"

        path = join(dirname(__file__), "data/seq.csv")
        scenery = SeqParser.parse(path, args)

        mean = sum(scenery.fitness) / len(scenery.fitness)
        self.assertAlmostEqual(mean, 0.5111, places=3)

    def test_load_mutation(self):
        """test load a mutation format csv"""
        args = MutArgs()
        args.fitness_label = "score"
        args.mutation_label = "variant"
        args.vt_offset = 0
        args.wile_type = "AAA"

        path = join(dirname(__file__), "data/mut.csv")
        scenery = MutParser.parse(path, args)

        mean = sum(scenery.fitness) / len(scenery.fitness)
        self.assertAlmostEqual(mean, 0.4625, places=3)

    def test_calculate_rug(self):
        """test calculate ruggness"""
        chars = list("AT")

        scenery = Scenery()
        scenery.sequence = ["AAA", "AAT", "ATA",
                            "TAA", "ATT", "TAT", "TTA", "TTT"]
        scenery.fitness = [0.1, 0.2, 0.4, 0.3, 0.3, 0.6, 0.8, 1.0]

        meta = MetaData(scenery, chars)
        meta.get_neighbour()

        calculator = Ruggness(meta)
        rug = calculator.calculate()

        self.assertAlmostEqual(rug, 0.0095, places=3)

    def test_calculate_epi(self):
        """test calculate epistasis"""
        chars = list("AT")

        scenery = Scenery()
        scenery.sequence = ["AAA", "AAT", "ATA",
                            "TAA", "ATT", "TAT", "TTA", "TTT"]
        scenery.fitness = [0.1, 0.2, 0.4, 0.3, 0.3, 0.6, 0.8, 1.0]

        calculator = Epistasis(scenery, 3, chars)
        epi = calculator.calculate()

        percent_1 = epi[(1,)][('A',)] / epi[(0,)][('A',)]
        percent_2 = epi[(2,)][('A',)] / epi[(0,)][('A',)]

        self.assertAlmostEqual(percent_1, 0.7647, places=3)
        self.assertAlmostEqual(percent_2, 0.2941, places=3)
