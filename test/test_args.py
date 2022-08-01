"""do the unit test of the argument client."""

import unittest
from os.path import join, dirname

from click.testing import CliRunner
from cliff.client import rug_mut, rug_seq, epi_mut, epi_seq


class TestArgCall(unittest.TestCase):
    """do the unit test of the argument client."""

    def test_rug_mut(self):
        """test calculate a ruggness on mutation format dataset"""
        path = join(dirname(__file__), "data/mut.csv")
        wile_type = "AAA"

        runner = CliRunner()
        result = runner.invoke(
            rug_mut, [path, '-w', wile_type, '-s', 'variant', '-f', 'score', '-c', 'AT'])

        self.assertEqual(result.exception, None)
        self.assertEqual(result.exit_code, 0)

    def test_rug_seq(self):
        """test calculate a ruggness on sequence format dataset"""
        path = join(dirname(__file__), "data/seq.csv")

        runner = CliRunner()
        result = runner.invoke(
            rug_seq, [path, '-s', 'Sequence', '-f', 'Fitness', '-c', 'ABCDEFGHIKL'])

        self.assertEqual(result.exception, None)
        self.assertEqual(result.exit_code, 0)

    def test_epi_mut(self):
        """test calculate a epistasis on mutation format dataset"""
        path = join(dirname(__file__), "data/mut.csv")
        wile_type = "AAA"

        runner = CliRunner()
        result = runner.invoke(
            epi_mut, [path, '-w', wile_type, '-s', 'variant', '-f', 'score', '-o', '1', '-c', 'AT'])

        self.assertEqual(result.exception, None)
        self.assertEqual(result.exit_code, 0)

    def test_epi_seq(self):
        """test calculate a epistasis on sequence format dataset"""
        path = join(dirname(__file__), "data/seq.csv")

        runner = CliRunner()
        result = runner.invoke(
            epi_seq, [path, '-s', 'Sequence', '-f', 'Fitness', '-o', '1', '-c', 'ABCDEFGHIKL'])

        self.assertEqual(result.exception, None)
        self.assertEqual(result.exit_code, 0)
