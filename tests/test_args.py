import unittest
from os.path import join, dirname, abspath

from click.testing import CliRunner
from cliff.cli import rug_mut, rug_seq, epi_mut, epi_seq


class TestArgCall(unittest.TestCase):
    def test_rug_mut(self):
        path = join(dirname(__file__), "data/mut.csv")
        wile_type = "AAA"

        runner = CliRunner()
        result = runner.invoke(
            rug_mut, [path, '-w', wile_type, '-s', 'variant', '-f', 'score', '-c', 'AT'])
        self.maxDiff = None

        self.assertEqual(result.exception, None)
        self.assertEqual(result.exit_code, 0)

    def test_rug_seq(self):
        path = join(dirname(__file__), "data/seq.csv")

        runner = CliRunner()
        result = runner.invoke(
            rug_seq, [path, '-s', 'Sequence', '-f', 'Fitness', '-c', 'ABCDEFGHIKL'])

        self.assertEqual(result.exception, None)
        self.assertEqual(result.exit_code, 0)

    def test_epi_mut(self):
        path = join(dirname(__file__), "data/mut.csv")
        wile_type = "AAA"

        runner = CliRunner()
        result = runner.invoke(
            epi_mut, [path, '-w', wile_type, '-s', 'variant', '-f', 'score', '-o', '1', '-c', 'AT'])

        self.assertEqual(result.exception, None)
        self.assertEqual(result.exit_code, 0)

    def test_epi_seq(self):
        path = join(dirname(__file__), "data/seq.csv")

        runner = CliRunner()
        result = runner.invoke(
            epi_seq, [path, '-s', 'Sequence', '-f', 'Fitness', '-o', '1', '-c', 'ABCDEFGHIKL'])

        self.assertEqual(result.exception, None)
        self.assertEqual(result.exit_code, 0)
