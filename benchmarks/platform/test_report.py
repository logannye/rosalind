import io
import tempfile
import unittest
from pathlib import Path

from bcftools_features import write_rows
from report import compare


class PlatformEvidenceTests(unittest.TestCase):
    def fixture(self, left, right):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        a, b = Path(temporary.name) / 'a.tsv', Path(temporary.name) / 'b.tsv'
        a.write_text(left)
        b.write_text(right)
        return a, b

    def test_merge_join_counts_missing_and_different_rows_in_reference_order(self):
        a, b = self.fixture('chr2\t1\tA\nchr2\t3\tT\nchr1\t1\tA\n', 'chr2\t2\tA\nchr2\t3\tC\nchr1\t1\tA\n')
        result = compare(a, b, [2], {'chr2': 0, 'chr1': 1})
        self.assertEqual(result['left_rows'], 3)
        self.assertEqual(result['right_rows'], 3)
        self.assertEqual(result['left_only'], 1)
        self.assertEqual(result['right_only'], 1)
        self.assertEqual(result['common_rows'], 2)
        self.assertEqual(result['different_common_rows'], 1)

    def test_duplicate_or_unsorted_rows_are_not_silently_overwritten(self):
        for contents in ('chr1\t1\tA\nchr1\t1\tT\n', 'chr1\t2\tA\nchr1\t1\tA\n'):
            a, b = self.fixture(contents, 'chr1\t1\tA\n')
            with self.assertRaises(ValueError):
                compare(a, b, [2], {'chr1': 0})

    def test_bcftools_normalization_consumes_only_one_row_at_a_time(self):
        output = io.StringIO()
        def source():
            for position in range(1, 100):
                if position > 1:
                    self.assertIn(f'chr1\t{position - 1}\t', output.getvalue())
                yield f'chr1\t{position}\tA\t5\t5,0\n'
        write_rows(source(), output)
        self.assertEqual(len(output.getvalue().splitlines()), 100)

    def test_empty_stream_and_numeric_position_order(self):
        a, b = self.fixture('', 'chr1\t2\tA\nchr1\t10\tC\n')
        result = compare(a, b, [2], {'chr1': 0})
        self.assertEqual(result['right_only'], 2)
        self.assertEqual(result['left_rows'], 0)


if __name__ == '__main__':
    unittest.main()
