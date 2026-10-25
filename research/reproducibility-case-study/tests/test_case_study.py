"""Contract tests; authored expectations do not invoke the native implementation."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import case_common as c


class RequestTests(unittest.TestCase):
    def test_fixed_cardinalities_changed_alt_and_mask_only_control(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contigs = [{'name':'chr20', 'length':64000000}]
            rows = []
            for i in range(1000):
                ref = 'ACGT'[i%4]
                rows.append(dict(contig='chr20',pos=10000001+i*1000,ref=ref,alt=next(b for b in 'ACGT' if b!=ref)))
            c.write_sites(root/'probes.vcf', rows, contigs)
            result = c.freeze_queries(root/'probes.vcf', root/'requests', contigs)
            self.assertEqual([len(r['rows']) for r in result], [250,500,1000])
            self.assertEqual(result[0]['rows'], rows[:250])
            self.assertEqual([r['pos'] for r in result[1]['rows']], [r['pos'] for r in rows[::2]])
            self.assertTrue(all(r['alt'] != old['alt'] and r['alt'] != r['ref'] for r,old in zip(result[1]['rows'],rows[::2])))
            self.assertEqual((root/'requests/absent.bed').read_text(),'chr20\t10000001\t10000002\n')
            wrong = c.probes(root/'requests/wrong-ref.vcf')[0]
            self.assertNotEqual(wrong['ref'], rows[0]['ref'])
            self.assertNotEqual(wrong['alt'], wrong['ref'])
            c.verify_inventory(root/'requests', c.read_json(root/'requests/freeze.json')['files'])
            with self.assertRaises(FileExistsError):
                c.freeze_queries(root/'probes.vcf',root/'requests',contigs)
            self.assertEqual(len(c.probes(root/'requests/Q1.vcf')),250)

    def test_wrong_scope_duplicate_or_unsorted_refused(self):
        rows = [dict(contig='chr20',pos=10000001+i*1000,ref='A',alt='C') for i in range(1000)]
        contigs=[dict(name='chr20',length=64000000)]
        for changed in (rows[:-1], list(reversed(rows)), [rows[0]]+rows[:-1], [dict(rows[0],pos=1)]+rows[1:]):
            with self.subTest(kind=changed[0]), self.assertRaises(ValueError):
                c.request_definitions(changed, contigs)


class ReducerTests(unittest.TestCase):
    def write_fixture(self, root):
        header=['#contig','pos','ref','requested_alts','callable_depth','a','c','g','t',
                'a_fwd','a_rev','c_fwd','c_rev','g_fwd','g_rev','t_fwd','t_rev']
        queries=[]
        rows=[]
        for i,depth in enumerate([0,4,5,9,10,19,20,2**64-1],1):
            alt='C' if i%2 else 'G'
            query=dict(contig='chr20',pos=i,ref='A',alt=alt)
            queries.append(query)
            row=dict.fromkeys(header,'0')
            row.update({'#contig':'chr20','pos':str(i),'ref':'A','requested_alts':alt,'callable_depth':str(depth),
                        alt.lower():str(depth), alt.lower()+'_fwd':str(depth)})
            rows.append(row)
        with (root/'native.tsv').open('w') as output:
            writer=csv.DictWriter(output,fieldnames=header,delimiter='\t',lineterminator='\n')
            writer.writeheader();writer.writerows(rows)
        return dict(id='Q3',rows=queries,thresholds=[5,10,20]),rows,header

    def test_hand_authored_thresholds_zero_and_uint64(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            request,_,_=self.write_fixture(root)
            actual=c.reduce_native(root/'native.tsv',request)
            oracle=c.oracle_answers(root/'native.tsv',request)
            self.assertEqual(actual,oracle)
            self.assertEqual([r['depth_ge_5'] for r in actual],[False,False,True,True,True,True,True,True])
            self.assertEqual([r['depth_ge_10'] for r in actual],[False,False,False,False,True,True,True,True])
            self.assertEqual([r['depth_ge_20'] for r in actual],[False,False,False,False,False,False,True,True])
            self.assertEqual(actual[0]['alt_count'],0)
            self.assertEqual(actual[-1]['alt_count'],2**64-1)

    def test_ref_alt_duplicate_missing_and_integer_overflow_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            request,rows,header=self.write_fixture(root)
            variants=[rows[:-1],rows+[rows[-1]], [dict(rows[0],ref='T')]+rows[1:],
                      [dict(rows[0],requested_alts='G')]+rows[1:],
                      [dict(rows[0],callable_depth=str(2**64))]+rows[1:]]
            for bad in variants:
                with (root/'bad.tsv').open('w') as output:
                    writer=csv.DictWriter(output,fieldnames=header,delimiter='\t',lineterminator='\n')
                    writer.writeheader();writer.writerows(bad)
                with self.assertRaises(ValueError):
                    c.reduce_native(root/'bad.tsv',request)


class FailureRetentionTests(unittest.TestCase):
    def test_unexpected_subprocess_failure_retained_and_nonempty_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)/'failure'
            with self.assertRaises(RuntimeError):
                c.command([sys.executable,'-c','import sys; print("retained"); sys.exit(7)'],root)
            self.assertEqual(c.read_json(root/'command.json')['exit_code'],7)
            self.assertEqual((root/'stdout.log').read_text(),'retained\n')
            before=c.inventory(root)
            with self.assertRaises(FileExistsError):
                c.command([sys.executable,'-c','pass'],root)
            self.assertEqual(c.inventory(root),before)

    def test_changed_protocol_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'protocol.json'
            value=c.load_protocol(HERE/'protocol-v1.json')
            value['native']['workers']=2
            c.write_json(path,value)
            with self.assertRaises(ValueError):
                c.load_protocol(path)

    def test_inventory_relative_paths_and_missing_required_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'payload').write_text('evidence')
            before=c.inventory(root)
            c.verify_inventory(root,before)
            (root/'payload').unlink()
            with self.assertRaises(ValueError):
                c.verify_inventory(root,before)
            for name in ['../escape','/absolute','a/../b','a//b','a\\b']:
                with self.subTest(name=name),self.assertRaises(ValueError):
                    c.safe_path(root,name)


if __name__=='__main__':
    unittest.main()
