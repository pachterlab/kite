import kite
import unittest
from collections import OrderedDict



class TestMakeMismatchMap(unittest.TestCase):
    def test_even_tag_length_error(self):
        len8_HTO_tags_return = kite.make_mismatch_map('/home/munfred/kite/tests/len8_HTO_tags.fasta')
        self.assertEqual(len8_HTO_tags_return, None)

    def test_odd_tag_length(self):
        len3_microtgs_return = kite.make_mismatch_map('/home/munfred/kite/tests/len3_microtags.fasta')

        expected_microtags_return = OrderedDict(
            [('microtagA', 'GAC'), ('microtagA+-+0-1', 'TAC'), ('microtagA+-+0-2', 'AAC'), ('microtagA+-+0-3', 'CAC'),
             ('microtagA+-+1-1', 'GTC'), ('microtagA+-+1-2', 'GGC'), ('microtagA+-+1-3', 'GCC'),
             ('microtagA+-+2-1', 'GAT'), ('microtagA+-+2-2', 'GAG'), ('microtagA+-+2-3', 'GAA'), ('microtagB', 'ATG'),
             ('microtagB+-+0-1', 'TTG'), ('microtagB+-+0-2', 'GTG'), ('microtagB+-+0-3', 'CTG'),
             ('microtagB+-+1-1', 'AAG'), ('microtagB+-+1-2', 'AGG'), ('microtagB+-+1-3', 'ACG'),
             ('microtagB+-+2-1', 'ATT'), ('microtagB+-+2-2', 'ATA'), ('microtagB+-+2-3', 'ATC'), ('microtagC', 'CTT'),
             ('microtagC+-+0-1', 'TTT'), ('microtagC+-+0-2', 'GTT'), ('microtagC+-+0-3', 'ATT'),
             ('microtagC+-+1-1', 'CAT'), ('microtagC+-+1-2', 'CGT'), ('microtagC+-+1-3', 'CCT'),
             ('microtagC+-+2-1', 'CTA'), ('microtagC+-+2-2', 'CTG'), ('microtagC+-+2-3', 'CTC')])

        self.assertEqual(len3_microtgs_return, expected_microtags_return)

class TestStringMethods(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

if __name__ == '__main__':
    unittest.main()