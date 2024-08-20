import unittest
from src.scripts.standardise_output_tools_wg import compare_sequences, check_and_update_if_haplotype_exists
import os
import pandas as pd

def make_haplotype_file():
    haplotype_file = 'test_file.tsv'

    f = open(haplotype_file, 'w')
    f.write("haplotype_id\tregion\trel_abundance\tsequence\n")
    f.write("1\tregion1\t0.2\tATCN\n")
    f.close()   
    return

def cleanup():
    haplotype_file = 'test_file.tsv'
    os.remove(haplotype_file)
    return

class TestSequenceComparison(unittest.TestCase):

    def test_compare_sequences_only_differs_in_one_base(self):
        sequence1 = 'ATCG'
        sequence2 = 'ATCA'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertFalse(bool_value)
        # assert that the corrected sequence is None
        self.assertIsNone(corrected_seq)
        

    def test_compare_sequences_same_shorter_sequence(self):
        sequence1 = 'ATCG'
        sequence2 = 'AT'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertTrue(bool_value)
        # assert that the corrected sequence is
        self.assertEqual(corrected_seq, 'ATCG')
    
    def test_compare_sequences_differs_longer_sequence(self):
        sequence1 = 'AC'
        sequence2 = 'ATCG'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertFalse(bool_value)
        # assert that the corrected sequence is
        self.assertIsNone(corrected_seq)
    
    def test_compare_sequences_same_sequence(self):
        sequence1 = 'ATCG'
        sequence2 = 'ATCG'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertTrue(bool_value)
        # assert that the corrected sequence is
        self.assertEqual(corrected_seq, 'ATCG')
    
    def test_compare_sequences_empty_sequence(self):
        sequence1 = ''
        sequence2 = 'ATCG'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertTrue(bool_value)
        # assert that the corrected sequence is
        self.assertEqual(corrected_seq, 'ATCG')
    
    def test_differs_only_by_Ns(self):
        sequence1 = 'ATCG'
        sequence2 = 'ATNG'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertTrue(bool_value)
        # assert that the corrected sequence is
        self.assertEqual(corrected_seq, 'ATCG')
    
    def test_differs_only_by_Ns_and_shorter_sequence(self):
        sequence1 = 'ATCG'
        sequence2 = 'ATN'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertTrue(bool_value)
        # assert that the corrected sequence is
        self.assertEqual(corrected_seq, 'ATCG')
    
    def test_differs_by_Ns_and_one_base(self):
        sequence1 = 'ATCG'
        sequence2 = 'ATNT'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertFalse(bool_value)
        # assert that the corrected sequence is None
        self.assertIsNone(corrected_seq)

    def test_differs_by_Ns_and_one_base_and_shorter_sequence(self):
        sequence1 = 'ATCG'
        sequence2 = 'ACN'
        bool_value, corrected_seq = compare_sequences(sequence1, sequence2)
        self.assertFalse(bool_value)
        # assert that the corrected sequence is None
        self.assertIsNone(corrected_seq)

    def test_same_seq_check_and_update_if_haplotype_exists_haplotype_exists(self):
        make_haplotype_file()
        haplotype_file = 'test_file.tsv'
        exists_bool = check_and_update_if_haplotype_exists(haplotype_file, 'region1', 'ATCG', 0.1)
        self.assertTrue(exists_bool)
        # open the file and check if the sequence exists
        f = pd.read_csv(haplotype_file, sep='\t', header=0)
        # get first row
        row = f.iloc[0]
        self.assertEqual(row['rel_abundance'], 0.3)
        self.assertEqual(row['sequence'], 'ATCG')
        cleanup()
    
    def test_check_and_update_if_haplotype_exists_haplotype_does_not_exist(self):
        make_haplotype_file()
        haplotype_file = 'test_file.tsv'
        exists_bool = check_and_update_if_haplotype_exists(haplotype_file, 'region1', 'AAC', 0.1)
        self.assertFalse(exists_bool)
        # open the file and check if the sequence exists
        f = pd.read_csv(haplotype_file, sep='\t', header=0)
        # get first row
        row = f.iloc[0]
        self.assertEqual(row['rel_abundance'], 0.2)
        cleanup()


if __name__ == '__main__':
    unittest.main()

        

   