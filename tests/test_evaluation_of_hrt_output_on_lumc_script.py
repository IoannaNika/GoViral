import unittest
import os
import pandas as pd
from src.scripts.evaluate_hrt_output_lumc import closest_haplotype, calculate_recall, calculate_duplication_ratio, calculate_relative_absolute_abundance_error, calculate_average_number_of_haplotypes, calculate_average_edit_distance


class TestEvaluationOfHRTOutputOnLUMCScript(unittest.TestCase):

    def test_closest_haplotype(self):
        seq = "ACGT"
        omicron_hap = "ATTC"
        wuhan_hap = "ACGT"
        sample_name = "test"
        region = "test_region"
        true_abs_per_sample_per_region  = {"test": {"test_region":{"Wuhan": 0.25, "Omicron": 0.75}}}
        ab_per_region = {"test_region": {"Wuhan": 0.0, "Omicron": 0.75}}
        rel_ab = 0.25

        closest_hap, distance = closest_haplotype(seq, rel_ab, omicron_hap, wuhan_hap, true_abs_per_sample_per_region, region, ab_per_region, sample_name)

        self.assertEqual(closest_hap, "Wuhan")
        self.assertEqual(distance, 0)

    def test_closest_haplotype_shared_region(self):
        seq = "ACGT"
        omicron_hap = "ACGT"
        wuhan_hap = "ACGT"
        sample_name = "test"
        region = "test_region"
        true_abs_per_sample_per_region  = {"test": {"test_region":{"Wuhan": 1, "Omicron": 0}}}
        ab_per_region = {"test_region": {"Wuhan": 0.5, "Omicron": 0}}
        rel_ab = 0.25

        closest_hap, distance = closest_haplotype(seq, rel_ab, omicron_hap, wuhan_hap, true_abs_per_sample_per_region, region, ab_per_region, sample_name)

        self.assertEqual(closest_hap, "Wuhan")
        self.assertEqual(distance, 0)

    
    def test_calculate_recall_and_duplication_ratio(self): 
        num_haps_per_region = {"test_region":{"Wuhan": 1, "Omicron": 1}}
        sample_name = "test"
        true_num_of_haps_per_sample_per_region = {"test": {"test_region": {"Wuhan": 1, "Omicron": 0}}}
        recall_wuhan, recall_omicron = calculate_recall(num_haps_per_region, true_num_of_haps_per_sample_per_region, sample_name)

        self.assertEqual(recall_wuhan, 1)
        self.assertEqual(recall_omicron, 1)

        dup_ratio = calculate_duplication_ratio(num_haps_per_region, true_num_of_haps_per_sample_per_region, sample_name)
        self.assertEqual(dup_ratio, 2)

    def test_calculate_absolute_relative_abundance_error(self): 
        true_abs_per_sample_per_region  = {"test": {"test_region":{"Wuhan": 0.25, "Omicron": 0.75}}}
        ab_per_region = {"test_region": {"Wuhan": 0.0, "Omicron": 0.75}}
        sample_name = "test"

        abs_error = calculate_relative_absolute_abundance_error(ab_per_region, true_abs_per_sample_per_region, sample_name)

        self.assertEqual(abs_error, 0.5)

    
    def test_calculate_average_number_of_haplotypes(self):
        num_haps_per_region = {"test_region_1":{"Wuhan": 1, "Omicron": 1}, "test_region_2":{"Wuhan": 2, "Omicron": 0}, "test_region_3":{"Wuhan": 0, "Omicron": 1}}
        avg_num_of_haps = calculate_average_number_of_haplotypes(num_haps_per_region)

        self.assertEqual(avg_num_of_haps, 1.67)

    def test_calculate_average_edit_distance(self):
        edit_distances = {"test_region_1":{ "Wuhan": 1, 'Omicron': 0}, "test_region_2":{ "Wuhan": 2, 'Omicron': 1}}
   
        avg_edit_distance = calculate_average_edit_distance(edit_distances)

        self.assertEqual(avg_edit_distance, 1)

if __name__ == '__main__':
    unittest.main()

    



