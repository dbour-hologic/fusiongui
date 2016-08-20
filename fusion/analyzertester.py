import unittest
import pandas as pd
from .fanalyzer import FusionAnalysis
from .fcombiner import FusionCombiner

class FileCombinerTest(unittest.TestCase):
    
    def setUp(self):
        
        import os
        
        lis_files = os.path.join(os.getcwd(), 'fusion', 'data', 'LIS')
        pcr_files = os.path.join(os.getcwd(), 'fusion', 'data', 'PCR Data')
        lis_list = [files for files in os.listdir(lis_files)]
        pcr_list = [files for files in os.listdir(pcr_files)]
        lis_file_read = os.path.join(lis_files, lis_list[0])
        pcr_file_read = os.path.join(pcr_files, pcr_list[0])
        self.dataframe = FusionAnalysis(pcr_file_read, lis_file_read, "P 1/2/3/4")

    def test_valid_datatype(self):
        self.assertTrue(self.dataframe.valid_data)

    def test_file_available(self):
        self.assertEqual(self.dataframe.pcr_file['Specimen Barcode'][0], '1011113418731122501243', "PCR File did not load properly")
        self.assertEqual(self.dataframe.lis_file['Specimen Barcode'][0], '1011113418731122501243', "LIS File did not load properly")

    def test_column_name_changes(self):
        self.assertRaises(KeyError, lambda: self.dataframe.pcr_file['RFU Range'])
        self.assertRaises(KeyError, lambda: self.dataframe.lis_file['Interpretation 2'])
        self.assertTrue('Unrounded RFU Range' in self.dataframe.pcr_file)
        self.assertTrue('IC Rounded Ct' in self.dataframe.lis_file)

    def test_trimmer(self):

        id_map = {
            "CapAndVialTrayID":"444383",
            "ElutionBufferRFID":"165629",
            "ReconstitutionBufferRFID":"144364"
        }

        for columns, values in id_map.items():
            self.assertEqual(self.dataframe.pcr_file[columns][0], values, "ID does not match")

    def test_combination_file(self):
        combination = self.dataframe.combine_files("P 1/2/3/4")
        self.assertTrue('WellID' in combination)
        self.assertTrue('Specimen Barcode' in combination)

    def test_combination_file_filter(self):
        combination = self.dataframe.combine_files("P 1/2/3/4", name_space=['Specimen Barcode'])
        self.assertTrue('Specimen Barcode' in combination)
        self.assertRaises(KeyError, lambda: combination['WellID'])

class FusionCombinerTest(unittest.TestCase):

    def setUp(self):

        import os

        lis_files = os.path.join(os.getcwd(), 'fusion', 'data', 'LIS')
        pcr_files = os.path.join(os.getcwd(), 'fusion', 'data', 'PCR Data')
        lis_list = [files for files in os.listdir(lis_files)]
        pcr_list = [files for files in os.listdir(pcr_files)]

        giant_list = []

        for x in lis_list:
            giant_list.append(os.path.join(lis_files, x))
        for y in pcr_list:
            giant_list.append(os.path.join(pcr_files,y))

        self.combined = FusionCombiner(giant_list, "P 1/2/3/4")

    def test_has_pairs(self):

        print(self.combined.all_items)
        self.assertTrue(len(self.combined.all_items['paired']) > 0)

if __name__ == '__main__':
    unittest.main()