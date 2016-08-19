import unittest
import pandas as pd
from .fanalyzer import FusionAnalysis

class FileLoadsCorrectly(unittest.TestCase):
    
    def setUp(self):
        
        import os
        
        lis_files = os.path.join(os.getcwd(), 'fusion', 'data', 'LIS')
        pcr_files = os.path.join(os.getcwd(), 'fusion', 'data', 'PCR Data')
        lis_list = [files for files in os.listdir(lis_files)]
        pcr_list = [files for files in os.listdir(pcr_files)]
        lis_file_read = os.path.join(lis_files, lis_list[0])
        pcr_file_read = os.path.join(pcr_files, pcr_list[0])
        self.dataframe = FusionAnalysis(pcr_file_read, lis_file_read, "P 1/2/3/4")
        
    def test_file_available(self):
        self.assertEqual(self.dataframe.pcr_file['Specimen Barcode'][0], '1011113418731122501243', "PCR File did not load properly")
        self.assertEqual(self.dataframe.lis_file['Specimen Barcode'][0], '1011113418731122501243', "LIS File did not load properly")
        
if __name__ == '__main__':
    unittest.main()