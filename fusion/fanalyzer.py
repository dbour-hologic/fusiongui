"""
Fanalyzer --
Using the PANDAS dataframe to manipulate,
consolidate, and perform data analysis on the
FUSION Paraflu and Flu assays
"""

import pandas as pd
import numpy as np

class FusionAnalysis():

	""" FusionAnalysis Class
	Purpose is to create a pandas dataframe from where
	data manipulation and data analysis of the LIS and 
	PCR files can occur.

	Features:
		1) Combine PCR & LIS files
		2) Statistical Analysis (future update)
		3) PQ (Performance Qualification) (future update)

	"""

	def __init__(self, lis_path, pcr_path, assay_type):

		self.assay_type = assay_type

		if (assay_type == 'Paraflu'):

			self.pcr_file = pd.read_csv(pcr_path,
										delimiter=',',
										encoding='utf-8-sig',
										dtype={'CapAndVialTrayID': object,
											   'ElutionBufferRFID': object,
											   'OilRFID': object,
											   'ReconstitutionBufferRFID': object,
											   'Test order #': object
											   }
									 	)

			self.lis_file = pd.read_csv(lis_path,
										delimiter='\t',
										encoding='utf-8-sig',
										dtype={'Test order #': object}
										)

	def combine_files(self, assay_profile, save_to):

		""" Combines the PCR & LIS Files 
		Args:
			assay_profile - the type of assay to manipulate (str)
			save_to - destination to save the file (str)
		Returns:
			None
		Output:
			Combined LIS & PCR file in *.csv format
		"""

		# Universal Settings - Can move to a JSON configuration later
		CHANGE_PCR_COLUMN_NAMES = True
		CHANGE_PCR_COLUMN_DICT = {'RFU Range':'Unrounded RFU Range','LR_Ct_NonNormalized':'Unrounded Ct'}
		CHANGE_LIS_COLUMN_NAMES = False
		CHANGE_LIS_COLUMN_DICT = {}
		PCR_COLUMNS_KEEP = ['Specimen Barcode', 'Analyte', 'Run ID', 
                            'Channel', 'Unrounded RFU Range', 'EstimatedBaseline',
                        	'Unrounded Ct', 'LR_TSlope_NonNormalized','Cartridge Lot #',
                         	'CapAndVialTrayID','Test order #', 'FCRBarcode', 'FERBarcode',
                         	'ElutionBufferRFID','ReconstitutionBufferRFID', 'OilRFID',
                         	'WellID','FusionTestOrder']

		if CHANGE_PCR_COLUMN_NAMES:
			self.pcr_file.rename(columns=CHANGE_PCR_COLUMN_DICT, inplace=True)
		if CHANGE_LIS_COLUMN_NAMES:
			self.lis_file.rename(columnns=CHANGE_LIS_COLUMN_DICT, inplace=True)

		# --- STARTING HERE IS PARAFLU SPECIFIC

		# Partition the barcode numbers
		self.pcr_file['CapAndVialTrayID'] = self.pcr_file['CapAndVialTrayID'].apply(lambda num: self.trimmer(num, trim_front=2, trim_back=10)) 
		self.pcr_file['FCRBarcode'] = self.pcr_file['FCRBarcode'].apply(lambda num: self.trimmer(num, trim_front=4, trim_back=11))
		self.pcr_file['FERBarcode'] = self.pcr_file['FERBarcode'].apply(lambda num: self.trimmer(num, trim_front=4, trim_back=11))
		self.pcr_file['ElutionBufferRFID'] = self.pcr_file['ElutionBufferRFID'].apply(lambda num: self.trimmer(num, trim_front=4, trim_back=11))
		self.pcr_file['ReconstitutionBufferRFID'] = self.pcr_file['ReconstitutionBufferRFID'].apply(lambda num: self.trimmer(num,trim_front=4, trim_back=11))
		self.pcr_file['OilRFID'] = self.pcr_file['OilRFID'].apply(lambda num: self.trimmer(num, trim_front=4, trim_back=11))

		pcr_file_filtered_columns = self.pcr_file[PCR_COLUMNS_KEEP]

		# Add an extra column to help group the subjects. Will later to be used for combining with another dataframe.
		pcr_file_filtered_columns = pcr_file_filtered_columns.assign(UniqueID = pcr_file_filtered_columns['Specimen Barcode'] + "_" + 
                                                             					pcr_file_filtered_columns['Run ID'] + "_" + 
                                                             				    pcr_file_filtered_columns['Test order #'])

		# Remove the '[end]' from 'Specimen Barcode', a designation for end of file that came from the automated Panther Software
		pcr_file_filtered_columns = pcr_file_filtered_columns[~pcr_file_filtered_columns['Specimen Barcode'].str.contains('[end]')]

		# Pivot the PCR File - the long to wide conversion
		pivot_pcr_file = pcr_file_filtered_columns.pivot(index='UniqueID', columns='Channel')

		# Merge columns from multilayer to single layer headers
		pivot_pcr_file.columns = [name[1]+"-"+name[0] for name in pivot_pcr_file.columns.values]

		# --- END PARAFLU SPECIFIC

		print(pivot_pcr_file[:2])

	def trimmer(self, number, trim_front=0, trim_back=0):

		"""
		Trimmer takes in an ID number, (obj) representation
		and parses out the selected segment depending on the
		parameters.

		Args:
		    number - the number to parse (obj)
		    trim_front - how many characters to remove from the beginning of number (int)
		    trim_back - how many characters to remove starting from end to beginning of number (int)
		Returns:
		     the result of the parsed number string (str)
		"""
		# Convert to an actual object representation
		object_to_string = str(number)

		if trim_front > 0 and trim_back > 0:
		    # Convert to trim backwards
		    trim_back *= -1
		    after_trim_beginning = object_to_string[trim_front:]
		    after_trim_ending = after_trim_beginning[:trim_back]
		    return after_trim_ending
		elif trim_front > 0:
		    return object_to_string[trim_front:]
		elif trim_back > 0:
		    trim_back *= -1
		    return object_to_string[:trim_back]

		return object_to_string

if __name__ == '__main__':
	f = FusionAnalysis("@Pt22090000574-20160701-183343-000574-20160701-04.lis", "@DI2090000574-20160701-215050-000574-20160701-04.csv", "Paraflu")
	f.combine_files('Paraflu', 'here')