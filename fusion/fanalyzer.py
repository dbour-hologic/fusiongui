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

	def __init__(self, pcr_path, lis_path, assay_type):

		self.assay_type = assay_type

		if (assay_type == 'P 1/2/3/4'):

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



	def check_assay_types(self):

		""" Checks if the assay types designated
		on the file is the same as the one that is 
		specified """

		if (self.pcr_file['Analyte'][0] != self.assay_type):
			return False
		return True


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
		CHANGE_PCR_COLUMN_DICT = {'RFU Range':'Unrounded RFU Range',
								  'LR_Ct_NonNormalized':'Unrounded Ct'}
		CHANGE_LIS_COLUMN_NAMES = True
		CHANGE_LIS_COLUMN_DICT = {'Interpretation 1':'FAM Rounded Ct',
                           		  'Interpretation 2':'HEX Rounded Ct',
		                          'Interpretation 3':'ROX Rounded Ct',
		                          'Interpretation 4':'RED647 Rounded Ct',
		                          'Interpretation 5':'IC Rounded Ct',
		                          'Interpretation 6':'POS/NEG/Invalid for HPIV-1',
		                          'Interpretation 7':'POS/NEG/Invalid for HPIV-2',
		                          'Interpretation 8':'POS/NEG/Invalid for HPIV-3',
		                          'Interpretation 9':'POS/NEG/Invalid for HPIV-4',
		                          'Interpretation 10':'Valid/Invalid for IC',
		                          'OtherData 1':'FAM Rounded RFU Range (HPIV-1)',
		                          'OtherData 2':'HEX Rounded RFU Range (HPIV-2)',
		                          'OtherData 3':'IC Rounded RFU Range',
		                          'OtherData 4':'RED647 Rounded RFU Range (HPIV-4)',
		                          'OtherData 5':'ROX Rounded RFU Range (HPIV-3)',
                          		 }
		PCR_COLUMNS_KEEP = ['Specimen Barcode', 'Analyte', 'Run ID', 
                            'Channel', 'Unrounded RFU Range', 'EstimatedBaseline',
                        	'Unrounded Ct', 'LR_TSlope_NonNormalized','Cartridge Lot #',
                         	'CapAndVialTrayID','Test order #', 'FCRBarcode', 'FERBarcode',
                         	'ElutionBufferRFID','ReconstitutionBufferRFID', 'OilRFID',
                         	'WellID','FusionTestOrder']

		LIS_COLUMNS_KEEP = ['Specimen Barcode','Analyte','Run ID','Instrument Flags',
		                      'FAM Rounded Ct', 'HEX Rounded Ct', 'ROX Rounded Ct', 'RED647 Rounded Ct', 
		                      'IC Rounded Ct', 'POS/NEG/Invalid for HPIV-1', 'POS/NEG/Invalid for HPIV-2',
		                      'POS/NEG/Invalid for HPIV-3', 'POS/NEG/Invalid for HPIV-4',
		                      'Valid/Invalid for IC','Overall_Validity', 'Serial Number', 'Sample Type', 'Sample Name',
		                      'Test order #', 'FAM Rounded RFU Range (HPIV-1)', 'HEX Rounded RFU Range (HPIV-2)',
		                      'IC Rounded RFU Range', 'RED647 Rounded RFU Range (HPIV-4)', 
		                      'ROX Rounded RFU Range (HPIV-3)'
		                    ]

		if CHANGE_PCR_COLUMN_NAMES:
			self.pcr_file.rename(columns=CHANGE_PCR_COLUMN_DICT, inplace=True)
		if CHANGE_LIS_COLUMN_NAMES:
			self.lis_file.rename(columns=CHANGE_LIS_COLUMN_DICT, inplace=True)

		# -----------> STARTING HERE IS PARAFLU SPECIFIC

		# --- START PCR MODIFICATIONS

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

		# --- END PCR MODIFICATIONS

		# --- START LIS MODIFICATIONS

		# Remove the '[end]' from 'Specimen Barcode', a designation for end of file that came from the automated Panther Software
		self.lis_file = self.lis_file[~self.lis_file['Specimen Barcode'].str.contains('[end]')]

		# Check for overall validity
		self.lis_file['Overall_Validity'] = "Valid"
		self.lis_file['Overall_Validity'][(
		                                (self.lis_file['POS/NEG/Invalid for HPIV-1']).str.contains('neg') &
		                                (self.lis_file['POS/NEG/Invalid for HPIV-2']).str.contains('neg') &
		                                (self.lis_file['POS/NEG/Invalid for HPIV-3']).str.contains('neg') &
		                                (self.lis_file['POS/NEG/Invalid for HPIV-4']).str.contains('neg') &
		                                (self.lis_file['Valid/Invalid for IC']).str.contains('Invalid')
		                             ) |
		                                (self.lis_file['POS/NEG/Invalid for HPIV-1']).str.contains('Invalid') |
		                                (self.lis_file['POS/NEG/Invalid for HPIV-2']).str.contains('Invalid') |
		                                (self.lis_file['POS/NEG/Invalid for HPIV-3']).str.contains('Invalid') |
		                                (self.lis_file['POS/NEG/Invalid for HPIV-4']).str.contains('Invalid')] = "Invalid"

		# Filter Columns
		lis_file_filtered_columns = self.lis_file[LIS_COLUMNS_KEEP]

		# Now let's give the table a unique id like we did for the PCR data so we can re-combine them!  
		# Add an extra column to help group the subjects. Will later to be used for combining with another dataframe.
		lis_file_filtered_columns = lis_file_filtered_columns.assign(UniqueID = lis_file_filtered_columns['Specimen Barcode'] + "_" + 
		                                                            			lis_file_filtered_columns['Run ID'] + "_" + 
		                                                             			lis_file_filtered_columns['Test order #'])
		lis_file_filtered_columns = lis_file_filtered_columns.set_index(['UniqueID'])

		# --- END LIS MODIFICATIONS
		# --- START COMBINING

		# Combine Files on Unique Key
		pcr_and_lis = lis_file_filtered_columns.join(pivot_pcr_file)

		# Logic calls to see if a positive hit was found, if not, mark the RFU Range channel with a "-".
		pcr_and_lis.loc[pcr_and_lis['POS/NEG/Invalid for HPIV-1'].str.contains("neg", case=False), 'FAM Rounded RFU Range (HPIV-1)'] = "-"
		pcr_and_lis.loc[pcr_and_lis['POS/NEG/Invalid for HPIV-2'].str.contains("neg", case=False), 'HEX Rounded RFU Range (HPIV-2)'] = "-"
		pcr_and_lis.loc[pcr_and_lis['POS/NEG/Invalid for HPIV-3'].str.contains("neg", case=False), 'ROX Rounded RFU Range (HPIV-3)'] = "-"
		pcr_and_lis.loc[pcr_and_lis['POS/NEG/Invalid for HPIV-4'].str.contains("neg", case=False), 'RED647 Rounded RFU Range (HPIV-4)'] = "-"
		pcr_and_lis.loc[pcr_and_lis['Valid/Invalid for IC'].str.contains("Invalid", case=False), 'IC Rounded RFU Range'] = "-"

		# Consolidate these columns since they more or less have the same information per channel
		pcr_and_lis.loc[:,'WellID'] = pcr_and_lis['FAM-WellID']
		pcr_and_lis.loc[:,'CapAndVialTrayID'] = pcr_and_lis['FAM-CapAndVialTrayID']
		pcr_and_lis.loc[:,'Cartridge Lot #'] = pcr_and_lis['FAM-Cartridge Lot #']
		pcr_and_lis.loc[:,'FCRBarcode'] = pcr_and_lis['FAM-FCRBarcode']
		pcr_and_lis.loc[:,'FERBarcode'] = pcr_and_lis['FAM-FERBarcode']
		pcr_and_lis.loc[:,'ElutionBufferRFID'] = pcr_and_lis['FAM-ElutionBufferRFID']
		pcr_and_lis.loc[:,'ReconstitutionBufferRFID'] = pcr_and_lis['FAM-ReconstitutionBufferRFID']
		pcr_and_lis.loc[:,'OilRFID'] = pcr_and_lis['FAM-OilRFID']
		pcr_and_lis.loc[:,'FusionTestOrder'] = pcr_and_lis['FAM-FusionTestOrder']

		# Final columns to keep
		try:
			save_as = pcr_and_lis[['Specimen Barcode', 
			             'Sample Type', 'Analyte', 
			             'Run ID', 
			             'FAM-WellID', 
			             'ROX-WellID',
			             'HEX-WellID',
			             'RED647-WellID', 
			             'IC-WellID',
			             'Test order #', 
			             'Instrument Flags',
			             'FAM Rounded Ct',
			             'FAM Rounded RFU Range (HPIV-1)',
			             'HEX Rounded Ct',
			             'HEX Rounded RFU Range (HPIV-2)',
			             'ROX Rounded Ct',
			             'ROX Rounded RFU Range (HPIV-3)',
			             'RED647 Rounded Ct',
			             'RED647 Rounded RFU Range (HPIV-4)',
			             'IC Rounded Ct',
			             'IC Rounded RFU Range',
			             'POS/NEG/Invalid for HPIV-1',
			             'POS/NEG/Invalid for HPIV-2',
			             'POS/NEG/Invalid for HPIV-3',
			             'POS/NEG/Invalid for HPIV-4',
			             'Valid/Invalid for IC',
			             'Overall_Validity',
			             'Serial Number',
			             'FAM-CapAndVialTrayID',
			             'ROX-CapAndVialTrayID',
			             'HEX-CapAndVialTrayID',
			             'RED647-CapAndVialTrayID',
			             'IC-CapAndVialTrayID',
			             'FAM-Cartridge Lot #',
			             'ROX-Cartridge Lot #',
			             'HEX-Cartridge Lot #',
			             'RED647-Cartridge Lot #',
			             'IC-Cartridge Lot #',
			             'FAM-FCRBarcode',
			             'ROX-FCRBarcode',
			             'HEX-FCRBarcode',
			             'RED647-FCRBarcode',
			             'IC-FCRBarcode',
			             'FAM-ElutionBufferRFID',
			             'ROX-ElutionBufferRFID',
			             'HEX-ElutionBufferRFID',
			             'RED647-ElutionBufferRFID',
			             'IC-ElutionBufferRFID',
			             'FAM-ReconstitutionBufferRFID',
			             'ROX-ReconstitutionBufferRFID',
			             'HEX-ReconstitutionBufferRFID',
			             'RED647-ReconstitutionBufferRFID',
			             'IC-ReconstitutionBufferRFID',
			             'FAM-OilRFID',
			             'ROX-OilRFID',
			             'HEX-OilRFID',
			             'RED647-OilRFID',
			             'IC-OilRFID',
			             'FAM-FusionTestOrder',
			             'ROX-FusionTestOrder',
			             'HEX-FusionTestOrder',
			             'RED647-FusionTestOrder',
			             'IC-FusionTestOrder',
			             'FAM-EstimatedBaseline',
			             'ROX-EstimatedBaseline',
			             'HEX-EstimatedBaseline',
			             'RED647-EstimatedBaseline',
			             'IC-EstimatedBaseline',
			             'FAM-LR_TSlope_NonNormalized',
			             'ROX-LR_TSlope_NonNormalized',
			             'HEX-LR_TSlope_NonNormalized',
			             'RED647-LR_TSlope_NonNormalized',
			             'IC-LR_TSlope_NonNormalized',
			             'FAM-Unrounded RFU Range',
			             'ROX-Unrounded RFU Range',
			             'HEX-Unrounded RFU Range',
			             'RED647-Unrounded RFU Range',
			             'IC-Unrounded RFU Range']]
		except KeyError:
			print("Error:", pcr_and_lis)

		# Save destination
		save_as.to_excel(save_to)

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
	pass