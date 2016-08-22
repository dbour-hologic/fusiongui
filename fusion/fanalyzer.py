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
	"""

	def __init__(self, pcr_path, lis_path, assay_type, *args, **kwargs):

		# Holds the columns to rename
		self.pcr_columns_rename = {}
		self.lis_columns_rename = {}

		# Holds the trim parameters
		self.trim_map = {}

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
										dtype={'Test order #': object,
											'Serial Number':object
											}
										)

		self.__set_variables(assay_type)
		self.__change_column_names(assay_type, self.lis_file, self.lis_columns_rename)
		self.__change_column_names(assay_type, self.pcr_file, self.pcr_columns_rename)
		self.__trim_column_values(assay_type, self.pcr_file, self.trim_map)

		self.valid_data = self.__is_valid_dataset(assay_type)


	def __is_valid_dataset(self, assay_type):
		""" (PRIVATE) Checks if the pair files are of a matching 
		assay type or empty.

		Args:
			assay_type - assay_type of data (str)
		Returns:
			true if pass, false if doesn't pass (bool)
		"""

		try:
			if self.lis_file['Analyte'][0] == assay_type and self.pcr_file['Analyte'][0] == assay_type:
				return True
		except KeyError:
			return False
		return False

	def __set_variables(self, assay_type):
		""" (PRIVATE) Method is currently used as a substitution for
		a JSON file to set various variable settings for the different assay
		types. 

		Args:
			assay_type - assay type to run <flu, paraflu, adeno, etc> (str)
		"""
		if assay_type == "P 1/2/3/4":
			self.pcr_columns_rename =  {'RFU Range':'Unrounded RFU Range','LR_Ct_NonNormalized':'Unrounded Ct'}
			self.lis_columns_rename =  {'Interpretation 1':'FAM Rounded Ct',
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
							'OtherData 5':'ROX Rounded RFU Range (HPIV-3)'}	

			self.trim_map = {
				"CapAndVialTrayID":(2,10),
				"FCRBarcode":(4,11),
				"FERBarcode":(4,11),
				"ElutionBufferRFID":(4,11),
				"ReconstitutionBufferRFID":(4,11),
				"OilRFID":(4,11)
			}

			self.mapping_source = {
				"2090000330":"F39",
				"2090000186":"F39",
				"2090000323":"F41",
				"2090000327":"F42",
				"2090000775":"F45",
				"2090000735":"F51",
				"2090001197":"F62",
				"2090000574":"F64"
			}

			self.software_source = {
				"0.91.4":"GAP 6.0.6.15",
				"0.93.4":"GAP 6.0.6.17",
				"0.95.4":"GAP 6.9.6.17",
			}

	def __change_column_names(self, assay_type, dframe, name_map):
		""" (PRIVATE) Change column name is an internal method used
		for manipulating the PCR FILE column names or the LIS FILE 
		column names. Ideally all of the the settings for specific assays will
		be loaded through JSON, but are hard coded for the time being

		Assay Types:
		(1) Paraflu =  ' P 1/2/3/4'
		(2) Flu = (future)
		(3) Adeno = (future)

		Args:
			assay_type - assay type to run <flu, paraflu, adeno, etc> (str)
			dframe - the dataframe to manipulate <LIS or PCR> (obj)
			name_map - mapping of column and their name changes (dict)
		Returns:
			None

		"""

		dframe.rename(columns=name_map, inplace=True)

	def __trim_column_values(self, assay_type, dframe, name_map):
		""" (PRIVATE) Trim column values according to a user
	 	specified settings file. Ideally all of the settings for specific
	 	assays will be loaded through JSON, but are hard coded for now.

	 	Main purpose: Some 'barcodes' have extra numbers appended
	 	as a suffix or prefix. This will help extract the actual barcode number
	 	from the string of numbers

	 	Args:
			assay_type - assay type to run <flu, paraflu, adeno, etc> (str)
			dframe - the dataframe to manipulate <LIS or PCR> (obj)
			name_map - mapping of column and coordinates (dict)
		Returns:
			None
		"""

		if assay_type == "P 1/2/3/4":
			for column, coord in name_map.items():
				dframe[column] = dframe[column].apply(lambda num: self.__trimmer(num, trim_front=coord[0], trim_back=coord[1]))

	def __trimmer(self, number, trim_front=0, trim_back=0):

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

	def combine_files(self, assay_type, instrument_mapping=None, software_mapping=None, *args, **kwargs):

		""" Combines the PCR & LIS Files 
		Args:
			assay_type - the type of assay to manipulate (str)
			save_to - destination to save the file (str)
		Optional:
			*args: None
			**kwargs:
				name_space - the columns you want to keep (list)
		Returns:
			the combined file dataframe object (obj)
		"""

		name_space = kwargs.get("name_space", None)


		if assay_type == "P 1/2/3/4":


			pcr_unique_id = self.pcr_file['Specimen Barcode'] + '_' + self.pcr_file['Run ID'] + "_" + self.pcr_file['Test order #']
			self.pcr_file = self.pcr_file.assign(UniqueID = pcr_unique_id)
			# Remove end of the file desgination set by Panther (Messes up 'pivot')
			# The '~' is a special PANDAS feature, it acts as a not equal ( i.e. != )
			self.pcr_file = self.pcr_file[~self.pcr_file['Specimen Barcode'].str.contains('[end]', regex=False)]
			self.pivot_pcr_file = self.pcr_file.pivot(index='UniqueID', columns='Channel')
			# This essentially appends "Channel name" + "Column name"
			self.pivot_pcr_file.columns = [name[1] + '-' + name[0] for name in self.pivot_pcr_file.columns.values]

			lis_unique_id = self.lis_file['Specimen Barcode'] + '_' + self.lis_file['Run ID'] + "_" + self.lis_file['Test order #']
			self.lis_file = self.lis_file[~self.lis_file['Specimen Barcode'].str.contains('[end]',  regex=False)]
			self.lis_file = self.lis_file.assign(UniqueID = lis_unique_id)
			self.lis_file = self.lis_file.set_index(['UniqueID'])

			# Magic pivot ! Like EXCEL
			pcr_and_lis = self.lis_file.join(self.pivot_pcr_file)
		
			# Reduces the amount of columns; we use only one channel to make a general channel name
			# Since they all effectly have the values across each channel. It should never be different! hah...
			self.consolidation_map = {
				"WellID": "FAM-WellID",
				"CapAndVialTrayID":"FAM-CapAndVialTrayID",
				"Cartridge Lot #":"FAM-Cartridge Lot #",
				"FCRBarcode":"FAM-FCRBarcode",
				"FERBarcode":"FAM-FERBarcode",
				"ElutionBufferRFID":"FAM-ElutionBufferRFID",
				"ReconstitutionBufferRFID":"FAM-ReconstitutionBufferRFID",
				"OilRFID":"FAM-OilRFID",
				"FusionTestOrder":"FAM-FusionTestOrder",
			}

			for new_col, old_col in self.consolidation_map.items():
				pcr_and_lis.loc[:,new_col] = pcr_and_lis[old_col]


			if instrument_mapping:
				pcr_and_lis = self.__apply_instrument_mapping(instrument_mapping, pcr_and_lis)
			if software_mapping:
				pcr_and_lis = self.__apply_software_mapping(software_mapping, pcr_and_lis)


			if name_space:
				return pcr_and_lis[name_space]
			return pcr_and_lis

	def __apply_instrument_mapping(self, mapping, dframe):
		""" Give a mapping for instruments """

		try:
			dframe['Serial Number'] = dframe['Serial Number'].apply(lambda element: mapping[element])
		except KeyError:
			print("No serial number mapping available.")
		
		return dframe

	def __apply_software_mapping(self, mapping, dframe):
		try:
			dframe['Software Revision'] = dframe['Software Revision'].apply(lambda element: mapping[element])
		except KeyError as e:
			print("No software version matching!", e)

		return dframe

	def save_combined(self, df_combined, save_to):
		"""
		Saves dataframe as an Excel file. This is done on purpose
		to prevent the long string of numbers from being automatically
		converted in Excel when opening it as a .csv/.tsv

		Args:
			df_combined - a combined dataframe
			save_to - path to save combined file
		Returns:
			None
		Output:
			Combined LIS & PCR file in *.xlsx format
		"""

		df_combined.to_excel(save_to)

	def _write_out(self, fields_to_output, save_dir, file_id, is_valid):

		import os

		validity_tag = "PASS"

		if not is_valid:
			validity_tag = "FAIL"
		
		file_prex = "%s_PQ-" % (validity_tag)
		file_name = file_prex + file_id + '-pq_results.tsv'

		file_directory = os.path.join(save_dir, file_name)

		try:
			with open(file_directory, 'w') as writeout:

				if len(fields_to_output) > 0:

					for header, data in fields_to_output.items():

						if header == "is_failed_positive":

							writeout.write("FAIL REASON:\t" + "DID NOT MEET PQ THRESHOLD\n\n")
							writeout.write("SPECIMEN BARCODE"+"\t"+"RUN ID"+"\t"+"TEST ORDER #\n")
							for channel_type, sample_list in data.items():
								if len(sample_list) > 0:
									writeout.write("CHANNEL:\t" + channel_type + "\n")
									for samples in sample_list:
										writeout.write('\t'.join(samples) + "\n")
						else:
							writeout.write("FAIL REASON:\t" + header + '\n')
							writeout.write("SPECIMEN BARCODE"+"\t"+"RUN ID"+"\t"+"TEST ORDER #\n")
							for row in data:
								join_lines = '\t'.join(row)
								writeout.write(join_lines + "\n")
				else:

					writeout.write("PQ passed.")
		except IOError:
			print("Error writing PQ analysis output file.")


if __name__ == '__main__':

	import os

	lis_files = os.path.join(os.getcwd(), 'data', 'LIS')
	pcr_files = os.path.join(os.getcwd(), 'data', 'PCR Data')

	lis_list = [files for files in os.listdir(lis_files)]
	pcr_list = [files for files in os.listdir(pcr_files)]

	lis_file_read = os.path.join(lis_files, lis_list[0])
	pcr_file_read = os.path.join(pcr_files, pcr_list[0])

	lis_file_read2 = os.path.join(lis_files, lis_list[1])
	pcr_file_read2 = os.path.join(pcr_files, pcr_list[1])

	p = FusionAnalysis(pcr_file_read, lis_file_read, "P 1/2/3/4")
	pp = p.combine_files('P 1/2/3/4', os.path.join(os.getcwd(),'sampletest.xlsx'))
	p2 = FusionAnalysis(pcr_file_read2, lis_file_read2, "P 1/2/3/4")
	pp2 = p2.combine_files('P 1/2/3/4', os.path.join(os.getcwd(),'sampletest2.xlsx'))
	p.pq_analysis(p.combine_files('P 1/2/3/4', os.path.join(os.getcwd(),'samplepq.xlsx')), os.getcwd(), 'lol')

	ff = FusionCombiner(pp, pp2)
	ff.aggr()
	

	# p.pq_analysis(p.combine_files('P 1/2/3/4', os.path.join(os.getcwd(),'sampletest.xlsx')), os.getcwd(), '@DI209000001-20160811-165658-000001-20160810-01')