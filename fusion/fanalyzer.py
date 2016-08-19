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
										dtype={'Test order #': object}
										)

		self.__set_variables(assay_type)
		self.__change_column_names(assay_type, self.lis_file, self.lis_columns_rename)
		self.__change_column_names(assay_type, self.pcr_file, self.pcr_columns_rename)
		self.__trim_column_values(assay_type, self.pcr_file, self.trim_map)

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

	def combine_files(self, assay_type, *args, **kwargs):

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

			if name_space:
				return pcr_and_lis[name_space]
			return pcr_and_lis

	def save_combed(self, df_combined, save_to):
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
			Combined LIS & PCR file in *.csv format
		"""

		df_combined.to_excel(save_to)


	def pq_analysis(self, combined_file, save_to, file_id):

		""" Performs PQ Analysis on the combined file
		Args:
			combined_file - the combined file from 'combine_files' (obj)
			save_to - the directory to save the summary results to (str)
			file_id - the unique file id (str)
		Returns:

			(primary) if run passes every PQ criteria, return boolean True, else false
			(secondary) file with results of the PQ analysis 

		# LAB CONTROL SETTINGS
		# These positive control usually begin with these prefix'd numbers
		# Flu - 101101..####
		# AMR - 102101..####
		# Paraflu - 103101..#####
		# Neg Control - 101111..#####
		"""

		# Assume valid
		is_valid = True

		# LAB CONTROL SETTINGS
		flu_ctrl = r'101101.*'
		amr_ctrl = r'102101.*'
		par_ctrl = r'103101.*'
		neg_ctrl = r'101111.*'

		# The channel to compare and the new column to hold the PQ results
		channel_comparisons = {
			'FAM Rounded RFU Range (HPIV-1)':['PQ-FAM-RFU','FAM'],
			'HEX Rounded RFU Range (HPIV-2)':['PQ-HEX-RFU','HEX'],
			'ROX Rounded RFU Range (HPIV-3)':['PQ-ROX-RFU','ROX'],
			'RED647 Rounded RFU Range (HPIV-4)':['PQ-RED647-RFU','RED647']
		}

		for channel, pq_result in channel_comparisons.items():
			combined_file[pq_result[0]] = combined_file.loc[
				combined_file['Specimen Barcode'].str.contains('PANEL', case=False) | 
				combined_file['Specimen Barcode'].str.contains(par_ctrl)
			][channel].apply(lambda rfu: self._pq_threshold(rfu, pq_result[1]))

		# Check PQs
		failed_positive_pqs = self._pq_getposfailedsamples(combined_file)

		# Check for validity
		invalid_positives = self._check_invalid_positives(combined_file, par_ctrl)
		invalid_negatives = self._check_invalid_negatives(combined_file, neg_ctrl)

		# Check for false positives
		false_pos = self._check_false_positives(combined_file, neg_ctrl)

		overall_result = self._overall_validity(invalid_positives,invalid_negatives,false_pos,failed_positive_pqs,combined_file)

		key_map = {
			'invalid_positives':invalid_positives,
			'invalid_negatives':invalid_negatives,
			'is_false_positive':false_pos,
			'is_failed_positive':failed_positive_pqs
		}

		fields_to_write = {}

		for validity_cat, status in overall_result.items():
			if status == True:
				fields_to_write[validity_cat] = key_map[validity_cat]
				is_valid = False
		
		self._write_out(fields_to_write, save_to, file_id, is_valid)

		return is_valid

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

	def _overall_validity(self, invalid_pos, invalid_neg, false_positives, failed_positives, combined_file):

		""" Overall Validity is determined on the following criterias
		(1) No invalids
		(2) No false positive samples
		(3) Controls have to work
		(4) PQ thresholds have to be met
			(4a) RED647 can have up to 2 failures

		Args:
			invalid_pos - a list of 'positive' samples that were marked as invalid (includes controls)
			invalid_neg -  a list of 'negative' samples that were marked as invalid (includes controls)
			false_positives - a list of false positive samples (includes controls)
			failed_positives - a dictionary with channels and samples that failed in the specific channel  (includes controls)
			combined_file - modified dataframe from _pq_threshold (obj)
		Returns:
			validity_checks - a dictionary containing all the categories(key) and if they passed the validity checks(values)
		"""
		
		# Assume everything is correct at first
		validity_checks = {
			'invalid_positives':False,
			'invalid_negatives':False,
			'is_false_positive':False,
			'is_failed_positive':False
		}

		if len(invalid_pos) > 1:
			validity_checks['invalid_positives'] = True
	
		if len(invalid_neg) > 1:
			validity_checks['invalid_negatives'] = True
	
		if len(false_positives) > 1:
			validity_checks['is_false_positive'] = True

		# Can only have one 'failed' to meet PQ threshold
		# Exception:
		# We can have up to 2 'failed to meet PQ threshold' for RED647

		for channel_type, channel_column in failed_positives.items():

			if channel_type == 'HEX':
				if len(channel_column) > 2:
					validity_checks['is_failed_positive'] = True
					break
			if len(channel_column) > 1:
				validity_checks['is_failed_positive'] = True
				break

		return validity_checks

	def _check_invalid_negatives(self, combined_file, neg_ctrl_id):

		""" Checks if any of the 'negative' samples are invalid
		(1) If 'negative' samples are invalid and IC is invalid - invalid run

		A 'negative' sample here is something that has no organisms in it 
		at all. Completely blank and we expect it to be blank.

		Args:
			combined_file - modified dataframe from _pq_threshold (obj) 
			neg_ctrl_id - the negative control id (regex str) 
		Returns:
			a list of invalid results. The list contains a list of each row
			which correspondings to one sample. (Specimen Barcode) - (Run ID) - (Test order #) (list)
		"""

		neg_column = combined_file.loc[combined_file['Specimen Barcode'].str.contains('neg', case=False)]

		neg_aggregated = neg_column[
			~(neg_column['POS/NEG/Invalid for HPIV-1']).str.contains('neg', case=False) | 
			~(neg_column['POS/NEG/Invalid for HPIV-2']).str.contains('neg', case=False) |
			~(neg_column['POS/NEG/Invalid for HPIV-3']).str.contains('neg', case=False) |
			~(neg_column['POS/NEG/Invalid for HPIV-4']).str.contains('neg', case=False) |
			~(neg_column['Valid/Invalid for IC']).str.contains('valid', case=False)
		]

		# Adds a check for the negative control

		neg_ctrl_check = combined_file.loc[combined_file['Specimen Barcode'].str.contains(neg_ctrl_id)]
		invalid_neg_ctrl = neg_ctrl_check[
			~(neg_ctrl_check['POS/NEG/Invalid for HPIV-1']).str.contains('neg', case=False) | 
			~(neg_ctrl_check['POS/NEG/Invalid for HPIV-2']).str.contains('neg', case=False) |
			~(neg_ctrl_check['POS/NEG/Invalid for HPIV-3']).str.contains('neg', case=False) |
			~(neg_ctrl_check['POS/NEG/Invalid for HPIV-4']).str.contains('neg', case=False) |
			~(neg_ctrl_check['Valid/Invalid for IC']).str.contains('valid', case=False)
		]

		neg_agg_list = neg_aggregated[['Specimen Barcode','Run ID','Test order #']].values.tolist()
		neg_ctrl_agg_list = invalid_neg_ctrl[['Specimen Barcode','Run ID','Test order #']].values.tolist()
		final_result = neg_agg_list + neg_ctrl_agg_list

		return final_result

	def _check_invalid_positives(self, combined_file, pos_ctrl_id):

		""" Checks if any of the non-negative samples 'positive' are invalid.
		(1) If non-negative samples are invalid - it is invalid NO MATTER WHAT IC IS
		(2) If non-negative samples are valid - it is valid NO MATTER WHAT IC IS
		(3) If non-negative samples are negative - it is valid if IC IS POSITIVE, INVALID otherwise	

		A non-negative sample here is something we know that isn't used as a negative control. 
		It has some sort of product or we think it does.

		Args:
			combined_file - modified dataframe from _pq_threshold (obj)
			pos_ctrl_id - the positive control id (regex str) 
		Returns:
			a list of invalid results. The list contains a list of each row
			which correspondings to one sample. (Specimen Barcode) - (Run ID) - (Test order #) (list)
		"""

		pos_column = combined_file.loc[combined_file['Specimen Barcode'].str.contains('PANEL', case=False)]

		# What about 'Invalid positives' BUT 'Valid IC' --- ASK THE TEAM
		pos_aggregated = pos_column[
		    (pos_column['POS/NEG/Invalid for HPIV-1'].str.contains('Invalid', case=False) & pos_column['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_column['POS/NEG/Invalid for HPIV-2'].str.contains('Invalid', case=False) & pos_column['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_column['POS/NEG/Invalid for HPIV-3'].str.contains('Invalid', case=False) & pos_column['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_column['POS/NEG/Invalid for HPIV-4'].str.contains('Invalid', case=False) & pos_column['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_column['POS/NEG/Invalid for HPIV-1'].str.contains('neg', case=False) & pos_column['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_column['POS/NEG/Invalid for HPIV-2'].str.contains('neg', case=False) & pos_column['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_column['POS/NEG/Invalid for HPIV-3'].str.contains('neg', case=False) & pos_column['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_column['POS/NEG/Invalid for HPIV-4'].str.contains('neg', case=False) & pos_column['Valid/Invalid for IC'].str.contains('Invalid', case=False))
		]

		# Adds a check for the positive control

		pos_ctrl_check = combined_file.loc[combined_file['Specimen Barcode'].str.contains(pos_ctrl_id)]

		invalid_pos_ctrl = pos_ctrl_check[
		    (pos_ctrl_check['POS/NEG/Invalid for HPIV-1'].str.contains('Invalid', case=False) & pos_ctrl_check['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_ctrl_check['POS/NEG/Invalid for HPIV-2'].str.contains('Invalid', case=False) & pos_ctrl_check['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_ctrl_check['POS/NEG/Invalid for HPIV-3'].str.contains('Invalid', case=False) & pos_ctrl_check['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_ctrl_check['POS/NEG/Invalid for HPIV-4'].str.contains('Invalid', case=False) & pos_ctrl_check['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_ctrl_check['POS/NEG/Invalid for HPIV-1'].str.contains('neg', case=False) & pos_ctrl_check['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_ctrl_check['POS/NEG/Invalid for HPIV-2'].str.contains('neg', case=False) & pos_ctrl_check['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_ctrl_check['POS/NEG/Invalid for HPIV-3'].str.contains('neg', case=False) & pos_ctrl_check['Valid/Invalid for IC'].str.contains('Invalid', case=False)) |
		    (pos_ctrl_check['POS/NEG/Invalid for HPIV-4'].str.contains('neg', case=False) & pos_ctrl_check['Valid/Invalid for IC'].str.contains('Invalid', case=False))
		]

		pos_agg_list = pos_aggregated[['Specimen Barcode','Run ID','Test order #']].values.tolist()
		pos_ctrl_agg_list = invalid_pos_ctrl[['Specimen Barcode','Run ID','Test order #']].values.tolist()
		final_result = pos_agg_list + pos_ctrl_agg_list

		return final_result

	def _check_false_positives(self, combined_file, neg_ctrl_id):

		""" Checks if any of the 'negative' samples are false positives
		Args:
			combined_file - modified dataframe from _pq_threshold (obj)
			neg_ctrl_id - the negative control id (regex str) 
		Returns:
			a list of false positive results. The list contains a list of each row
			which correspondings to one sample. (Specimen Barcode) - (Run ID) - (Test order #) (list)
		"""

		negative_columns = combined_file.loc[combined_file['Specimen Barcode'].str.contains('neg', case=False)]

		neg_aggregated_false_pos = negative_columns[
			(negative_columns['POS/NEG/Invalid for HPIV-1']).str.contains('pos', case=False) | 
			(negative_columns['POS/NEG/Invalid for HPIV-2']).str.contains('pos', case=False) |
			(negative_columns['POS/NEG/Invalid for HPIV-3']).str.contains('pos', case=False) |
			(negative_columns['POS/NEG/Invalid for HPIV-4']).str.contains('pos', case=False) 
		]

		# Adds a check for the negative control
		neg_ctrl_column = combined_file.loc[combined_file['Specimen Barcode'].str.contains(neg_ctrl_id)]
		neg_ctrl_false_pos = neg_ctrl_column[
			(neg_ctrl_column['POS/NEG/Invalid for HPIV-1']).str.contains('pos', case=False) | 
			(neg_ctrl_column['POS/NEG/Invalid for HPIV-2']).str.contains('pos', case=False) |
			(neg_ctrl_column['POS/NEG/Invalid for HPIV-3']).str.contains('pos', case=False) |
			(neg_ctrl_column['POS/NEG/Invalid for HPIV-4']).str.contains('pos', case=False) 
		]

		neg_false_pos = neg_aggregated_false_pos[['Specimen Barcode','Run ID','Test order #']].values.tolist()
		neg_ctrl_false_pos = neg_ctrl_false_pos[['Specimen Barcode','Run ID','Test order #']].values.tolist()
		final_result = neg_false_pos + neg_ctrl_false_pos

		return final_result

	def _pq_threshold(self, rfu_range, channel_type):
		""" Checks if the RFU range meets the specifications set in the 
		PQ (performance quality) documentation.

		Args:
			rfu_range - the RFU range of the specified channel (str)
			channel_type - the type of channel (str)
		Returns
			(str) containing the values 'pass' or 'fail' depending on the criteria
		"""

		# PQ THRESHOLD SETTINGS
		FAM_RFU_THRESHOLD_MIN = 1200
		FAM_RFU_THRESHOLD_MAX = None
		HEX_RFU_THRESHOLD_MIN = 2000
		HEX_RFU_THRESHOLD_MAX = None
		ROX_RFU_THRESHOLD_MIN = 1500
		ROX_RFU_THRESHOLD_MAX = None
		RED647_RFU_THRESHOLD_MIN = 400
		RED647_RFU_THRESHOLD_MAX = None
		IC_RFU_THRESHOLD_MIN = None
		IC_RFU_THRESHOLD_MAX = None

		THRESHOLD_SETTINGS = {
		    'FAM': [FAM_RFU_THRESHOLD_MIN, FAM_RFU_THRESHOLD_MAX],
		    'HEX': [HEX_RFU_THRESHOLD_MIN, HEX_RFU_THRESHOLD_MAX],
		    'ROX': [ROX_RFU_THRESHOLD_MIN, ROX_RFU_THRESHOLD_MAX],
		    'RED647': [RED647_RFU_THRESHOLD_MIN, RED647_RFU_THRESHOLD_MAX],
		    'IC': [IC_RFU_THRESHOLD_MIN, IC_RFU_THRESHOLD_MAX],
		}

		# Determine the comparison cases

		# Comparison Modes
		# 0 = checks if RFU range is greater than min
		# 1 = checks if RFU range is less than max
		# 2 = checks if RFU range is between min and max
		comparison_mode = 2
		channel_settings = THRESHOLD_SETTINGS[channel_type]

		# Quick hacky fix for postive samples with no data
		if (rfu_range == '-'):
			rfu_range = 0

		rfu_range = float(rfu_range)


		# If both are 'None', mark as 'pass' by default 
		if channel_settings[0] == None and channel_settings[1] == None:
		    return 'pass'

		for number, check_none in enumerate(channel_settings):
		    if check_none is not None:
		        comparison_mode = number

		if comparison_mode == 0:
		    if rfu_range > channel_settings[0]:
		        return 'pass'
		elif comparison_mode == 1:
		    if rfu_range < channel_settings[1]:
		        return 'pass'
		else:
		    if channel_settings[0] < rfu_range < channel_settings[1]:
		        return 'pass'

		return 'fail'		

	def _pq_getposfailedsamples(self, combined_file):
		""" Get all of the samples marked as 'fail' (False)
		from the _pq_threshold method.

		Args:
			combined_file - modified dataframe from _pq_threshold (obj)
		Returns:
			a dict of invalid results. The dict contains a channel(key) & list(value) of each row
			which correspondings to one sample. (Specimen Barcode) - (Run ID) - (Test order #) (list)
		"""

		channels = {
			'FAM':'PQ-FAM-RFU',
			'HEX':'PQ-HEX-RFU',
			'ROX':'PQ-ROX-RFU',
			'RED647':'PQ-RED647-RFU'
		}

		failed_results = {}

		for channel_type, channel_column in channels.items():
			pq_group = combined_file.loc[
				combined_file['Specimen Barcode'].str.contains('PANEL', case=False) &
				combined_file[channel_column].str.contains('fail')
			]
			failed_results[channel_type] = pq_group[['Specimen Barcode','Run ID','Test order #']].values.tolist()

		return failed_results

class FusionCombiner(FusionAnalysis):

	def __init__(self,  reference_dataframe, *args, **kwargs):

		self.list_of_dataframes = [df for df in args]
		self.list_of_dataframes.append(reference_dataframe)
		self.combination = pd.concat(self.list_of_dataframes)

	def printout(self):
		pd.set_option('display.max_rows', None)
		pd.set_option('display.max_columns', None)
		pd.set_option('display.width', 1000)
		pd.set_option('display.height', 1000)
		pd.set_option('display.expand_frame_repr', False)
		print(self.combination)
		# self.combination.to_excel(os.path.join(os.getcwd(),'sampletest_combineds.xlsx'))

	def aggr(self):
		g = self.combination[self.combination['Specimen Barcode'].str.contains('Panel', case=False)]
		h = g[~g['FAM Rounded Ct'].str.contains('Invalid', case=False)]
		ok = h['FAM Rounded Ct'].astype(float)
		print(ok.mean())


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