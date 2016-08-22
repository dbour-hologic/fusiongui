""" Leverages the power of the single file combiner to create a 
big fat combined file for FUSION data. 
"""

from .fanalyzer import FusionAnalysis
import pandas as pd

class FusionCombiner():

	def __init__(self, file_list, assay_type, instrument_mapping=None, software_mapping=None, *args, **kwargs):
		""" Constructor is used to hold the big fat list of files.
		It will automatically handle matching pairs of PCR & LIS.
		It will log files that did not match or does not have data.

		Args:
			file_list - list of  full length file paths (list)
			assay_type - assay type to look out for (str)
		"""

		# Contains all of the files to combine
		self.file_list = file_list

		# Contains all of the matched files, unpaired files, and unknown files
		self.all_items = self.__match_files()

		# Contains all of the individual valid combined files
		self.all_combined_items = self.__generate_file_combiner_obj("P 1/2/3/4")

		# Contains the MEGA MERGED FILE
		self.mega_combination = pd.concat(self.all_combined_items['valid_fusion'])

	def write_combined_multiple(self, save_to):
		""" Saves the mega merged dataframes as an Excel file. This is done on purpose
		to prevent the long string of numbers from being automatically converted
		in Excel when opening it as a .csv/.tsv.

		Args:
			save_to - path to save combined file
		Returns:
			None
		Output:
			mega merged file in *.xlsx format
		"""

		self.mega_combination.to_excel(save_to)

	def __generate_file_combiner_obj(self, assay_type, instrument_mapping=None, software_mapping=None):
		""" (PRIVATE) Creates all valid objects by checking if the
		data <assay_type> is the proper assay type. 

		i.e. If looking for Paraflu,
		LIS file must be analyte "P 1/2/3/4"
		PCR file must be analyte "P 1/2/3/4"

		Args:
			assay_type - assay type to look out for (str)
		Returns:
			combined_list - a mapping of invalid and valid fusion data; based on
			assay type specified. (dict)
		"""

		combined_list = {"invalid_fusion":[], "valid_fusion":[]}

		for unique_id, list_of_pairs in self.all_items.get("paired").items():

			pcr_file = ""
			lis_file = ""

			for file in list_of_pairs:

				if file.find("@DI") > -1:
					pcr_file = file
			
				if file.find("@Pt2") > -1:
					lis_file = file

			fusion_combined = FusionAnalysis(pcr_file, lis_file, "P 1/2/3/4")

			if fusion_combined.valid_data:
				combined_list['valid_fusion'].append(fusion_combined.combine_files("P 1/2/3/4", instrument_mapping=fusion_combined.mapping_source, software_mapping=fusion_combined.software_source))
			else:
				combined_list['invalid_fusion'].append(fusion_combined)

		return combined_list

	def __match_files(self):
		""" (PRIVATE) Match files will pair up the LIS & PCR files. Unpaired
		files will be separated out from the pack. Matches are dependent
		on how the files are named; does not check the content.

		Args:
			None
		Returns:
			match_results - a dictionary containing  (dict)
			(1) no_pairs : <list of files>
			(2) paired: <list of {mapped} files>
			(3) random: <random list of files that don't match PCR or LIS>
		"""

		# Categories to return
		no_pairs = []
		paired = {}
		random = []

		# Temporary list to separate out each type of file
		pcr_file_list = []
		lis_file_list = []
		random_file_list = []

		# ____ This section categorizes the files into PCR, LIS, or Unknown
		for files in self.file_list:
			if files.find("@DI") > -1:
				pcr_file_list.append(files)
			elif files.find("@Pt2") > -1:
				lis_file_list.append(files)
			else:
				random_file_list.append(files)

		# ____ This section does the mapping to find pairs
		for pcr_file in pcr_file_list:
			get_only_pcr_fname = pcr_file[pcr_file.find("@DI"):]
			partition_pcr_file = get_only_pcr_fname.split("-")
			pcr_unique_id = partition_pcr_file[0].replace("@DI","") \
			+ partition_pcr_file[3] + "_" \
			+ partition_pcr_file[4] + "_" \
			+ partition_pcr_file[5].replace(".csv","")

			paired[pcr_unique_id] = [pcr_file]

		for lis_file in lis_file_list:
			get_only_lis_fname = lis_file[lis_file.find("@Pt2"):]
			partition_lis_file = get_only_lis_fname.split("-")
			lis_unique_id = partition_lis_file[0].replace("@Pt2","") \
			+ partition_lis_file[3] + "_" \
			+ partition_lis_file[4] + "_" \
			+ partition_lis_file[5].replace(".lis","")

			has_pair = paired.get(lis_unique_id, None)
			if has_pair:
				paired[lis_unique_id].append(lis_file)
			else:
				no_pairs.append(lis_file)

		# ____ PCR Keys to remove that have no pair
		pcr_keys_to_remove = []

		for uniqueID, pairs in paired.items():
			if (len(pairs) < 2):
				pcr_keys_to_remove.append(uniqueID)
				no_pairs.append(pairs)

		for pcr_keys in pcr_keys_to_remove:
			del paired[pcr_keys]

		final_result = {"paired":paired, "no_pairs":no_pairs, "random":random}

		return final_result



	def get_logs(self):
		pass