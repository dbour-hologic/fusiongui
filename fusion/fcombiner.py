""" Leverages the power of the single file combiner to create a 
big fat combined file for FUSION data. 
"""

from .fanalyzer import FusionAnalysis

class FusionCombiner():

	def __init__(self, file_list, assay_type, *args, **kwargs):
		""" Constructor is used to hold the big fat list of files.
		It will automatically handle matching pairs of PCR & LIS.
		It will log files that did not match or does not have data.

		Args:
			file_list - list of  full length file paths (list)
			assay_type - assay type to look out for (str)
		"""
		self.file_list = file_list
		self.all_items = self.__match_files()

	def __match_files(self):
		""" Match files will pair up the LIS & PCR files. Unpaired
		files will be moved separated out from the pack. 

		Args:
			None
		Returns:
			match_results - a dictionary containing  (dict)
			(1) no_pairs : <list of files>
			(2) paired: <list of {mapped} files>
		"""

		no_pairs = []
		paired = {}

		# LOGICAL ORDERING IS INCORRECT IN THIS
		# HAVE TO LOOP THROUGH ONE FIRST TO 
		# ESTBALISH THE DICT, THEN CHECK WITH OTHER

		# for files in self.file_list:

		# 	if files.find("@DI") > -1:
		# 		get_only_pcr_fname = files[files.find("@DI"):]
		# 		partition_pcr_file = get_only_pcr_fname.split("-")
		# 		pcr_unique_id = partition_pcr_file[0].replace("@DI","") \
		# 		+ partition_pcr_file[3] + "_" \
		# 		+ partition_pcr_file[4] + "_" \
		# 		+ partition_pcr_file[5].replace(".csv","")

		# 		paired[pcr_unique_id] = [files]

		# 	elif files.find("@Pt2") > -1:
		# 		get_only_lis_fname = files[files.find("@Pt2"):]
		# 		partition_lis_file = get_only_lis_fname.split("-")
		# 		lis_unique_id = partition_lis_file[0].replace("@Pt2","") \
		# 		+ partition_lis_file[3] + "_" \
		# 		+ partition_lis_file[4] + "_" \
		# 		+ partition_lis_file[5].replace(".lis","")

		# 		# paired[lis_unique_id]
		# 		print(lis_unique_id)
		# 		print("PP", paired)

		# 		has_pair = paired.get(lis_unique_id, None)
		# 		if has_pair:
		# 			print("Hola")
		# 			paired[lis_unique_id].append(files)
		# 		else:
		# 			no_pairs.append(files)
		# 	else:
		# 		no_pairs.append(files)




		# PCR Keys to remove that have no pair
		# pcr_keys_to_remove = []

		# for uniqueID, pairs in paired.items():
		# 	if (len(pairs) < 2):
		# 		pcr_keys_to_remove.append(uniqueID)
		# 		no_pairs.append(pairs)

		# for pcr_keys in pcr_keys_to_remove:
		# 	del paired[pcr_keys]

		# final_result = {"paired":paired, "no_pairs":no_pairs}

		return final_result

	def get_logs(self):
		pass