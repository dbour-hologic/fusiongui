""" 
Fusion PQ (Performance Quality) analysis program. The following
program takes in a modified combined LIS & PCR file and performs
checks on the data to see if it meets specifications set by 
a guideline
"""

class FusionPQ():

	def __init__(self, dframe, assay_type, pos_label, neg_label, *args, **kwargs):

			self.dframe = dframe

			if assay_type == "P 1/2/3/4":
				self.settings = {
					"assay_type":assay_type,
					"pos_ctrl":'103101.*',
					"neg_ctrl":'101111.*',
					"pos_label":pos_label.lower(),
					"neg_label":neg_label.lower()
				}

	def run_pq(self):

		POS_CTRL = self.settings['pos_ctrl']
		NEG_CTRL = self.settings['neg_ctrl']
		POS_LBL = self.settings['pos_label']
		NEG_LBL = self.settings['neg_label']

		self.set_labels(POS_LBL, NEG_LBL, POS_CTRL, NEG_CTRL)
		self.check_validity()


	def set_labels(self, pos_label, neg_label, pos_ctrl, neg_ctrl):
		""" 
		Creates a temporary column that allows us to categorize the
		samples as a positive sample, a negative sample, or controls.

		Args:
			pos_label - a string pattern that is shared with all positive samples (str)
			neg_label - a string pattern that is shared with all negative samples (str)
			pos_ctrl - a regex pattern that corresponds to positive controls (str)
			neg_ctrl - a regex pattern that corresponds to negative controls (str)
		Returns:
			None
		"""
		self.dframe['SAMPLE Category'] = self.dframe['Specimen Barcode'].apply(lambda row: self.__set_labels_helper(row, pos_label, neg_label, pos_ctrl, neg_ctrl))

	def __set_labels_helper(self, row, pos_label, neg_label, positive_ctrl, negative_ctrl):
		""" The column marker logic for set_labels method """

		import re

		pos_pattern = re.compile(r'%s' % positive_ctrl)
		neg_pattern = re.compile(r'%s' % negative_ctrl)

		if row.lower().find(pos_label) > -1 or pos_pattern.search(row.lower()):
			return "POS"
		elif row.lower().find(neg_label) > -1 or neg_pattern.search(row.lower()):
			return "NEG"
		else:
			return ""

	def check_validity(self):
		""" Checks the validity status of each channel.
		This only checks if the samples are considered 'valid' which is
		different from the samples passing PQ criterias. The validity status
		will be stored in a temporary column.

		Samples are classified as follows:
		(1) If the interpretation of a channel is invalid, sample is INVALID.
		(2) If a 'positive' sample is positive and IC is valid, sample is VALID.
		(3) If a 'positive' sample is negative and IC is valid, sample is VALID.
		(4) If a 'positive' sample is positive and IC is invalid, sample is VALID. *
		(5) If a 'positive' sample is negative and IC is invalid, sample is INVALID.
		(6) If a 'negative' sample is positive and IC is invalid, sample is INVALID.
		(7) If a 'negative' sample is positive and IC is valid, sample is VALID.
		(8) If a 'negative' sample is negative and IC is valid, sample is VALID.
		(9) If a 'negative' sample is negative and IC is invalid, sample is INVALID.

		* Not a typo, invalid IC's can be overshadowed by a high positive value.
		"""

		columns_to_validate = {
			"FAM": "POS/NEG/Invalid for HPIV-1",
			"HEX": "POS/NEG/Invalid for HPIV-2",
			"RED647": "POS/NEG/Invalid for HPIV-4",
			"ROX": "POS/NEG/Invalid for HPIV-3"
		}

		for channel, column_name in columns_to_validate.items():
			results_column = "Validity for {0}".format(column_name)
			self.dframe[results_column] = self.dframe.apply(self.__check_validity_helper, axis=1, column_type=column_name)

	def __check_validity_helper(self, row, column_type):
		""" The column validity helper for check validity method 

		Args:
			row - the row to modify (pandas obj)
			column_type - the specific channel that's being interpreted. (str)
		Returns:
			'Valid' if the data is valid. (str)
			'Invalid' if the data is invalid. (str)
		"""

		import re

		channel_examine = row[column_type]

		# Helps differentiate between positive or negative samples to apply validity logic.
		sample_category = row['SAMPLE Category']

		ic_column = row['Valid/Invalid for IC']

		label_for_validity = re.compile(r'\bvalid\b')

		# Sample type 'positive'
		if sample_category.lower().find("pos") > -1:
			if channel_examine.lower().find("pos") > -1:
				return "VALID"
			elif channel_examine.lower().find("neg") > -1 and label_for_validity.search(ic_column.lower()):
				return "VALID"
			else:
				return "INVALID"
		# Sample type 'negative'
		elif sample_category.lower().find("neg") > -1:
			if channel_examine.lower().find("neg") > -1 and label_for_validity.search(ic_column.lower()):
				return "VALID"
			elif channel_examine.lower().find("pos") > -1 and label_for_validity.search(ic_column.lower()):
				return "VALID"
			else:
				return "INVALID"
		else:
				return ""

	def check_false_calls(self):
		""" 
		Has to be called after labels are set.

		Checks for False Positives or False Negatives. Creates a temporary
		column to hold the results for false positives or false negatives. If the
		data does not contain false pos/false neg, it will have no status marked in this column.

		Args:
			None
		Returns:
			None
		"""

		columns_to_examine_false = {
			"FAM": "POS/NEG/Invalid for HPIV-1",
			"HEX": "POS/NEG/Invalid for HPIV-2",
			"RED647": "POS/NEG/Invalid for HPIV-4",
			"ROX": "POS/NEG/Invalid for HPIV-3"    		
		}

		for channel, column_name in columns_to_examine_false.items():
			results_column = "TRUTHS table for {0}".format(column_name)
			self.dframe[results_column] = self.dframe.apply(self.__check_false_calls_helper, axis=1, column_type=column_name)

	def __check_false_calls_helper(self, row, column_type):
		""" Helper method for check false calls
		Args:
			row - the specific row to analyze (pandas obj)
			column_type - the specific channel that's being interpreted. (str)
		Returns:
			'False Negative' if result is false negative (str)
			'False Positive' if result is false positive (str)
		"""

		sample_barcode = row['SAMPLE Category']
		channel_column = row[column_type]

		if sample_barcode.lower().find("POS") > -1:
			if channel_column.lower().find("neg") > -1:
				return "False Negative"
			return ""
		elif sample_barcode.lower().find("NEG") > -1:
			if channel_column.lower().find("pos") > -1:
				return "False Positive"
			return ""
		else:
			return ""
