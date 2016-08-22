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


	def set_labels(self, pos_label, neg_label, pos_ctrl, neg_ctrl):
		""" 
		Creates a temporary column that allows us to categorize the
		samples as a positive sample, a negative sample, or controls.

		Args:
			pos_label - a string pattern that is shared with all positive controls (str)
			neg_label - a string pattern that is shared with all negative controls (str)
			pos_ctrl - a regex pattern that corresponds to positive controls (str)
			neg_ctrl - a regex pattern that corresponds to negative controls (str)
		Returns:
			None
		"""
		self.dframe['SAMPLE Category'] = self.dframe['Specimen Barcode'].apply(lambda row: self.__set_labels_helper(row, pos_label, neg_label, pos_ctrl, neg_ctrl))

	def __set_labels_helper(self, row, pos_label, neg_label, positive_ctrl, negative_ctrl):
		""" The column modifier for set_labels method """

		import re

		pos_pattern = re.compile(r'%s' % positive_ctrl)
		neg_pattern = re.compile(r'%s' % negative_ctrl)

		if row.lower().find(pos_label) > -1 or pos_pattern.search(row.lower()):
			return "POS"
		elif row.lower().find(neg_label) > -1 or neg_pattern.search(row.lower()):
			return "NEG"
		else:
			return ""