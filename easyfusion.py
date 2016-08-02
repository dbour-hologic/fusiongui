import sys, os
from PyQt4 import QtGui

class FusionGui(QtGui.QWidget):

	def __init__(self):
		super(FusionGui, self).__init__()

		self.initUI()

		self.pcr_files = []
		self.lis_files = []

	def initUI(self):

		grid = QtGui.QGridLayout()

		lisUploadBox = QtGui.QGroupBox("LIS Uploads", self)
		pcrUploadBox = QtGui.QGroupBox("PCR Uploads", self)
		assayTypeBox = QtGui.QGroupBox("Assay Types", self)
		executeBox = QtGui.QGroupBox(self)

		upload_lis_button = QtGui.QPushButton("Upload LIS Files") 
		upload_pcr_button = QtGui.QPushButton("Upload PCR Files")

		clear_upload_lis_button = QtGui.QPushButton("Clear LIS Uploads")
		clear_upload_pcr_button =  QtGui.QPushButton("Clear PCR Uploads")


		self.lisFileList = QtGui.QListWidget()
		self.pcrFileList = QtGui.QListWidget()

		vboxLISFiles = QtGui.QVBoxLayout()
		vboxLISFiles.addWidget(upload_lis_button)
		vboxLISFiles.addWidget(self.lisFileList)
		vboxLISFiles.addWidget(clear_upload_lis_button)
		vboxLISFiles.addStretch(1)

		vboxPCRFiles = QtGui.QVBoxLayout()
		vboxPCRFiles.addWidget(upload_pcr_button)
		vboxPCRFiles.addWidget(self.pcrFileList)
		vboxPCRFiles.addWidget(clear_upload_pcr_button)
		vboxPCRFiles.addStretch(1)

		assay_type_paraflu = QtGui.QRadioButton("Paraflu")
		assay_type_paraflu.setChecked(True)
		

		vboxAssay = QtGui.QVBoxLayout()
		vboxAssay.addWidget(assay_type_paraflu)
		vboxAssay.addStretch(1)

		self.status_msg = QtGui.QTextEdit()
		self.status_msg.setReadOnly(True)
		self.status_msg.setStyleSheet("background-color: #EEF3F9;")
		execute_run = QtGui.QPushButton("Save && Combine Files")

		vboxRun = QtGui.QVBoxLayout()
		vboxRun.addWidget(self.status_msg)
		vboxRun.addWidget(execute_run)

		lisUploadBox.setLayout(vboxLISFiles)
		pcrUploadBox.setLayout(vboxPCRFiles)
		assayTypeBox.setLayout(vboxAssay)
		executeBox.setLayout(vboxRun)
		
		grid.addWidget(lisUploadBox, 1, 1)
		grid.addWidget(pcrUploadBox, 1, 2)
		grid.addWidget(assayTypeBox, 1, 3)
		grid.addWidget(executeBox, 2, 1, 1, 3)
		self.setLayout(grid)

		upload_lis_button.clicked.connect(self.populate_fields)
		upload_pcr_button.clicked.connect(self.populate_fields)
		clear_upload_lis_button.clicked.connect(self.clear_upload_fields)
		clear_upload_pcr_button.clicked.connect(self.clear_upload_fields)
		execute_run.clicked.connect(self.run_program)

		self.setGeometry(300,300,700,200)
		self.setWindowTitle('Fusion LIS & PCR Combiner')
		self.show()

	def clear_upload_fields(self):
		"""
		Clears the upload list of PCR uploads or LIS uploads
		"""

		if (self.sender().text() == 'Clear LIS Uploads'):
			self.lisFileList.clear()
		elif (self.sender().text() == 'Clear PCR Uploads'):
			self.pcrFileList.clear()
		else:
			print("Error", self.sender().text())

	def populate_fields(self):
		"""
		Populates a list widget with the files selected
		"""
		
		# Sets the type of list to populate (LIS/PCR)
		list_to_populate = None
		# Set the file type
		file_type = None

		if (self.sender().text() == 'Upload LIS Files'):
			list_to_populate = self.lisFileList
			file_type = "*.lis"
		elif (self.sender().text() == 'Upload PCR Files'):
			list_to_populate = self.pcrFileList
			file_type = "*.csv"
		else:
			print("Error", self.sender().text())

		file_names = QtGui.QFileDialog.getOpenFileNames(self, "Select Files", '/home', "%s" % file_type)

		if file_names:
			for file_name in file_names:
				list_to_populate.addItem(file_name)

	def run_program(self):
		""" Executes the program """

		from fusion.fanalyzer import FusionAnalysis

		save_directory = QtGui.QFileDialog.getExistingDirectory(self,'Select Save Directory')

		lisFileTextList = [str(self.lisFileList.item(i).text()) for i in range(self.lisFileList.count())]
		pcrFileTextList = [str(self.pcrFileList.item(x).text()) for x in range(self.pcrFileList.count())]
		
		match_database = self.match_files(lisFileTextList, pcrFileTextList)

		self.status_msg.clear()

		for identifier, keypairs in match_database.items():
			if ((identifier == 'missing_lis') or (identifier == 'missing_pcr')):
				if (identifier == 'missing_lis'):
					for lis in match_database['missing_lis']:
						self.status_msg.insertHtml('<b>MISSING FILE</b>: The LIS file %s is missing a PCR pair.<br>' % lis)
				if (identifier == match_database['missing_pcr']):
					for pcr in missing_pcr:
						self.status_msg.insertHtml('<b>MISSING FILE</b>: The PCR file %s is missing a LIS pair.<br>' % pcr)
			else:
				successful_run = FusionAnalysis(keypairs[0], keypairs[1], "P 1/2/3/4")
				if successful_run.check_assay_types():
					new_file_name = identifier + ".xlsx"
					save_file_path = os.path.join(str(save_directory), new_file_name)
					successful_run.combine_files("Paraflu", save_file_path)
				else:
					self.status_msg.insertHtml('<b>ASSAY TYPE WARNING</b>: The files %s are not of the specified assay type<br>' % keypairs)

		self.status_msg.insertHtml('<b>RUN IS COMPLETE. Results saved at %s</b>' % save_directory)
		


	def match_files(self, list_of_LIS_files, list_of_PCR_files):

		""" Main purpose of the program is to match
		the PCR file to its appropriate LIS file. The matching 
		is done on the worklist-ID and date.

		Args:
			list_of_LIS_files - all of the LIS files in a list (list)
			list_of_PCR_files - all of the PCR files in a list (list)
		Returns:
			the unique_id as a (key) and list of [lis_file, pcr_file], also
			notfies of missing pcr or lis files. (dict)
		"""

		# LIS FILES WITH NO PCR PAIR
		missing_lis_files = []
		# PCR FILES WITH NO LIS PAIR
		missing_pcr_files = []
		# MAIN INDEX
		fusion_file_index = {}

		for pcr_files in list_of_PCR_files:
			partition_pcr_file = pcr_files.split("-")
			get_pcr_header = partition_pcr_file[0][partition_pcr_file[0].find("@DI"):]
			unique_id = get_pcr_header.replace("@DI",'') + "_" + \
						partition_pcr_file[3] + "_" + \
						partition_pcr_file[4] + "_" + \
						partition_pcr_file[5].replace('.csv','')

			fusion_file_index[unique_id] = [pcr_files]

		for lis_files in list_of_LIS_files:
			partition_lis_file = lis_files.split("-")
			get_lis_header = partition_lis_file[0][partition_lis_file[0].find("@Pt2"):]
			lis_unique_id = get_lis_header.replace('@Pt2','') + "_" + \
							partition_lis_file[3] + "_" + \
							partition_lis_file[4] + "_" + \
							partition_lis_file[5].replace('.lis','')

		# Checks if the LIS file is paired with an available PCR file
			try:
				fusion_file_index[lis_unique_id].append(lis_files)
			except KeyError:
				missing_lis_files.append(lis_files)
				print("No matching PCR file found for %s" % lis_files)

	 	# Remove PCR Keys
	 	pcr_keys_to_remove = []

		# Checks if the PCR file is paired with a LIS file
		for uniqueID, pairs in fusion_file_index.items():
			if (len(pairs) < 2):
				print("No matching LIS file found for %s" % pairs[0])
				missing_pcr_files.append(pairs[0])
				# Remove the unique_id from analysis
				pcr_keys_to_remove.append(uniqueID)

		for pcr_keys in pcr_keys_to_remove:
			del fusion_file_index[pcr_keys]

		fusion_file_index['missing_lis'] = missing_lis_files
		fusion_file_index['missing_pcr'] = missing_pcr_files

		return fusion_file_index




def main():

	app = QtGui.QApplication(sys.argv)
	fusion_widget = FusionGui()
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()