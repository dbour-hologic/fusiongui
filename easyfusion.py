import sys, os
from PyQt4 import QtGui


class Fusion(QtGui.QWidget):

	def __init__(self):
		super(Fusion, self).__init__()

		self.initUI()

		self.pcr_files = []
		self.lis_files = []

	def initUI(self):

		grid = QtGui.QGridLayout()

		lisUploadBox = QtGui.QGroupBox("LIS Uploads", self)
		pcrUploadBox = QtGui.QGroupBox("PCR Uploads", self)
		assayTypeBox = QtGui.QGroupBox("Assay Types", self)
		executeBox = QtGui.QGroupBox(self)

		upload_lis_button = QtGui.QPushButton("Upload LIS File(s)")
		upload_pcr_button = QtGui.QPushButton("Upload PCR File(s)")


		self.lisFileList = QtGui.QListWidget()
		self.pcrFileList = QtGui.QListWidget()

		vboxLISFiles = QtGui.QVBoxLayout()
		vboxLISFiles.addWidget(upload_lis_button)
		vboxLISFiles.addWidget(self.lisFileList)
		vboxLISFiles.addStretch(1)

		vboxPCRFiles = QtGui.QVBoxLayout()
		vboxPCRFiles.addWidget(upload_pcr_button)
		vboxPCRFiles.addWidget(self.pcrFileList)
		vboxPCRFiles.addStretch(1)

		assay_type_paraflu = QtGui.QRadioButton("Paraflu")
		

		vboxAssay = QtGui.QVBoxLayout()
		vboxAssay.addWidget(assay_type_paraflu)
		vboxAssay.addStretch(1)

		status_msg = QtGui.QTextEdit()
		status_msg.setReadOnly(True)
		status_msg.setStyleSheet("background-color: #EEF3F9;")
		execute_run = QtGui.QPushButton("Combine Files")

		vboxRun = QtGui.QVBoxLayout()
		vboxRun.addWidget(status_msg)
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

		upload_lis_button.clicked.connect(self.populateFields)
		upload_pcr_button.clicked.connect(self.populateFields)

		self.setGeometry(300,300,500,250)
		self.setWindowTitle('Fusion LIS & PCR Combiner')
		self.show()

	def populateFields(self):
		"""
		Populates a list widget with the files selected
		"""
		
		# Sets the type of list to populate (LIS/PCR)
		list_to_populate = None

		if (self.sender().text() == 'Upload LIS File(s)'):
			list_to_populate = self.lisFileList
		elif (self.sender().text() == 'Upload PCR File(s)'):
			list_to_populate = self.pcrFileList
		else:
			print("Error", self.sender().text())

		file_names = QtGui.QFileDialog.getOpenFileNames(self, "Select Files", '/home', "*.lis *.csv")

		if file_names:
			for file_name in file_names:
				list_to_populate.addItem(file_name)

	



def main():

	app = QtGui.QApplication(sys.argv)
	fusion_widget = Fusion()
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()