import sys
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
		assayTypeBox = QtGui.QGroupBox("Assay Type", self)
		executeBox = QtGui.QGroupBox(self)

		upload_lis_button = QtGui.QPushButton("Upload LIS File")
		upload_pcr_button = QtGui.QPushButton("Upload PCR File")


		lisFileList = QtGui.QListWidget()
		pcrFileList = QtGui.QListWidget()

		vboxLISFiles = QtGui.QVBoxLayout()
		vboxLISFiles.addWidget(upload_lis_button)
		vboxLISFiles.addWidget(lisFileList)
		vboxLISFiles.addStretch(1)

		vboxPCRFiles = QtGui.QVBoxLayout()
		vboxPCRFiles.addWidget(upload_pcr_button)
		vboxPCRFiles.addWidget(pcrFileList)
		vboxPCRFiles.addStretch(1)

		assay_type_paraflu = QtGui.QRadioButton("Paraflu")
		assay_type_flu = QtGui.QRadioButton("Flu")

		vboxAssay = QtGui.QVBoxLayout()
		vboxAssay.addWidget(assay_type_paraflu)
		vboxAssay.addWidget(assay_type_flu)
		vboxAssay.addStretch(1)

		execute_run = QtGui.QPushButton("Run")

		vboxRun = QtGui.QVBoxLayout()
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

		self.setGeometry(300,300,500,250)
		self.setWindowTitle('Fusion LIS & PCR Combiner')
		self.show()

	def populateFields(self):
		"""
		Populates a list widget with the files selected
		"""
		pass

def main():

	app = QtGui.QApplication(sys.argv)
	fusion_widget = Fusion()
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()