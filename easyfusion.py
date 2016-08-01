import sys
from PyQt4 import QtGui


class Fusion(QtGui.QWidget):

	def __init__(self):
		super(Fusion, self).__init__()

		self.initUI()

	def initUI(self):

		self.setGeometry(300,300,500,300)
		self.setWindowTitle('Fusion LIS & PCR Analysis')
		self.show()


def main():

	app = QtGui.QApplication(sys.argv)
	fusion_widget = Fusion()
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()