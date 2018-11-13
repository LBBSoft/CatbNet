from App import Application
from PyQt4 import QtGui

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = Application()
    ui.show()
    sys.exit(app.exec_())
