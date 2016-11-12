import sys
import fipy
import meshio as mio

# Qt imports
from PyQt5 import QtCore, QtWidgets

# Matplotlib-Qt imports
import matplotlib
matplotlib.use('Qt5Agg')  # Make sure that we are using QT5
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from scarce import (constant, plot, geometry)

PROJECT_NAME = 'Scarce'
WIDTH = 1024
HEIGHT = 768

class SensorCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        points, cells = geometry.mesh_3D_sensor(x=250,
                                    y=50,
                                    n_pixel_x=1, 
                                    n_pixel_y=1,
                                    radius=6,
                                    nD=2,
                                    resolution=50)
                                     
        mio.write('sensor.msh', points, cells)
        mesh = fipy.GmshImporter2D('sensor.msh')
        
        fig = Figure(figsize=(width, height), dpi=dpi)
        plot.get_mesh_plot(fig, mesh)
        
        fig.get_axes()[0].set_aspect('equal')

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setWindowTitle(PROJECT_NAME)
        self.init_UI()
        
    def init_UI(self):
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        self.menu()

        self.main_widget = QtWidgets.QWidget(self)

        l = QtWidgets.QVBoxLayout(self.main_widget)
        sc = SensorCanvas(self.main_widget, width=50, height=40, dpi=100)
        
        # Add zoom / safe navigation bar
        NavigationToolbar(sc, sc)
        
        l.addWidget(sc)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.statusBar().showMessage("Hello and welcome!", 2000)
        
    def menu(self):
        self.file_menu = QtWidgets.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtWidgets.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.help_menu.addAction('&About', self.about)

    def fileQuit(self):
        self.close()

    def closeEvent(self, _):
        self.fileQuit()

    def about(self):
        QtWidgets.QMessageBox.about(self, "About",
                                    """Scarce version %s""" % constant.VERSION)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    aw = ApplicationWindow()
    aw.show()
    sys.exit(app.exec_())