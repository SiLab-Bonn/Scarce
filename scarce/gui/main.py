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
        # Sensor figure
        fig = Figure(figsize=(width, height), dpi=dpi)

        # Plot 3D sensor array for now
        width_x, width_y = 250., 50.
        n_pixel_x, n_pixel_y = 3, 3
        radius, nD = 6., 2
        resolution = 10
        mesh = geometry.mesh_3D_sensor(width_x=width_x,
                                       width_y=width_y,
                                       n_pixel_x=n_pixel_x,
                                       n_pixel_y=n_pixel_y,
                                       radius=radius,
                                       nD=nD,
                                       resolution=resolution)

        plot.get_3D_sensor_plot(fig=fig,
                                width_x=width_x,
                                width_y=width_y,
                                radius=radius,
                                nD=nD,
                                n_pixel_x=n_pixel_x,
                                n_pixel_y=n_pixel_y,
                                V_bias=0,
                                V_readout=1,
                                potential_function=None,
                                field_function=None,
                                mesh=mesh,
                                title=None)

        # Init figure
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
