import unittest
import os
import sys

from PyQt5 import QtWidgets

from scarce.gui import main


class TestGui(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Linux CI travis runs headless, thus virtual x server is needed for gui testing
        if os.getenv('TRAVIS', False):
            from xvfbwrapper import Xvfb
            cls.vdisplay = Xvfb()
            cls.vdisplay.start()

        # Create Gui
        cls.app = QtWidgets.QApplication(sys.argv)
        cls.gui = main.ApplicationWindow()

    @classmethod
    def tearDownClass(cls):
        cls.gui.close()

    def test_device_tab(self):
        ''' Not implemented
        '''

        pass


if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    unittest.main()
