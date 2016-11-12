import unittest
import os

from scarce.scripts import plot_properties


class Test(unittest.TestCase):

    def setUp(self):
        if os.getenv('TRAVIS', False):
            from xvfbwrapper import Xvfb
            self.vdisplay = Xvfb()
            self.vdisplay.start()

    def tearDown(self):
        filelist = [ f for f in os.listdir(".") if f.endswith(".pdf") ]
        for f in filelist:
            os.remove(f)

    def test_plot_properties(self):
        ''' Check script to plot all define silicon properties.
        '''
        plot_properties.create_plots()
        
if __name__ == "__main__":
    unittest.main()