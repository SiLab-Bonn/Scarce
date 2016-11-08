''' Checks the semiempirical functions by comparing to test data. 
'''

import unittest
import numpy as np

from siliconproperties.python_files.getDepletionDepth import get_depletion_depth
from siliconproperties.python_files.getDepletionVoltage import get_depletion_voltage
from siliconproperties.python_files.getDiffusionPotential import get_diffusion_potential
from siliconproperties.python_files.getEffAcceptorConcentration import get_eff_acceptor_concentration
from siliconproperties.python_files.getFreePath import get_free_path
from siliconproperties.python_files.getMobility import get_mobility
from siliconproperties.python_files.getResistivity import get_resistivity
from siliconproperties.python_files.getTrapping import get_trapping

from siliconproperties.testing import tools


class Test(unittest.TestCase):


    def test_depletion_depth(self):
        self.assertTrue(tools.check_with_fixture(get_depletion_depth,
                                 V_bias=np.linspace(0, 100., 1000.),
                                 n_eff=[1, 10, 100], 
                                 temperature=[200, 250, 300]))
        
    def test_get_depletion_voltage(self):
        self.assertTrue(tools.check_with_fixture(get_depletion_voltage,
                             n_eff = np.linspace(0.1, 100., 1000.),
                             distance = range(10, 300, 10)))
        
    def test_diffusion_potential(self):
        self.assertTrue(tools.check_with_fixture(get_diffusion_potential,
                                 n_eff=np.linspace(0.1, 100., 1000.),
                                 temperature=range(200, 350, 20)))
        
    def test_eff_acceptor_concentration(self):
        self.assertTrue(tools.check_with_fixture(get_eff_acceptor_concentration,
                             fluence = np.logspace(0., 4., 1000.),
                             n_eff_0=[1, 1.5, 1.8, 2.0],
                             is_ntype=[True, False], 
                             is_oxygenated=[True, False]))
        
    def test_free_path(self):
        self.assertTrue(tools.check_with_fixture(get_free_path,
                             fluence = np.logspace(12., 15., 1000.),
                             e_field=[1e4, 1e5, 1e6],
                             temperature=range(150, 350, 50), 
                             is_electron=[False, True]))
         
    def test_get_mobility(self):
        self.assertTrue(tools.check_with_fixture(get_mobility,
                             e_field = np.logspace(3., 5., 1000.),
                             temperature=range(150, 350, 50), 
                             is_electron=[False, True]))
        
    def test_resistivity(self):
        self.assertTrue(tools.check_with_fixture(get_resistivity,
                             n_eff = np.logspace(11., 15., 1000.),
                             is_n_type=[False, True], 
                             temperature=range(150, 350, 50), 
                             e_field=[1e3, 1e4, 1e5, 1e6]))


    def test_trapping(self):
        self.assertTrue(tools.check_with_fixture(get_trapping,
                             fluence = np.logspace(12., 15., 1000.),
                             is_electron=[False, True], 
                             paper=[1, 2]))



if __name__ == "__main__":
    unittest.main()