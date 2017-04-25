''' Minimizes the chi2 between simulation and measurement by finding better trapping
time constants.
'''

import os
import argparse

# import objgraph
from scipy.interpolate import UnivariateSpline
import numpy as np
import tables as tb
from collections import defaultdict
from gc import get_objects
from scipy.interpolate import interp1d

from scipy import optimize
from scarce import plot, solver, geometry, silicon, fields, tools, analysis


TEMP = 300
NPIXEL = 9
WIDTH = 50.
PITCH = 30.
THICKNESS_1 = 250.
THICKNESS_2 = 200.
RESOLUTION_1 = 287
RESOLUTION_2 = 251
SMOOTHING = 0.05
FLUENCE_1 = 1000.
FLUENCE_2 = 7000.
BIASES_1 = [-250, -350, -550, -800]
BIASES_1_2 = [-900, -1400]
BIASES_2 = [-250, -350, -550, -800]

NEFF0 = 1.475e12  # before irradiation


def chi2_to_spline(spl, x, y):
    ''' Calculated the chi2 for the spline description of points '''
    return np.sum(np.square(spl(x) - y))


def plot_cce_mpv(data_files, data_files_iv, start_i, s=1):
    def get_voltage(vbias, v, i):
        interpolation = interp1d(x=v, y=i, kind='slinear', bounds_error=True)
        # 100 kOhm * uA = 100 Ohm * mA
        return vbias - 100. * interpolation(vbias) * 0.001
    data = np.loadtxt(fname=data_files[0], dtype=np.float, skiprows=1)
    data = data[data[:, 0].argsort()[::-1]]
    data = data[1:]
    if data_files_iv:
        with tb.open_file(data_files_iv[0]) as in_file_h5:
            v = in_file_h5.root.IV_data[:]['voltage']
            i = in_file_h5.root.IV_data[:]['current']
            v[-1] = 800  # Fix last entry
            bias = get_voltage(np.abs(data[:, 0]), v, i)
    else:
        v = np.array([700, 900, 1000, 1100, 1300, 1400, 1500])
        i = np.array([26, 38, 41, 47, 70, 80, 95])
        bias = get_voltage(np.abs(data[:, 0]), v, i)
    mpv = data[:, 1] / unirrad_charge * 100.
    # Add Fit error
    mpv_low = np.abs(data[:, 3]) / unirrad_charge * 100.
    mpv_high = data[:, 2] / unirrad_charge * 100.
    # Add PlsrDAC calib error of data
    mpv_low = np.sqrt(mpv_low**2 + (data[:, 7] / unirrad_charge * 100.)**2)
    mpv_high = np.sqrt(mpv_high**2 + (data[:, 7] / unirrad_charge * 100.)**2)
    # Add PlsrDAC calib error of reference measurement
    mpv_low = np.sqrt(mpv_low**2 + (16. / unirrad_charge * 100.)**2)
    mpv_high = np.sqrt(mpv_high**2 + (16 / unirrad_charge * 100.)**2)

    spl = UnivariateSpline(bias[start_i:], mpv[start_i:],
                           w=np.max([mpv_low[start_i:],
                                     mpv_high[start_i:]],
                                    axis=0),
                           s=s)
#     plt.clf()
#     plt.errorbar(x=bias, y=mpv, yerr=[mpv_low, mpv_high],
#                  fmt='.', label='$1\cdot10^{15}\ \mathrm{N_{eq}/cm^2}$',
#                  linewidth=2, color='green', ecolor='black')
#     x = np.arange(bias[0], bias[-1], 1)
#     plt.plot(x, spl(x), color='red', alpha=0.5)
#
#     plt.plot(plt.xlim(), [unirrad_charge / unirrad_charge * 100.,
#                           unirrad_charge / unirrad_charge * 100.],
#              '-', label='Unirradiated', linewidth=2, color='black')
#     plt.show()
    return spl

def get_mesh(n_pixel, width, thickness, resolution):
    try:
        return tools.load(os.path.join(DATAFOLDER, 'mesh_%d_%d_%d_%d' % (n_pixel, int(width), int(thickness), int(resolution))))
    except IOError:
        mesh = geometry.mesh_planar_sensor(
                n_pixel=n_pixel,
                width=width,
                thickness=thickness,
                resolution=resolution)
        tools.save(mesh, os.path.join(DATAFOLDER, 'mesh_%d_%d_%d_%d' % (n_pixel, int(width), int(thickness), int(resolution))))
        return mesh

def get_potential(V_bias, n_eff, thickness):
    try:
        return tools.load(os.path.join(DATAFOLDER, 'pot_descr_%d_%d_%d' % (V_bias, n_eff, thickness)))
    except IOError:
        V_readout = 0.
        V_bi = -silicon.get_diffusion_potential(n_eff / 1e12, temperature=TEMP)

#         mesh = MESH_1 if thickness == THICKNESS_1 else MESH_2
        mesh = get_mesh(n_pixel=NPIXEL, width=WIDTH, thickness=thickness, resolution=RESOLUTION_1 if thickness == THICKNESS_1 else RESOLUTION_2)

        potential = fields.calculate_planar_sensor_potential(
            mesh=mesh,
            width=WIDTH,
            pitch=PITCH,
            n_pixel=NPIXEL,
            thickness=thickness,
            n_eff=n_eff,
            V_bias=V_bias,
            V_readout=V_readout,
            V_bi=V_bi)

        min_x = float(mesh.getFaceCenters()[0, :].min())
        max_x = float(mesh.getFaceCenters()[0, :].max())
        nx = WIDTH * NPIXEL
        ny = thickness

        pot_descr = fields.Description(potential,
                                       min_x=min_x,
                                       max_x=max_x,
                                       min_y=0,
                                       max_y=thickness,
                                       nx=nx,
                                       ny=ny,
                                       smoothing=SMOOTHING)

        tools.save(pot_descr, os.path.join(
            DATAFOLDER, 'pot_descr_%d_%d_%d' % (V_bias, n_eff, thickness)))

    return pot_descr

def get_w_potential(thickness):
    try:
        return tools.load(os.path.join(DATAFOLDER, 'pot_w_descr_%d' % thickness))
    except IOError:
        mesh = get_mesh(n_pixel=NPIXEL, width=WIDTH, thickness=thickness, resolution=RESOLUTION_1 if thickness == THICKNESS_1 else RESOLUTION_2)

        w_potential = fields.calculate_planar_sensor_w_potential(
                        mesh=mesh,
                        width=WIDTH,
                        pitch=PITCH,
                        n_pixel=NPIXEL,
                        thickness=thickness)

        pot_w_descr = fields.Description(w_potential,
                                     min_x=float(
                                         mesh.getFaceCenters()[0, :].min()),
                                     max_x=float(
                                         mesh.getFaceCenters()[0, :].max()),
                                     min_y=0,
                                     max_y=thickness,
                                     nx=WIDTH * NPIXEL,
                                     ny=thickness,
                                     smoothing=SMOOTHING)

        tools.save(pot_w_descr, os.path.join(DATAFOLDER, 'pot_w_descr_%d' % thickness))

    return pot_w_descr

def charge_collected(t_e_trapping, t_h_trapping,
                     t_e_t1, t_h_t1, thickness,
                     pot_descr):
    potw_descr = get_w_potential(thickness) 
    return analysis.get_charge_planar(WIDTH, thickness, pot_descr, pot_w_descr=potw_descr,
                                      t_e_trapping=t_e_trapping, t_h_trapping=t_h_trapping,
                                      t_e_t1=t_e_t1, t_h_t1=t_h_t1, t_r=None,
                                      grid_x=5, grid_y=25, n_pairs=10, dt=0.001,
                                      n_steps=20000, temperature=TEMP, multicore=False)

def ccs(t_trappings, t_t1, n_eff, thickness, biases):
    cc = []
    for v in biases:
        pot_descr=get_potential(V_bias=v,
                                n_eff=n_eff,
                                thickness=thickness)
        _, _, charge = charge_collected(t_e_trapping=t_trappings,
                                        t_h_trapping=t_trappings,
                                        t_e_t1=t_t1,
                                        t_h_t1=t_t1,
                                        thickness=thickness,
                                        pot_descr=pot_descr)
        cc.append(charge)
#         objgraph.show_growth()

    return np.array(cc)


def find_trappig():
    # Calculated unirradiated case
    _, _, cc_0_1 = charge_collected(t_e_trapping=300000,
                                    t_h_trapping=300000,
                                    t_e_t1=0.,
                                    t_h_t1=0.,
                                    thickness=THICKNESS_1,
                                    pot_descr=get_potential(V_bias=-80,
                                                            n_eff=NEFF0,
                                                            thickness=THICKNESS_1))
#     _, _, cc_0_2 = charge_collected(t_e_trapping=300000,
#                                     t_h_trapping=300000,
#                                     t_e_t1=0.,
#                                     t_h_t1=0.,
#                                     thickness=THICKNESS_2,
#                                     pot_descr=get_potential(V_bias=-80,
#                                                             n_eff=NEFF0,
#                                                             thickness=THICKNESS_2))

    with open(RESULT, "a+") as myfile:
        myfile.write('Biases 1 %s\n' % str(BIASES_1))
#         myfile.write('Data 1 %s\n' % str(SPLINE_1(np.abs(BIASES_1))))
#         myfile.write('Biases 12 %s\n' % str(BIASES_1_2))
        myfile.write('Data 12 %s\n' % str(SPLINE_1_2(np.abs(BIASES_1_2))))
#         myfile.write('Biases 2 %s\n' % str(BIASES_2))
#         myfile.write('Data 2 %s\n' % str(SPLINE_2(np.abs(BIASES_2))))
        myfile.write(
            'Fluence\tN_eff\tt_e\tt_h\tt1_e\tt1_h\tCCE12\tChi2\n')

    def minimize_me(args):
        tr = args[0]
        tr1 = args[1]
#         # Sensor 1 1e15
#         cc_1 = ccs(t_trappings=tr, t_t1=tr1, n_eff=NEFF,
#                    thickness=THICKNESS_1, biases=BIASES_1)
#         cce_1 = cc_1.sum(axis=(1, 2)) / cc_0_1.sum() * 100.
#         chi2 = chi2_to_spline(SPLINE_1, np.abs(BIASES_1), cce_1)
        # Sensor 1 7e15
        cc_1_2 = ccs(t_trappings=tr, t_t1=tr1, n_eff=NEFF2, thickness=THICKNESS_1, biases=BIASES_1_2)
        cce_1_2 = cc_1_2.sum(axis=(1, 2)) / cc_0_1.sum() * 100.
        chi2 = chi2_to_spline(SPLINE_1_2, np.abs(BIASES_1_2), cce_1_2)
#         # Sensor 2
#         cc_2 = ccs(t_trappings=tr, t_t1=tr1, n_eff=NEFF,
#                    thickness=THICKNESS_2, biases=BIASES_2)
#         cce_2 = cc_2.sum(axis=(1, 2)) / cc_0_2.sum() * 100.
#         chi2 = np.sqrt(
#             chi2**2 + chi2_to_spline(SPLINE_2, np.abs(BIASES_2), cce_2)**2)
        try:
            with open(RESULT, "a+") as myfile:
                myfile.write('%d\t%d\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%s\t%1.3f\n' % (
                    FLUENCE_1, NEFF, tr, tr, tr1, tr1, str(cce_1_2), chi2))
        except:
            print 'FAILED WRITING'
        print 'chi2', chi2
                
        return chi2

    if FIT:
        optimize.minimize(fun=minimize_me,
                          x0=[TR0, TR1],
                          args=(),
                          method='SLSQP',
                          jac=None,
                          bounds=[(2., 6.), (0., 0.5)],
                          constraints=(),
                          tol=None,
                          callback=None,
                          options={'disp': False,
                                   'iprint': 1,
                                   'eps': 0.01,
                                   'maxiter': 5,
                                   'ftol': 1e-04})
    else:
        minimize_me(args=[TR0, TR1])

#     optimize.basinhopping(func=minimize_me,
#                           x0=[5., 0.],
#                           niter=100,
#                           T=1.0,
#                           stepsize=0.5,
#                           minimizer_kwargs={'method': 'SLSQP',
#                                             'bounds': [(2., 6.), (0., 1.)],
#                                             'options': {'disp': False,
#                                                         'iprint': 1,
#                                                         'eps': 0.01,
#                                                         'maxiter': 100,
#                                                         'ftol': 1e-04}
#                                             },
#                           take_step=None,
#                           accept_test=None,
#                           callback=None,
#                           interval=50,
#                           disp=False,
#                           niter_success=None)


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    parser = argparse.ArgumentParser(
        description='Minimize to best trapping constants with start value')
    parser.add_argument(
        '--neff', type=float, nargs=1, help='Effective doping concentration [1e12 / cm2]', required=True)
    parser.add_argument(
        '--tr0', type=float, nargs=1, help='Start trapping parameter ns', required=True)
    parser.add_argument('--tr1', type=float, nargs=1,
                        help='Start trapping parameter per E-Field ns / V/um', required=True)
    parser.add_argument('--datafolder', type=str, nargs=1,
                        help='Folder with input and output data', required=True)
    parser.add_argument(
        '--output', type=str, nargs=1, help='Filename of output', required=True)
     
    parser.add_argument('-o', '--optimize', action='store_true')
  
    args = parser.parse_args()
    DATAFOLDER = args.datafolder[0]
    RESULT = args.output[0]
    TR0 = args.tr0[0]
    TR1 = args.tr1[0]
    NEFF = args.neff[0] * 1e12
    NEFF2 = NEFF
    FIT = args.optimize

#     TR0=4.
#     TR1=0.
#     DATAFOLDER=r'/media/documents/tmp/'
#     RESULT=r'/media/documents/tmp/result.txt'
#     NEFF = silicon.get_eff_acceptor_concentration(FLUENCE_1, NEFF0 / 1e12, is_ntype=True, is_oxygenated=True) * 1e12
#     NEFF2 = silicon.get_eff_acceptor_concentration(FLUENCE_2, NEFF0 / 1e12, is_ntype=True, is_oxygenated=True) * 1e12
    
    logging.info('Data folder: %s', DATAFOLDER)
    logging.info('Result file: %s', RESULT)
    logging.info('Neff: %1.2e', NEFF)
    logging.info('Neff2: %1.2e', NEFF2)
    logging.info('Start values: %1.2f ns, %1.2f ns*um/V', TR0, TR1)
    
    unirrad_charge = 17600.
#     data_files_1 = [
#         os.path.normpath(os.path.join(DATAFOLDER, 'mpv_bias_irrad1_scc30.txt'))]
#     data_files_iv_1 = [
#         os.path.normpath(os.path.join(DATAFOLDER, 'SCC_30_10e15_-30C.h5'))]
# 
#     SPLINE_1 = plot_cce_mpv(data_files_1, data_files_iv_1, start_i=7)

    data_files_1_2 = [
        os.path.normpath(os.path.join(DATAFOLDER, 'mpv_bias_irrad2_scc30.txt'))]
    data_files_iv_1_2 = None
    SPLINE_1_2 = plot_cce_mpv(
        data_files_1_2, data_files_iv_1_2, start_i=0, s=100)

#     unirrad_charge = 14570.
#     data_files_2 = [
#         os.path.normpath(os.path.join(DATAFOLDER, 'mpv_bias_irrad1_proto7.txt'))]
#     data_files_iv_2 = [
#         os.path.normpath(os.path.join(DATAFOLDER, r'Proto_7_10e15_-30C.h5'))]
# 
#     SPLINE_2 = plot_cce_mpv(data_files_2, data_files_iv_2, start_i=4)

    find_trappig()
