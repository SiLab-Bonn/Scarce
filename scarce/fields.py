import fipy
import numpy as np
import meshio as mio

from scarce import silicon
from scarce import geometry
from scarce import plot


def get_weighting_potential(x, y, D, S, W=None, is_planar=True):
    """ Planar sensor:
        From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D with:

        x [um] is the offset from the middle of the electrode
        y [um] the position in the sensor
        D [um] the sensor thickness
        S [um] the pixel pitch
        W [um] the electrode width

        3D sensor:
        Weighting potential for two cylinders with:
        D [um] distance between columns
        S [um] is the radius
        W number of readout columns
    """

    # Wheighting potential for one pixel
    if is_planar:
        if not W:  # Special case: 100% fill factor, amnalytic solution exists
            W = S

        xbar = np.pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D
        return -1. / np.pi * (np.arctan(np.tan(ybar / 2) * np.tanh((xbar + wbar / 2.) / 2.)) -
                              np.arctan(np.tan(ybar / 2) * np.tanh((xbar - wbar / 2.) / 2.)))
    else:
        R = S
        D = D / 2.  # D is the total distance between the columns
        a = np.sqrt(D * D - R * R)
        Phi_w = 1. / (4 * np.arccosh(D / R)) * \
            np.log(((x - a) ** 2 + y ** 2) / ((x + a) ** 2 + y ** 2)) + 0.5

        # Stability
        Phi_w = np.ma.masked_where(
            np.sqrt((x + D) * (x + D) + y * y) < R, Phi_w)
        Phi_w = np.ma.masked_where(
            np.sqrt((x - D) * (x - D) + y * y) < R, Phi_w)
        Phi_w = np.ma.masked_where(Phi_w < 0., Phi_w)
        Phi_w = np.ma.masked_where(Phi_w > 1., Phi_w)

        return Phi_w


def get_weighting_field(x, y, D, S, is_planar=True):
    """ From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D
        with x [um] is the position in the sensor [0:thickness], y [um] the offset from the
        middle of the electrode, D [um] the sensor thickness and S [um] the
        eletrode width. The field is calculated from the drivation of the
        potential in x and y.
    """

    if is_planar:
        xbar = np.pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D

        # Not easy to find a more simple form
        denom = (np.cosh(1. / 2. * (wbar - 2. * xbar)) + np.cos(ybar)) * \
            (np.cosh(1. / 2. * (wbar + 2. * xbar)) + np.cos(ybar)) * D

        E_x = - np.sin(ybar) * np.sinh(wbar / 2.) * np.sinh(xbar) / denom

        E_y = np.sinh(
            wbar / 2.) * (np.cosh(wbar / 2.) + np.cos(ybar) * np.cosh(xbar)) / denom

        return E_x, E_y
    else:
        # 3D sensor:
        # From the analytical derivation of the get_weighting_potential function
        # Weighting potential for two cylinders with:
        # S [um] is the radius
        # D [um] distance between columns

        R = S
        D = D / 2.
        a = np.sqrt(D * D - R * R)

        E_x = a / (np.arccosh(D / R)) * (a ** 2 - x ** 2 + y ** 2) / \
            (((a - x) ** 2 + y ** 2) * ((a + x) ** 2 + y ** 2))
        E_y = -2 * a / (np.arccosh(D / R)) * (x * y) / \
            (((a - x) ** 2 + y ** 2) * ((a + x) ** 2 + y ** 2))

        E_x = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, E_x)
        E_x = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, E_x)
        E_y = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, E_y)
        E_y = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, E_y)

        return -E_x, -E_y


def get_electric_field(x, y, V_bias, n_eff, D, S=None, is_planar=True):
    """ Calculates the 2D electric field E_x, E_y [V/um]

    Planar sensor:
        Calculates the field E_y[V/um], E_x = 0 in a planar sensor as a
        function of the position x between the electrodes [um],
        the bias voltage V_bias [V], the effective doping
        concentration n_eff [cm^-3] and the sensor Width D [um].
        The analytical function from the detector book p. 93 is used.

    3D sensor:
        Calculates the field E_x/E_y [V/um] in a 3d sensor as a function of the position
        x,y between the electrodes [um], the bias Voltage V_bias [V], the effective
        doping concentration n_eff [cm^-3], the electrode distance D [um] and radius R [um].
        So far the same field like the weighting field is used --> space charge is ignored.
    """

    if is_planar:
        if S:
            raise NotImplementedError(
                'The electrode width cannot be set, only full fill factor supported!')
        V_dep = silicon.get_depletion_voltage(n_eff, D)  # Depletion voltage
        a = (V_bias - V_dep) / D
        b = -2. * V_dep / (D ** 2)
        E_y = a - b * y
        return np.zeros_like(E_y), E_y
    else:
        E_x, E_y = get_weighting_field(x, y, D, S, is_planar=False)
        E_x = E_x * V_bias
        E_y = E_y * V_bias

        return E_x, E_y


def calculate_3D_sensor_potential(pitch_x, pitch_y, n_pixel_x, n_pixel_y, radius, resolution, V_readout, V_bias, nD=2):
    points, cells = geometry.mesh_3D_sensor(x=pitch_x,
                                    y=pitch_y,
                                    n_pixel_x=n_pixel_x, 
                                    n_pixel_y=n_pixel_y,
                                    radius=radius,
                                    nD=nD,
                                    resolution=resolution)
                                     
    mio.write('sensor.msh', points, cells)
    mesh = fipy.GmshImporter2D('sensor.msh')
    
    plot.plot_mesh(mesh)
    
#     potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
#     permittivity = 1.
#     potential.equation = (fipy.DiffusionTerm(coeff=permittivity) == 0.)
#     
#     bcs = []
#     allfaces = mesh.getExteriorFaces()
#     X,Y =  mesh.getFaceCenters()
#     
#     # Readout pillars
#     for pillar in range(nD):
#         position = pitch_x / nD * (pillar + 1. / 2.) - pitch_x / 2.
#         ring = allfaces & ( (X-position)**2+(Y)**2 < (radius)**2) 
#         bcs.append(fipy.FixedValue(value=V_readout,faces=ring))
#         
#     # Bias pillars
#     # Edges
#     positions = [(- pitch_x / 2., - pitch_y / 2.),
#                  (+ pitch_x / 2., - pitch_y / 2.),
#                  (+ pitch_x / 2., + pitch_y / 2.),
#                  (- pitch_x / 2., + pitch_y / 2.)]
#     # Sides
#     positions += [(0, - pitch_y / 2.),
#                  (0, + pitch_y / 2.)]
# 
#     for pos_x, pos_y in positions:
#         ring = allfaces & ( (X-pos_x)**2+(Y-pos_y)**2 < (radius)**2) 
#         bcs.append(fipy.FixedValue(value=V_bias, faces=ring))

#     # Calculate boundaries
#     p_pillars = mesh.getFaces()
#     n_pillars = mesh.getFacesTop()
# 
#     electrodes = readout_plane
#     bcs = [fipy.FixedValue(value=V_backplane, faces=backplane)]
#     
#     for pixel in range(n_pixel):
#         pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
#         bcs.append(fipy.FixedValue(value=V_readout,
#                                    faces=electrodes &
#                                    (X > pixel_position - pitch / 2.) &
#                                    (X < pixel_position + pitch / 2.)))

#     potential.equation.solve(var=potential, boundaryConditions=bcs)
#     return potential

if __name__ == '__main__':
    pitch_x = 250.
    pitch_y = 50.
    n_pixel_x, n_pixel_y = 1, 1
    radius = 6.
    resolution = 50.
    V_readout, V_bias,  = 0, -1
       
    potential = calculate_3D_sensor_potential(pitch_x, pitch_y, n_pixel_x, n_pixel_y, radius, resolution, V_readout, V_bias)
#     plot.plot_mesh(potential.mesh)
#     viewer = fipy.viewers.Viewer(vars=(potential, ))
#     viewer.plot("3D.png")
  
#     min_x, max_x = np.min(np.array(potential.mesh.getFaceCenters()[0])), np.max(np.array(potential.mesh.getFaceCenters()[0]))
#     min_y, max_y = np.min(np.array(potential.mesh.getFaceCenters()[1])), np.max(np.array(potential.mesh.getFaceCenters()[1]))
#      
#     print 'Interpolate'
#   
#     xnew = np.linspace(min_x, max_x, 1000)
#     ynew = np.linspace(min_y, max_y, 1000)
#     xnew_plot, ynew_plot = np.meshgrid(xnew, ynew)
#        
#     potential_function = interpolate_potential_2(potential)
#     print 'Done'
#        
#     plot.plot_3D_sensor(potential_function,
#                         pitch_x, 
#                         pitch_y, 
#                         n_pixel, 
#                         radius,
#                         V_bias,
#                         V_readout, 
#                         min_x, 
#                         max_x, 
#                         min_y,
#                         max_y
#                         )
