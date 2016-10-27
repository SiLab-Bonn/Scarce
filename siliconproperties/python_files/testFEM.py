''' gmsh has to be installed and put into an environment variable (e.g. PATH)
    in a way that the command gmsh from the terminal/console starts it.
'''

import pygmsh as pg
import numpy as np
from _hotshot import resolution


def generate_3D_sensor_mesh(x=250., y=50., radius=6., nD=2, resolution=5.):
    geom = pg.Geometry()

    pillars = []

    # Create readout pillars
    for pillar in range(nD):
        position = x / nD * (pillar + 1. / 2.) - x/2.
#         print position
#         continue
        circle = geom.add_circle(x0=[position, 0.0, 0.0],
                                      radius=radius,
                                      lcar=resolution / 2.,
                                      num_sections=4,
                                      # If compound==False, the section borders have to be points of the
                                      # discretization. If using a compound circle, they don't; gmsh can
                                      # choose by itself where to point the circle points.
                                      compound=False
                                      )
        pillars.append(geom.add_line_loop(circle))

#     # Create bias pillars
#     
#     circle = geom.add_circle(x0=[-x/2. + radius, -y/2., 0.0],
#                                   radius=radius,
#                                   lcar=resolution / 2.,
#                                   num_sections=2,
#                                   # If compound==False, the section borders have to be points of the
#                                   # discretization. If using a compound circle, they don't; gmsh can
#                                   # choose by itself where to point the circle points.
#                                   compound=False
#                                   )
#     print circle
#     pillars.append(geom.add_line_loop(circle))
    
#     raise
#     circle_right = geom.add_circle(x0=[1.5, 0.0, 0.0],
#                                    radius=radius,
#                                    lcar=resolution,
#                                    num_sections=4,
#                                    # If compound==False, the section borders have to be points of the
#                                    # discretization. If using a compound circle, they don't; gmsh can
#                                    # choose by itself where to point the circle points.
#                                    compound=False
#                                    )
#     
#     circle_right_ll = geom.add_line_loop(circle_right)
    
#     circle_top_left = geom.add_circle(x0=[1.5, 0.0, 0.0],
#                                radius=0.5,
#                                lcar=resolution,
#                                num_sections=3,
#                                # If compound==False, the section borders have to be points of the
#                                # discretization. If using a compound circle, they don't; gmsh can
#                                # choose by itself where to point the circle points.
#                                compound=False)
#            
    print -y/2
    circle_points = []
    circle_points.append(geom.add_point([-x/4, -y/2, 0], lcar=resolution))
    circle_points.append(geom.add_point([0, y, 0], lcar=resolution))
    circle_points.append(geom.add_point([x/4, -y/2, 0], lcar=resolution))
    circle_half = geom.add_circle_sector(circle_points)
    
    t = [circle_points[0] + circle_points[2]]
    #circle_half_ll = geom.add_compound_line(circle_points)

    points = []
    
    points_xyz = [
        #[-x/4, -y/4, 0],
        #[x/4, -y/4, 0],
        [x/4, -y/2, 0],
        [x/2, -y/2, 0],
        [x/2, y/2, 0],
        [-x/2, y/2, 0],
        [-x/2, -y/2, 0],
        [-x/4, -y/2, 0]
        ]
    for point in points_xyz:
        points.append(geom.add_point(point, lcar=resolution))
        
    #points = circle_points + points
    
    # Create lines
    lines = [geom.add_line(points[i], points[i+1]) for i in range(len(points)-1)]
    #lines.append(geom.add_line(points[-1], points[0]))
    
    line_loop = geom.add_line_loop(lines)
    geom.add_plane_surface([line_loop])
    
#     geom.add_polygon_loop(X, lcar=resolution)
# 
#     geom.add_rectangle(xmin=-x/2,
#                        xmax=x/2,
#                        ymin=-y/2,
#                        ymax=y/2,
#                        z=0,
#                        lcar=resolution,
#                        holes=pillars)

    return geom


if __name__ == '__main__':
    points, cells = pg.generate_mesh(generate_3D_sensor_mesh())
    import matplotlib.pyplot as plt

    # Plot triangle cells
    plt.plot(points[cells['triangle'].T, 0], points[cells['triangle'].T, 1], '-', color='black')
    plt.plot(points[cells['line'].T, 0], points[cells['line'].T, 1], '-', color='black')
#     plt.xlim((-2.1, 2.1))
#     plt.ylim((-2.1, 2.1))
    plt.show()

    # Save mesh
    import meshio as mio
    mio.write('circle.msh', points, cells)
