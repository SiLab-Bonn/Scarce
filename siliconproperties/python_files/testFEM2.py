import fipy

import pygmsh as pg
import meshio as mio


def generate(resolution):
    geom = pg.Geometry()

    # Two readout pillars
    circle = geom.add_circle(x0=[10., 0.0, 0.0],
                                   radius=5,
                                   lcar=resolution/2.,
                                   num_sections=4,
                                   # If compound==False, the section borders have to be points of the
                                   # discretization. If using a compound circle, they don't; gmsh can
                                   # choose by itself where to point the circle points.
                                   compound=False
                                   )
    circle_ll = geom.add_line_loop(circle)

    geom.add_rectangle(xmin=0,
                       xmax=20,
                       ymin=-10,
                       ymax=10,
                       z=0,
                       lcar=resolution,
                       holes=[circle_ll])

    return geom


if __name__ == '__main__':
    # Save mesh
    points, cells = pg.generate_mesh(generate(resolution=0.3))
    mio.write('test0.msh', points, cells)
#     import matplotlib.pyplot as plt
#     plt.plot(points[cells['triangle'].T, 0], points[cells['triangle'].T, 1], '-', color='black')
#     plt.plot(points[cells['line'].T, 0], points[cells['line'].T, 1], '-', color='black')
#     plt.savefig('test.svg')
    #plt.show()
    #raise
    mesh = fipy.GmshImporter2D('test0.msh')
    X,Y =  mesh.getFaceCenters()
    
    x0,y0 = 10, 0
    r = 5
    
    potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
    permittivity = 1.
    potential.equation = (fipy.DiffusionTerm(coeff = permittivity) == 0.)
    
    allfaces = mesh.getExteriorFaces()
    ring = allfaces & ( (X-x0)**2+(Y-y0)**2 < (r+1.5)**2) 
    
    bcs = (
        fipy.FixedValue(value=5.,faces=ring),
        fipy.FixedValue(value=0.,faces=mesh.getFacesLeft()),
    )
    
    potential.equation.solve(var=potential, boundaryConditions=bcs)
    
    viewer = fipy.viewers.Viewer(vars=(potential, ))
    viewer.plot("Example.png")
