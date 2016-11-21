''' Example that creates a 3D pixel array with a given geometry (number of pixels in x/y, height, width, and readout columns per pixel)
    and calculates the weighting potential and fields.
    
    .. NOTE::
       This example is not complete yet.
'''

from scarce import fields

if __name__ == '__main__':
    pitch_x = 250.
    pitch_y = 50.
    n_pixel_x, n_pixel_y = 3, 3
    radius = 6.
    resolution = 50.
    nD = 2  # Number of columns per pixel

    potential = fields.calculate_3D_sensor_w_potential(pitch_x, 
                                                     pitch_y, 
                                                     n_pixel_x, 
                                                     n_pixel_y, 
                                                     radius, 
                                                     resolution, 
                                                     nD=nD)
