import numpy as np
import pylab as py
from astropy.io import fits

class PSF_grid(object):
    """
    Container for a grid of PSFs. This function hosts
    all of the interpolation routines such that you can "get_psf"
    anywhere in the focal plane.
    """


    def __init__(self, psf,
                 wave_array=[487, 625, 770, 870, 1020, 1250, 1650, 2120],
                 grid_shape=[11,11]):
        """
        Load up a FITS file that contains at least 1, or possibly a grid
        of PSFs. These are stored and you can interpolate between them
        if necessary.
        """
        self.img_scale = 0.10  # arcseconds per pixel
        
        # Fix wave_array to be a float
        wave_array = np.array(wave_array, dtype=float)
        wave_shape = psf.shape[0]

        if wave_shape != len(wave_array):
            print('Problem with PSF shape and wave_array shape')
            
        # Reshape the array to get the X and Y positions
        psf = psf.reshape((wave_shape, grid_shape[0], grid_shape[1],
                           psf.shape[2], psf.shape[3]))
        psf = np.swapaxes(psf, 1, 2)

        # scale array = lambda / 2D (Nyquist sampled)
        tel_diam = 2.235  # meters
        psf_scale = wave_array * (206264.8 * 1e-9) / (2.0 * tel_diam) # arcsec / pixel
        
        # Calculate the positions of all these PSFs. We assume that the
        # outermost PSFs are at the corners such that all observed stars
        # are internal to these corners.
        x_pos = np.mgrid[0:grid_shape[1]]
        y_pos = np.mgrid[0:grid_shape[0]]

        # Need to multiply by some of the array size properties.
        # Note that this assumes a pixel scale.
        fov = 10.    # arcmin
        fov *= 60.   # arcsec
        fov /= self.img_scale # pixels

        x_pos *= fov / x_pos[-1]
        y_pos *= fov / y_pos[-1]

        
        self.psf = psf
        self.psf_x = x_pos  # 1D array
        self.psf_y = y_pos  # 1D array
        self.psf_wave = wave_array
        self.wave_shape = wave_shape
        self.grid_shape = grid_shape
        self.psf_scale = psf_scale

        return

    @classmethod
    def from_file(cls, psf_file,
                  wave_array=[487, 625, 770, 870, 1020, 1250, 1650, 2120],
                  grid_shape=[11,11]):
        # 4D array with [wave, grid_idx, flux_x, flux_y]
        psf = fits.getdata(psf_file)

        return cls(psf, wave_array=wave_array, grid_shape=grid_shape)

        
    def get_local_psf(self, x, y, wave_idx, method='bilinear'):
        """
        Return an interpolated PSF at the requested [x, y] location.
        Interpolation method is nearest neighbor ('neighbor', default) or
        bilinear interpolation ('bilinear').
        """
        psf_x = self.psf_x
        psf_y = self.psf_y

        # We can't handle stars that are outside our PSF grid.
        # But fail gracefully with an exception. 
        if (x < psf_x[0]) or (x > psf_x[-1]):
            raise ValueError('x is outside the valid PSF grid region')
        if (y < psf_y[0]) or (y > psf_y[-1]):
            raise ValueError('y is outside the valid PSF grid region')
        
        # Find the PSF
        if method == 'neighbor':
            xidx = np.argmin(abs(psf_x - x))
            yidx = np.argmin(abs(psf_y - y))

            psf_loc = self.psf[wave_idx, yidx, xidx]

        elif method == 'bilinear':
            xidx_lo = np.where(psf_x <= x)[0][-1]
            yidx_lo = np.where(psf_y <= y)[0][-1]
            xidx_hi = xidx_lo + 1
            yidx_hi = yidx_lo + 1

            psf_xlo_ylo = self.psf[wave_idx, yidx_lo, xidx_lo]
            psf_xhi_ylo = self.psf[wave_idx, yidx_lo, xidx_hi]
            psf_xlo_yhi = self.psf[wave_idx, yidx_hi, xidx_lo]
            psf_xhi_yhi = self.psf[wave_idx, yidx_hi, xidx_hi]

            dx = 1. * (x - psf_x[xidx_lo]) / (psf_x[xidx_hi] - psf_x[xidx_lo])
            dy = 1. * (y - psf_y[yidx_lo]) / (psf_y[yidx_hi] - psf_y[yidx_lo])

            psf_loc = ((1 - dx) * (1 - dy) * psf_xlo_ylo +
                       (1 - dx) * (  dy  ) * psf_xlo_yhi +
                       (  dx  ) * (1 - dy) * psf_xhi_ylo +
                       (  dx  ) * (  dy  ) * psf_xhi_yhi)
            
        return psf_loc

    def plot_psf_grid(self, wave_idx, psf_size=[50,50]):
        # Chop down the PSFs to the plotting region
        # and at the wavelength requested.
        psf_shape = self.psf.shape[-2:]
        psf_x_lo = int((psf_shape[1] / 2.0) - (psf_size[1] / 2.0))
        psf_y_lo = int((psf_shape[0] / 2.0) - (psf_size[0] / 2.0))
        psf_x_hi = psf_x_lo + psf_size[1]
        psf_y_hi = psf_y_lo + psf_size[0]

        psf = self.psf[wave_idx, :, :, psf_y_lo:psf_y_hi, psf_x_lo:psf_x_hi]
        psf_shape = psf.shape[-2:]
        grid_shape = psf.shape[0:2]

        img = np.zeros((psf_shape[0] * grid_shape[0],
                        psf_shape[1] * grid_shape[1]), dtype=float)

        for xx in range(grid_shape[1]):
            for yy in range(grid_shape[0]):
                xlo = 0 + (xx * psf_shape[1])
                xhi = xlo + psf_shape[1]
                ylo = 0 + (yy * psf_shape[0])
                yhi = ylo + psf_shape[0]
                
                img[ylo:yhi, xlo:xhi] = psf[yy, xx, :, :]

        py.clf()
        py.imshow(img)

        return



    
