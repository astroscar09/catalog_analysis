import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from matplotlib.patches import Circle
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from astropy import units as u

import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.family'] = 'serif'
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13

zscale = ZScaleInterval()

class AstronomyImage:
    def __init__(self, fits_file, extension):
        
        self.fits_file = fits_file
        self.hdul = fits.open(fits_file)
        self.image_data = self.hdul[extension].data
        self.wcs = WCS(self.hdul[extension].header)
        
        
    
    def add_marker(self, ra, dec, marker='o', color='red', size=100):
        
        props = {'edgecolor':'white', 'facecolor':color, 
                 'marker': marker, 's':size, 'zorder':99}
        
        props['transform'] = self.ax.get_transform('fk5')
        self.ax.scatter(ra, dec, **props)
    
    def add_circular_region(self, ra, dec, radius, edgecolor='blue', facecolor='none', lw=2):
        coord = SkyCoord(ra, dec, unit='deg')
        
        spherical_props['transform'] =  self.ax.get_transform('fk5')
        circle = SphericalCircle((ra * u.degree, dec* u.degree), 
                                    radius * u.arcsecond.to(u.degree) * u.degree, **spherical_props)
        self.ax.add_patch(circle)
        
    def create_cutout(self, ra, dec, size):
        
        coord = SkyCoord(ra, dec, unit='deg')
        
        cutout = Cutout2D(self.image_data, position=coord, size=size*u.arcsec, wcs=self.wcs)
        self.cutout_data = cutout.data
        self.cutout_wcs = cutout.wcs
        self.fig, self.ax = plt.subplots(subplot_kw={'projection': self.cutout_wcs})
        #self.norm = simple_norm(self.cutout_data, 'sqrt')
        
        interval = zscale.get_limits(cutout.data)

        vmin = interval[0]
        vmax = interval[1]
        
        self.ax.imshow(self.cutout_data, vmin = vmin, vmax = vmax, origin='lower', cmap='gray_r')
        self.add_marker(ra, dec, marker = '*', color = 'orange')
        
        lon = self.ax.coords[0]
        lat = self.ax.coords[1]
        
        lon.set_axislabel('RA', fontsize = 15)
        lat.set_axislabel('DEC', fontsize = 15)
        
        lon.set_ticklabel(exclude_overlapping=True)
        lat.set_ticklabel(exclude_overlapping=True)

    def plot_image(self):
        
        self.fig, self.ax = plt.subplots(subplot_kw={'projection': self.wcs})
        
        interval = zscale.get_limits(self.image_data)

        vmin = interval[0]
        vmax = interval[1]
        
        self.ax.imshow(self.image_data, vmin = vmin, vmax = vmax, origin='lower', cmap='gray_r')

        lon = self.ax.coords[0]
        lat = self.ax.coords[1]
        
        lon.set_axislabel('RA', fontsize = 15)
        lat.set_axislabel('DEC', fontsize = 15)
        
        lon.set_ticklabel(exclude_overlapping=True)
        lat.set_ticklabel(exclude_overlapping=True)
    
    def show(self):
        plt.show()
        
    def get_ax(self):
        return self.ax
    
    def save(self, output_file):
        self.fig.savefig(output_file)