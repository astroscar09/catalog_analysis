from astropy.table import Table
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from plotting_euclid_image_class import AstronomyImage

g_band = AstronomyImage('../../EUCLID_IMAGES/H20_EDFN_v2.3_HSC-G_bkgsub.fits', 1)

matplotlib.use('TkAgg')

euclid_cat = Table.read('edfn_dawn_catalog_sept2024.fits')


flux_cols = [    'cfht_u_flux_ujy',
                 'hsc_g_flux_ujy',
                 'hsc_r_flux_ujy',
                 'hsc_i_flux_ujy',
                 'hsc_z_flux_ujy',
                 'hsc_y_flux_ujy',
                 'nisp_y_flux_ujy',
                 'nisp_j_flux_ujy',
                 'nisp_h_flux_ujy',
                 'irac_ch1_flux_ujy',
                 'irac_ch2_flux_ujy',
                 'tf_vis_03_ujy',
                 'tf_vis_05_ujy',
                 'tf_vis_10_ujy',
                 'tf_vis_15_ujy']

flux_err_cols = [    'cfht_u_flux_ujy_err',
                     'hsc_g_flux_ujy_err',
                     'hsc_r_flux_ujy_err',
                     'hsc_i_flux_ujy_err',
                     'hsc_z_flux_ujy_err',
                     'hsc_y_flux_ujy_err',
                     'nisp_y_flux_ujy_err',
                     'nisp_j_flux_ujy_err',
                     'nisp_h_flux_ujy_err',
                     'irac_ch1_flux_ujy_err',
                     'irac_ch2_flux_ujy_err',
                     'te_vis_03_ujy',
                     'te_vis_05_ujy',
                     'te_vis_10_ujy',
                     'te_vis_15_ujy']



flux_only_cat = euclid_cat[flux_cols].to_pandas()
flux_only_cat = flux_only_cat.to_numpy()
flux_err_only_cat = euclid_cat[flux_err_cols].to_pandas()
flux_err_only_cat = flux_err_only_cat.to_numpy()
#flux_only_cat/flux_err_only_cat

error_boosts = np.array([1, 2, 4, 3, 2, 1, 3, 3, 3, 4, 2, 1, 1, 1, 1])

boosted_errors = (flux_err_only_cat * error_boosts)

snr = flux_only_cat/flux_err_only_cat
snr_boosted = flux_only_cat/boosted_errors

plot_cols = ['CFHT-u', 'HSC-g', 'HSC-r', 'HSC-i', 'HSC-z', 'HSC-y', 
             'NISP-y', 'NISP-j', 'NISP-h', 'IRAC Ch1', 'IRAC Ch2', 'VIS 03', 
             'VIS 05', 'VIS 10', 'VIS 15']

cols = np.arange(0, len(plot_cols), 1)

fig, ax = plt.subplots(3, 5, figsize = (20, 15), constrained_layout = True)
ax = ax.flatten()

for title, a, idx_col in zip(plot_cols, ax, cols):
    
    a.hist(snr[:, idx_col], range = (-100, 1000), bins = 100, log = True)
    a.set_title(title, fontsize = 15)
    a.set_ylabel('Counts', fontsize = 15)

fig.savefig('EUCLID_SNR_per_Band_No_ErrorBoost.png', dpi = 150)


fig, ax = plt.subplots(3, 5, figsize = (20, 15), constrained_layout = True)
ax = ax.flatten()

for title, a, idx_col in zip(plot_cols, ax, cols):
    
    a.hist(snr_boosted[:, idx_col], range = (-100, 1000), bins = 100, log = True)
    a.set_title(title, fontsize = 15)
    a.set_ylabel('Counts', fontsize = 15)

fig.savefig('EUCLID_SNR_per_Band_w_ErrorBoost.png', dpi = 150)


#Limiting magnitude for the bands
# band :: pc/area[deg] for 50, 80, and 95th
# cfht_u :: 26.6/10.6, 26.9/7.7, 27.1/2.7
# hsc_g :: 27.3/9.4, 27.6/5.3, 27.9/2.6
# hsc_r :: 27.3/9.6, 27.6/7.3, 27.7/3.5
# hsc_i :: 26.7/9.9, 27.0/6.8, 27.2/3.8
# hsc_z :: 26.1/9.4, 26.4/6.3, 26.6/2.3
# hsc_y :: 24.4/10.1, 24.9/5.3, 25.2/2.3
# nisp_y :: 26.0/11.6, 26.3/11.6, 26.5/10.4
# nisp_j :: 26.1/11.6, 26.4/11.6, 26.6/10.2
# nisp_h :: 26.0/11.6, 26.2/11.6, 26.4/10.6
# irac_ch1 :: 25.6/11.0, 25.9/10.9, 26.4/10.8
# irac_ch2 :: 24.9/9.9, 25.1/9.3, 25.4/6.5

def ABmag_to_Jy(M):
    
    return 10**((M - 8.9)/-2.5)

def ABmag_to_uJy(M):
    
    Jy = ABmag_to_Jy(M)
    uJy = Jy*1e6
    
    return uJy


def uJy_to_ABmag(uJy):

    Jy = uJy*1e-6
    
    return -2.5*np.log10(Jy) + 8.9




limit_u_Mag = np.array([26.6, 26.9, 27.1])
limit_g_Mag = np.array([27.3, 27.6, 27.9])
limit_r_Mag = np.array([27.3, 27.6, 27.7])
limit_i_Mag = np.array([26.7, 27.0, 27.2])
limit_z_Mag = np.array([26.1, 26.4, 26.6])
limit_y_Mag = np.array([24.4, 24.9, 25.2])
limit_nisp_y_Mag = np.array([26.0, 26.3, 26.5])
limit_nisp_j_Mag = np.array([ 26.1, 26.4, 26.6])
limit_nisp_h_Mag = np.array([26.0, 26.2, 26.4])
limit_ch1_Mag = np.array([25.6, 25.9, 26.4])
limit_ch2_Mag = np.array([24.9, 25.1, 25.4])

stacked_Mag_limits = (limit_u_Mag,
                       limit_g_Mag, 
                       limit_r_Mag,
                       limit_i_Mag,
                       limit_z_Mag,
                       limit_y_Mag, 
                       limit_nisp_y_Mag,
                       limit_nisp_j_Mag,
                       limit_nisp_h_Mag,
                       limit_ch1_Mag,
                       limit_ch2_Mag,
                       )
    
limit_u = ABmag_to_Jy(limit_u_Mag)
limit_g = ABmag_to_Jy(limit_g_Mag)
limit_r = ABmag_to_Jy(limit_r_Mag)
limit_i = ABmag_to_Jy(limit_i_Mag)
limit_z = ABmag_to_Jy(limit_z_Mag)
limit_y = ABmag_to_Jy(limit_y_Mag)
limit_nisp_y = ABmag_to_Jy(limit_nisp_y_Mag)
limit_nisp_j = ABmag_to_Jy(limit_nisp_j_Mag)
limit_nisp_h = ABmag_to_Jy(limit_nisp_h_Mag)
limit_ch1 = ABmag_to_Jy(limit_ch1_Mag)
limit_ch2 = ABmag_to_Jy(limit_ch2_Mag)

stacked_flux_limits = (limit_u,
                       limit_g, 
                       limit_r,
                       limit_i,
                       limit_z,
                       limit_y, 
                       limit_nisp_y,
                       limit_nisp_j,
                       limit_nisp_h,
                       limit_ch1,
                       limit_ch2,
                       )



def error_prop_mag(flux_err):
    
    a = 2.5
    sigma_a = flux_err
    A  = 1
    ln10 = np.log(10)

    err_mag = (a*sigma_a)/(A*ln10)

    return err_mag
    


flux_limits = np.vstack(stacked_flux_limits)

flim_df = pd.DataFrame(flux_limits, columns = ['flim_l16', 'flim_50', 'flim_u84'])
flim_df.to_csv('Euclid_flim.txt', sep = ' ')


eff_waves = np.loadtxt('filter_lambda_eff.dat')

filter_band_name = ['CFHT-u', 'HSC-g', 
                    'HSC-r', 'HSC-i', 
                    'HSC-z', 'HSC-Y', 'VIS', 'NISP Y', 
                    'NISP J', 'NISP H', 
                    'IRAC Ch1', 'IRAC Ch2']

map_filter_band = dict(zip(filter_band_name, eff_waves))

flux_cols_bands = ['CFHT-u', 'HSC-g', 
                    'HSC-r', 'HSC-i', 
                    'HSC-z', 'HSC-Y', 'NISP Y', 
                    'NISP J', 'NISP H', 
                    'IRAC Ch1', 'IRAC Ch2', 'VIS', 'VIS', 'VIS', 'VIS']

photom_waves = [map_filter_band[band] for band in flux_cols_bands]

def plot_SED_photom(photom_waves, photom, photom_err):

    fig, ax = plt.subplots(1, 1, figsize = (10, 5))
    ax.errorbar(photom_waves, photom, yerr = photom_err, fmt = 'o', color = 'black')
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Wavelength (Angstroms)', fontsize = 15)
    ax.set_ylabel('Flux (uJy)', fontsize = 15)
    ax.set_title('SED', fontsize = 15)
    ax.grid()

    return fig, ax

