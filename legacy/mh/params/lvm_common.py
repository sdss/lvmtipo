import numpy as np

r_outer = 44.5 / 2     # mm, outer radius (constrained by IFU/chip separation, defines the telescope FoV)
                       # taken from SDSS-V_0129_LVMi_PDR.pdf Table 7
image_scale = 8.92     # 1arcsec in microns from SDSS-V_0129_LVMi_PDR.pdf Table 7
pix_scale = 1.01       # arcsec per pixel
                    # taken from SDSS-V_0129_LVMi_PDR.pdf Table 13
a_telescope = np.pi * (16.2 / 2) ** 2

chip_size_pix = [1600, 1100]
chip_size_mm = [14.4, 10.2] # mm, guide chip height from SDSS-V_0129_LVMi_PDR.pdf Table 13

bias = 100
gain = 5
dark_current = 15
readout_noise = 5

inner_search_radius = 0 # degrees
outer_search_radius = 1000 * r_outer / image_scale / 3600

# guiding limits
mag_lim_lower = 17      # mag: 16.44 mag would be the limit, taken from LVM-0059 Table 8 (optimistic)
mag_lim_upper = 0.                # mag, currently no limit
min_neighbour_distance = -1.      #arcsec,

flux_of_vega = a_telescope * 1.6e6 # e-/sec/cm2 
# This is for the optimistic case, in the pessimistic case the number would be 8.71e5
# flux_of_vega = a_telescope * 8.71e5
zp = -2.5 * np.log10(flux_of_vega)
    
    
focus_stage_best_position = 0.0
