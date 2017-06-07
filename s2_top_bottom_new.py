from hax.minitrees import TreeMaker
from pax.InterpolatingMap import InterpolatingMap
import numpy as np

class S2TopBottom(TreeMaker):
    """Extra information for analyses using the S2 just from the top or bottom PMT array.

    Provides:
    - s2_top: The uncorrected area in pe of the main interaction's S2 from the top array.
    - s2_bottom: The uncorrected area in pe of the main interaction's S2 from the bottom array.
    - cs2_top: The corrected area in pe of the main interaction's S2 from the top array.
    - cs2_bottom: The corrected area in pe of the main interaction's S2 from the bottom array.
    - cs2_new: The corrected area of the S2 using the new x, y correction map from 24th Feb 2017.
    - cs2_top_new: The corrected area of the S2 using the new x, y correction map from 24th Feb 2017.
    - cs2_bottom_new: The corrected area of the S2 using the new x, y correction map from 24th Feb 2017.
    - x_observed: The reconstructed observed interaction x coordinate, before the r, z correction.
    - y_observed: The reconstructed observed interaction y coordinate, before the r, z correction.
    - z_observed: The reconstructed observed interaction z coordinate, before the r, z correction.

    Notes:
    - 'corrected' refers to applying all available position- and/or saturation corrections
      (see https://github.com/XENON1T/pax/blob/master/pax/plugins/interaction_processing/BuildInteractions.py#L105)
    """
    __version__ = '0.4'
    extra_branches = ['peaks.s2_saturation_correction',
                      'peaks.s2_bottom_spatial_correction',
                      'peaks.s2_top_spatial_correction',
                      'interactions.s2_lifetime_correction',
                      'peaks.area_fraction_top',
                      'peaks.area',
                      'interactions.x',
                      'interactions.y',
                      'interactions.z',
                      'interactions.r_correction',
                      'interactions.z_correction']
    
    xymap_new = InterpolatingMap('/home/amb243/s2_xy_XENON1T_24Feb2017.json')

    def extract_data(self, event):
        result = dict()

        if not len(event.interactions):
            return result

        interaction = event.interactions[0]
        s2 = event.peaks[interaction.s2]
        
        r_corrected = np.sqrt(interaction.x ** 2 + interaction.y ** 2)
        r_observed = r_corrected - interaction.r_correction
        z_observed = interaction.z - interaction.z_correction
        x_observed = (r_observed / r_corrected) * interaction.x
        y_observed = (r_observed / r_corrected) * interaction.y
        
        s2_spatial_correction_new = 1.0 / S2TopBottom.xymap_new.get_value(x_observed, y_observed)
        s2_top_spatial_correction_new = 1.0 / S2TopBottom.xymap_new.get_value(x_observed, y_observed, map='map_top')
        s2_bottom_spatial_correction_new = 1.0 / S2TopBottom.xymap_new.get_value(x_observed, y_observed, map='map_bottom')

        s2_non_spatial_correction = interaction.s2_lifetime_correction
        s2_top_correction = s2_non_spatial_correction * s2.s2_top_spatial_correction 
        s2_bottom_correction = s2_non_spatial_correction * s2.s2_bottom_spatial_correction
        s2_correction_new = s2_non_spatial_correction * s2_spatial_correction_new
        s2_top_correction_new = s2_non_spatial_correction * s2_top_spatial_correction_new
        s2_bottom_correction_new = s2_non_spatial_correction * s2_bottom_spatial_correction_new

        s2_top = s2.area * s2.area_fraction_top
        s2_bottom = s2.area * (1.0 - s2.area_fraction_top)

        result['s2_top'] = s2_top
        result['s2_bottom'] = s2_bottom
        result['cs2_top'] = s2_top * s2_top_correction
        result['cs2_bottom'] = s2_bottom * s2_bottom_correction
        result['cs2_new'] = s2.area * s2_correction_new
        result['cs2_top_new'] = s2_top * s2_top_correction_new
        result['cs2_bottom_new'] = s2_bottom * s2_bottom_correction_new
        result['z_observed'] = z_observed
        result['x_observed'] = x_observed
        result['y_observed'] = y_observed

        return result
