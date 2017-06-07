from hax.minitrees import TreeMaker
import pax.utils
from pax.InterpolatingMap import InterpolatingMap


class S2CorrectInHax(TreeMaker):
    """Extra information for analyses using the S2 just from the top or bottom PMT array.

    Provides:
    - s2_top: The uncorrected area in pe of the main interaction's S2 from the top array.
    - s2_bottom: The uncorrected area in pe of the main interaction's S2 from the bottom array.
    - cs2_hax: The corrected area in pe of the main interaction's S2
    - cs2_top: The corrected area in pe of the main interaction's S2 from the top array.
    - cs2_bottom: The corrected area in pe of the main interaction's S2 from the bottom array.

    Notes:
    - The S2 correction is applied at hax level using the map harcoded into this TreeMaker class
    - 'corrected' refers to applying all available position- and/or saturation corrections
      (see https://github.com/XENON1T/pax/blob/master/pax/plugins/interaction_processing/BuildInteractions.py#L105)
    """
    
    __version__ = '0.3.5'
    extra_branches = ['peaks.s2_saturation_correction',
                      'peaks.s2_spatial_correction',
                      'peaks.s2_bottom_spatial_correction',
                      'peaks.s2_top_spatial_correction',
                      'interactions.x',
                      'interactions.y',
                      'interactions.s2_lifetime_correction',
                      'peaks.area_fraction_top',
                      'peaks.area']

    # Which correction map to use
    filename = 's2_xy_map_v2.1.json'
    filename = pax.utils.data_file_name(filename)
    xymap = InterpolatingMap(filename)
    
    def extract_data(self, event):
        result = dict()

        if not len(event.interactions):
            return result

        interaction = event.interactions[0]
        s2 = event.peaks[interaction.s2]

        s2_top = s2.area * s2.area_fraction_top
        s2_bottom = s2.area * (1.0 - s2.area_fraction_top)
        s2_non_spatial_correction = s2.s2_saturation_correction * interaction.s2_lifetime_correction
        s2_total_pos_correction = 1.0 / S2CorrectInHax.xymap.get_value(interaction.x, interaction.y)
        s2_top_pos_correction = 1.0 / S2CorrectInHax.xymap.get_value(interaction.x, interaction.y, map='map_top')
        s2_bottom_pos_correction = 1.0 / S2CorrectInHax.xymap.get_value(interaction.x, interaction.y, map='map_bottom')
        s2_total_correction = s2_non_spatial_correction * s2_total_pos_correction
        s2_top_correction = s2_non_spatial_correction  * s2_top_pos_correction
        s2_bottom_correction = s2_non_spatial_correction * s2_bottom_pos_correction
        result['s2_top'] = s2_top
        result['s2_bottom'] = s2_bottom
        result['cs2_hax'] = s2.area * s2_total_correction
        result['cs2_top'] = s2_top * s2_top_correction
        result['cs2_bottom'] = s2_bottom * s2_bottom_correction
        result['s2_lifetime_correction'] = s2.s2_saturation_correction
        result['s2_position_correction'] = s2_total_pos_correction
        result['s2_correction_pax'] = s2.s2_spatial_correction
        result['s2_area_correction'] = interaction.s2_area_correction
        
        return result
