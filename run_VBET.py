import classVBET


class RunVBET:
    def __init__(self):
        self.params = {
            'network': '/path/to/stream/network.shp',
            'dem': '/path/to/dem.tif',
            'out': '/path/to/vbet/output.shp',
            'scratch': 'path/to/scratch/workspace',
            'lg_da': 250,
            'med_da': 25,
            'lg_slope': 3,
            'med_slope': 4,
            'sm_slope': 5,
            'lg_buf': 500,
            'med_buf': 200,
            'sm_buf': 80,
            'min_buf': 10,
            'dr_area': None,
            'da_field': 'TotDASqKm',
            'lg_depth': 3,
            'med_depth': 2,
            'sm_depth': 1.5
            }

    def run(self):
        vb = classVBET.VBET(**self.params)
        if self.params['da_field'] is None:
            vb.add_da()
        vb.valley_bottom()


# vbrun = RunVBET()
# vbrun.run()
