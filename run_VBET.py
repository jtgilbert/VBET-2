import classVBET
from datetime import datetime


class RunVBET:
    def __init__(self):
        self.params = {
            'network': '/path/to/stream/network.shp',
            'dem': '/path/to/dem.tif',
            'out': '/path/to/save/output.shp',
            'scratch': '/path/of/scratch/workspace',
            'lg_da': 300,
            'med_da': 30,
            'lg_slope': 2,
            'med_slope': 3,
            'sm_slope': 4,
            'lg_buf': 500,
            'med_buf': 200,
            'sm_buf': 80,
            'min_buf': 10,
            'dr_area': '/path/to/drainage/area/raster.tif',
            'lg_depth': None,
            'med_depth': None,
            'sm_depth': None
            }

    def run(self):
        print('started: ', datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        vb = classVBET.VBET(**self.params)
        vb.add_da()
        vb.valley_bottom()
        print('ended: ', datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        metatxt = '{out}_metadata.txt'.format(out=self.params['out'])
        L = ['network: {} \n'.format(self.params['network']),
             'dem: {} \n'.format(self.params['dem']),
             'output: {} \n'.format(self.params['out']),
             'scratch workspace: {} \n'.format(self.params['scratch']),
             'large drainage area threshold: {} \n'.format(self.params['lg_da']),
             'medium drainage area threshold: {} \n'.format(self.params['med_da']),
             'large slope threshold: {} \n'.format(self.params['lg_slope']),
             'medium slope threshold: {} \n'.format(self.params['med_slope']),
             'small slope threshold: {} \n'.format(self.params['sm_slope']),
             'large buffer: {} \n'.format(self.params['lg_buf']),
             'medium buffer: {} \n'.format(self.params['med_buf']),
             'small buffer: {} \n'.format(self.params['sm_buf']),
             'minimum buffer: {} \n'.format(self.params['min_buf']),
             'large depth: {} \n'.format(self.params['lg_depth']),
             'medium depth: {} \n'.format(self.params['med_depth']),
             'small depth: {} \n'.format(self.params['sm_depth'])
             ]
        md = open(metatxt, 'w+')
        md.writelines(L)
        md.close()
