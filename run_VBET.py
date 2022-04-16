import classVBET
from datetime import datetime


class RunVBET:
    def __init__(self):
        self.params = {
            'network': '/home/jordan/Documents/Riverscapes/scratch/10190004_network.shp',
            'dem': '/home/jordan/Documents/Riverscapes/scratch/dem_10190004.tif',
            'out': '/home/jordan/Documents/Riverscapes/scratch/vbet2_out2.shp',
            'scratch': '/home/jordan/Documents/Riverscapes/scratch/scratch',
            'lg_da': 300,
            'med_da': 30,
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
        print('started: ', datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        vb = classVBET.VBET(**self.params)
        if self.params['da_field'] is None:
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


# vbrun = RunVBET()
# vbrun.run()
