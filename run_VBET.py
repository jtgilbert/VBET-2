import VBET


params = {
    'network': '/media/jordan/Seagate Backup Plus Drive/sidework/Rio_Grande/2mDEM_Projected/hrvbet/NHD_700m.shp',
    'dem': '/media/jordan/Seagate Backup Plus Drive/sidework/Rio_Grande/2mDEM_Projected/hrvbet/2mDEM.tif',
    'out': '/media/jordan/Seagate Backup Plus Drive/sidework/Rio_Grande/2mDEM_Projected/hrvbet/valley.shp',
    'scratch': '/media/jordan/Seagate Backup Plus Drive/sidework/Rio_Grande/2mDEM_Projected/hrvbet/scratch',
    'lg_da': 250,
    'med_da': 25,
    'lg_slope': 3,
    'med_slope': 6,
    'sm_slope': 10,
    'lg_buf': 1000,
    'med_buf': 250,
    'sm_buf': 10,
    'min_buf': 8,
    'dr_area': '/media/jordan/Seagate Backup Plus Drive/sidework/Rio_Grande/VBET Data/topo/Rio_Grande_DrArea.tif',
    'lg_depth': 2.5,
    'med_depth': 1.5
}

vb = VBET.VBET(**params)
vb.add_da()
vb.add_elev()
vb.valley_bottom()