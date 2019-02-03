import VBET


params = {
    'network': '',
    'dem': '',
    'scratch': '',
    'lg_da': None,
    'med_da': None,
    'log_slope': None,
    'med_slope': None,
    'sm_slope': None,
    'lg_buf': None,
    'med_buf': None,
    'sm_buf': None,
    'min_buf': None,
    'dr_area': None,
    'lg_depth': None,
    'med_depth': None
}

vb = VBET.VBET(**params)
vb.add_da()
vb.add_elev()
vb.valley_bottom()