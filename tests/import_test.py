# -*- coding: utf-8 -*-
import pytest

def test_import():
    import csu_radartools
def test_import_fhc():
    from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain,
                            csu_dsd, csu_kdp, csu_misc, fundamentals)

