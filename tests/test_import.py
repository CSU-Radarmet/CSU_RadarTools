# -*- coding: utf-8 -*-
import pytest


def test_import_package():
    import csu_radartools


def test_import_csu_fhc():
    from csu_radartools import csu_fhc


def test_import_csu_liquid_ice_mass():
    from csu_radartools import csu_liquid_ice_mass


def test_import_csu_blended_rain():
    from csu_radartools import csu_blended_rain


def test_import_csu_blended_rain_tropical():
    from csu_radartools import csu_blended_rain_tropical


def test_import_csu_dsd():
    from csu_radartools import csu_dsd


def test_import_csu_kdp():
    from csu_radartools import csu_kdp


def test_import_csu_misc():
    from csu_radartools import csu_misc


def test_import_csu_t_traps():
    from csu_radartools import csu_t_traps


def test_import_fundamentals():
    from csu_radartools import fundamentals


def test_import_csu_fhc_melt():
    from csu_radartools import csu_fhc_melt

def test_import_csu_fhc_ml1():
    from csu_radartools import csu_fhc_ml1


def test_import_csu_fhc_winter():
    from csu_radartools import csu_fhc_winter


def test_import_common():
    from csu_radartools import common


def test_import_beta_functions():
    from csu_radartools import beta_functions


def test_import_version():
    from csu_radartools import _version


def test_import_calc_kdp_ray_fir():
    try:
        from csu_radartools import calc_kdp_ray_fir
    except ImportError:
        pytest.fail("Failed to import calc_kdp_ray_fir. Needs Compilation.")
