import pytest
import montu

def test_time_zone_default():
    t = montu.Time('2024-05-01 12:00:00')
    assert t.readable.datepro.startswith('2024-05-01 12:00:00')

def test_time_zone_string_utc():
    # Example from user prompt:
    t = montu.Time('139-07-20 06:00:00', zone='UTC-5', calendar='mixed')
    # Without zone correction, time is 06:00:00.
    # With zone correction (UTC-5):
    # self - zonedt = self - (-5 * montu.HOUR) = self + 5 * montu.HOUR
    # So 06:00:00 + 5 hours = 11:00:00
    assert t.readable.datemix.startswith('139-07-20 11:00:00')

def test_time_zone_numeric_integer():
    t = montu.Time('2024-05-01 12:00:00', zone=-5)
    # 12:00:00 - (-5) = 17:00:00
    assert t.readable.datepro.startswith('2024-05-01 17:00:00')

def test_time_zone_numeric_float():
    t = montu.Time('2024-05-01 12:00:00', zone=5.5)
    # 12:00:00 - 5.5 = 06:30:00
    assert t.readable.datepro.startswith('2024-05-01 06:30:00')

def test_time_zone_string_colon():
    t = montu.Time('2024-05-01 12:00:00', zone='UTC+5:30')
    # 12:00:00 - 5.5 = 06:30:00
    assert t.readable.datepro.startswith('2024-05-01 06:30:00')

def test_time_zone_observer():
    # Create an observer with longitude = 75 degrees.
    # Timezone offset: lon/15 = 75/15 = 5.0
    obs = montu.Observer(lon=75.0, lat=0)
    t = montu.Time('2024-05-01 12:00:00', zone=obs)
    # 12:00:00 - 5 = 07:00:00
    assert t.readable.datepro.startswith('2024-05-01 07:00:00')
