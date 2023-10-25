# Developer guide

## Third packages dependencies

MontuPython depends on different astronomical packages:

- `pymeeus`: routines of Meeus. Dependencies are:

    - Calculation of deltat-t: routine `pymeeus_tt2ut` in class `MonTime`
    - Calculation of date in mixed calendar: routine `pymeeus_Epoch` in class `MonTime`
    - Determine if date is julian: routine `is_julian` of `pymeeus_Epoch`.

- `astropy`: 
    - Class `Time` in class `MonTime`
    - 

- `pyephem`:
    - Class `Date` in class `MonTime`

- `pyplanets`:
    - Class `Epoch` in class `MonTime`
    - Routine `pyplanets_true_obliquity` in class `Montime`
    - Routine `pyplanets_nutation_longitude` in class `Montime``
    - Routine `apparent_sidereal_time` in class `Montime``

- `spiceypy`:
    - Routine `furnsh` in load_kernels.
    - Routine `pxform` in class `MonTime`
    