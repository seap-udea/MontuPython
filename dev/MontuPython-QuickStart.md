# Montu Python 
## Astronomical ephemerides for the Ancient World

## Quickstart

In this section you should provide the most simple instructions to use
your package.

You may import the package using:


```python
from montu import *
```

or for a safe import:


```python
import montu
%load_ext autoreload
%autoreload 2
```

    The autoreload extension is already loaded. To reload it, use:
      %reload_ext autoreload


It is important that before using the most interesting commands of the package, load relevant data:


```python
# Useful aliases
from montu.util import D2H, PRINTDF, TABLEDF
# Load kerneks
montu.Util.load_kernels(montu.PRECISION_KERNELS,dir='montmp/')
# Load stars
allstars = montu.Stars()
```

> **NOTE**: 
 > - *Precision kernels* are binary files used by the NASA NAIF SPICE information system. 
 > - These files are used by `MontuPython` to calculate precise position of the planets. 
 > - Some of the required kernels are large files (in the order of Giga bytes). You need enough space to download them and patience the first time. 
 > - Once you have download them, you may save the files in a safe place and load from them `montu.Util.load_kernels(dir='my_safe_directory')`.

Two very basic codes that do something non-trivial albeit simple calculations in `MontuPython are:`

1. Compute the position of mars at a given time and while observing from a given site on Earth:


```python
mtime = montu.Time('-2500-01-01 12:00:00.00')
tebas = montu.Observer(lon=33,lat=24,planet=montu.Planet('Earth'))
mars = montu.Planet('Mars')

# Get the position (do some magic! - 'heka' in ancient egyptian)
position = montu.Heka.where_in_sky(mars,at=mtime,site=tebas)

print(f"Position in the celestial sphere (right ascension, declination) at JED {float(position.epoch_jed)}:",
      D2H(position.RAEpoch),D2H(position.DecEpoch))
print(f"Position above the horizon (azimuth, elevation) at JED {float(position.epoch_jed)}:",
      D2H(position.az),D2H(position.el))
```

    Position in the celestial sphere (right ascension, declination) at JED 807954.0: 08:32:9.796 24:06:28.555
    Position above the horizon (azimuth, elevation) at JED 807954.0: 06:11:24.275 -41:38:31.094


2. Obtain the information about a star from the stellar catalogue and, as in the case of the planet, obtain the position of the star in the sky.


```python
mtime = montu.Time('-2500-01-01 12:00:00.00')
tebas = montu.Observer(lon=33,lat=24,planet=montu.Planet('Earth'))
aldebaran = allstars.get_stars(ProperName='Aldebaran')
position = montu.Heka.where_in_sky(aldebaran,at=mtime,site=tebas)
print(f"Position in the celestial sphere (right ascension, declination) at JED {float(position.epoch_jed)}:",
      D2H(position.RAEpoch),D2H(position.DecEpoch))
print(f"Position above the horizon (azimuth, elevation) at JED {float(position.epoch_jed)}:",
      D2H(position.az),D2H(position.el))
```

    Position in the celestial sphere (right ascension, declination) at JED 807954.0: 00:36:36.821 -2:15:10.754
    Position above the horizon (azimuth, elevation) at JED 807954.0: 107:35:52.939 29:35:33.617


### Working with stars

The stellar catalogue included with `MontuPython` contains more than 119 000 stars, including almost 9 000 visible to the naked eye.  Information about stars is stored in a `pandas Data Frame` whose columns are:


```python
allstars.data.columns
```




    Index(['MN', 'HD', 'HR', 'HIP', 'Gl', 'Name', 'OtherDesignations',
           'ProperName', 'Bayer', 'Flamsteed', 'Constellation', 'RAJ2000',
           'DecJ2000', 'pmRA', 'pmDec', 'RadVel', 'Distance', 'Vmag', 'Vmag_min',
           'Vmag_max', 'B-V', 'SpType', 'Luminosity', 'XJ2000', 'YJ2000', 'ZJ2000',
           'VXJ2000', 'VYJ2000', 'VZJ2000', 'Primary', 'MultipleID', 'IsMultiple',
           'IsVariable'],
          dtype='object')



Altough you may manipulate this DataFrame using the conventional commands in pandas, we have designed several useful methods to obtain subsets of the catalogue. For instance if you want to extract the stars visible to naked eye, the command would be:


```python
stars = allstars.get_stars(Vmag=[-2,6.5])
print(f"There is {stars.number} visible to the naked eye in the catalogue")
```

    There is 8920 visible to the naked eye in the catalogue


You can use any of the properties of the stars to filter them. A common filter is to look for single stars:


```python
aldebaran = stars.get_stars(ProperName='Aldebaran')
TABLEDF(aldebaran.data)
```

    |    |   MN |    HD |   HR |   HIP | Gl        | Name      | OtherDesignations                                                                     | ProperName   | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |   B-V | SpType   |   Luminosity |   XJ2000 |   YJ2000 |   ZJ2000 |   VXJ2000 |   VYJ2000 |   VZJ2000 |   Primary | MultipleID   |   IsMultiple |   IsVariable |
    |----|------|-------|------|-------|-----------|-----------|---------------------------------------------------------------------------------------|--------------|---------|-------------|-----------------|-----------|------------|--------|---------|----------|------------|--------|------------|------------|-------|----------|--------------|----------|----------|----------|-----------|-----------|-----------|-----------|--------------|--------------|--------------|
    | 14 |   15 | 29139 | 1457 | 21421 | Gl 171.1A | Aldebaran | 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau | Aldebaran    | α Tau   | 87 Tau      | Tau             |   4.59868 |    16.5093 |  62.78 | -189.36 |     54.5 |    20.4332 |   0.87 |      0.888 |      0.858 | 1.538 | K5III    |       163.23 |  7.02722 |  18.2876 |  5.80666 | -2.14e-06 | 5.709e-05 | 1.528e-05 |     21368 | Gl 171.1     |            1 |            1 |


Another useful method included with the class `Stars` is that of filtering the getting the stars close to a given point in the sky. For illustrare, below is the command to obtain all stars in the sky with magnitudes less than 5 and that are at 5.5 degrees or less than Aldebaran:


```python
hyades = stars.get_stars_around(center=[aldebaran.data.RAJ2000,aldebaran.data.DecJ2000],radius=5.5,Vmag=[-1,5])
TABLEDF(hyades.data)
```

    |      |   MN |    HD |   HR |   HIP | Gl        | Name           | OtherDesignations                                                                     | ProperName     | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |   B-V | SpType   |   Luminosity |   XJ2000 |   YJ2000 |   ZJ2000 |    VXJ2000 |    VYJ2000 |    VZJ2000 |   Primary | MultipleID   |   IsMultiple |   IsVariable |
    |------|------|-------|------|-------|-----------|----------------|---------------------------------------------------------------------------------------|----------------|---------|-------------|-----------------|-----------|------------|--------|---------|----------|------------|--------|------------|------------|-------|----------|--------------|----------|----------|----------|------------|------------|------------|-----------|--------------|--------------|--------------|
    |   14 |   15 | 29139 | 1457 | 21421 | Gl 171.1A | Aldebaran      | 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau | Aldebaran      | α Tau   | 87 Tau      | Tau             |   4.59868 |    16.5093 |  62.78 | -189.36 |     54.5 |    20.4332 |   0.87 |      0.888 |      0.858 | 1.538 | K5III    |    163.23    |  7.02722 |  18.2876 |  5.80666 | -2.14e-06  |  5.709e-05 |  1.528e-05 |     21368 | Gl 171.1     |            1 |            1 |
    |  259 |  260 | 28319 | 1412 | 20894 | nan       | Chamukuy       | 78 Tau/78The2Tau/Chamukuy/HD 28319/HIP 20894/HR 1412/HYG 20842/MN 260/θ2 Tau          | Chamukuy       | θ2 Tau  | 78 Tau      | Tau             |   4.47771 |    15.8709 | 108.66 |  -26.39 |     40   |    46.1042 |   3.4  |      3.411 |      3.391 | 0.179 | A7III    |     80.8351  | 17.2097  |  40.8716 | 12.6082  |  5.51e-06  |  4.718e-05 | -6.49e-06  |     20842 | nan          |            0 |            1 |
    |  301 |  302 | 28305 | 1409 | 20889 | nan       | Ain            | 74 Tau/74Eps Tau/Ain/HD 28305/HIP 20889/HR 1409/HYG 20837/MN 302/ε Tau                | Ain            | ε Tau   | 74 Tau      | Tau             |   4.47694 |    19.1804 | 107.23 |  -36.77 |     39   |    44.964  |   3.53 |    nan     |    nan     | 1.014 | K0III    |     68.1711  | 16.4884  |  39.1368 | 14.7728  |  5.53e-06  |  4.622e-05 | -5.89e-06  |     20837 | nan          |            0 |            0 |
    |  344 |  345 | 27371 | 1346 | 20205 | nan       | Prima Hyadum   | 54 Tau/54Gam Tau/HD 27371/HIP 20205/HR 1346/HYG 20155/MN 345/Prima Hyadum/γ Tau       | Prima Hyadum   | γ Tau   | 54 Tau      | Tau             |   4.32989 |    15.6276 | 115.29 |  -23.86 |     39   |    49.5295 |   3.65 |    nan     |    nan     | 0.981 | G8III    |     74.0628  | 20.1974  |  43.2117 | 13.3426  |  5.23e-06  |  4.792e-05 | -8.16e-06  |     20155 | nan          |            0 |            0 |
    |  394 |  395 | 27697 | 1373 | 20455 | nan       | Secunda Hyadum | 61 Tau/61Del1Tau/HD 27697/HIP 20455/HR 1373/HYG 20405/MN 395/Secunda Hyadum/δ1 Tau    | Secunda Hyadum | δ1 Tau  | 61 Tau      | Tau             |   4.38225 |    17.5425 | 107.75 |  -28.84 |     39   |    47.7099 |   3.77 |    nan     |    nan     | 0.983 | G8III    |     61.546   | 18.696   |  41.472  | 14.3805  |  5.66e-06  |  4.675e-05 | -6.26e-06  |     20405 | nan          |            0 |            0 |
    |  429 |  430 | 28307 | 1411 | 20885 | nan       | θ1 Tau         | 77 Tau/77The1Tau/HD 28307/HIP 20885/HR 1411/HYG 20833/MN 430/θ1 Tau                   | nan            | θ1 Tau  | 77 Tau      | Tau             |   4.47625 |    15.9622 | 104.76 |  -15.01 |     40   |    47.3261 |   3.84 |    nan     |    nan     | 0.952 | G7III    |     56.8068  | 17.6738  |  41.929  | 13.0149  |  7.94e-06  |  4.645e-05 | -6.5e-06   |     20833 | nan          |            0 |            0 |
    |  556 |  557 | 31421 | 1580 | 22957 | nan       | ο2 Ori         | 9 Ori/9Omi2Ori/HD 31421/HIP 22957/HR 1580/HYG 22904/MN 557/ο2 Ori                     | nan            | ο2 Ori  | 9 Ori       | Ori             |   4.93952 |    13.5145 | -77.77 |  -45.97 |      1   |    57.0125 |   4.06 |    nan     |    nan     | 1.158 | K2III    |     67.2977  | 15.1933  |  53.3112 | 13.3233  | -1.211e-05 | -2.08e-06  |  2.176e-05 |     22904 | nan          |            0 |            0 |
    |  693 |  694 | 29388 | 1473 | 21589 | nan       | 90 Tau         | 90 Tau/HD 29388/HIP 21589/HR 1473/HYG 21536/MN 694                                    | nan            | nan     | 90 Tau      | Tau             |   4.63596 |    12.5108 | 101.73 |  -14.9  |     45   |    47.081  |   4.27 |    nan     |    nan     | 0.122 | A6V      |     37.8094  | 16.0671  |  43.0638 | 10.199   |  6.65e-06  |  5.09e-05  | -5.79e-06  |     21536 | nan          |            0 |            0 |
    |  721 |  722 | 27962 | 1389 | 20648 | nan       | δ3 Tau         | 68 Tau/68Del3Tau/HD 27962/HIP 20648/HR 1389/HYG 20597/MN 722/δ3 Tau                   | nan            | δ3 Tau  | 68 Tau      | Tau             |   4.42483 |    17.9279 | 108.26 |  -32.47 |     35   |    45.5373 |   4.3  |      4.309 |      4.299 | 0.049 | A2IV     |     34.4191  | 17.3648  |  39.6945 | 14.0174  |  4.2e-06   |  4.28e-05  | -7.36e-06  |     20597 | nan          |            0 |            1 |
    |  894 |  895 | 28052 | 1394 | 20713 | nan       | 71 Tau         | 71 Tau/HD 28052/HIP 20713/HR 1394/HYG 20661/MN 895                                    | nan            | nan     | 71 Tau      | Tau             |   4.43909 |    15.6183 | 114.66 |  -33.3  |     38   |    49.0918 |   4.48 |      4.488 |      4.468 | 0.262 | F0V...   |     33.8844  | 18.7872  |  43.3865 | 13.2169  |  2.83e-06  |  4.715e-05 | -9.32e-06  |     20661 | nan          |            0 |            1 |
    | 1072 | 1073 | 28910 | 1444 | 21273 | nan       | ρ Tau          | 86 Tau/86Rho Tau/HD 28910/HIP 21273/HR 1444/HYG 21220/MN 1073/ρ Tau                   | nan            | ρ Tau   | 86 Tau      | Tau             |   4.56414 |    14.8444 | 103.69 |  -25.94 |     40   |    48.5201 |   4.65 |      4.662 |      4.642 | 0.255 | A8V      |     28.3139  | 17.2181  |  43.6262 | 12.4307  |  4.58e-06  |  4.719e-05 | -7.59e-06  |     21220 | nan          |            0 |            1 |
    | 1100 | 1101 | 29488 | 1479 | 21683 | nan       | σ2 Tau         | 92 Tau/92Sig2Tau/HD 29488/HIP 21683/HR 1479/HYG 21630/MN 1101/σ2 Tau                  | nan            | σ2 Tau  | 92 Tau      | Tau             |   4.65458 |    15.918  |  82.4  |  -19.53 |     36   |    47.6872 |   4.67 |    nan     |    nan     | 0.147 | A5Vn     |     26.8411  | 15.8209  |  43.0435 | 13.0788  |  5.76e-06  |  4.097e-05 | -5.24e-06  |     21630 | nan          |            0 |            0 |
    | 1123 | 1124 | 28100 | 1396 | 20732 | nan       | π Tau          | 73 Tau/73Pi Tau/HD 28100/HIP 20732/HR 1396/HYG 20680/MN 1124/π Tau                    | nan            | π Tau   | 73 Tau      | Tau             |   4.44344 |    14.7138 |  -7.57 |  -31.15 |     32   |   127.714  |   4.69 |    nan     |    nan     | 0.979 | G8III    |    188.973   | 48.9561  | 113.411  | 32.4382  | -1.034e-05 |  3.17e-05  |  1.879e-05 |     20680 | nan          |            0 |            0 |
    | 1144 | 1145 | 30959 | 1556 | 22667 | nan       | ο1 Ori         | 4 Ori/4Omi1Ori/HD 30959/HIP 22667/HR 1556/HYG 22614/MN 1145/ο1 Ori                    | nan            | ο1 Ori  | 4 Ori       | Ori             |   4.87554 |    14.2506 |  -2.62 |  -56.13 |     -8   |   199.601  |   4.71 |      4.758 |      4.668 | 1.773 | M3Sv     |    453.315   | 56.1317  | 185.136  | 49.1345  | -5.465e-05 |  4.47e-06  |  4e-06     |     22614 | nan          |            0 |            1 |
    | 1242 | 1243 | 28527 | 1427 | 21029 | Gl 170.1  | HD 28527       | Gl 170.1/HD 28527/HIP 21029/HR 1427/HYG 20977/MN 1243                                 | nan            | nan     | nan         | Tau             |   4.50934 |    16.194  | 104.98 |  -25.14 |     41   |    43.1965 |   4.78 |    nan     |    nan     | 0.17  | A6IV     |     19.8976  | 15.781   |  38.364  | 12.0472  |  6.64e-06  |  4.696e-05 | -4.45e-06  |     20977 | nan          |            0 |            0 |
    | 1275 | 1276 | 27819 | 1380 | 20542 | nan       | δ2 Tau         | 64 Tau/64Del2Tau/HD 27819/HIP 20542/HR 1380/HYG 20490/MN 1276/δ2 Tau                  | nan            | δ2 Tau  | 64 Tau      | Tau             |   4.4016  |    17.4441 | 109.99 |  -33.47 |     39   |    49.4805 |   4.8  |    nan     |    nan     | 0.154 | A7V      |     25.633   | 19.182   |  43.1321 | 14.8332  |  4.3e-06   |  4.769e-05 | -7.67e-06  |     20490 | nan          |            0 |            0 |
    | 1485 | 1486 | 27045 | 1329 | 19990 | nan       | ω2 Tau         | 50 Tau/50Ome2Tau/HD 27045/HIP 19990/HR 1329/HYG 19940/MN 1486/ω2 Tau                  | nan            | ω2 Tau  | 50 Tau      | Tau             |   4.28768 |    20.5786 | -40.88 |  -61.45 |     16   |    28.9436 |   4.93 |    nan     |    nan     | 0.259 | A3m      |      7.78395 | 11.7443  |  24.4195 | 10.1735  | -2.32e-06  |  1.405e-05 |  1.312e-05 |     19940 | nan          |            0 |            0 |
    | 1544 | 1545 | 28292 | 1407 | 20877 | nan       | 75 Tau         | 75 Tau/HD 28292/HIP 20877/HR 1407/HYG 20825/MN 1545                                   | nan            | nan     | 75 Tau      | Tau             |   4.47399 |    16.3597 |   8.12 |   17.79 |     18   |    57.241  |   4.96 |    nan     |    nan     | 1.137 | K2IIIvar |     29.621   | 21.3634  |  50.5985 | 16.1229  |  9.92e-06  |  1.587e-05 |  4.25e-06  |     20825 | nan          |            0 |            0 |


We can map the stars:


```python
fig,axs = hyades.plot_stars()
fig.savefig('gallery/hyades.png')
```


    
![png](MontuPython-QuickStart_files/MontuPython-QuickStart_24_0.png)
    


Now you can precess the stars:


```python
mtime = montu.Time('-2500-01-01 12:00:00.00')
montu.Heka.precess_to_epoch(hyades,at=mtime,inplace=True)
TABLEDF(hyades.data)
```

    |      |   MN |    HD |   HR |   HIP | Gl        | Name           | OtherDesignations                                                                     | ProperName     | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |   B-V | SpType   |   Luminosity |   XJ2000 |   YJ2000 |   ZJ2000 |    VXJ2000 |    VYJ2000 |    VZJ2000 |   Primary | MultipleID   |   IsMultiple |   IsVariable |   RAJ2000t |   DecJ2000t |   RAEpoch |   DecEpoch |     epoch_tt |   epoch_jed |
    |------|------|-------|------|-------|-----------|----------------|---------------------------------------------------------------------------------------|----------------|---------|-------------|-----------------|-----------|------------|--------|---------|----------|------------|--------|------------|------------|-------|----------|--------------|----------|----------|----------|------------|------------|------------|-----------|--------------|--------------|--------------|------------|-------------|-----------|------------|--------------|-------------|
    |   14 |   15 | 29139 | 1457 | 21421 | Gl 171.1A | Aldebaran      | 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau | Aldebaran      | α Tau   | 87 Tau      | Tau             |   4.59868 |    16.5093 |  62.78 | -189.36 |     54.5 |    20.4332 |   0.87 |      0.888 |      0.858 | 1.538 | K5III    |    163.23    |  7.02722 |  18.2876 |  5.80666 | -2.14e-06  |  5.709e-05 |  1.528e-05 |     21368 | Gl 171.1     |            1 |            1 |    4.59345 |     16.746  |  0.610228 |  -2.25299  | -1.42006e+11 |      807954 |
    |  259 |  260 | 28319 | 1412 | 20894 | nan       | Chamukuy       | 78 Tau/78The2Tau/Chamukuy/HD 28319/HIP 20894/HR 1412/HYG 20842/MN 260/θ2 Tau          | Chamukuy       | θ2 Tau  | 78 Tau      | Tau             |   4.47771 |    15.8709 | 108.66 |  -26.39 |     40   |    46.1042 |   3.4  |      3.411 |      3.391 | 0.179 | A7III    |     80.8351  | 17.2097  |  40.8716 | 12.6082  |  5.51e-06  |  4.718e-05 | -6.49e-06  |     20842 | nan          |            0 |            1 |    4.46865 |     15.9039 |  0.50886  |  -3.52984  | -1.42006e+11 |      807954 |
    |  301 |  302 | 28305 | 1409 | 20889 | nan       | Ain            | 74 Tau/74Eps Tau/Ain/HD 28305/HIP 20889/HR 1409/HYG 20837/MN 302/ε Tau                | Ain            | ε Tau   | 74 Tau      | Tau             |   4.47694 |    19.1804 | 107.23 |  -36.77 |     39   |    44.964  |   3.53 |    nan     |    nan     | 1.014 | K0III    |     68.1711  | 16.4884  |  39.1368 | 14.7728  |  5.53e-06  |  4.622e-05 | -5.89e-06  |     20837 | nan          |            0 |            0 |    4.46801 |     19.2264 |  0.452316 |  -0.31723  | -1.42006e+11 |      807954 |
    |  344 |  345 | 27371 | 1346 | 20205 | nan       | Prima Hyadum   | 54 Tau/54Gam Tau/HD 27371/HIP 20205/HR 1346/HYG 20155/MN 345/Prima Hyadum/γ Tau       | Prima Hyadum   | γ Tau   | 54 Tau      | Tau             |   4.32989 |    15.6276 | 115.29 |  -23.86 |     39   |    49.5295 |   3.65 |    nan     |    nan     | 0.981 | G8III    |     74.0628  | 20.1974  |  43.2117 | 13.3426  |  5.23e-06  |  4.792e-05 | -8.16e-06  |     20155 | nan          |            0 |            0 |    4.32028 |     15.6575 |  0.374225 |  -4.29655  | -1.42006e+11 |      807954 |
    |  394 |  395 | 27697 | 1373 | 20455 | nan       | Secunda Hyadum | 61 Tau/61Del1Tau/HD 27697/HIP 20455/HR 1373/HYG 20405/MN 395/Secunda Hyadum/δ1 Tau    | Secunda Hyadum | δ1 Tau  | 61 Tau      | Tau             |   4.38225 |    17.5425 | 107.75 |  -28.84 |     39   |    47.7099 |   3.77 |    nan     |    nan     | 0.983 | G8III    |     61.546   | 18.696   |  41.472  | 14.3805  |  5.66e-06  |  4.675e-05 | -6.26e-06  |     20405 | nan          |            0 |            0 |    4.37327 |     17.5786 |  0.392493 |  -2.2482   | -1.42006e+11 |      807954 |
    |  429 |  430 | 28307 | 1411 | 20885 | nan       | θ1 Tau         | 77 Tau/77The1Tau/HD 28307/HIP 20885/HR 1411/HYG 20833/MN 430/θ1 Tau                   | nan            | θ1 Tau  | 77 Tau      | Tau             |   4.47625 |    15.9622 | 104.76 |  -15.01 |     40   |    47.3261 |   3.84 |    nan     |    nan     | 0.952 | G7III    |     56.8068  | 17.6738  |  41.929  | 13.0149  |  7.94e-06  |  4.645e-05 | -6.5e-06   |     20833 | nan          |            0 |            0 |    4.46752 |     15.9809 |  0.506504 |  -3.45939  | -1.42006e+11 |      807954 |
    |  556 |  557 | 31421 | 1580 | 22957 | nan       | ο2 Ori         | 9 Ori/9Omi2Ori/HD 31421/HIP 22957/HR 1580/HYG 22904/MN 557/ο2 Ori                     | nan            | ο2 Ori  | 9 Ori       | Ori             |   4.93952 |    13.5145 | -77.77 |  -45.97 |      1   |    57.0125 |   4.06 |    nan     |    nan     | 1.158 | K2III    |     67.2977  | 15.1933  |  53.3112 | 13.3233  | -1.211e-05 | -2.08e-06  |  2.176e-05 |     22904 | nan          |            0 |            0 |    4.946   |     13.5719 |  0.996366 |  -3.887    | -1.42006e+11 |      807954 |
    |  693 |  694 | 29388 | 1473 | 21589 | nan       | 90 Tau         | 90 Tau/HD 29388/HIP 21589/HR 1473/HYG 21536/MN 694                                    | nan            | nan     | 90 Tau      | Tau             |   4.63596 |    12.5108 | 101.73 |  -14.9  |     45   |    47.081  |   4.27 |    nan     |    nan     | 0.122 | A6V      |     37.8094  | 16.0671  |  43.0638 | 10.199   |  6.65e-06  |  5.09e-05  | -5.79e-06  |     21536 | nan          |            0 |            0 |    4.62748 |     12.5295 |  0.71675  |  -6.18807  | -1.42006e+11 |      807954 |
    |  721 |  722 | 27962 | 1389 | 20648 | nan       | δ3 Tau         | 68 Tau/68Del3Tau/HD 27962/HIP 20648/HR 1389/HYG 20597/MN 722/δ3 Tau                   | nan            | δ3 Tau  | 68 Tau      | Tau             |   4.42483 |    17.9279 | 108.26 |  -32.47 |     35   |    45.5373 |   4.3  |      4.309 |      4.299 | 0.049 | A2IV     |     34.4191  | 17.3648  |  39.6945 | 14.0174  |  4.2e-06   |  4.28e-05  | -7.36e-06  |     20597 | nan          |            0 |            1 |    4.41581 |     17.9685 |  0.425378 |  -1.72072  | -1.42006e+11 |      807954 |
    |  894 |  895 | 28052 | 1394 | 20713 | nan       | 71 Tau         | 71 Tau/HD 28052/HIP 20713/HR 1394/HYG 20661/MN 895                                    | nan            | nan     | 71 Tau      | Tau             |   4.43909 |    15.6183 | 114.66 |  -33.3  |     38   |    49.0918 |   4.48 |      4.488 |      4.468 | 0.262 | F0V...   |     33.8844  | 18.7872  |  43.3865 | 13.2169  |  2.83e-06  |  4.715e-05 | -9.32e-06  |     20661 | nan          |            0 |            1 |    4.42954 |     15.6599 |  0.476447 |  -3.90778  | -1.42006e+11 |      807954 |
    | 1072 | 1073 | 28910 | 1444 | 21273 | nan       | ρ Tau          | 86 Tau/86Rho Tau/HD 28910/HIP 21273/HR 1444/HYG 21220/MN 1073/ρ Tau                   | nan            | ρ Tau   | 86 Tau      | Tau             |   4.56414 |    14.8444 | 103.69 |  -25.94 |     40   |    48.5201 |   4.65 |      4.662 |      4.642 | 0.255 | A8V      |     28.3139  | 17.2181  |  43.6262 | 12.4307  |  4.58e-06  |  4.719e-05 | -7.59e-06  |     21220 | nan          |            0 |            1 |    4.5555  |     14.8768 |  0.607583 |  -4.20029  | -1.42006e+11 |      807954 |
    | 1100 | 1101 | 29488 | 1479 | 21683 | nan       | σ2 Tau         | 92 Tau/92Sig2Tau/HD 29488/HIP 21683/HR 1479/HYG 21630/MN 1101/σ2 Tau                  | nan            | σ2 Tau  | 92 Tau      | Tau             |   4.65458 |    15.918  |  82.4  |  -19.53 |     36   |    47.6872 |   4.67 |    nan     |    nan     | 0.147 | A5Vn     |     26.8411  | 15.8209  |  43.0435 | 13.0788  |  5.76e-06  |  4.097e-05 | -5.24e-06  |     21630 | nan          |            0 |            0 |    4.64772 |     15.9424 |  0.674712 |  -2.82064  | -1.42006e+11 |      807954 |
    | 1123 | 1124 | 28100 | 1396 | 20732 | nan       | π Tau          | 73 Tau/73Pi Tau/HD 28100/HIP 20732/HR 1396/HYG 20680/MN 1124/π Tau                    | nan            | π Tau   | 73 Tau      | Tau             |   4.44344 |    14.7138 |  -7.57 |  -31.15 |     32   |   127.714  |   4.69 |    nan     |    nan     | 0.979 | G8III    |    188.973   | 48.9561  | 113.411  | 32.4382  | -1.034e-05 |  3.17e-05  |  1.879e-05 |     20680 | nan          |            0 |            0 |    4.44407 |     14.7527 |  0.505239 |  -4.73345  | -1.42006e+11 |      807954 |
    | 1144 | 1145 | 30959 | 1556 | 22667 | nan       | ο1 Ori         | 4 Ori/4Omi1Ori/HD 30959/HIP 22667/HR 1556/HYG 22614/MN 1145/ο1 Ori                    | nan            | ο1 Ori  | 4 Ori       | Ori             |   4.87554 |    14.2506 |  -2.62 |  -56.13 |     -8   |   199.601  |   4.71 |      4.758 |      4.668 | 1.773 | M3Sv     |    453.315   | 56.1317  | 185.136  | 49.1345  | -5.465e-05 |  4.47e-06  |  4e-06     |     22614 | nan          |            0 |            1 |    4.87576 |     14.3208 |  0.916527 |  -3.46506  | -1.42006e+11 |      807954 |
    | 1242 | 1243 | 28527 | 1427 | 21029 | Gl 170.1  | HD 28527       | Gl 170.1/HD 28527/HIP 21029/HR 1427/HYG 20977/MN 1243                                 | nan            | nan     | nan         | Tau             |   4.50934 |    16.194  | 104.98 |  -25.14 |     41   |    43.1965 |   4.78 |    nan     |    nan     | 0.17  | A6IV     |     19.8976  | 15.781   |  38.364  | 12.0472  |  6.64e-06  |  4.696e-05 | -4.45e-06  |     20977 | nan          |            0 |            0 |    4.50059 |     16.2254 |  0.533144 |  -3.10185  | -1.42006e+11 |      807954 |
    | 1275 | 1276 | 27819 | 1380 | 20542 | nan       | δ2 Tau         | 64 Tau/64Del2Tau/HD 27819/HIP 20542/HR 1380/HYG 20490/MN 1276/δ2 Tau                  | nan            | δ2 Tau  | 64 Tau      | Tau             |   4.4016  |    17.4441 | 109.99 |  -33.47 |     39   |    49.4805 |   4.8  |    nan     |    nan     | 0.154 | A7V      |     25.633   | 19.182   |  43.1321 | 14.8332  |  4.3e-06   |  4.769e-05 | -7.67e-06  |     20490 | nan          |            0 |            0 |    4.39243 |     17.486  |  0.411738 |  -2.27077  | -1.42006e+11 |      807954 |
    | 1485 | 1486 | 27045 | 1329 | 19990 | nan       | ω2 Tau         | 50 Tau/50Ome2Tau/HD 27045/HIP 19990/HR 1329/HYG 19940/MN 1486/ω2 Tau                  | nan            | ω2 Tau  | 50 Tau      | Tau             |   4.28768 |    20.5786 | -40.88 |  -61.45 |     16   |    28.9436 |   4.93 |    nan     |    nan     | 0.259 | A3m      |      7.78395 | 11.7443  |  24.4195 | 10.1735  | -2.32e-06  |  1.405e-05 |  1.312e-05 |     19940 | nan          |            0 |            0 |    4.29109 |     20.6554 |  0.267759 |   0.458064 | -1.42006e+11 |      807954 |
    | 1544 | 1545 | 28292 | 1407 | 20877 | nan       | 75 Tau         | 75 Tau/HD 28292/HIP 20877/HR 1407/HYG 20825/MN 1545                                   | nan            | nan     | 75 Tau      | Tau             |   4.47399 |    16.3597 |   8.12 |   17.79 |     18   |    57.241  |   4.96 |    nan     |    nan     | 1.137 | K2IIIvar |     29.621   | 21.3634  |  50.5985 | 16.1229  |  9.92e-06  |  1.587e-05 |  4.25e-06  |     20825 | nan          |            0 |            0 |    4.47332 |     16.3374 |  0.505882 |  -3.09337  | -1.42006e+11 |      807954 |


Or compute their horizontal positions:


```python
montu.Heka.where_in_sky(hyades,at=mtime,site=tebas,inplace=True)
TABLEDF(hyades.data)
```

    |      |   MN |    HD |   HR |   HIP | Gl        | Name           | OtherDesignations                                                                     | ProperName     | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |   B-V | SpType   |   Luminosity |   XJ2000 |   YJ2000 |   ZJ2000 |    VXJ2000 |    VYJ2000 |    VZJ2000 |   Primary | MultipleID   |   IsMultiple |   IsVariable |   RAJ2000t |   DecJ2000t |   RAEpoch |   DecEpoch |     epoch_tt |   epoch_jed |      HA |      az |      el |     zen |
    |------|------|-------|------|-------|-----------|----------------|---------------------------------------------------------------------------------------|----------------|---------|-------------|-----------------|-----------|------------|--------|---------|----------|------------|--------|------------|------------|-------|----------|--------------|----------|----------|----------|------------|------------|------------|-----------|--------------|--------------|--------------|------------|-------------|-----------|------------|--------------|-------------|---------|---------|---------|---------|
    |   14 |   15 | 29139 | 1457 | 21421 | Gl 171.1A | Aldebaran      | 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau | Aldebaran      | α Tau   | 87 Tau      | Tau             |   4.59868 |    16.5093 |  62.78 | -189.36 |     54.5 |    20.4332 |   0.87 |      0.888 |      0.858 | 1.538 | K5III    |    163.23    |  7.02722 |  18.2876 |  5.80666 | -2.14e-06  |  5.709e-05 |  1.528e-05 |     21368 | Gl 171.1     |            1 |            1 |    4.59345 |     16.746  |  0.610228 |  -2.25299  | -1.42006e+11 |      807954 | 20.2635 | 107.598 | 29.5927 | 60.4073 |
    |  259 |  260 | 28319 | 1412 | 20894 | nan       | Chamukuy       | 78 Tau/78The2Tau/Chamukuy/HD 28319/HIP 20894/HR 1412/HYG 20842/MN 260/θ2 Tau          | Chamukuy       | θ2 Tau  | 78 Tau      | Tau             |   4.47771 |    15.8709 | 108.66 |  -26.39 |     40   |    46.1042 |   3.4  |      3.411 |      3.391 | 0.179 | A7III    |     80.8351  | 17.2097  |  40.8716 | 12.6082  |  5.51e-06  |  4.718e-05 | -6.49e-06  |     20842 | nan          |            0 |            1 |    4.46865 |     15.9039 |  0.50886  |  -3.52984  | -1.42006e+11 |      807954 | 20.3648 | 109.749 | 30.2713 | 59.7287 |
    |  301 |  302 | 28305 | 1409 | 20889 | nan       | Ain            | 74 Tau/74Eps Tau/Ain/HD 28305/HIP 20889/HR 1409/HYG 20837/MN 302/ε Tau                | Ain            | ε Tau   | 74 Tau      | Tau             |   4.47694 |    19.1804 | 107.23 |  -36.77 |     39   |    44.964  |   3.53 |    nan     |    nan     | 1.014 | K0III    |     68.1711  | 16.4884  |  39.1368 | 14.7728  |  5.53e-06  |  4.622e-05 | -5.89e-06  |     20837 | nan          |            0 |            0 |    4.46801 |     19.2264 |  0.452316 |  -0.31723  | -1.42006e+11 |      807954 | 20.4214 | 106.977 | 32.6045 | 57.3955 |
    |  344 |  345 | 27371 | 1346 | 20205 | nan       | Prima Hyadum   | 54 Tau/54Gam Tau/HD 27371/HIP 20205/HR 1346/HYG 20155/MN 345/Prima Hyadum/γ Tau       | Prima Hyadum   | γ Tau   | 54 Tau      | Tau             |   4.32989 |    15.6276 | 115.29 |  -23.86 |     39   |    49.5295 |   3.65 |    nan     |    nan     | 0.981 | G8III    |     74.0628  | 20.1974  |  43.2117 | 13.3426  |  5.23e-06  |  4.792e-05 | -8.16e-06  |     20155 | nan          |            0 |            0 |    4.32028 |     15.6575 |  0.374225 |  -4.29655  | -1.42006e+11 |      807954 | 20.4995 | 111.728 | 31.6008 | 58.3992 |
    |  394 |  395 | 27697 | 1373 | 20455 | nan       | Secunda Hyadum | 61 Tau/61Del1Tau/HD 27697/HIP 20455/HR 1373/HYG 20405/MN 395/Secunda Hyadum/δ1 Tau    | Secunda Hyadum | δ1 Tau  | 61 Tau      | Tau             |   4.38225 |    17.5425 | 107.75 |  -28.84 |     39   |    47.7099 |   3.77 |    nan     |    nan     | 0.983 | G8III    |     61.546   | 18.696   |  41.472  | 14.3805  |  5.66e-06  |  4.675e-05 | -6.26e-06  |     20405 | nan          |            0 |            0 |    4.37327 |     17.5786 |  0.392493 |  -2.2482   | -1.42006e+11 |      807954 | 20.4812 | 109.492 | 32.4238 | 57.5762 |
    |  429 |  430 | 28307 | 1411 | 20885 | nan       | θ1 Tau         | 77 Tau/77The1Tau/HD 28307/HIP 20885/HR 1411/HYG 20833/MN 430/θ1 Tau                   | nan            | θ1 Tau  | 77 Tau      | Tau             |   4.47625 |    15.9622 | 104.76 |  -15.01 |     40   |    47.3261 |   3.84 |    nan     |    nan     | 0.952 | G7III    |     56.8068  | 17.6738  |  41.929  | 13.0149  |  7.94e-06  |  4.645e-05 | -6.5e-06   |     20833 | nan          |            0 |            0 |    4.46752 |     15.9809 |  0.506504 |  -3.45939  | -1.42006e+11 |      807954 | 20.3672 | 109.7   | 30.3374 | 59.6626 |
    |  556 |  557 | 31421 | 1580 | 22957 | nan       | ο2 Ori         | 9 Ori/9Omi2Ori/HD 31421/HIP 22957/HR 1580/HYG 22904/MN 557/ο2 Ori                     | nan            | ο2 Ori  | 9 Ori       | Ori             |   4.93952 |    13.5145 | -77.77 |  -45.97 |      1   |    57.0125 |   4.06 |    nan     |    nan     | 1.158 | K2III    |     67.2977  | 15.1933  |  53.3112 | 13.3233  | -1.211e-05 | -2.08e-06  |  2.176e-05 |     22904 | nan          |            0 |            0 |    4.946   |     13.5719 |  0.996366 |  -3.887    | -1.42006e+11 |      807954 | 19.8773 | 106.073 | 23.739  | 66.261  |
    |  693 |  694 | 29388 | 1473 | 21589 | nan       | 90 Tau         | 90 Tau/HD 29388/HIP 21589/HR 1473/HYG 21536/MN 694                                    | nan            | nan     | 90 Tau      | Tau             |   4.63596 |    12.5108 | 101.73 |  -14.9  |     45   |    47.081  |   4.27 |    nan     |    nan     | 0.122 | A6V      |     37.8094  | 16.0671  |  43.0638 | 10.199   |  6.65e-06  |  5.09e-05  | -5.79e-06  |     21536 | nan          |            0 |            0 |    4.62748 |     12.5295 |  0.71675  |  -6.18807  | -1.42006e+11 |      807954 | 20.157  | 110.552 | 26.244  | 63.756  |
    |  721 |  722 | 27962 | 1389 | 20648 | nan       | δ3 Tau         | 68 Tau/68Del3Tau/HD 27962/HIP 20648/HR 1389/HYG 20597/MN 722/δ3 Tau                   | nan            | δ3 Tau  | 68 Tau      | Tau             |   4.42483 |    17.9279 | 108.26 |  -32.47 |     35   |    45.5373 |   4.3  |      4.309 |      4.299 | 0.049 | A2IV     |     34.4191  | 17.3648  |  39.6945 | 14.0174  |  4.2e-06   |  4.28e-05  | -7.36e-06  |     20597 | nan          |            0 |            1 |    4.41581 |     17.9685 |  0.425378 |  -1.72072  | -1.42006e+11 |      807954 | 20.4483 | 108.659 | 32.2636 | 57.7364 |
    |  894 |  895 | 28052 | 1394 | 20713 | nan       | 71 Tau         | 71 Tau/HD 28052/HIP 20713/HR 1394/HYG 20661/MN 895                                    | nan            | nan     | 71 Tau      | Tau             |   4.43909 |    15.6183 | 114.66 |  -33.3  |     38   |    49.0918 |   4.48 |      4.488 |      4.468 | 0.262 | F0V...   |     33.8844  | 18.7872  |  43.3865 | 13.2169  |  2.83e-06  |  4.715e-05 | -9.32e-06  |     20661 | nan          |            0 |            1 |    4.42954 |     15.6599 |  0.476447 |  -3.90778  | -1.42006e+11 |      807954 | 20.3973 | 110.413 | 30.4955 | 59.5045 |
    | 1072 | 1073 | 28910 | 1444 | 21273 | nan       | ρ Tau          | 86 Tau/86Rho Tau/HD 28910/HIP 21273/HR 1444/HYG 21220/MN 1073/ρ Tau                   | nan            | ρ Tau   | 86 Tau      | Tau             |   4.56414 |    14.8444 | 103.69 |  -25.94 |     40   |    48.5201 |   4.65 |      4.662 |      4.642 | 0.255 | A8V      |     28.3139  | 17.2181  |  43.6262 | 12.4307  |  4.58e-06  |  4.719e-05 | -7.59e-06  |     21220 | nan          |            0 |            1 |    4.5555  |     14.8768 |  0.607583 |  -4.20029  | -1.42006e+11 |      807954 | 20.2661 | 109.554 | 28.6578 | 61.3422 |
    | 1100 | 1101 | 29488 | 1479 | 21683 | nan       | σ2 Tau         | 92 Tau/92Sig2Tau/HD 29488/HIP 21683/HR 1479/HYG 21630/MN 1101/σ2 Tau                  | nan            | σ2 Tau  | 92 Tau      | Tau             |   4.65458 |    15.918  |  82.4  |  -19.53 |     36   |    47.6872 |   4.67 |    nan     |    nan     | 0.147 | A5Vn     |     26.8411  | 15.8209  |  43.0435 | 13.0788  |  5.76e-06  |  4.097e-05 | -5.24e-06  |     21630 | nan          |            0 |            0 |    4.64772 |     15.9424 |  0.674712 |  -2.82064  | -1.42006e+11 |      807954 | 20.199  | 107.622 | 28.4722 | 61.5278 |
    | 1123 | 1124 | 28100 | 1396 | 20732 | nan       | π Tau          | 73 Tau/73Pi Tau/HD 28100/HIP 20732/HR 1396/HYG 20680/MN 1124/π Tau                    | nan            | π Tau   | 73 Tau      | Tau             |   4.44344 |    14.7138 |  -7.57 |  -31.15 |     32   |   127.714  |   4.69 |    nan     |    nan     | 0.979 | G8III    |    188.973   | 48.9561  | 113.411  | 32.4382  | -1.034e-05 |  3.17e-05  |  1.879e-05 |     20680 | nan          |            0 |            0 |    4.44407 |     14.7527 |  0.505239 |  -4.73345  | -1.42006e+11 |      807954 | 20.3685 | 110.975 | 29.701  | 60.299  |
    | 1144 | 1145 | 30959 | 1556 | 22667 | nan       | ο1 Ori         | 4 Ori/4Omi1Ori/HD 30959/HIP 22667/HR 1556/HYG 22614/MN 1145/ο1 Ori                    | nan            | ο1 Ori  | 4 Ori       | Ori             |   4.87554 |    14.2506 |  -2.62 |  -56.13 |     -8   |   199.601  |   4.71 |      4.758 |      4.668 | 1.773 | M3Sv     |    453.315   | 56.1317  | 185.136  | 49.1345  | -5.465e-05 |  4.47e-06  |  4e-06     |     22614 | nan          |            0 |            1 |    4.87576 |     14.3208 |  0.916527 |  -3.46506  | -1.42006e+11 |      807954 | 19.9572 | 106.291 | 24.9907 | 65.0093 |
    | 1242 | 1243 | 28527 | 1427 | 21029 | Gl 170.1  | HD 28527       | Gl 170.1/HD 28527/HIP 21029/HR 1427/HYG 20977/MN 1243                                 | nan            | nan     | nan         | Tau             |   4.50934 |    16.194  | 104.98 |  -25.14 |     41   |    43.1965 |   4.78 |    nan     |    nan     | 0.17  | A6IV     |     19.8976  | 15.781   |  38.364  | 12.0472  |  6.64e-06  |  4.696e-05 | -4.45e-06  |     20977 | nan          |            0 |            0 |    4.50059 |     16.2254 |  0.533144 |  -3.10185  | -1.42006e+11 |      807954 | 20.3406 | 109.109 | 30.1737 | 59.8263 |
    | 1275 | 1276 | 27819 | 1380 | 20542 | nan       | δ2 Tau         | 64 Tau/64Del2Tau/HD 27819/HIP 20542/HR 1380/HYG 20490/MN 1276/δ2 Tau                  | nan            | δ2 Tau  | 64 Tau      | Tau             |   4.4016  |    17.4441 | 109.99 |  -33.47 |     39   |    49.4805 |   4.8  |    nan     |    nan     | 0.154 | A7V      |     25.633   | 19.182   |  43.1321 | 14.8332  |  4.3e-06   |  4.769e-05 | -7.67e-06  |     20490 | nan          |            0 |            0 |    4.39243 |     17.486  |  0.411738 |  -2.27077  | -1.42006e+11 |      807954 | 20.462  | 109.343 | 32.1637 | 57.8363 |
    | 1485 | 1486 | 27045 | 1329 | 19990 | nan       | ω2 Tau         | 50 Tau/50Ome2Tau/HD 27045/HIP 19990/HR 1329/HYG 19940/MN 1486/ω2 Tau                  | nan            | ω2 Tau  | 50 Tau      | Tau             |   4.28768 |    20.5786 | -40.88 |  -61.45 |     16   |    28.9436 |   4.93 |    nan     |    nan     | 0.259 | A3m      |      7.78395 | 11.7443  |  24.4195 | 10.1735  | -2.32e-06  |  1.405e-05 |  1.312e-05 |     19940 | nan          |            0 |            0 |    4.29109 |     20.6554 |  0.267759 |   0.458064 | -1.42006e+11 |      807954 | 20.6059 | 107.797 | 35.3978 | 54.6022 |
    | 1544 | 1545 | 28292 | 1407 | 20877 | nan       | 75 Tau         | 75 Tau/HD 28292/HIP 20877/HR 1407/HYG 20825/MN 1545                                   | nan            | nan     | 75 Tau      | Tau             |   4.47399 |    16.3597 |   8.12 |   17.79 |     18   |    57.241  |   4.96 |    nan     |    nan     | 1.137 | K2IIIvar |     29.621   | 21.3634  |  50.5985 | 16.1229  |  9.92e-06  |  1.587e-05 |  4.25e-06  |     20825 | nan          |            0 |            0 |    4.47332 |     16.3374 |  0.505882 |  -3.09337  | -1.42006e+11 |      807954 | 20.3678 | 109.339 | 30.5308 | 59.4692 |


And plot them again:


```python
fig,axs = hyades.plot_stars(coords=['RAEpoch','DecEpoch'])
fig.savefig('gallery/hyades-precessed.png')
```


    
![png](MontuPython-QuickStart_files/MontuPython-QuickStart_30_0.png)
    


### Working with time

One of the most interesting and basic functionalities of MontuPython is to convert date among 
different type of calendars and astronomical scales.  You may taste these functionalities using:


```python
mtime = montu.Time('bce2501-01-01 12:00:00')
```

other alternative formats for the same date are:


```python
mtime = montu.Time('2501 b.c.e. 01-01 12:00:00')
mtime = montu.Time('-2500-01-01 12:00:00')
```

If you print this time object you will get:


```python
print(mtime)
```

    Montu Time Object:
    --------------------------
    Date in proleptic UTC: -2500-01-01 12:00:00.0002
    Date in mixed UTC: -2500-01-22 12:00:00
    Date in SPICE format: 2501 B.C. 01-01 12:00:00.200
    General:
        Components: [-1, 2500, 1, 1, 12, 0, 0, 200]
        Is bce: True
        Is Julian: True
    Uniform scales:
        Terrestrial time:
            tt: -142006202700.3199
            jtd: 807954.69096852
        UTC time:
            et: -142006262399.99988
            jed: 807954.0
        Delta-t = TT - UTC = 59699.68000000001
    Objects:
        Date in datetime64 format: -2500-01-01T12:00:00.000200
        Date in PyPlanet Epoch: 807954.0
        Date in PyEphem Epoch: -2501/1/22 12:00:00
        Date in AstroPy Time: 807954.69096852
    Astronomical properties at Epoch:
        True obliquity of ecliptic: 23:58:33.587
        True nutation longitude: 00:00:10.214
        Greenwhich Meridian Sidereal Time: 18:40:25.323
    


Notice that the date in Gregorian proleptic will be 2501 b.c.e 01-01 but in the mixed calendar that uses Julian calendar before its adoption at 1582-10-04, will be 2501 bce 01-22.

You may add or substract time to a given date. This is done by adding or substracting seconds to the reference time:


```python
mtime = montu.Time(0,format='tt')
(mtime,
 mtime - 12*montu.HOUR, 
 mtime + 1*montu.DAY, 
 mtime - 3*montu.CALYEAR, 
 mtime + 20*montu.JULYEAR)
```




    (Time('2000-01-01 12:00:00.0000'/'2000-01-01 12:00:00'),
     Time('2000-01-01 00:00:00.0000'/'2000-01-01 00:00:00'),
     Time('2000-01-02 12:00:00.0000'/'2000-01-02 12:00:00'),
     Time('1997-01-01 12:00:01.5319'/'1997-01-01 12:00:00'),
     Time('2020-01-01 11:59:52.2521'/'2020-01-01 11:59:59'))



As you may notice, adding or substracting and integer number of seconds not necesarily correspond to adding or sutracting days or years to the calendar. This is because of different factors. One of them are the leap seconds. Normally leapseconds are included every once in a while. However to calculate ephemerides in the ancient world, `MontuPython` uses a continuous model of deltat that produces these discrepancies.

### Working with planets and observing sites

`MontuPython` allows calculate the position of all planets in the solar system, including the moon:


```python
earth = montu.Planet('Earth')
mars = montu.Planet('Mars')
neptune = montu.Planet('neptune')
jupiter = montu.Planet('JUPITER')
moon = montu.Planet('Moon')
```

You may create an observing site in the surface of any planet of the solar system:


```python
tebas = montu.Observer(planet=earth,lon=33,lat=24,height=0)
```

The straight routine method to calculate position of the planet in the sky at any date is:


```python
mtime = montu.Time('-2500-01-01 12:00:00')
mars.where_among_stars(mtime,tebas)
TABLEDF(mars.position)
```

    |    |    jed |           tt |   RAJ2000 |   DecJ2000 |         pmRA |         paRA |       pmDec |       paDec |   LonJ2000 |   LatJ2000 |        pmLon |        paLon |       pmLat |        paLat |   site_distance |   sun_distance |   elongation |   phase |   mag |      XJ2000 |       YJ2000 |      ZJ2000 |   VXJ2000 |   VYJ2000 |   VZJ2000 |
    |----|--------|--------------|-----------|------------|--------------|--------------|-------------|-------------|------------|------------|--------------|--------------|-------------|--------------|-----------------|----------------|--------------|---------|-------|-------------|--------------|-------------|-----------|-----------|-----------|
    |  0 | 807954 | -1.42006e+11 |   12.5302 |    1.62005 | -2.06254e+07 | -3.39321e+08 | 1.72852e+08 | 1.61982e+09 |    186.663 |    4.64342 | -3.53556e+08 | -5.33309e+09 | 3.65425e+07 | -5.54062e+08 |         0.66045 |        1.62611 |      157.822 | 13.3589 |  -1.1 | -9.7214e+07 | -1.31655e+07 | 3.00274e+06 |   6.77554 |   5.83293 |   2.44421 |


You may get also the position at different times and store it into a data frame:


```python
import numpy as np

# Reset the data stored in the planet
mars.reset_store()

# Loop on different times
for deltat in np.arange(0,100*montu.YEAR,10*montu.YEAR):
    mars.where_among_stars(mtime + deltat,tebas,store=True)

# Show results
TABLEDF(mars.data)
```

    |    |    jed |           tt |   RAJ2000 |   DecJ2000 |         pmRA |         paRA |        pmDec |        paDec |   LonJ2000 |   LatJ2000 |        pmLon |        paLon |        pmLat |        paLat |   site_distance |   sun_distance |   elongation |   phase |   mag |       XJ2000 |       YJ2000 |       ZJ2000 |    VXJ2000 |   VYJ2000 |   VZJ2000 |
    |----|--------|--------------|-----------|------------|--------------|--------------|--------------|--------------|------------|------------|--------------|--------------|--------------|--------------|-----------------|----------------|--------------|---------|-------|--------------|--------------|--------------|------------|-----------|-----------|
    |  0 | 807954 | -1.42006e+11 | 12.5302   |    1.62005 | -2.06254e+07 | -3.39321e+08 |  1.72852e+08 |  1.61982e+09 |  186.663   |   4.64342  | -3.53556e+08 | -5.33309e+09 |  3.65425e+07 | -5.54062e+08 |        0.66045  |        1.62611 |     157.822  | 13.3589 |  -1.1 | -9.7214e+07  | -1.31655e+07 |  3.00274e+06 |   6.77554  |   5.83293 |   2.44421 |
    |  1 | 811607 | -1.41691e+11 | 21.0132   |  -18.2028  |  6.89538e+07 | -6.19843e+07 |  2.79034e+08 |  1.48896e+09 |  312.39    |  -1.16497  |  1.02157e+09 | -2.09215e+07 | -9.67462e+06 |  5.73308e+07 |        2.13311  |        1.39928 |      32.4582 | 22.4411 |   1.1 |  2.1788e+08  | -2.10529e+08 | -9.83786e+07 |  32.0002   |  35.8586  |  15.1836  |
    |  2 | 815259 | -1.41375e+11 |  1.90954  |   11.529   |  5.66079e+07 |  2.49856e+07 |  3.37651e+08 | -8.20501e+08 |   30.6923  |  -0.198064 |  8.97521e+08 | -2.24781e+08 |  2.59682e+07 | -4.81603e+07 |        2.05806  |        1.54147 |      45.9176 | 27.6332 |   1.2 |  2.63799e+08 |  1.47972e+08 |  6.31053e+07 | -11.33     |  38.9387  |  18.1761  |
    |  3 | 818912 | -1.41059e+11 | 10.1451   |   16.7684  | -3.01404e+07 |  2.04569e+08 |  1.54548e+08 | -2.24682e+09 |  148.212   |   4.98701  | -4.61319e+08 |  3.64542e+09 | -8.27678e+06 | -8.28787e+08 |        0.670374 |        1.64769 |     162.849  | 10.26   |  -1.1 | -8.48593e+07 |  4.54329e+07 |  2.92148e+07 |   0.550876 |   7.14267 |   3.27063 |
    |  4 | 822564 | -1.40744e+11 | 20.0935   |  -21.1787  |  7.01507e+07 | -4.4581e+07  |  1.92042e+08 |  1.77227e+09 |  299.072   |  -0.852742 |  9.99803e+08 |  7.72502e+07 | -1.55048e+07 |  1.93085e+07 |        1.93104  |        1.42446 |      45.5776 | 29.9319 |   0.9 |  1.4328e+08  | -2.27509e+08 | -1.03374e+08 |  33.6395   |  28.0052  |  11.5538  |
    |  5 | 826217 | -1.40428e+11 |  1.15967  |    6.70872 |  5.78953e+07 |  5.74816e+06 |  3.80035e+08 | -6.09309e+08 |   18.5964  |  -0.625838 |  9.42242e+08 | -3.39326e+08 |  2.4192e+07  | -1.48386e+07 |        2.21636  |        1.49875 |      34.0177 | 21.8064 |   1.3 |  3.13711e+08 |  1.02268e+08 |  4.0497e+07  |  -6.38497  |  44.2337  |  20.4051  |
    |  6 | 829869 | -1.40113e+11 |  7.95318  |   24.87    | -3.22157e+06 |  4.3325e+08  | -1.50498e+07 | -1.38206e+09 |  116.429   |   4.07891  | -4.01637e+07 |  6.06472e+09 | -2.33179e+07 | -2.0166e+08  |        0.819955 |        1.65841 |     131.799  | 26.5728 |  -0.6 | -5.48691e+07 |  9.79048e+07 |  5.19942e+07 |  -4.92466  |   9.86752 |   4.71003 |
    |  7 | 833522 | -1.39797e+11 | 19.0946   |  -22.9612  |  6.89809e+07 | -1.73085e+07 |  8.90434e+07 |  1.82221e+09 |  285.087   |  -0.377725 |  9.56733e+08 |  1.05139e+08 | -1.86458e+07 | -2.50144e+07 |        1.68988  |        1.46031 |      59.3723 | 35.8998 |   0.7 |  6.86005e+07 | -2.21528e+08 | -9.79208e+07 |  32.2347   |  20.4358  |   8.16299 |
    |  8 | 837174 | -1.39482e+11 |  0.424627 |    1.70368 |  5.93218e+07 | -1.43502e+07 |  4.04973e+08 | -3.12415e+08 |    6.52179 |  -0.96445  |  9.77213e+08 | -3.74749e+08 |  2.03235e+07 |  2.1864e+07  |        2.32982  |        1.45715 |      22.1397 | 14.9142 |   1.4 |  3.46306e+08 |  4.28028e+07 |  1.22545e+07 |   0.500054 |  48.0946  |  21.9142  |
    |  9 | 840827 | -1.39166e+11 |  6.24893  |   26.3567  |  2.50937e+07 |  2.7097e+08  | -1.66469e+07 | -4.05003e+08 |   93.3498  |   2.96077  |  3.38048e+08 |  3.67477e+09 | -7.88252e+06 | -6.50616e+07 |        1.04125  |        1.65814 |     109.013  | 34.5667 |  -0.1 | -9.88316e+06 |  1.40473e+08 |  6.97311e+07 |  -9.37521  |  13.8042  |   6.66527 |


Once calculated, you can compute the horizontal coordinates as we show before, but, since each table already has a date, you don't need to provide an epoch:


```python
montu.Heka.where_in_sky(mars,site=tebas,inplace=True)
TABLEDF(mars.data)
```

    |    |    jed |           tt |   RAJ2000 |   DecJ2000 |         pmRA |         paRA |        pmDec |        paDec |   LonJ2000 |   LatJ2000 |        pmLon |        paLon |        pmLat |        paLat |   site_distance |   sun_distance |   elongation |   phase |   mag |       XJ2000 |       YJ2000 |       ZJ2000 |    VXJ2000 |   VYJ2000 |   VZJ2000 |   RAJ2000t |   DecJ2000t |   RAEpoch |   DecEpoch |     epoch_tt |   epoch_jed |        HA |        az |        el |      zen |
    |----|--------|--------------|-----------|------------|--------------|--------------|--------------|--------------|------------|------------|--------------|--------------|--------------|--------------|-----------------|----------------|--------------|---------|-------|--------------|--------------|--------------|------------|-----------|-----------|------------|-------------|-----------|------------|--------------|-------------|-----------|-----------|-----------|----------|
    |  0 | 807954 | -1.42006e+11 | 12.5302   |    1.62005 | -2.06254e+07 | -3.39321e+08 |  1.72852e+08 |  1.61982e+09 |  186.663   |   4.64342  | -3.53556e+08 | -5.33309e+09 |  3.65425e+07 | -5.54062e+08 |        0.66045  |        1.62611 |     157.822  | 13.3589 |  -1.1 | -9.7214e+07  | -1.31655e+07 |  3.00274e+06 |   6.77554  |   5.83293 |   2.44421 |  12.5302   |     1.62005 |   8.53605 |    24.1079 | -1.42006e+11 |      807954 | 12.3376   |   6.19008 | -41.642   | 131.642  |
    |  1 | 811607 | -1.41691e+11 | 21.0132   |  -18.2028  |  6.89538e+07 | -6.19843e+07 |  2.79034e+08 |  1.48896e+09 |  312.39    |  -1.16497  |  1.02157e+09 | -2.09215e+07 | -9.67462e+06 |  5.73308e+07 |        2.13311  |        1.39928 |      32.4582 | 22.4411 |   1.1 |  2.1788e+08  | -2.10529e+08 | -9.83786e+07 |  32.0002   |  35.8586  |  15.1836  |  21.0132   |   -18.2028  |  16.564   |   -23.1532 | -1.41691e+11 |      811607 | 16.3919   | 103.816   | -30.2122  | 120.212  |
    |  2 | 815259 | -1.41375e+11 |  1.90954  |   11.529   |  5.66079e+07 |  2.49856e+07 |  3.37651e+08 | -8.20501e+08 |   30.6923  |  -0.198064 |  8.97521e+08 | -2.24781e+08 |  2.59682e+07 | -4.81603e+07 |        2.05806  |        1.54147 |      45.9176 | 27.6332 |   1.2 |  2.63799e+08 |  1.47972e+08 |  6.31053e+07 | -11.33     |  38.9387  |  18.1761  |   1.90954  |    11.529   |  22.0745  |   -12.6072 | -1.41375e+11 |      815259 | 22.9628   | 155.776   |  50.3606  |  39.6394 |
    |  3 | 818912 | -1.41059e+11 | 10.1451   |   16.7684  | -3.01404e+07 |  2.04569e+08 |  1.54548e+08 | -2.24682e+09 |  148.212   |   4.98701  | -4.61319e+08 |  3.64542e+09 | -8.27678e+06 | -8.28787e+08 |        0.670374 |        1.64769 |     162.849  | 10.26   |  -1.1 | -8.48593e+07 |  4.54329e+07 |  2.92148e+07 |   0.550876 |   7.14267 |   3.27063 |  10.1451   |    16.7684  |   5.72227 |    28.5327 | -1.41059e+11 |      818912 |  3.39657  | 287.203   |  44.4213  |  45.5787 |
    |  4 | 822564 | -1.40744e+11 | 20.0935   |  -21.1787  |  7.01507e+07 | -4.4581e+07  |  1.92042e+08 |  1.77227e+09 |  299.072   |  -0.852742 |  9.99803e+08 |  7.72502e+07 | -1.55048e+07 |  1.93085e+07 |        1.93104  |        1.42446 |      45.5776 | 29.9319 |   0.9 |  1.4328e+08  | -2.27509e+08 | -1.03374e+08 |  33.6395   |  28.0052  |  11.5538  |  20.0935   |   -21.1787  |  15.66    |   -20.2988 | -1.40744e+11 |      822564 |  5.54013  | 248.715   |  -2.18906 |  92.1891 |
    |  5 | 826217 | -1.40428e+11 |  1.15967  |    6.70872 |  5.78953e+07 |  5.74816e+06 |  3.80035e+08 | -6.09309e+08 |   18.5964  |  -0.625838 |  9.42242e+08 | -3.39326e+08 |  2.4192e+07  | -1.48386e+07 |        2.21636  |        1.49875 |      34.0177 | 21.8064 |   1.3 |  3.13711e+08 |  1.02268e+08 |  4.0497e+07  |  -6.38497  |  44.2337  |  20.4051  |   1.15967  |     6.70872 |  21.3219  |   -16.7997 | -1.40428e+11 |      826217 | 11.9592   | 355.331   | -82.7769  | 172.777  |
    |  6 | 829869 | -1.40113e+11 |  7.95318  |   24.87    | -3.22157e+06 |  4.3325e+08  | -1.50498e+07 | -1.38206e+09 |  116.429   |   4.07891  | -4.01637e+07 |  6.06472e+09 | -2.33179e+07 | -2.0166e+08  |        0.819955 |        1.65841 |     131.799  | 26.5728 |  -0.6 | -5.48691e+07 |  9.79048e+07 |  5.19942e+07 |  -4.92466  |   9.86752 |   4.71003 |   7.95318  |    24.87    |   3.43869 |    22.8408 | -1.40113e+11 |      829869 | 17.9236   |  68.5454  |   8.10871 |  81.8913 |
    |  7 | 833522 | -1.39797e+11 | 19.0946   |  -22.9612  |  6.89809e+07 | -1.73085e+07 |  8.90434e+07 |  1.82221e+09 |  285.087   |  -0.377725 |  9.56733e+08 |  1.05139e+08 | -1.86458e+07 | -2.50144e+07 |        1.68988  |        1.46031 |      59.3723 | 35.8998 |   0.7 |  6.86005e+07 | -2.21528e+08 | -9.79208e+07 |  32.2347   |  20.4358  |   8.16299 |  19.0946   |   -22.9612  |  14.7523  |   -16.122  | -1.39797e+11 |      833522 | 18.6904   | 108.92    |   2.56904 |  87.431  |
    |  8 | 837174 | -1.39482e+11 |  0.424627 |    1.70368 |  5.93218e+07 | -1.43502e+07 |  4.04973e+08 | -3.12415e+08 |    6.52179 |  -0.96445  |  9.77213e+08 | -3.74749e+08 |  2.03235e+07 |  2.1864e+07  |        2.32982  |        1.45715 |      22.1397 | 14.9142 |   1.4 |  3.46306e+08 |  4.28028e+07 |  1.22545e+07 |   0.500054 |  48.0946  |  21.9142  |   0.424627 |     1.70368 |  20.5393  |   -20.2896 | -1.39482e+11 |      837174 |  0.984469 | 199.222   |  43.434   |  46.566  |
    |  9 | 840827 | -1.39166e+11 |  6.24893  |   26.3567  |  2.50937e+07 |  2.7097e+08  | -1.66469e+07 | -4.05003e+08 |   93.3498  |   2.96077  |  3.38048e+08 |  3.67477e+09 | -7.88252e+06 | -6.50616e+07 |        1.04125  |        1.65814 |     109.013  | 34.5667 |  -0.1 | -9.88316e+06 |  1.40473e+08 |  6.97311e+07 |  -9.37521  |  13.8042  |   6.66527 |   6.24893  |    26.3567  |   1.94529 |    14.7574 | -1.39166e+11 |      840827 |  7.65831  | 294.41    | -15.5445  | 105.544  |


#### Daily motion of planet and stars

We want to study how the elevation and azimuth of a star changes as time passes at a given epoch.  Let's take the case for instance of Aldebaran in the latest example.


```python
aldebaran_path = montu.Heka.move_over_nut(aldebaran.data,at=mtime,site=tebas,
                                          during=6*montu.HOUR,each=1*montu.HOUR)
TABLEDF(aldebaran_path)
```

    |    |           tt |        HA |      az |      el |     zen | isvisible   |
    |----|--------------|-----------|---------|---------|---------|-------------|
    |  0 | -1.42006e+11 | 20.2635   | 107.598 | 29.5927 | 60.4073 | True        |
    |  1 | -1.42006e+11 | 21.2662   | 117.59  | 42.2888 | 47.7112 | True        |
    |  2 | -1.42006e+11 | 22.2689   | 132.53  | 53.5817 | 36.4183 | True        |
    |  3 | -1.42006e+11 | 23.2717   | 156.473 | 61.6793 | 28.3207 | True        |
    |  4 | -1.42006e+11 |  0.274424 | 189.232 | 63.4436 | 26.5564 | True        |
    |  5 | -1.42006e+11 |  1.27716  | 217.992 | 57.8112 | 32.1888 | True        |
    |  6 | -1.42006e+11 |  2.2799   | 236.477 | 47.6476 | 42.3524 | True        |


We can plot the path on the sky of the star:


```python
fig,axs = montu.Heka.plot_over_nut(aldebaran_path)
```


    
![png](MontuPython-QuickStart_files/MontuPython-QuickStart_57_0.png)
    


We may make the same with planets:


```python
mars.reset_store()
mars.where_among_stars(mtime,tebas,store=True)
mars_path = montu.Heka.move_over_nut(mars.data.iloc[0],at=mtime,site=tebas,
                                     during=12*montu.HOUR,each=1*montu.HOUR)
```


```python
fig,axs = montu.Heka.plot_over_nut([aldebaran_path,mars_path])
```


    
![png](MontuPython-QuickStart_files/MontuPython-QuickStart_60_0.png)
    


## An indepth example: evolution of pole stars

Choose from database all bright stars that according to [wikipedia](https://en.wikipedia.org/wiki/Pole_star#Precession_of_the_equinoxes) were or will be close to the celestial North pole:


```python
star_names = ('Polaris','Vega','Thuban','Deneb','Alderamin','Kochab')
stars = stars.get_stars(ProperName=star_names)
TABLEDF(stars.data)
```

    |     |   MN |     HD |   HR |    HIP | Gl     | Name      | OtherDesignations                                                                   | ProperName   | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |    B-V | SpType       |   Luminosity |     XJ2000 |     YJ2000 |    ZJ2000 |    VXJ2000 |   VYJ2000 |    VZJ2000 |   Primary |   MultipleID |   IsMultiple |   IsVariable |
    |-----|------|--------|------|--------|--------|-----------|-------------------------------------------------------------------------------------|--------------|---------|-------------|-----------------|-----------|------------|--------|---------|----------|------------|--------|------------|------------|--------|--------------|--------------|------------|------------|-----------|------------|-----------|------------|-----------|--------------|--------------|--------------|
    |   5 |    6 | 172167 | 7001 |  91262 | Gl 721 | Vega      | 3 Lyr/3Alp Lyr/Gl 721/HD 172167/HIP 91262/HR 7001/HYG 90978/MN 6/Vega/α Lyr         | Vega         | α Lyr   | 3 Lyr       | Lyr             |  18.6156  |    38.7837 | 201.02 |  287.46 |    -12.1 |     7.6787 |   0.03 |    nan     |    nan     | -0.001 | A0Vvar       |      49.9344 |   0.960565 |   -5.90801 |   4.80973 |  5.9e-07   | 1.734e-05 |  4.76e-06  |     90979 |          nan |            0 |            0 |
    |  21 |   22 | 197345 | 7924 | 102098 | nan    | Deneb     | 50 Cyg/50Alp Cyg/Deneb/HD 197345/HIP 102098/HR 7924/HYG 101766/MN 22/α Cyg          | Deneb        | α Cyg   | 50 Cyg      | Cyg             |  20.6905  |    45.2803 |   1.56 |    1.55 |     -5   |   432.9    |   1.25 |      1.294 |      1.224 |  0.092 | A2Ia         |   51617.9    | 197.251    | -232.113   | 307.601   | -1.34e-06  | 6.62e-06  | -1.33e-06  |    101767 |          nan |            0 |            1 |
    |  48 |   49 |   8890 |  424 |  11767 | nan    | Polaris   | 1 UMi/1Alp UMi/HD 8890/HIP 11767/HR 424/HYG 11734/MN 49/Polaris/α UMi               | Polaris      | α UMi   | 1 UMi       | UMi             |   2.52975 |    89.2641 |  44.22 |  -11.74 |    -17   |   132.626  |   1.97 |      1.993 |      1.953 |  0.636 | F7:Ib-IIv SB |    2495.74   |   1.3431   |    1.04763 | 132.615   | -1.748e-05 | 2.692e-05 | -1.171e-05 |     11734 |          nan |            0 |            1 |
    |  59 |   60 | 131873 | 5563 |  72607 | nan    | Kochab    | 7 UMi/7Bet UMi/HD 131873/HIP 72607/HR 5563/HYG 72379/Kochab/MN 60/β UMi             | Kochab       | β UMi   | 7 UMi       | UMi             |  14.8451  |    74.1555 | -32.29 |   11.91 |     17   |    40.1445 |   2.07 |    nan     |    nan     |  1.465 | K4IIIvar     |     208.545  |  -8.05816  |   -7.42971 |  38.6194  |  1.736e-05 | 2.91e-06  | -6.11e-06  |     72380 |          nan |            0 |            0 |
    |  91 |   92 | 203280 | 8162 | 105199 | Gl 826 | Alderamin | 5 Cep/5Alp Cep/Alderamin/Gl 826/HD 203280/HIP 105199/HR 8162/HYG 104860/MN 92/α Cep | Alderamin    | α Cep   | 5 Cep       | Cep             |  21.3096  |    62.5856 | 149.91 |   48.27 |    -11.5 |    15.0376 |   2.45 |    nan     |    nan     |  0.257 | A7IV-V       |      20.6253 |   5.27611  |   -4.4832  |  13.3488  | -8.82e-06  | 1.386e-05 |  5.7e-07   |    104861 |          nan |            0 |            0 |
    | 351 |  352 | 123299 | 5291 |  68756 | nan    | Thuban    | 11 Dra/11Alp Dra/HD 123299/HIP 68756/HR 5291/HYG 68536/MN 352/Thuban/α Dra          | Thuban       | α Dra   | 11 Dra      | Dra             |  14.0732  |    64.3758 | -56.52 |   17.19 |    -13   |    92.9368 |   3.67 |    nan     |    nan     | -0.049 | A0III SB     |     256.094  | -34.416    |  -20.7588  |  83.7964  | -8.64e-06  | 2.838e-05 | -2.25e-06  |     68537 |          nan |            0 |            0 |


Now precess the position of all stars from -20 000 to 20 000 years from 2000:


```python
import pandas as pd
import tqdm

now = montu.Time()
df = pd.DataFrame()
for dt in tqdm.tqdm(np.linspace(-20000*montu.YEAR,20000*montu.YEAR,1000)):
    past = now + dt
    pstars = montu.Heka.precess_to_epoch(stars.data,at=past)
    row = dict(tt = past.tt)
    for star in star_names:
        row.update({star:float(pstars[pstars.ProperName == star].DecEpoch)})
    df = pd.concat([df,pd.DataFrame([row])])
```

      0%|          | 0/1000 [00:00<?, ?it/s]100%|██████████| 1000/1000 [00:40<00:00, 24.73it/s]


Now plot declinations as a function of time:


```python
import matplotlib.pyplot as plt

fig,ax = plt.subplots(figsize=(10,6))
for star in star_names:
    ax.plot(df['tt'],df[star],label=star)

ax.legend(loc='lower center',ncol=len(star_names))
ax.set_xlabel("Time [year]")
ax.set_ylabel("Declination [deg]")
ax.axvline(montu.Time().tt,color='k',lw=3)
ax.text(0.5,1.01,'Now',ha='center',transform=ax.transAxes)
ax.margins(0)
ax.set_xticks(np.linspace(df['tt'].min(),df['tt'].max(),14))
ax.grid()
montu.Time.set_time_ticks(ax)
montu.Util.montu_mark(ax)
fig.savefig('gallery/pole-stars.png')
```


    
![png](MontuPython-QuickStart_files/MontuPython-QuickStart_67_0.png)
    


Check date when star is close to the pole:


```python
for star in star_names:
    imax = df[star].argmax()
    mtime = montu.Time(df.iloc[imax].tt)
    print(f"Star {star} will be the closest to the pole at {mtime.datespice} (declination {montu.D2H(df.iloc[imax][star])})")
```

    Star Polaris will be the closest to the pole at 2083-11-17 10:17:16.2830 (declination 89:31:50.395)
    Star Vega will be the closest to the pole at 11612 B.C. 11-20 03:05:46.95900 (declination 86:21:54.919)
    Star Thuban will be the closest to the pole at 2803 B.C. 11-22 10:26:57.411400 (declination 89:55:30.839)
    Star Deneb will be the closest to the pole at 14735 B.C. 09-08 22:33:57.715900 (declination 86:57:7.337)
    Star Alderamin will be the closest to the pole at 7529-06-08 04:08:34.4940 (declination 87:58:44.091)
    Star Kochab will be the closest to the pole at 1041 B.C. 09-09 09:46:16.893700 (declination 83:29:35.893)
