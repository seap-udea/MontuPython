###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import pandas as pd

###############################################################
# Module constants
###############################################################
PLANETARY_DATAFILE = "planets-jpl.csv"


def load_planets():
    """Load the planetary data table shipped with the package.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by planet name, containing orbital and physical
        parameters. A derived column ``SynodicOrbit`` (synodic period in
        years) is added automatically.

    Examples
    --------
    >>> import montu
    >>> planets = montu.load_planets()
    >>> planets.loc["Mars", "SiderealOrbit"]
    1.8808...
    """
    planets = pd.read_csv(
        montu.Util._data_path(PLANETARY_DATAFILE, check=True),
        sep=";",
    )
    planets.set_index("Planet", inplace=True)

    planets["SynodicOrbit"] = abs(
        1 / (1 / planets.loc["Earth", "SiderealOrbit"] - 1 / planets["SiderealOrbit"])
    )

    return planets
