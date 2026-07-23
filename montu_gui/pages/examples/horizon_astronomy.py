"""
Horizon Astronomy with MontuPython

This script demonstrates how to compute and plot the topographic 
horizon for a historical observer site, both as a pure landscape 
profile and combined with the celestial sphere (stars, planets, etc.).
"""
import montu

# 1. Define the observer location and date
print("Setting up observer at Giza...")
giza = montu.Observer(site="giza")

# 2. Calculate the horizon
giza.horizon_profile(
    max_dist=30,       # km
    az_step=1,         # degrees
    coarse_step=2,     # km
    verbose=1
)

# 3. Plotting the bare horizon profile
# Let's focus on the Eastern horizon (az_center = 90)
print("\nPlotting bare topographical horizon...")
fig1 = giza.horizon.plot_horizon(
    az_center=90,
    az_delta=90
)

# 4. Plotting the horizon with the celestial sphere at a specific time
at = montu.Time("-2550-01-22 06:10:00", calendar="mixed")
print("Plotting horizon with stars, constellation boundaries, and the Sun...")
fig2 = giza.horizon.plot_horizon(
    at=at,
    az_center=110,   # Sunrise azimuth in winter
    az_delta=60,
    elev_view=30,
    show_boundaries=False,
    show_asterism=True,
    show_starnames=True,
    show_planets=["Sun", "Moon"]
)