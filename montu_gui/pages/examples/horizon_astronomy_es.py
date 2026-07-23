"""
Astronomía en el Horizonte con MontuPython

Este script demuestra cómo calcular y graficar el horizonte 
topográfico para un sitio de observación histórico, tanto como 
un perfil de paisaje puro como combinado con la esfera celeste 
(estrellas, planetas, etc.).
"""
import montu

# 1. Definir la ubicación del observador y la fecha
print("Configurando el observador en Giza...")
giza = montu.Observer(site="giza")

# 2. Calcular el horizonte
giza.horizon_profile(
    max_dist=30,       # km
    az_step=1,         # grados
    coarse_step=2,     # km
    verbose=1
)

# 3. Graficar el perfil del horizonte desnudo
# Nos enfocamos en el horizonte Este (az_center = 90)
print("\nGraficando el horizonte topográfico desnudo...")
fig1 = giza.horizon.plot_horizon(
    az_center=90,
    az_delta=90
)

# 4. Graficar el horizonte con la esfera celeste en un momento específico
at = montu.Time("-2550-01-22 06:10:00", calendar="mixed")
print("Graficando el horizonte con estrellas, constelaciones y el Sol...")
fig2 = giza.horizon.plot_horizon(
    at=at,
    az_center=110,   # Azimut del amanecer en invierno
    az_delta=60,
    elev_view=30,
    show_boundaries=False,
    show_asterism=True,
    show_starnames=True,
    show_planets=["Sun", "Moon"]
)

# Limpia el cache
montu.Horizon.clean_cache(verbose=True)