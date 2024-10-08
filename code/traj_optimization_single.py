
##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
################################# Single Run #################################
##############################################################################

# ----------------------------------------------------------------------------

"""
NOTE:
    - be aware of assumptions!

    -> assumptions:
        - For solving the lambert problem, start and end position of the spacecraft is the center of Earth and of Jupiter respectively.
        - The flyby is assumed to be an instantaneuos maneuver, so no time pasts between arriving and leaving Jupiter.
"""

# ----------------------------------------------------------------------------

""" IMPORTS """

import numpy as np
from astropy.time import Time
from astropy import units as u
from poliastro.bodies import Sun, Earth, Jupiter, Saturn
from poliastro.twobody import Orbit
from poliastro.plotting.static import StaticOrbitPlotter
from poliastro.ephem import Ephem
import matplotlib.pyplot as plt
import pipeline
import time

# ----------------------------------------------------------------------------

""" PARAMETERS """

"""
The following parameters need to be specified:
    * departure_time:   Time() object of the astropy.time package
    * flyby_time:       Time() object of the astropy.time package
    * height_flyby:     float in [km] but has to be multiplied by u.km (astropy.units); corresponds to the closest distance of the flyby to Jupiter wrt. the surface
    * flyby_angle:      float in [rad]; corresponds to the entry angle of the s/c at the flyby of Jupiter 
"""

# Template
# departure_time = Time("YYYY-MM-DD", scale="tdb")
# flyby_time = Time("YYYY-MM-DD", scale="tdb")
# height_flyby = XX * u.km
# flyby_angle = XX

# Voyager 1 parameters (unknown entry angle)
# departure_time = Time("1977-09-05", scale="tdb")
# flyby_time = Time("1979-03-05", scale="tdb")
# height_flyby = 349000 * u.km
# flyby_angle = 0.

# Voyager 2 parameters (unknown entry angle)
# departure_time = Time("1977-08-20", scale="tdb")
# flyby_time = Time("1979-07-09", scale="tdb")
# height_flyby = 570000 * u.km
# flyby_angle = 0.

# Possible parameters for our mission
# departure_time = Time("2038-02-01", scale="tdb")
# flyby_time = Time("2040-01-18", scale="tdb")
# height_flyby = 391000 * u.km
# flyby_angle = 0.06860000000000004

# Selected parameters for the project mission
departure_time = Time("2037-09-25", scale="tdb")
flyby_time = Time("2039-07-29", scale="tdb")
height_flyby = 1159500 * u.km
flyby_angle = -0.10265

# ----------------------------------------------------------------------------

""" CODE """

# Start timer for runtime analysis
start_time = time.time()

# ---------- Execute pipeline ----------

v_sc_departure, delta_v_leo, v_sc_flybyDeparture, distance_scToSaturn, crossing_time, flag_saturn_crossed, distance_scToSaturn_vector = pipeline.simulate_run(departure_time, flyby_time, height_flyby, flyby_angle)


# ---------- Print results to console ----------

print(f"\nSpeed of s/c required at Earth [km/s]:\t\t\t\t {np.linalg.norm(v_sc_departure.value)}")
print(f"Delta-V that needs to be applied [km/s]:\t\t\t {delta_v_leo.value}")

if flag_saturn_crossed: print("\nSatellite meets Saturn!\n")    
else: print("\nSatellite does NOT meet Saturn!\n")

print(f"Closest distance s/c to Saturn [km]:\t\t\t\t {distance_scToSaturn}")
print(f"Distance vector s/c to Saturn [km]:\t\t\t\t {distance_scToSaturn_vector}")
print(f"Distance vector s/c to Saturn [normalized]:\t\t\t {distance_scToSaturn_vector / distance_scToSaturn}")
print(f"Time of closest distance of spacecraft to saturn orbit:\t\t {crossing_time}\n")

# End timer for runtime analysis
end_time = time.time()
print(f"\nSimulation finished in {end_time-start_time} seconds!\n")


# ---------- Plots ----------

# Get Earth, Jupiter, and Saturn's ephemeris and create their orbits
earth_ephem = Ephem.from_body(Earth, departure_time)
jupiter_ephem = Ephem.from_body(Jupiter, flyby_time)
saturn_ephem = Ephem.from_body(Saturn, crossing_time)

earth = Orbit.from_ephem(Sun, earth_ephem, departure_time)
jupiter = Orbit.from_ephem(Sun, jupiter_ephem, flyby_time)
saturn = Orbit.from_ephem(Sun, saturn_ephem, crossing_time)

# Create spacecraft orbit after departure and after flyby
initial_orbit = Orbit.from_vectors(Sun, earth.r, v_sc_departure, epoch=departure_time)

post_flyby_orbit = Orbit.from_vectors(Sun, jupiter.r, v_sc_flybyDeparture, epoch=flyby_time)

# Create figure
fig, ax = plt.subplots(figsize=(10, 10))
plotter = StaticOrbitPlotter(ax)

# Plot spacecraft's initial trajectory (Earth to Jupiter)
plotter.plot(initial_orbit, label="Spacecraft (Earth to Jupiter)", color='purple')

# Plot post-flyby trajectory (Jupiter to Saturn)
plotter.plot(post_flyby_orbit, label="Spacecraft (Jupiter to Saturn)", color='red')

# Plot Earth, Jupiter, Saturn orbits
plotter.plot(earth, label="Earth's Orbit", color='blue', trail="true")
plotter.plot(jupiter, label="Jupiter's Orbit", color='orange', trail="true")
plotter.plot(saturn, label="Saturn's Orbit", color='green', trail="true")

# Show plot with legend
ax.legend()
plt.show()