
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
TO-DO:
    - EXPLAIN ALL PARAMETERS IN THE BEGINNING 
"""

# Define time of departure and flyby maneuver
# departure_time = Time("2030-01-01", scale="tdb")
# flyby_time = Time("2032-01-01", scale="tdb")
# height_flyby = 500000 * u.km                   # height of flyby wrt. surface of Jupiter

# Voyager 2 parameters
# departure_time = Time("1977-08-20", scale="tdb")
# flyby_time = Time("1979-07-09", scale="tdb")
# height_flyby = 570000 * u.km

# Voyager 1 parameters
# departure_time = Time("1977-09-05", scale="tdb")
# flyby_time = Time("1979-03-05", scale="tdb")
# height_flyby = 349000 * u.km

# Possible parameters for our mission
# departure_time = Time("2038-02-01", scale="tdb")
# flyby_time = Time("2040-01-18", scale="tdb")
# height_flyby = 391000 * u.km
# flyby_angle = 0.06860000000000004

# Possible parameters for our mission
# departure_time = Time("2038-02-01", scale="tdb")
# flyby_time = Time("2040-01-18", scale="tdb")
# height_flyby = 570000 * u.km

# Possible parameters for our mission
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

print(f"\nSpeed of s/c required at departure [km/s]:\t\t\t {np.linalg.norm(v_sc_departure.value)}\n")

if flag_saturn_crossed: print("\nSatellite crosses Saturn's orbit!\n")    
else: print("\nSatellite does NOT cross Saturn's orbit!\n")

print(f"Distance s/c to Saturn [km]:\t\t\t\t\t {distance_scToSaturn}")
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