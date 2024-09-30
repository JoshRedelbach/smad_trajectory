
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
        - position of s/c at departure from Earth: center of Earth
        - positions of s/c before flyby (used for Lambert problem): center of Jupiter
        - positions of s/c after flyby (used for orbit determination/propagation): center of Jupiter
        - the time of crossing is determined by comparing the distance of the spacecraft to the sun and the semi-major axis of the Saturn
            -> so assumed that the s/c is in the right plane and that trajectory of spacecraft is in the orbit plane of Saturn

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
departure_time = Time("2038-02-01", scale="tdb")
flyby_time = Time("2040-01-18", scale="tdb")
height_flyby = 391000 * u.km


# ----------------------------------------------------------------------------


""" CODE """

# Start timer for runtime analysis
start_time = time.time()

# ---------- Execute pipeline ----------

delta_v, v_sc_flybyDeparture, distance_scToSaturn, crossing_time, flag_saturn_crossed = pipeline.simulate_run(departure_time, flyby_time, height_flyby)



# ---------- Print results to console ----------

print(f"\nSpeed of s/c required at departure [km/s]:\t\t\t {np.linalg.norm(delta_v.value)}\n")

if flag_saturn_crossed: print("\nSatellite crosses Saturn's orbit!\n")    
else: print("\nSatellite does NOT cross Saturn's orbit!\n")

print(f"Distance s/c to Saturn [km]:\t\t\t\t\t {distance_scToSaturn}")
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
initial_orbit = Orbit.from_vectors(Sun, earth.r, delta_v, epoch=departure_time)

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