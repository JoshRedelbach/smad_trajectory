
##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
################################# Single Run #################################
##############################################################################

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OLD VERSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!! DIFFERENT COMPUTATION OF CROSSING !!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
from poliastro.iod import lambert
from poliastro.plotting.static import StaticOrbitPlotter
from poliastro.ephem import Ephem
from poliastro.threebody.flybys import compute_flyby
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------

""" PARAMETERS """

# Define time of departure and flyby maneuver
departure_time = Time("2030-01-01", scale="tdb")
flyby_time = Time("2032-01-01", scale="tdb")
height_flyby = 500000 * u.km                   # height of flyby wrt. surface of Jupiter

# Voyager 2 parameters
# departure_time = Time("1977-08-20", scale="tdb")
# flyby_time = Time("1979-07-09", scale="tdb")
# height_flyby = 570000 * u.km

# Voyager 1 parameters
# departure_time = Time("1977-09-05", scale="tdb")
# flyby_time = Time("1979-03-05", scale="tdb")
# height_flyby = 349000 * u.km

# # Possible parameters for our mission
# departure_time = Time("2038-02-01", scale="tdb")
# flyby_time = Time("2040-02-01", scale="tdb")
# height_flyby = 410000 * u.km

# ----------------------------------------------------------------------------

""" CODE """

# ---------- Get Earth, Jupiter, and Saturn's ephemeris and create their orbits ---------- 
earth_ephem = Ephem.from_body(Earth, departure_time)
jupiter_ephem = Ephem.from_body(Jupiter, flyby_time)
saturn_ephem = Ephem.from_body(Saturn, flyby_time)

earth = Orbit.from_ephem(Sun, earth_ephem, departure_time)
jupiter = Orbit.from_ephem(Sun, jupiter_ephem, flyby_time)
saturn = Orbit.from_ephem(Sun, saturn_ephem, flyby_time)


# ---------- Solve Lambert Problem (Earth to Jupiter transfer) ---------- 
v_sc_departure, v_sc_flybyArrival = lambert(Sun.k, earth.r, jupiter.r, (flyby_time - departure_time).to(u.s))
print(f"\nSpeed of s/c required at departure [km/s]:\t\t\t {np.linalg.norm(v_sc_departure.value)}\n")


# ---------- Create spacecraft orbit after departure ----------
initial_orbit = Orbit.from_vectors(Sun, earth.r, v_sc_departure, epoch=departure_time)


# ---------- Compute spacecraft trajectory post-flyby ----------
_, v_jupiter_flyby = jupiter_ephem.rv()
v_jupiter_flyby = v_jupiter_flyby[0]
v_jupiter_flyby.to(u.km / u.day)
radius_flyby = Jupiter.R + height_flyby                 # radius of flyby wrt. Jupiter center

v_sc_flybyDeparture, delta_ = compute_flyby(v_sc_flybyArrival, v_jupiter_flyby, Jupiter.k, radius_flyby)

#  Create the orbit after the flyby
post_flyby_orbit = Orbit.from_vectors(Sun, jupiter.r, v_sc_flybyDeparture, epoch=flyby_time)


# ---------- Calculate when the spacecraft crosses Saturn's orbit ----------
saturn_semimajoraxis = saturn.a         # Saturn's semi-major axis (orbital radius)
step_size = 0.25 * u.day                # Time step for propagation
max_time = 5 * u.year                   # Maximum time to propagate after flyby

current_time = flyby_time

r_spacecraft_list = list()
distance_spacecraft_saturnorbit_list = list()
timestamps = list()

best_distance = 1e20
best_timing = None

# Search for timestamp when the spacecraft is crossing the orbit of Saturn
while current_time < flyby_time + max_time:
    r_spacecraft = np.linalg.norm(post_flyby_orbit.propagate(current_time).r.to(u.km).value)        # Distance from the Sun to the spacecraft
    distance_spacecraft_saturnorbit = np.abs(r_spacecraft - saturn_semimajoraxis.value)             # Distance from spacecraft to orbit of Saturn
    
    # Check if current distance is closer than the best one so far, if so -> save timestamp and new best distance
    if distance_spacecraft_saturnorbit < best_distance:
        best_distance = distance_spacecraft_saturnorbit
        best_timing = current_time
    
    # Collect data for further analysis
    r_spacecraft_list.append(r_spacecraft)
    distance_spacecraft_saturnorbit_list.append(distance_spacecraft_saturnorbit)
    timestamps.append(current_time.value)

    # Increase timestamps
    current_time += step_size

# Check if spacecraft crosses the orbit of Saturn
if best_distance < 1e6:
    crossing_time = best_timing
    print("\nSatellite crosses Saturn's orbit!\n")
    saturn_ephem = Ephem.from_body(Saturn, crossing_time)
    saturn = Orbit.from_ephem(Sun, saturn_ephem, crossing_time)
else:
    print("\nSatellite does NOT cross Saturn's orbit!\n")

print(f"Closest distance of spacecraft to saturn orbit [km]:\t\t {best_distance}")
print(f"Time of closest distance of spacecraft to saturn orbit:\t\t {best_timing}\n")


# ---------- Determine distance of spacecraft to Saturn when the spacecraft crosses the orbit of Saturn ----------
pos_sc_crossing = post_flyby_orbit.propagate(crossing_time).r.to(u.km).value
print(post_flyby_orbit.propagate(crossing_time).v.value)
pos_saturn_crossing = saturn.r.to(u.km).value
distance_scToSaturn_crossing = np.linalg.norm(pos_sc_crossing - pos_saturn_crossing)

print(f"Position s/c when crossing [km]:\t\t\t\t {pos_sc_crossing}")
print(f"Position Saturn when crossing [km]:\t\t\t\t {pos_saturn_crossing}")
print(f"Distance s/c to Saturn [km]:\t\t\t\t\t {distance_scToSaturn_crossing}")

print("\n")


# ---------- Plots ----------
# plt.plot(timestamps, distance_spacecraft_saturnorbit_list)
# plt.show()

fig, ax = plt.subplots(figsize=(10, 10))
plotter = StaticOrbitPlotter(ax)

# Plot Earth, Jupiter, Saturn orbits
plotter.plot(earth, label="Earth's Orbit", color='blue', trail="true")
plotter.plot(jupiter, label="Jupiter's Orbit", color='orange', trail="true")
plotter.plot(saturn, label="Saturn's Orbit", color='green', trail="true")

# Plot spacecraft's initial trajectory (Earth to Jupiter)
plotter.plot(initial_orbit, label="Spacecraft (Earth to Jupiter)", color='purple')

# Plot post-flyby trajectory (Jupiter to Saturn)
# plotter.plot(post_flyby_orbit, label="Spacecraft (Jupiter to Saturn)", color='red')
post_flyby_orbit = Orbit.from_vectors(Sun, jupiter.r, v_sc_flybyDeparture, epoch=crossing_time)
plotter.plot(post_flyby_orbit, label="Spacecraft (Jupiter to Saturn)", color='red')


# Show plot with legends
ax.legend()
plt.show()

print(post_flyby_orbit.inc * 57.3)
print(saturn.inc * 57.3)





# #  Create the orbit after the flyby
# post_flyby_orbit_2 = Orbit.from_vectors(Sun, jupiter.r, v_sc_flybyDeparture, epoch=flyby_time)
# post_flyby_orbit_2.inc = saturn.inc



# # ---------- Calculate when the spacecraft crosses Saturn's orbit ----------
# saturn_semimajoraxis = saturn.a         # Saturn's semi-major axis (orbital radius)
# step_size = (1./24.) * u.day                # Time step for propagation
# max_time = 5 * u.year                   # Maximum time to propagate after flyby

# current_time = flyby_time

# r_spacecraft_list = list()
# distance_spacecraft_saturnorbit_list = list()
# timestamps = list()

# best_distance = 1e20
# best_timing = None

# # Search for timestamp when the spacecraft is crossing the orbit of Saturn
# while current_time < flyby_time + max_time:
#     r_spacecraft = np.linalg.norm(post_flyby_orbit_2.propagate(current_time).r.to(u.km).value)        # Distance from the Sun to the spacecraft
#     distance_spacecraft_saturnorbit = np.abs(r_spacecraft - saturn_semimajoraxis.value)             # Distance from spacecraft to orbit of Saturn
    
#     # Check if current distance is closer than the best one so far, if so -> save timestamp and new best distance
#     if distance_spacecraft_saturnorbit < best_distance:
#         best_distance = distance_spacecraft_saturnorbit
#         best_timing = current_time
    
#     # Collect data for further analysis
#     r_spacecraft_list.append(r_spacecraft)
#     distance_spacecraft_saturnorbit_list.append(distance_spacecraft_saturnorbit)
#     timestamps.append(current_time.value)

#     # Increase timestamps
#     current_time += step_size

# # Check if spacecraft crosses the orbit of Saturn
# if best_distance < 1e6:
#     crossing_time = best_timing
#     print("\nSatellite crosses Saturn's orbit!\n")
#     saturn_ephem = Ephem.from_body(Saturn, crossing_time)
#     saturn = Orbit.from_ephem(Sun, saturn_ephem, crossing_time)
# else:
#     print("\nSatellite does NOT cross Saturn's orbit!\n")

# print(f"Closest distance of spacecraft to saturn orbit [km]:\t\t {best_distance}")
# print(f"Time of closest distance of spacecraft to saturn orbit:\t\t {best_timing}\n")


# # ---------- Determine distance of spacecraft to Saturn when the spacecraft crosses the orbit of Saturn ----------
# pos_sc_crossing = post_flyby_orbit_2.propagate(crossing_time).r.to(u.km).value
# pos_saturn_crossing = saturn.r.to(u.km).value
# distance_scToSaturn_crossing = np.linalg.norm(pos_sc_crossing - pos_saturn_crossing)

# print(f"Position s/c when crossing [km]:\t\t\t\t {pos_sc_crossing}")
# print(f"Position Saturn when crossing [km]:\t\t\t\t {pos_saturn_crossing}")
# print(f"Distance s/c to Saturn [km]:\t\t\t\t\t {distance_scToSaturn_crossing}")

# print("\n")
# print(post_flyby_orbit_2.inc * 57.3)
# print(saturn.inc * 57.3)

# # ---------- Plots ----------
# # plt.plot(timestamps, distance_spacecraft_saturnorbit_list)
# # plt.show()

# fig, ax = plt.subplots(figsize=(10, 10))
# plotter = StaticOrbitPlotter(ax)

# # Plot Earth, Jupiter, Saturn orbits
# plotter.plot(earth, label="Earth's Orbit", color='blue', trail="true")
# plotter.plot(jupiter, label="Jupiter's Orbit", color='orange', trail="true")
# plotter.plot(saturn, label="Saturn's Orbit", color='green', trail="true")

# # Plot spacecraft's initial trajectory (Earth to Jupiter)
# plotter.plot(initial_orbit, label="Spacecraft (Earth to Jupiter)", color='purple')

# # Plot post-flyby trajectory (Jupiter to Saturn)
# plotter.plot(post_flyby_orbit_2, label="Spacecraft (Jupiter to Saturn)", color='red')


# # Show plot with legends
# ax.legend()
# plt.show()

