##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
################################# Pipeline ###################################
##############################################################################

import numpy as np
from astropy.time import Time
from astropy import units as u
from poliastro.bodies import Sun, Earth, Jupiter, Saturn
from poliastro.twobody import Orbit
from poliastro.iod import lambert
from poliastro.ephem import Ephem
from poliastro.threebody.flybys import compute_flyby


def simulate_run(departure_time, flyby_time, height_flyby):

    # ---------- Get Earth, Jupiter, and Saturn's ephemeris and create their orbits ---------- 
    earth_ephem = Ephem.from_body(Earth, departure_time)
    jupiter_ephem = Ephem.from_body(Jupiter, flyby_time)

    earth = Orbit.from_ephem(Sun, earth_ephem, departure_time)
    jupiter = Orbit.from_ephem(Sun, jupiter_ephem, flyby_time)


    # ---------- Solve Lambert Problem (Earth to Jupiter transfer) ---------- 
    v_sc_departure, v_sc_flybyArrival = lambert(Sun.k, earth.r, jupiter.r, (flyby_time - departure_time).to(u.s))


    # ---------- Compute spacecraft trajectory post-flyby ----------
    _, v_jupiter_flyby = jupiter_ephem.rv()
    v_jupiter_flyby = v_jupiter_flyby[0]
    v_jupiter_flyby.to(u.km / u.day)
    radius_flyby = Jupiter.R + height_flyby                 # radius of flyby wrt. Jupiter center

    v_sc_flybyDeparture, delta_ = compute_flyby(v_sc_flybyArrival, v_jupiter_flyby, Jupiter.k, radius_flyby)

    #  Create the orbit after the flyby
    post_flyby_orbit = Orbit.from_vectors(Sun, jupiter.r, v_sc_flybyDeparture, epoch=flyby_time)


    # Optional: Adjust orbital plance of s/c orbit, to the same as Saturn
    # TODO

    # ---------- Propagate S/C orbit after flyby and calculate the distance to Saturn ----------
    step_size = 0.5 * u.day                # Time step for propagation
    max_time = 3 * u.year                   # Maximum time to propagate after flyby

    current_time = flyby_time

    distance_scToSaturn_list = list()
    timestamps = list()

    best_distance = 1e20
    best_timing = None
    
    saturn_ephem = Ephem.from_body(Saturn, flyby_time)

    # Search for timestamp when the spacecraft has closest distance to Saturn
    while current_time < flyby_time + max_time:
        r_spacecraft = post_flyby_orbit.propagate(current_time).r.to(u.km).value                # Position of s/c

        saturn_ephem = Ephem.from_body(Saturn, current_time)                                    
        saturn = Orbit.from_ephem(Sun, saturn_ephem, current_time)
        r_saturn = saturn.r.to(u.km).value                                                      # Position of Saturn

        distance_scToSaturn = np.linalg.norm(r_spacecraft - r_saturn)                           # Distance of s/c to Saturn
        
        # Check if current distance is closer than the best one so far, if so -> save timestamp and new best distance
        if distance_scToSaturn < best_distance:
            best_distance = distance_scToSaturn
            best_timing = current_time
        
        # Collect data for further analysis
        distance_scToSaturn_list.append(distance_scToSaturn)
        timestamps.append(current_time.value)

        # Increase timestamps
        current_time += step_size

    crossing_time = best_timing
    # Check if s/c crosses Saturn
    flag_saturn_crossed = False
    if best_distance < 1e8: flag_saturn_crossed = True

    return v_sc_departure, v_sc_flybyDeparture, best_distance, crossing_time, flag_saturn_crossed
