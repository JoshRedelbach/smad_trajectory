##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
################################# Pipeline ###################################
##############################################################################

import numpy as np
from astropy import units as u
from poliastro.bodies import Sun, Earth, Jupiter, Saturn
from poliastro.twobody import Orbit
from poliastro.iod import lambert
from poliastro.ephem import Ephem
from poliastro.threebody.flybys import compute_flyby


def simulate_run(departure_time, flyby_time, height_flyby, flyby_entry_angle=0):

    # ---------- Get Earth, Jupiter, and Saturn's ephemeris and create their orbits ---------- 
    earth_ephem = Ephem.from_body(Earth, departure_time)
    jupiter_ephem = Ephem.from_body(Jupiter, flyby_time)

    earth = Orbit.from_ephem(Sun, earth_ephem, departure_time)
    jupiter = Orbit.from_ephem(Sun, jupiter_ephem, flyby_time)


    # ---------- Solve Lambert Problem (Earth to Jupiter transfer) ---------- 
    v_sc_departure, v_sc_flybyArrival = lambert(Sun.k, earth.r, jupiter.r, (flyby_time - departure_time).to(u.s))


    # ---------- Solve for delta v required in LEO ----------
    # v_inf_vec = v_sc - v_earth
    _, v_earth = earth.rv()
    v_inf_vec = v_sc_departure.to(u.km / u.s) - v_earth.to(u.km / u.s)
    
    # print(f"v_sc_departure: {v_sc_departure}")
    # print(f"norm(v_sc_departure): {np.linalg.norm(v_sc_departure.to(u.km / u.s))}")
    # print(f"v_earth: {v_earth.to(u.km / u.s)}")
    # print(f"norm(v_earth): {np.linalg.norm(v_earth.to(u.km / u.s))}")
    # print(f"v_inf_vec: {v_inf_vec}")
    # print(f"norm(v_inf_vec): {np.linalg.norm(v_inf_vec.to(u.km / u.s))}")

    # Get the semi-major axis of the hyperbolic 
    # v_inf = |v_inf_vec|
    # v_inf = sqrt( - u / a)   =>  a = - u / (v_inf)^2
    v_inf = np.linalg.norm(v_inf_vec).to(u.km / u.s)
    # print(f"V_inf: {v_inf}")
    # print(f"Earth.k: {Earth.k}")
    a_orbit_escape = - Earth.k.to((u.km**3  ) / (u.s**2)) / (v_inf.to(u.km / u.s) **2)
    # print(f"a_orbit_escape: {a_orbit_escape}")

    # Energy of the orbit
    energy_orbit_escape = - Earth.k.to((u.km**3  ) / (u.s**2)) / (2. * a_orbit_escape)
    # print(f"energy_orbit_escape: {energy_orbit_escape}")

    # Get the velocity we need at out LEO to escape
    # E = (v_p^2 / 2) - (u / r_p)
    # v_p = sqrt(2 * (E + u / r_p))
    r_perigee_orbit_escape = Earth.R.to(u.km) + 200 * u.km
    # print(f"r_perigee_orbit_escape: {r_perigee_orbit_escape}")
    v_perigee_orbit_escape = np.sqrt(2. * (energy_orbit_escape + (Earth.k.to((u.km**3  ) / (u.s**2)) / r_perigee_orbit_escape)) )
    # print(f"v_perigee_orbit_escape: {v_perigee_orbit_escape}")

    # Calculate the delta v that we need to apply
    v_leo = np.sqrt(Earth.k.to((u.km**3  ) / (u.s**2)) / r_perigee_orbit_escape)
    delta_v_leo = v_perigee_orbit_escape - v_leo
    # print(f"v_leo: {v_leo}")
    # print(f"\n\nDELTA-V LEO: {delta_v_leo}\n\n")
    

    # ---------- Compute spacecraft trajectory post-flyby ----------
    _, v_jupiter_flyby = jupiter_ephem.rv()
    v_jupiter_flyby = v_jupiter_flyby[0]
    radius_flyby = Jupiter.R + height_flyby                 # radius of flyby wrt. Jupiter center
    v_jupiter_flyby.to(u.km / u.day)
    
    v_sc_flybyDeparture, delta_ = compute_flyby(v_sc_flybyArrival, v_jupiter_flyby, Jupiter.k, radius_flyby, flyby_entry_angle * u.rad)

    #  Create the orbit after the flyby
    post_flyby_orbit = Orbit.from_vectors(Sun, jupiter.r, v_sc_flybyDeparture, epoch=flyby_time)


    # ---------- Propagate S/C orbit after flyby and calculate the distance to Saturn ----------
    step_size = 1. / 12. * u.day              # Time step for propagation
    max_time = 3 * u.year                   # Maximum time to propagate after flyby

    current_time = flyby_time

    distance_scToSaturn_list = list()
    timestamps = list()

    best_distance_vector = None
    best_distance = 1e20
    best_timing = None
    
    saturn_ephem = Ephem.from_body(Saturn, flyby_time)

    # Search for timestamp when the spacecraft has closest distance to Saturn
    while current_time < flyby_time + max_time:
        r_spacecraft = post_flyby_orbit.propagate(current_time).r.to(u.km).value                # Position of s/c

        saturn_ephem = Ephem.from_body(Saturn, current_time)                                    
        saturn = Orbit.from_ephem(Sun, saturn_ephem, current_time)
        r_saturn = saturn.r.to(u.km).value                                                      # Position of Saturn

        distance_scToSaturn_vector = r_spacecraft - r_saturn
        distance_scToSaturn = np.linalg.norm(distance_scToSaturn_vector)                        # Distance of s/c to Saturn
        
        # Check if current distance is closer than the best one so far, if so -> save timestamp and new best distance
        if distance_scToSaturn < best_distance:
            best_distance = distance_scToSaturn
            best_timing = current_time
            best_distance_vector = distance_scToSaturn_vector
        
        # Collect data for further analysis
        distance_scToSaturn_list.append(distance_scToSaturn)
        timestamps.append(current_time.value)

        # Increase timestamps
        current_time += step_size

    crossing_time = best_timing
    # Check if s/c crosses Saturn
    flag_saturn_crossed = False
    if best_distance < 1e6: flag_saturn_crossed = True

    # print(f"\nVelocity when arriving at Saturn: \t\t {post_flyby_orbit.propagate(crossing_time).v.to(u.km / u.s).value}\n")
    # print(f"\nVelocity when arriving at Saturn: \t\t {np.linalg.norm(post_flyby_orbit.propagate(crossing_time).v.to(u.km / u.s).value)}\n")

    return v_sc_departure, delta_v_leo, v_sc_flybyDeparture, best_distance, crossing_time, flag_saturn_crossed, best_distance_vector
