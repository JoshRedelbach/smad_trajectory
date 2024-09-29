
##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
#################### Multiple Runs - Variable Flyby Date #####################
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

import file_handling as fh
import time

import pipeline

# ----------------------------------------------------------------------------

""" PARAMETERS """

# Fixed parameters
departure_time = Time("2038-02-01", scale="tdb")        # Time of departure
height_flyby = 410000 * u.km                            # Height of flyby wrt. surface of Jupiter

# Define parameters for multi-run simulation
number_of_iterations = 10                                # Number of iterations performed
flyby_time = Time("2039-12-01", scale="tdb")             # Initial time of flyby at Jupiter
iteration_step_size_day = 10 * u.day                     # Time step for iterations

# ----------------------------------------------------------------------------

""" CODE """

# Start timer for runtime analysis
start_time = time.time()

# Initialize lists for collecting the data to be saved
departure_time_list = list()                    # Time of departure
flyby_time_list = list()                        # Time of flyby at Jupiter
distance_of_flyby_list = list()                 # Distance of flyby wrt. surface of Jupiter
delta_v_list = list()                           # Delta V
distance_scToSaturn_crossing_list = list()      # Distance of s/c to Saturn when crossing Saturn orbit
arrival_time_list = list()                      # Time of arrival at Saturn (when s/c crosses Saturn orbit)
flag_saturn_crossed_list = list()               # Flag that indicates if the trajectory of the spacecraft crosses the orbit of Saturn

iteration = 0


# ---------- Start iterations ----------
print("\nSimulation started!\n")
while iteration < number_of_iterations:

    # Update console
    if (iteration+1 % 50 == 0):
        print(f"Iteration\t\t\t{iteration+1} / {number_of_iterations}")

    # Execute pipeline
    delta_v, _, distance_scToSaturn, crossing_time, flag_saturn_crossed = pipeline.simulate_run(departure_time, flyby_time, height_flyby)

    # Collect results of current run
    departure_time_list.append(departure_time)
    flyby_time_list.append(flyby_time)
    distance_of_flyby_list.append(height_flyby)
    delta_v_list.append(np.linalg.norm(delta_v.value))
    distance_scToSaturn_crossing_list.append(distance_scToSaturn)
    arrival_time_list.append(crossing_time)
    flag_saturn_crossed_list.append(flag_saturn_crossed)

    # Adjust flyby time for next iteration
    iteration += 1
    flyby_time += iteration_step_size_day



# ---------- Save data ----------
fh.save_data_in_csv(departure_time_list, flyby_time_list, distance_of_flyby_list, delta_v_list, distance_scToSaturn_crossing_list, arrival_time_list, flag_saturn_crossed_list)


# ---------- Save parameters ----------
flag_type_of_simulation = 1
fh.save_parameters_in_txt(flag_type_of_simulation, departure_time, height_flyby, number_of_iterations, flyby_time, iteration_step_size_day)

# End timer for runtime analysis
end_time = time.time()
print(f"\n\nSimulation finished in {end_time-start_time} seconds!\n")