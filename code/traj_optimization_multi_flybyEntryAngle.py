
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
        - For solving the lambert problem, start and end position of the spacecraft is the center of Earth and of Jupiter respectively.
        - The flyby is assumed to be an instantaneuos maneuver, so no time pasts between arriving and leaving Jupiter.
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

"""
The following parameters need to be specified:
    * departure_time:           Time() object of the astropy.time package
    * flyby_time:               Time() object of the astropy.time package
    * height_flyby:             float in [u.km]; corresponds to the closest distance of the flyby to Jupiter wrt. the surface

    * initial_angle:            float in [rad]; corresponds to the initial entry angle of the s/c at the flyby that is tested
    * end_angle:                float in [rad]; corresponds to the end entry angle of the s/c at the flyby that is tested
    * iteration_step_size_rad:  float in [rad]: corresponds to the angle step size
"""

# Fixed parameters
departure_time = Time("2037-09-24", scale="tdb")
height_flyby = 391000 * u.km
flyby_time = Time("2040-01-18", scale="tdb")

# Parameters for multi-run simulation
initial_angle = 0.068
end_angle = 0.070
iteration_step_size_rad = 0.001

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
flyby_entry_angle_rad_list = list()             # Flyby entry angle at Jupiter


current_angle = initial_angle
iteration = 0
number_of_iterations = int((end_angle - initial_angle) / iteration_step_size_rad)

# ---------- Start iterations ----------
print("\nSimulation started!\n")
while current_angle < end_angle:
    
    # Update console
    if (iteration+1 % 10 == 0):
        print(f"Iteration\t\t\t{iteration+1} / {number_of_iterations}")
    # Execute pipeline
    _, delta_v, _, distance_scToSaturn, crossing_time, flag_saturn_crossed, best_distance_scToSaturn_vector = pipeline.simulate_run(departure_time, flyby_time, height_flyby, current_angle)

    # Collect results of current run
    departure_time_list.append(departure_time)
    flyby_time_list.append(flyby_time)
    distance_of_flyby_list.append(height_flyby)
    delta_v_list.append(np.linalg.norm(delta_v.value))
    distance_scToSaturn_crossing_list.append(distance_scToSaturn)
    arrival_time_list.append(crossing_time)
    flag_saturn_crossed_list.append(flag_saturn_crossed)
    flyby_entry_angle_rad_list.append(current_angle)

    # Adjust flyby angle for next iteration
    iteration += 1
    current_angle += iteration_step_size_rad


# ---------- Save data ----------
fh.save_data_in_csv(departure_time_list, flyby_time_list, distance_of_flyby_list, delta_v_list, distance_scToSaturn_crossing_list, arrival_time_list, flag_saturn_crossed_list, flyby_entry_angle_rad_list)


# ---------- Save parameters ----------
flag_type_of_simulation = 3
fh.save_parameters_in_txt(flag_type_of_simulation, departure_time, flyby_time, None, None, height_flyby, None, None, initial_angle, end_angle, iteration_step_size_rad)

# End timer for runtime analysis
end_time = time.time()
print(f"\n\nSimulation finished in {end_time-start_time} seconds!\n")