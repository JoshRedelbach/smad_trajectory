
##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
################## Multiple Runs - Variable Flyby Distance ###################
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
    * entry_angle_flyby:        float in [rad]; corresponds to the entry angle of the s/c at the flyby

    * initial_height_flyby:     float in [u.km]; corresponds to the initial closest flyby height of the s/c during the flyby wrt. the surface of Jupiter
    * iteration_step_size_km:   float in [u.km]: corresponds to the height step size
    * number_of_iterations:     integer; indicates how often the step size is added to the initial flyby height
"""

# Fixed parameters
departure_time = Time("2038-02-01", scale="tdb")
flyby_time = Time("2040-01-18", scale="tdb")
entry_angle_flyby = 0.

# Parameters for multi-run simulation
initial_height_flyby = 380000 * u.km
iteration_step_size_km = 500 * u.km
number_of_iterations = 2

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

iteration = 0

height_flyby = initial_height_flyby

# ---------- Start iterations ----------
print("\nSimulation started!\n")
while iteration < number_of_iterations:

    # Update console
    if (iteration+1 % 10 == 0):
        print(f"Iteration\t\t\t{iteration+1} / {number_of_iterations}")

    # Execute pipeline
    _, delta_v, _, distance_scToSaturn, crossing_time, flag_saturn_crossed, best_distance_scToSaturn_vector = pipeline.simulate_run(departure_time, flyby_time, height_flyby, entry_angle_flyby)

    # Collect results of current run
    departure_time_list.append(departure_time)
    flyby_time_list.append(flyby_time)
    distance_of_flyby_list.append(height_flyby.value)
    delta_v_list.append(np.linalg.norm(delta_v.value))
    distance_scToSaturn_crossing_list.append(distance_scToSaturn)
    arrival_time_list.append(crossing_time)
    flag_saturn_crossed_list.append(flag_saturn_crossed)
    flyby_entry_angle_rad_list.append(entry_angle_flyby)

    # Adjust flyby time for next iteration
    iteration += 1
    height_flyby += iteration_step_size_km


# ---------- Save data ----------
fh.save_data_in_csv(departure_time_list, flyby_time_list, distance_of_flyby_list, delta_v_list, distance_scToSaturn_crossing_list, arrival_time_list, flag_saturn_crossed_list, flyby_entry_angle_rad_list)


# ---------- Save parameters ----------
flag_type_of_simulation = 2
fh.save_parameters_in_txt(flag_type_of_simulation, departure_time, flyby_time, None, None, height_flyby, number_of_iterations, iteration_step_size_km, entry_angle_flyby, None, None)

# End timer for runtime analysis
end_time = time.time()
print(f"\n\nSimulation finished in {end_time-start_time} seconds!\n")