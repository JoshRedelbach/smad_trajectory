
##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
######## Multiple Runs - Variable Flyby Date, Distance & Entry Angle #########
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

# Departure date
departure_time_initial = Time("2037-09-25", scale="tdb")    # Time of departure
departure_time_stepsize = 1.0 * u.day                       # Step size of iterations over departure time
departure_time_number_of_iterations = 1                    # Number of iterations the departure_time_stepsize is added to flyby_time_initial

# Flyby date
flyby_time_initial = Time("2039-07-29", scale="tdb")        # Initial time of flyby that is tested (included)
flyby_time_stepsize = 1.0 * u.day                           # Step size of iterations over flyby time
flyby_time_number_of_iterations = 1                         # Number of iterations the flyby_time_stepsize is added to flyby_time_initial

# Flyby distance
# !! Minimum flyby distance to the center of Jupiter 643.428km from the center of the Jupiter (radius of Jupiter 69.911km) !!
flyby_height_initial = 1159000                               # Height of flyby wrt. surface of Jupiter [km](included)
flyby_height_end = 1161000                                   # Height of flyby wrt. surface of Jupiter [km](included)
flyby_height_stepsize = 100                               # Step size of iterations over flyby distance [km]

# Flyby entry angle
flyby_angle_initial = 0.102                                   # Angle of flyby [rad] (included)
flyby_angle_end = 0.104                                       # Angle of flyby [rad] (included)
flyby_angle_stepsize = 0.00005                                 # Step size of iterations over flyby angle [rad]

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

number_of_total_iterations = departure_time_number_of_iterations * flyby_time_number_of_iterations * ((flyby_angle_end-flyby_angle_initial) / flyby_angle_stepsize) * ((flyby_height_end-flyby_height_initial) / flyby_height_stepsize)
total_iterator = 0

print("\n")
# ---------- Start iteration I -> DEPARTURE DATE ----------
current_departure_time = departure_time_initial
iteration_departure_time = 0

while iteration_departure_time < departure_time_number_of_iterations:

    # ---------- Start iteration II -> FLYBY DATE ----------
    current_flyby_time = flyby_time_initial
    iteration_flyby_time = 0

    while iteration_flyby_time < flyby_time_number_of_iterations:

        # ---------- Start iteration III -> FLYBY HEIGHT ----------
        current_flyby_height = flyby_height_initial
        while current_flyby_height <= flyby_height_end:

            # ---------- Start iteration IV -> FLYBY ANGLE ----------
            current_flyby_angle = flyby_angle_initial
            while current_flyby_angle <= flyby_angle_end:
                
                total_iterator += 1
                if(total_iterator % 100 == 0):
                    print(f"Iteration: \t {total_iterator} / {number_of_total_iterations}")
                # Execute Pipeline
                _, delta_v, _, distance_scToSaturn, crossing_time, flag_saturn_crossed, best_distance_scToSaturn_vector = pipeline.simulate_run(current_departure_time, current_flyby_time, current_flyby_height * u.km, -1. * current_flyby_angle)

                # Collect results of current run
                departure_time_list.append(current_departure_time)
                flyby_time_list.append(current_flyby_time.value)
                distance_of_flyby_list.append(current_flyby_height)
                delta_v_list.append(np.linalg.norm(delta_v.value))
                distance_scToSaturn_crossing_list.append(distance_scToSaturn)
                arrival_time_list.append(crossing_time)
                flag_saturn_crossed_list.append(flag_saturn_crossed)
                flyby_entry_angle_rad_list.append(-1. * current_flyby_angle)

                current_flyby_angle += flyby_angle_stepsize
            
            current_flyby_height += flyby_height_stepsize
        
        iteration_flyby_time += 1
        current_flyby_time += flyby_time_stepsize
    
    iteration_departure_time += 1
    current_departure_time += departure_time_stepsize


# ---------- Save data ----------
fh.save_data_in_csv(departure_time_list, flyby_time_list, distance_of_flyby_list, delta_v_list, distance_scToSaturn_crossing_list, arrival_time_list, flag_saturn_crossed_list, flyby_entry_angle_rad_list)


# ---------- Save parameters ----------
fh.save_parameters_in_txt_multi_varying(departure_time_initial, departure_time_stepsize, departure_time_number_of_iterations, flyby_time_initial, flyby_time_stepsize, flyby_time_number_of_iterations, flyby_height_initial, flyby_height_end, flyby_height_stepsize, flyby_angle_initial, flyby_angle_end, flyby_angle_stepsize)

# End timer for runtime analysis
end_time = time.time()
print(f"\n\nSimulation finished in {end_time-start_time} seconds!\n")