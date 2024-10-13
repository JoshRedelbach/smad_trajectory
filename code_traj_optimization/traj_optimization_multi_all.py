
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

"""
The following parameters need to be specified:
    * departure_time_initial:               Time() object of the astropy.time package; initial time of departure that is tested
    * departure_time_stepsize:              float in [u.day]: corresponds to the time step size
    * departure_time_number_of_iterations:  integer; indicates how often the step size is added to the initial departure date

    * flyby_time_initial:                   Time() object of the astropy.time package; initial date of flyby that is tested
    * flyby_time_stepsize:                  float in [u.day]: corresponds to the time step size
    * flyby_time_number_of_iterations:      integer; indicates how often the step size is added to the initial flyby date

    * flyby_height_initial:                 float in [km]; corresponds to the initial closest flyby height of the s/c during the flyby wrt. the surface of Jupiter that is tested
    * flyby_height_end:                     float in [km]; corresponds to the last closest flyby height of the s/c during the flyby wrt. the surface of Jupiter that is tested
    * flyby_height_stepsize:                float in [km]; indicates the step size that is added to the initial flyby height
    
    * flyby_angle_initial:                  float in [rad]; corresponds to the initial entry angle of the flyby that is tested
    * flyby_angle_end:                      float in [rad]; corresponds to the last entry angle of the flyby that is tested
    * flyby_angle_stepsize:                 float in [rad]; indicates the step size that is added to the initial flyby entry angle
"""

# Departure date
departure_time_initial = Time("2037-09-25", scale="tdb")
departure_time_stepsize = 1.0 * u.day
departure_time_number_of_iterations = 1

# Flyby date
flyby_time_initial = Time("2039-07-29", scale="tdb")
flyby_time_stepsize = 1.0 * u.day
flyby_time_number_of_iterations = 1

# Flyby distance
# !! Minimum flyby distance to the center of Jupiter 643 428 km from the center of the Jupiter due to radiation (radius of Jupiter 69.911km) !!
flyby_height_initial = 1161000
flyby_height_end = 1161100
flyby_height_stepsize = 100

# Flyby entry angle
flyby_angle_initial = 0.1
flyby_angle_end = 0.2
flyby_angle_stepsize = 0.1

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