##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
############################### File Handling ################################
##############################################################################



def save_data_in_csv(departure_time_list, flyby_time_list, distance_of_flyby_list, delta_v_list, distance_scToSaturn_crossing_list, arrival_time_list, flag_saturn_crossed_list, flyby_entry_angle_rad_list):
    """ Save the collected data of the multiple runs in a .csv-file. """

    print("\n\nStarted saving results.\n")
    print("...")

    path = "Results/NEW_RESULTS/data.csv"
    file = open(path, 'w')
    
    # Define coloumns
    file.write('Departure Time,Flyby Time,Arrival Time,Distance of Flyby to Jupiter,Delta V,Status Saturn crossed,Crossing distance S/C-Saturn,Entry Angle of Flyby\n')

    # Insert data
    for i in range(len(departure_time_list)):
        file.write(str(departure_time_list[i]))
        file.write(',')
        file.write(str(flyby_time_list[i]))
        file.write(',')
        file.write(str(arrival_time_list[i]))
        file.write(',')
        file.write(str(distance_of_flyby_list[i]))
        file.write(',')
        file.write(str(delta_v_list[i]))
        file.write(',')
        file.write(str(flag_saturn_crossed_list[i]))
        file.write(',')
        file.write(str(distance_scToSaturn_crossing_list[i]))
        file.write(',')
        file.write(str(flyby_entry_angle_rad_list[i]))
        file.write('\n')
    
    file.close()

    print("\nFinished saving results.\n")


# ----------------------------------------------------------------------------


def save_parameters_in_txt(flag_type_of_simulation, departure_time, flyby_time_initial, flyby_time_iterations, flyby_time_step_size, flyby_height_initial, flyby_height_iterations, flyby_height_step_size, flyby_angle_initial, flyby_angle_end, flyby_angle_step_size):

    """ Creates txt-file of the used parameters. """

    print("\n\nStarted saving parameters.\n")
    print("...")

    path = "Results/NEW_RESULTS/parameters.txt"
    file = open(path, 'w')
    
    file.write('\n################################################################################################## \n')
    file.write('########################################### Parameters ########################################### \n')
    file.write('################################################################################################## \n')

    if flag_type_of_simulation == 1:
        type_of_simulation = "Variable Time of Flyby - Fixed Distance of Flyby - Fixed Entry Angle of Flyby"
    elif flag_type_of_simulation == 2:
        type_of_simulation = "Variable Distance of Flyby - Fixed Time of Flyby - Fixed Entry Angle of Flyby"    
    elif flag_type_of_simulation == 3:
        type_of_simulation = "Variable Entry Angle of Flyby - Fixed Distance of Flyby - Fixed Time of Flyby"

    file.write('\n\tType of Simulation: \t\t\t\t\t\t' + str(type_of_simulation))
    file.write('\n\n\tDeparture Time: \t\t\t\t\t\t\t' + str(departure_time))

    file.write('\n\tFlyby Time Initial: \t\t\t\t\t\t' + str(flyby_time_initial))
    file.write('\n\tFlyby Time Iterations: \t\t\t\t\t\t' + str(flyby_time_iterations))
    file.write('\n\tFlyby Time Step Size [days]: \t\t\t\t' + str(flyby_time_step_size))

    file.write('\n\tFlyby Height Initial [km]: \t\t\t\t\t' + str(flyby_height_initial))
    file.write('\n\tFlyby Height Iterations [km]: \t\t\t\t' + str(flyby_height_iterations))
    file.write('\n\tFlyby Height Step Size [km]: \t\t\t\t' + str(flyby_height_step_size))

    file.write('\n\tFlyby Entry Angle Initial [rad]: \t\t\t' + str(flyby_angle_initial))
    file.write('\n\tFlyby Entry Angle End [rad]: \t\t\t\t' + str(flyby_angle_end))
    file.write('\n\tFlyby Entry Angle Step Size [rad]: \t\t\t' + str(flyby_angle_step_size))

    print("\nFinished saving parameters.\n")


# ----------------------------------------------------------------------------


def save_parameters_in_txt_multi_varying(initial_departure_time, departure_time_stepsize, departure_time_number_of_iterations, initial_flyby_time, flyby_time_stepsize, flyby_time_number_of_iterations, initial_height_flyby, end_height_flyby, height_flyby_stepsize, initial_flyby_angle, end_flyby_angle, flyby_angle_stepsize):
    """ Creates txt-file of the used parameters. """

    print("\n\nStarted saving parameters.\n")
    print("...")

    path = "Results/NEW_RESULTS/parameters.txt"
    file = open(path, 'w')
    
    file.write('\n################################################################################################## \n')
    file.write('########################################### Parameters ########################################### \n')
    file.write('################################################################################################## \n')

    file.write('\nParameters of the Iterations:')

    file.write('\n\n\tInitial Departure Time: \t\t\t\t' + str(initial_departure_time))
    file.write('\n\tDeparture Time Step Size: \t\t\t\t' + str(departure_time_stepsize))
    file.write('\n\tDeparture Time Number of Iterations: \t' + str(departure_time_number_of_iterations))

    file.write('\n\n\tInitial Flyby Time: \t\t\t\t\t' + str(initial_flyby_time))
    file.write('\n\tFlyby Time Step Size: \t\t\t\t\t' + str(flyby_time_stepsize))
    file.write('\n\tFlyby Time Number of Iterations: \t\t' + str(flyby_time_number_of_iterations))

    file.write('\n\n\tInitial Flyby Height: \t\t\t\t\t' + str(initial_height_flyby))
    file.write('\n\tEnd Flyby Height: \t\t\t\t\t\t' + str(end_height_flyby))
    file.write('\n\tFlyby Height Step Size: \t\t\t\t' + str(height_flyby_stepsize))

    file.write('\n\n\tInitial Flyby Angle: \t\t\t\t\t' + str(initial_flyby_angle))
    file.write('\n\tEnd Flyby Angle: \t\t\t\t\t\t' + str(end_flyby_angle))
    file.write('\n\tFlyby Angle Step Size: \t\t\t\t\t' + str(flyby_angle_stepsize))

    print("\nFinished saving parameters.\n")