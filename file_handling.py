##############################################################################
########################## Trajectory Optimization ###########################
##############################################################################
############################### File Handling ################################
##############################################################################



def save_data_in_csv(departure_time_list, flyby_time_list, distance_of_flyby_list, delta_v_list, distance_scToSaturn_crossing_list, arrival_time_list, flag_saturn_crossed_list):
    """ Save the collected data of the multiple runs in a .csv-file. """

    print("\n\nStarted saving results.\n")
    print("...")

    path = "Results/NEW_RESULTS/data.csv"
    file = open(path, 'w')
    
    # Define coloumns
    file.write('Departure Time,Flyby Time,Arrival Time,Distance of Flyby to Jupiter,Delta V,Status Saturn crossed,Crossing distance S/C-Saturn\n')

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
        file.write('\n')
    
    file.close()

    print("\nFinished saving results.\n")


# ----------------------------------------------------------------------------


def save_parameters_in_txt(flag_type_of_simulation, departure_time, height_flyby, number_of_iterations, flyby_time, iterarion_step_size):
    """ Creates txt-file of the used parameters. """

    print("\n\nStarted saving parameters.\n")
    print("...")

    path = "Results/NEW_RESULTS/parameters.txt"
    file = open(path, 'w')
    
    file.write('\n################################################################################################## \n')
    file.write('########################################### Parameters ########################################### \n')
    file.write('################################################################################################## \n')

    if flag_type_of_simulation == 1:
        type_of_simulation = "Variable Time of Flyby - Fixed Distance of Flyby"

        file.write('\n\tType of Simulation: \t\t\t\t' + str(type_of_simulation))
        file.write('\n\n\tDeparture Time: \t\t\t\t\t' + str(departure_time))
        file.write('\n\tIniial Flyby Time: \t\t\t\t\t' + str(flyby_time))
        file.write('\n\tDistance of Flyby: \t\t\t\t\t' + str(height_flyby))
        file.write('\n\n\tNumber of Iterations: \t\t\t\t' + str(number_of_iterations))
        file.write('\n\n\tIteration Step Size [days]: \t\t\t' + str(iterarion_step_size))
    
    elif flag_type_of_simulation == 2:
        type_of_simulation = "Variable Distance of Flyby - Fixed Time of Flyby"

        file.write('\n\tType of Simulation: \t\t\t\t\t' + str(type_of_simulation))
        file.write('\n\n\tDeparture Time: \t\t\t\t\t\t' + str(departure_time))
        file.write('\n\tFlyby Time: \t\t\t\t\t\t\t' + str(flyby_time))
        file.write('\n\tInitial Distance of Flyby: \t\t\t\t' + str(height_flyby))
        file.write('\n\n\tNumber of Iterations: \t\t\t\t\t' + str(number_of_iterations))
        file.write('\n\n\tIteration Step Size [km]: \t\t\t\t' + str(iterarion_step_size))

    print("\nFinished saving parameters.\n")
