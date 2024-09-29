# Trajectory Optimization with Python

This simulation should help to find and analyze the _best_ trajectory of a spacecraft to Saturn.
During the flight to Saturn, the spacecraft will perform a flyby maneuver at Jupiter.


## Structure of the project
The project is structure as follows:
* **traj_optimization_single.py**: This file executes one single simulation run for the assigned parameters of departure date, flyby date and flyby distance. The main results are printed in the console and a plot of the trajectories is created and shown.
* **traj_optimization_multi_flybyDistance.py**: This file executes multiple simuation runs with different distances of the s/c to Jupiter at the flyby. It takes the assigned parameters of the departure date and flyby date and varies the flyby distance during each run according to the parameters _initial\_height\_flyby_, _iterations\_step\_size\_km_ and _number\_of\_iterations_. The results are saved in .csv-file and the used parameters are saved in a .txt-file in the folder _Results/NEW\_RESULTS_.
* **traj_optimization_single_flybyDate.py**: This file executes multiple simuation runs with different flyby dates at Jupiter. It takes the assigned parameters of the departure date and flyby distance and varies the flyby date during each run according to the parameters _flyby\_time_, _iterations\_step\_size\_day_ and _number\_of\_iterations_. The results are saved in .csv-file and the used parameters are saved in a .txt-file in the folder _Results/NEW\_RESULTS_.
* **pipeline.py**: This file contains the pipeline of a single run. This pipeline gets called during the single run simulation as well as during the multi run simulation.
* **file_handling.py**: This file contains the functions for saving the results and the used parameters of the multi run simulations.


## How to setup the simulation
### TO-DO


# How to execute the simulation
### TO-DO