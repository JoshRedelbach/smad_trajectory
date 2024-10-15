# Interplanetary Trajectory Optimization with Python

This project includes simulations on designing and optimizing ...
* ... an interplanetary trajectory to Saturn performing a flyby at Jupiter,
* ... the entry trajectory when arriving at Saturn,
* ... the transfer trajectory from the entry trajectory to the final working orbit,
* ... the exit trajectory after successfully achieving the mission's objectives,
* ... the propulsion system.


## Development Details
  - Name: _Interplanetary Trajectory Optimization performing Flybys_
  - Current version number: 2.0.0
  - Date of initial coding: 2024-09-25
  - Date of first release (verison 1.0.0): 2024-09-29
  - Name of author: Joshua Redelbach, contact: <joshua.redelbach@t-online.de>
  - A list of all required packages, dependencies, the tested versions and how to install them is given in the section _How to Setup the Simulation_
  - Inputs and outputs are described in the section _Structure of the Project_ as well as the expected range of validity and the tested range.
  - The algorithms used in this project are indicated in the section _Idea_
  - A detailed description of each simulation that can be executed is given in the section _Structure of the Project_

## Parts
This project contains three different simulations:
1. Finding an interplanetary trajectory to Saturn.
2. Finding an entry, transfer and exit trajectory at Saturn.
3. Estimating the required mass of the propellant.

For details on the simulations of Part 1, please continue reading this README.md file. It presents the idea, the structure, how to setup and execute the project as well as the results regarding the final interplanetary trajectory.

Part 2 and Part 3 are presented in the form of Jupyter notebooks and can be found in the folders [code_saturn_entry](/code_saturn_entry/) and [code_mass_calculations](/code_mass_calculations/) respectively. As the notebooks contain comments and explainations for the entire simulation, no details on the idea behind the simulation is given in this README.md file. Furthermore, all notebooks are already executed and the results are uploaded as well, it is not needed to run the notebooks locally. Both folders include two different notebooks, as two different concepts for the exit trajectory are developed. Therefore, the required $\Delta v$ changes and thus the mass calculations need to be adapted accordingly. The two notebooks for calculating the mass have the exact same structure, but just different values inserted corresponding to the different exit trajectories. The notebooks of Part 3 can be locally executed as soon as a working Python interpreter is installed. In order to locally execute the notebooks of Part 2, the same requirements as for Part 1 are required which are introduced later in this README.md file.


## Idea
This simulation was developed for the course _Space Mission Analysis and Design_ of the _Instituo Tecnico Superior_ of Lisbon, Portugal. The task was to design a interplanetary space mission. The spacecraft should depart from a Low Earth Orbit with an altitude of 200 km with the goal of reaching Saturn in a desired working orbit. An optimal trajectory should be found with respect to required delta-v and duration of the flight. Based on previous space missions, e.g. Voyager 1 and 2, it was decided to use a gravity-assist flyby around Jupiter.

The problem was simplified as following. The Lambert problem is solved for the a given time of departure and a given time of flyby at Jupiter. For those dates the positions of the center of Earth and Jupiter are used for the initial and end position of the spacecraft respectively as input for the Lambert problem. For solving the Lambert problem the implemented solver in the package `poliastro` is used whose documentation can be found [here](https://docs.poliastro.space/en/stable/autoapi/poliastro/iod/izzo/index.html#poliastro.iod.izzo.lambert). The algorithm is based on the paper of [Izzo (2014)](/references/Revisiting%20Lambert's%20Problem%20(Izzo,%202014).pdf). It assumes that the spacecraft is only in the gravitational influence of the Sun (patched-conic approach). This algorithm outputs the required initial velocity at Earth and the velocity with which the spacecraft arrives at Jupiter, both with respect to the sun.

Based on the initial velocity and the velocity of the Earth with respect to the sun at the time of departure, the $v_{\inf}$ is computed relative to the Earth. Based on $v_{\inf}$ and the velocity a spacecraft has in a LEO with an altitude of 200 km, delta-v is calculated which needs to be applied to the spacecraft in the LEO.

Based on the velocity with which the spacecraft arrives at Jupiter, a given entry angle of the flyby and a given closest distance of the spacecraft to Jupiter during the flyby, the velocity is calculated with which the spacecraft is leaving Jupiter. For solving this, another algorihtm of the `poliastro` package is used whose documentation can be found [here](https://docs.poliastro.space/en/stable/autoapi/poliastro/core/flybys/index.html#poliastro.core.flybys.compute_flyby). Again, the patched-conic method is used, so only the influence of Jupiter is considered. 

Based on the obtained velocity of the spacecraft with respect to the sun, the new orbit after the flyby is propagated for 3 years starting from an approximated position of the spacecraft at the center of the Jupiter. For this task, the implemented propagator of the `poliastro` package is used, which is based on the work of [Farnocchi (2013)](/references/Robust%20resolution%20of%20Kepler’s%20equation%20in%20all%20eccentricity%20regimes%20(Farnocchi,%202013).pdf). The documentation can be found [here](https://docs.poliastro.space/en/stable/autoapi/poliastro/twobody/propagation/index.html). During the entire propagation the distance of the spacecraft to the Saturn is calculated every hour. The minimum distance and the respective time stamp is saved. If the minimum distance of the spacecraft to Saturn is less than 100 000 km (requirement given by the professor), we consider our trajectory as successfull.

So in summary, the following assumptions are used:
  - For solving the lambert problem, start and end position of the spacecraft is the center of Earth and of Jupiter respectively.
  - The flyby is assumed to be an instantaneous maneuver, so no time pasts between arriving and leaving Jupiter.


## Structure of the Project
The project is structured as follows. There are in total 5 different simulations that can be executed. Each of them are named in the format `traj_optimization_XXX.py`. One of them can be used to simulate and analyze a single simulation run, so a specific trajectory, in detail:
* [`traj_optimization_single.py`](/code_traj_optimization/traj_optimization_single.py): This file executes one single simulation run for the assigned parameters of departure date, flyby date, flyby distance and flyby angle. The main results are printed in the console and a 2D plot of the trajectories is created and shown. No results are saved, so make sure to save results you need manually by your own.

The other four can be used to test multiple different trajectories by iterating over a range of the parameters:
* [`traj_optimization_multi_flybyDistance.py`](/code_traj_optimization/traj_optimization_multi_flybyDistance.py): This file executes multiple simuation runs with different distances of the spacecraft to Jupiter at the flyby. It takes the assigned parameters of the departure date, flyby date and entry angle of flyby and varies the flyby distance during each run according to the parameters _initial\_height\_flyby_, _iterations\_step\_size\_km_ and _number\_of\_iterations_.
* [`traj_optimization_multi_flybyDate.py`](/code_traj_optimization/traj_optimization_multi_flybyDate.py): This file executes multiple simuation runs with different flyby dates at Jupiter. It takes the assigned parameters of the departure date, flyby distance and entry angle of flyby and varies the flyby date during each run according to the parameters _flyby\_time_, _iterations\_step\_size\_day_ and _number\_of\_iterations_.
* [`traj_optimization_multi_flybyEntryAngle.py`](/code_traj_optimization/traj_optimization_multi_flybyEntryAngle.py): This file executes multiple simuation runs with different entry angles of the flyby at Jupiter. It takes the assigned parameters of the departure date, flyby date and flyby distance and varies the flyby entry angle during each run according to the parameters _initial\_angle_, _iterations\_step\_size\_rad_ and _end\_angle_.
* [`traj_optimization_multi_all.py`](/code_traj_optimization/traj_optimization_multi_all.py): This file executes multiple simulation runs with different departure dates, flyby dates, flyby distances and flyby entry angles.

Details on the input parameters (description, unit, range, ...) is given as a comment next to the assignement in the code.

For each of these multi-run simulations, a .csv-file is created containing the parameters and the results of each run. Additionally, a .txt-file is created summarizing the parameters selected for running this simulation, so you can recreate the results accordingly.

All of the mentioned executable simulations are based on the following two `python` files.
* [`pipeline.py`](/code_traj_optimization/pipeline.py): This file contains the pipeline of a single run. This pipeline gets called during the single run simulation as well as during the multi run simulation.
* [`file_handling.py`](/code_traj_optimization/file_handling.py): This file contains the functions for saving the results and the used parameters of the multi run simulations.

**IMPORTANT**: The runtime of each run can last multiple seconds (executed on a MacBook Pro with a M1 chip: ~ 20-30s). So running the multi-run simulations, make sure that the range of the parameters is not to large. Otherwise, the simulation will take a very long time.

## How to Setup the Simulation
This simulation is only based on `Python`. Thus, to be able to setup the simulation make sure that you have installed a working Python interpreter and the `pip` python package manager.

**IMPORTANT**: as the package _astropy_ is required and as this package only runs with `Python 3.8` - `Python 3.10`, make sure that you have the correct version of Python installed. This simulation was developed and tested with `Python 3.9.13` and `pip3`.

First, clone the repository into the desired working directory on your local machine:
```console
cd to/your/working/directory
git clone https://github.com/JoshRedelbach/smad_trajectory.git
```
If you are not familiar with _git_, you can also download the project files clicking on the button **<>Code** -> **Download ZIP**. Unpack the ZIP-folder in the desired local working directory.


Then, create a virtual environment (venv):
```console
python3 -m venv name_of_your_virtual_environment_folder
```
Make sure the venv is activated, so that all packages are installed in the venv and not on the entire machine:
```console
source name_of_your_virtual_environment_folder/bin/activate
````

Afterwards, install the following packages using `pip`:
| Package | Command for Installation | Tested Version |
| -------- | ------- | ------- |
| numpy | ```pip3 install numpy``` | 1.26.4
| astropy | ```pip3 install astropy``` | 5.3.4
| jplephem | ```pip3 install jplephem``` | 2.22
| matplotlib | ```pip3 install matplotlib``` | 3.9.2
| scipy | ```pip3 install scipy``` | 1.13.1
| numba | ```pip3 install numba``` | 0.60.0
| plotly | ```pip3 install plotly``` | 5.24.1
| poliastro: | ```pip3 install poliastro``` | 0.17.0

You can check the installed versions of your packages using ```pip3 list```.
As the poliastro and astropy packages have a few more dependencies that are downloaded and installed automatically, the list will be way longer than the table shown above. This fact can be ignored.

If all of these steps are performed successfully, the simulation is set up and ready for execution.

## How to Execute the Simulation
Each of the mentioned executable simulation files includes the following sections:
1. **Notes**: a list of aspects that have to be considered and that have to be kept in mind to correctly interprete the results.
2. **Imports**: all necessary imports of the python packages required for executing the file.
3. **Parameters**: definition and assignement of all parameters which can be chosen in order to find the optimal/desired trajectory.
4. **Code**: Includes the code. Nothing has to be done here. Adjust it only if you really know what you are doing.

In order to execute the file, assign the parameters that should be tested. For all multi-run simulations, make sure that the folder [results/NEW_RESULTS](/results/NEW_RESULTS/) is empty, so that no previous results are overwritten.


## Approach and Summary of the Results
For the project of determining the optimal orbit from the LEO to Saturn with respect to delta-v and duration of the flight, the following approach was performed. First, the earliest date of departure was decided to be 2035 due to the long required preparation of such an interplanetary space mission. Then, by manually trying out different dates of departure and of flyby using [`traj_optimization_single.py`](/code_traj_optimization/traj_optimization_single.py) a promising range of dates of departure and of flyby were selected. For those dates, a multi-run simulation was executed varying these dates within the range using [`traj_optimization_multi_all.py`](/code_traj_optimization/traj_optimization_multi_all.py) (results can be seen [here](/results/project_results/selecting_dates)). After, finding the most promising dates, promising ranges for the flyby distance and for the flyby entry angle were found by manually trying out using [`traj_optimization_single.py`](/code_traj_optimization/traj_optimization_single.py). Using [`traj_optimization_multi_all.py`](/code_traj_optimization/traj_optimization_multi_all.py), two multi-run simulations were executed to firstly roughly find the best values for the flyby distance and entry angle and then to secondly refine these values (results can be seen [here](/results/project_results/selecting_distance_and_angle_1/) and [here](/results/project_results/selecting_distance_and_angle_2/)).
This resulted in the following parameters:

* date of departure: 2037-09-25
* date of flyby: 2039-07-29
* date of arrival at Saturn: 2041-09-14
* height of flyby: 1 159 500 km
* entry angle of flyby: -0.10265 rad ~ 5,8814 deg
* required delta-v: 6.848 km/s

These parameters are already assigned in the single run simulation. So to get insights of this final trajectory, the file [`traj_optimization_single.py`](/code_traj_optimization/traj_optimization_single.py) can be executed and details are printed to the console and the 2D plot of the involved orbits is created.

![Plot of the orbits involved](/results/project_results/2D_plot_final_orbit.png "2D Orbit Plot")