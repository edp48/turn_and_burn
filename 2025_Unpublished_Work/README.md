This is work done with the intent of publication, though at this time that has not been achieved. 

This work seeks to build off of the work that was published in "Optimal Cislunar Trajectories With Continuous, High-Thrust Nuclear-Thermal Propulsion" by expanding to use a more robust collocation optimization program. Here, we used Matthew Peter Kelly's OptimTraj collocation toolkit found at: https://zenodo.org/badge/latestdoi/40544279, and the accompanying paper at https://doi.org/10.1137/16M1062569. 

Goals of this work:
- Examine the use cases for a turn-and-burn style maneuver in the context of raising an orbit
- Compare factors such as Delta V and Time of Flight to examine if there us utility in this method over standard Hohmann Transfer style maneuvers

Current State:
 - Optimization is working and produces results. This requires that you NAIF's spice toolkit installed, as well as OptimTraj and Chebfun downloaded (instructions for the Chebfun are in Matthew Kelly's documentation

- There are 4 scripts in the folder that are useful. The first is a Lambert Solver set to run for the same set of parameters as the turn-and-burn style case. The second "Orbit_Raising_Main" provides a one shot optimization of a single trajectory. The last two, "Orbit_Raising_Cont_set1" and "Orbit_raising_cont_set2" do continuation analysis on the system, and will increase the thrust level of the spacecraft to start to assemble solutions in the thrust phase space. The two files are for two different thrust levels

*Notes*
- This is a work in progress, and the Lambert solver output on the Delta V comparison graphs are not optimal at the moment, as that DV cost is higher than it should be. 

This work will be updated as Erik has time to work on it. 

