%% Set up code for running mice

addpath('/INSERT_PATH_TO_MICE_HERE/mice/src/mice/')
addpath('/INSERT_PATH_TO_MICE_HERE/mice/lib/' )
cspice_furnsh('/INSERT_PATH_TO_MICE_HERE/mice/naif0012.tls')
cspice_furnsh('/INSERT_PATH_TO_MICE_HERE/mice/de440.bsp')
disp("MICE Setup Complete")
