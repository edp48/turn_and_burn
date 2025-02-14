%% Set up code for running mice

if ismac
    addpath('/Users/erikpayton/Documents/mice/src/mice/')
    addpath('/Users/erikpayton/Documents/mice/lib/' )
    cspice_furnsh('/Users/erikpayton/Documents/mice/naif0012.tls')
    cspice_furnsh('/Users/erikpayton/Documents/mice/de440.bsp')
    disp("MICE Setup Complete for Macintosh")

elseif ispc
    addpath('/Users/geoep/Documents/mice/src/mice/')
    addpath('/Users/geoep/Documents/mice/lib/' )
    cspice_furnsh('/Users/geoep/Documents/mice/naif0012.tls')
    cspice_furnsh('/Users/geoep/Documents/mice/de440.bsp')
    disp("MICE Setup Complete for Windows")
end 

