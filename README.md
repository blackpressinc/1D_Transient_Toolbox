Test

# 1D_Transient_Toolbox
A python script for general purpose 1D transient thermal analysis including ablation and radiation

This script was designed to help me with my day-to-day work as a thermal engineer in Aerospace. It can calculate the transient and psuedo-transient response of any number of material layers in a 1D stackup. Current boundary conditions are

Steady/Transient
- Conduction
- Convection
- Radiation
- Fixed Heat Flux

Ablation can be toggled on/off

Static boundary conditions are currently handled within the program while Transient boundary conditions are read from an excel file. 

Material properties are stored as a function of temperature and pressure as polynomial equations which effectively makes them optional. Most materials are effected somewhat by temperature but aerogels are really affected by pressure. Material properties are re-calculated at a user defined interval to relieve some overhead on the solver

This solver is just an explicit forward Euler so it's second order in space and first order in time. It introduces much less error than the uncertainty of any given boundary condition in a real scenario.

I'll be finishing up comments and adding a much more in depth tutorial soon.