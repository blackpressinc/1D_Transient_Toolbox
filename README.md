# 1D_Transient_Toolbox
A python script for general purpose 1D transient thermal analysis including ablation and radiation

This script was designed to help me grow an intuition about transient thermal problems. It will also serve as a test bed for implimenting uncertainty calculations in an explicit forward euler solver. It can calculate the transient and psuedo-transient response of any number of material layers in a 1D stackup. Current boundary conditions are

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

*Note: This tool is not related to my work at Firefly Aerospace, was developed in my own free time, and serves as a test bed for ideas. I do not garentee any accuracy. If you are looking to use this to solve real world problems I would highly suggest NASA FIAT. Only a fool would use free, untested software on the internet if they don't know what they are doing!