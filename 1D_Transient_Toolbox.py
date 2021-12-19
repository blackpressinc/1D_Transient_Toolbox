# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:51:15 2021

@author: fraser.jones
"""
import time #time library is just used to track the time per iteration of the solver and provide an estimate of the run length to the user
import numpy as np #the bulk of this code is written in numpy to make it compatible with Uncertainties.py in the future
import pandas as pd #Pandas is just used to read excel files because trying to do it in numpy or default python made me want to die
from scipy import interpolate #Interpolate was used to convert transient data provided at discreet time intervals into functions that could
#could be called by the user to get a value. So you give it time and T_wall and it gives you a function f(time) --> T_wall
#I may have to ditch this in the future because its not compatible with uncertainties
import matplotlib.pyplot as plt 
#import CSV_to_Python #Maybe use in the future but I don't want to have more modules
#%%

#Global Conditions
# =============================================================================
# time_final is a essentially when the simulation ends. Inside the solver loop there's a variable called "sim_time" that just adds up the time_step
# each iteration. Once that variable exceeds the time_final, the simulation ends. You may want to add 1 second to the final time just because the 
# once the time_final is exceeded by sim_time the simulation ends /before/ it gets to the data saving loop. So it might not save the last data point
# =============================================================================
time_final = 41.0 #[s]
plotting_interval = 10.0 #[s] Data will ***saved*** be plotted at this interval. Note that this controls when data is saved
T_initial = 274.0 #[K] All cells in the simulation will be initialized with this value. Once the Temperature_array has been created, you can add a custom temperature function or something
# but I have not yet coded that
P_initial = 1.0 #[atm] #pressure needs to be re-worked so that it can be made a function of time. This just acts as a place holder
steffan_boltzman_constant = 5.67*10**-8 #used for radiation calculations
# =============================================================================
# This solver supports transient material properties where available but calculating them each iteration makes the solver incredibly slow
# So you can set a time to recalculate the material properties every X seconds
# =============================================================================
Material_Properties_Calculation_Interval = 0.1 #[s]

# =============================================================================
# This is the actual material stackup you are simulating. Materials are assigned left to right based on the thicknesses. 
# Material list is just strings that have to match with the some name in the first row of the Material_Properties.xlsx
# Thickness is the total thickness of each material in millimeters (I will never use english to define a property god help me)
# =============================================================================
Material_list = ["T792","A12","P50"]
Thickness = np.array([3.1,0.2,10]) #[mm] Array defining the thickness of the individual layers

# =============================================================================
# These parameters actually setup the simulation mesh. Whichever is the tightest requirement will be used
# You can make a layer as thin as you want but timestep gets smaller as a function of the smallest cell thickness squared dt ≈ 1/(dx^2)
# =============================================================================
minimum_cell_count_per_layer = 4 #Thinnest Layer will be used to calculate dx for the mesh
maximum_cell_size = 0.1 #[mm] Cells will not be made larger than this value

#%%
# =============================================================================
# This region sets up the boundary conditions for the simulation. I tried to make it as flexible as possible and I succeeded
# That said, flexible code makes for dumb physics. You can do so much non-physical stuff but it's really fun.
# I'll describe each section as we get there but the simulation is basically

# (Left Boundary Condition) |------ your mesh of materials -------| (Right Boundary Condition)

# And boundary conditions can either be static which you can define here in python or transient
# which are defined in excel documents like Transient_Boundary_Conditions_for_Alpha_Flight_S1.xlsx

# Weirdly enough this solver can't actually do steady state. For that you can just set each material thickness to 1 cell
# And then set the final time to be very long. At that point, the timestep will be large and it should get you close to steady state

#Special note: I used _left and _right to make coding easier
# =============================================================================

# =============================================================================
# You can turn Ablation on/off with this boolean
# =============================================================================
Ablation = True

# =============================================================================
# This section defines the static left boundary conditions

# T_wall_left is used as a reference for pure conduction. So if you left boundary condition was just a fixed temperature wall

# h_left and T_conv_left are used for convection

# emissivity_left defines the outgoing radiation effectiveness of the left surface. This is seperate from incoming radiation
# absorptivity_left defines the effectiveness of aborption of radiation from either of the two default sources
# View_Factor_1_left and View_Factor_2_left are coupled and sum to one as default but there's nothing that stops you from
# adding totally nonsensical values so have fun. This section should be expanded in the future to just use arrays so you can define
# complex radiative networks but I didn't feel like including bouncing radiation
# T_rad_1_left and T_rad_2_left set the temperatures of the radiation sources sending energy into the left boundary

# q_in_left just lets you set an arbitrary heat flux boundary condition. Useful for heaters or just messing around

# Note that transient boundary conditions is essentially a copy of these boundary conditions but as a function of time
# =============================================================================
#Left Boundary Conditions
T_wall_left = 300 #[K]
h_left = 10 #[W/m²·K]
T_conv_left = 273 #[K]
emissivity_left = 0.03
absorptivity_left = 0.03
View_Factor_1_left = 1.0
View_Factor_2_left = 1 - View_Factor_1_left
T_rad_1_left = 270 #[K]
T_rad_2_left = 270 #[K]
q_in_left = 250 #[W/m²]

# =============================================================================
# This section of Booleans turns the various sources on/off. You can mix and match whatever you want even if it makes no sense
# =============================================================================
Conduction_left = False
Convection_left = True
Radiation_to_Boundary_left = False
Radiation_from_Boundary_left = False
Fixed_Heatflux_left = False

Transient_Conduction_left = False
Transient_Convection_left = True
Transient_Radiation_to_Boundary_left = False
Transient_Radiation_from_Boundary_left = False
Transient_Fixed_Heatflux_left = False

# =============================================================================
# If you are loading in transient boundary conditions, you'll need to point them to an excel file somewhere in the directory you are running this script from
# It shouldn't run if all of the Transient_ booleans are set to false
# =============================================================================
#Left Transient Boundary Conditions
if (Transient_Conduction_left+Transient_Convection_left+Transient_Radiation_to_Boundary_left+Transient_Radiation_from_Boundary_left+Transient_Fixed_Heatflux_left > 0): #Test if any transient boundary conditions are set to True (1)
    
# =============================================================================
# This code was added to address a comment by FrancoisGLEYZON where an IDE would not recognize the transient boundary conditions and throw and error
# Normally the globals() code will create the necessary variables at run time which also allows for creating new boundary condition variables not originally included
# in the code, but until I can find a way to do this that doesn't anger Sypder this will serve as a good stand-in
# =============================================================================
    Transient_T_wall_left = 0 #[K]
    Transient_h_left = 0 #[W/m²·K]
    Transient_T_conv_left = 0 #[K]
    Transient_emissivity_left = 0
    Transient_absorptivity_left = 0
    Transient_View_Factor_1_left = 0
    Transient_View_Factor_2_left = 0
    Transient_T_rad_1_left = 0 #[K]
    Transient_T_rad_2_left = 0 #[K]
    Transient_q_in_left = 0 #[W/m²]
    
    Left_Transient_BC = pd.read_excel('Example_Transient_Boundary_Conditions.xlsx',index_col=0,skiprows=1,header=1) #If transient boundary conditions are being used, load
    # an excel sheet, skipping the first row and using the second row as the headers for a pandas database. The first column is used as the index
    # The example transient boundary conditions should be Transient_Boundary_Conditions_for_Alpha_Flight002_Hotfire_S1.xlsx

#This move is illegal in 9 countries
    for col in Left_Transient_BC.columns: #for the number of columns (basically each variable you want to create) do the following
        globals()["Transient_"+str(col)+"_left"]=interpolate.interp1d(Left_Transient_BC.index.to_numpy(dtype=float),Left_Transient_BC[col].to_numpy(dtype=float))
# =============================================================================
# This section of code is a little tricky but it works as follows
# globals() creates a new global variable the name "Transient_" + Whatever the title of the column is in the xlsx" and then "_left"
# Using this format, I was able to copy all of the code I had already written for the static boundary conditions
# The right side of this uses scipy's interpolate to create a function that can be called like T_wall(current_time)
# Left_Transient_BC.index.to_numpy(dtype=float) is the first column of the excel sheet so it's the time
# Left_Transient_BC[col].to_numpy(dtype=float) is the variable like T_wall, h, etc.
# The interpolate.interp1d just creates a function where you can add any time and get the value of the function
# =============================================================================

# =============================================================================
# The Right Boundary Conditions section is an exact copy of the Left_Boundary_Conditions section but with the words "right"
# =============================================================================
#Right Boundary Conditions

T_wall_right = 300 #[K]
h_right = 25 #[W/m²·K]
T_conv_right = 300 #[K]
emissivity_right = 0.8
absorptivity_right = 0.8
View_Factor_1_right = 0.5
View_Factor_2_right = 1 - View_Factor_1_right
T_rad_1_right = 2200 #[K]
T_rad_2_right = 300 #[K]
q_in_right = 17500 #[W/m²]

Conduction_right = False
Convection_right = True
Radiation_to_Boundary_right = True
Radiation_from_Boundary_right = True
Fixed_Heatflux_right = False

Transient_Conduction_right = False
Transient_Convection_right = False
Transient_Radiation_to_Boundary_right = False
Transient_Radiation_from_Boundary_right = False
Transient_Fixed_Heatflux_right = False

if (Transient_Conduction_right+Transient_Convection_right+Transient_Radiation_to_Boundary_right+Transient_Radiation_from_Boundary_right+Transient_Fixed_Heatflux_right > 0):
#Right Transient Boundary Conditions

# =============================================================================
# Copy of the code in the _left transient boundary condition to address an issue for FrancoisGLEYZON
# =============================================================================
    Transient_T_wall_right = 0 #[K]
    Transient_h_right = 0 #[W/m²·K]
    Transient_T_conv_right = 0 #[K]
    Transient_emissivity_right = 0
    Transient_absorptivity_right = 0
    Transient_View_Factor_1_right = 0
    Transient_View_Factor_2_right = 0
    Transient_T_rad_1_right = 0 #[K]
    Transient_T_rad_2_right = 0 #[K]
    Transient_q_in_right = 0 #[W/m²]
    
    Right_Transient_BC = pd.read_excel('Transient_Boundary_Conditions_for_Alpha_Flight_S1.xlsx',index_col=0,skiprows=1,header=1)

    for col in Right_Transient_BC.columns:
        globals()["Transient_"+str(col)+"_right"]=interpolate.interp1d(Right_Transient_BC.index.to_numpy(dtype=float),Right_Transient_BC[col].to_numpy(dtype=float))

#%%
# =============================================================================
#This section contains material data in the form of Temperature and pressure polynomials
#The Temperature polynomial sets the value of the material properties
#The Pressure polynomial just scales the value of the already calculated material properties
#
#[K(T),rho(T),Cp(T)] = a1+a2*T+a3*T^2+a4*T^3+a5*T^5'
#[K(T,P),rho(T,P),Cp(T,P)] =[K(T),rho(T),Cp(T)]*(b1+b2*P+b3*P^2+b4*P^3+b5*P^5)
#
#In the future this should just be loaded in as a default materials database stored in excel
# =============================================================================

Materials = pd.read_excel('Material_Properties.xlsx',index_col=(0),usecols=("C:ZZ"),skiprows=(0)) #Note that you can use "usecols" to just select a subset of 
# columns or rows in microsoft excel. So in this case, it's skipping the first row, using the first column as the index, but the first column
# is actually column C and it goes to column W. So columns A and B are just ignored

# Same as for the transient conditions, but in this case just starting at column C and going until all the materials are loaded
for i in range(0,Materials.shape[1]):
    globals()[Materials.iloc[:,i].name]=Materials.iloc[:,i].to_numpy(copy=True) #globals is the python equivalent of excel's "Indirect"
    #i don't know why but copy=true was necessary for some reason

# =============================================================================
# This code calculates the required number of cells to resolve the layers based 
# on the maximum_cell_size or the minimum_cell_count_per_layer
# =============================================================================

cell_count = (Thickness/maximum_cell_size).astype(int) # Check what the expected cell count would be for each material using the maximum cell size and thickness
if min(cell_count) < minimum_cell_count_per_layer: # If the previously calculated cell counts for each material is lower than the min cell count
    cell_count = Thickness/(min(Thickness)/minimum_cell_count_per_layer) # Then calculate a new thickness such that the minimum cell count for each material is larger than the minimum

# Now that cell count has met the minimum standards, use it and the total material thickness to calculate the cell thickness "x_step"
x_step = sum(Thickness)/sum(cell_count)/1000 #note that x-step is in meters. Makes math easier


# =============================================================================
# This code creates a list of materials equal to the total number of cells calculated previously
# and fills in those values based on the material stackup from left to right
# =============================================================================
cell_materials = [] # create an empty list that will store the material for each cell
for i in range(cell_count.size): # count form 0 to the total number of cells calculated previously
    for j in range(0,cell_count[i].astype(int)): # for each cell I just created
        cell_materials.append(Material_list[i]) # find the corresponding material from the materials list and 
#%%

# =============================================================================
# This section creates functions to quickly calculate the K, rho, and Cp using the materials name
# Temperature, and pressure. 
# 
# I also have functions that calculate the temperature gradient using second order finite forward euler differencing
# For more information I recommend checking out this reference
# http://hplgit.github.io/num-methods-for-PDEs/doc/pub/diffu/sphinx/._main_diffu001.html
#
# I need to add in reference pressure and temperature for materials in the future. Otherwise T and P are a bit abstract
# =============================================================================

# =============================================================================
# Functions for K, rho, and Cp follow the same pattern. The material name, temperature, and pressure is sent into the function
# Based on the "name" parameter, the materials which have already imported are referenced for their temperature and pressure
# coefficients. Those coefficients are used to calculate the new material property which is returned.
# There is probably a better way to do this than cell-by-cell but for now it works fine
# =============================================================================

def K(name, Temperature, Pressure): # define a function "K" that takes a material "name", "Temperature", and "Pressure"
    index = 0 # Create a variable named "index" and set it equal to zero. Since my material properties are stored as 5 equation polynomials, this will count from 0, 1, 2, 3, 4 for reading each polynomial constant
    K_as_a_function_of_Temperature = 0 # Create a variable named "K_as_a_function_of_Temperature" and set it equal to zero. This will store the value of K when we calculate it
    for i in range (0,5): # This will read the polynomial constants that define K as a function of temperature
        K_as_a_function_of_Temperature += globals()[name][i+index]*Temperature**i # Calculate the contribution to K as a function of temperature from the ith polynomial constant (note += just adds the calculation to the existing value in K_as_a_function_of_Temperature)
        
    # In this section we calculate the modification to the calculated K to account for pressure. So the previously calculated K is multiplied by this pressure correction term    
    Pressure_Modification = 0
    for i in range (0,5):
        Pressure_Modification += globals()[name][i+index+15]*Pressure**i
    
    return(K_as_a_function_of_Temperature*Pressure_Modification) # Return the value of K calculated as a function of temperature and modified by multiplying by a pressure correction

# rho and Cp are identical to the previously coded section for K. I wonder if there is a way to collapse this code into a single block?
# Probably, but I'm not that good yet    
def rho(name, Temperature, Pressure):
    index = 5
    rho_as_a_function_of_Temperature = 0
    for i in range (0,5):
        rho_as_a_function_of_Temperature += globals()[name][i+index]*Temperature**i
        
    Pressure_Modification = 0
    for i in range (0,5):
        Pressure_Modification += globals()[name][i+index+15]*Pressure**i
    
    return(rho_as_a_function_of_Temperature*Pressure_Modification)

def Cp(name, Temperature, Pressure):
    index = 10
    Cp_as_a_function_of_Temperature = 0
    for i in range (0,5):
        Cp_as_a_function_of_Temperature += globals()[name][i+index]*Temperature**i
        
    Pressure_Modification = 0
    for i in range (0,5):
        Pressure_Modification += globals()[name][i+index+15]*Pressure**i
    
    return(Cp_as_a_function_of_Temperature*Pressure_Modification)

# =============================================================================
# This function calculates a temperature gradient based on a passed Temperature array
# I actually tried using numpy.gradient before hand coding it. The gradient function returned a stencil that
# caused oscillation in an implicit solver I was trying to make so I just went back to hand coding the stencil
# As a benefit, you can basically mess with this all you want. Maybe try third and forth order solvers if you want
# =============================================================================

def Temperature_Gradient(Temperature_array): # define a function named "Temperature_Gradient" that gets passed a temperature array
    # Note that none of these gradients include the dx term. That can be applied later along with the diffusion and time step to calculate the actual heat flux from one cell to another. 
    # This function just needs to return the dT term so technically it's not the gradient :< and I lied to you. The math will work out though I promise
    Thermal_Gradient = np.zeros_like(Temperature_array) # First, create an empty array the size of the temperature array that was just passed to Temperature Gradient
    Thermal_Gradient[1:-1] = Temperature_array[:-2]+Temperature_array[2:]-2*Temperature_array[1:-1] # This is the primary stencil and its simply the forward Euler scheme for the diffusion equation
    # You can find more info on the stencil at the github comment I linked above the functions but wikipedia will do in a pinch
    Thermal_Gradient[0] = Temperature_array[1]-Temperature_array[0] # Left side dT. At the edges I just use a first order discritization because effectively it's just the equation for the central difference above
    # excep with one set of dT/dx removed. I'll have more information in a powerpoint presentation later
    Thermal_Gradient[-1] = Temperature_array[-2]-Temperature_array[-1] # Right side dT
    
    return(Thermal_Gradient) #Return the temperature differences for the array. 

# =============================================================================
# So the above function just calculated the temperature "difference" (I know I lied it's not a gradient).
# The following function actually calculates the heat fluxes for each cell broken out by each source
# So this section is basically shuttling around the energy into and out of the boundaries and through each cell
# This section is kind of like the energy balance though it's not technically conservative because it's just finite difference
# =============================================================================

def Boundary_Fluxes(Temperature_array, sim_time): # Define a function called "Boundary_Fluxes" that used a temperature array and a current time
# The current time gets passed off to functions that interpolate the transient boundary conditions
# Notice that all the terms begin with the Boolean statement that lets you turn them on/off, but there's not restriction to stop you from doing weird/dumb stuff so go wild
    
# =============================================================================
#     Static Boundary Fluxes from the python file
#     The conduction section closes out the non-transient part of the Forward Euler equation
#     Thermal conductivity is taken from the last cell so if you wanted to switch materials at the boundary (like metal plate sitting on concrete) you want to include one cell of the boundary material in the simulation
# =============================================================================
    
    Conduction_Heatflux_left = Conduction_left*(T_wall_left-Temperature_array[0])*(K_array[0]/x_step) # Calculate the heat flux from conduction at the left boundary node
    # Conduction_left is a boolean that turns this term on/off. T_wall_left is the boundary temperature. Temperature_array[0] is the left most cell in the simulation
    # And the K_array[0]/x_step is the scaling for thermal conduction. I think there's like, an old british term for K/x but I don't care at this point
    Conduction_Heatflux_right = Conduction_right*(T_wall_right-Temperature_array[-1])*(K_array[-1]/x_step) # Same as above but for the right boundary
    
    Convection_Heatflux_left = Convection_left*(h_left*(T_conv_left-Temperature_array[0])) # Convective boundary condition
    # Follows the same logic as above but for convection using h_left, T_conv_left and the left most temperature of the cell array
    Convection_Heatflux_right = Convection_right*(h_right*(T_conv_right-Temperature_array[-1])) # Same as above
    
    # Okay, so this one is a bit weird. I need a better way to store an array of radiative temperatures and view factors. Maybe I should include wavelength dependancy but that's long term. Right now the reason I don't use an array
    # and instead have two individually defined boundaries is because it reads the View Factors and Radiative Temperatures from columns in the excel file. I could probably add some string parsing but there's a lot of features I want to add
    
    # Raditive heat flux is broken out into incoming and outgoing radiation because that's way more realistic than this (T^4 - T^4) crap you find in undergrad text books. Emission and absorption are two totally different processes
    # Radiation_Heatflux_to_Boundary_XXX is energy going into the boundary cell. So this represents all the incoming energy into the simulation. The multiple view factors lets you simulate, for example, a hot source and a background
    # Otherwise it's just the generic radiative heating equation
    Radiation_Heatflux_to_Boundary_left = Radiation_to_Boundary_left*((absorptivity_left*View_Factor_1_left*steffan_boltzman_constant*T_rad_1_left**4)+(absorptivity_left*View_Factor_2_left*steffan_boltzman_constant*T_rad_2_left**4))
    Radiation_Heatflux_to_Boundary_right = Radiation_to_Boundary_right*((absorptivity_right*View_Factor_1_right*steffan_boltzman_constant*T_rad_1_right**4)+(absorptivity_right*View_Factor_2_right*steffan_boltzman_constant*T_rad_2_right**4))
    
    # Outgoing Radiative heatflux is a lot simpler since you don't care about view factor. Again this could be a function of wavelength or even directionality, but starting small.
    Radiation_Heatflux_from_Boundary_left = Radiation_from_Boundary_left*(-emissivity_left*steffan_boltzman_constant*Temperature_array[0]**4)
    Radiation_Heatflux_from_Boundary_right = Radiation_from_Boundary_right*(-emissivity_right*steffan_boltzman_constant*Temperature_array[-1]**4)
    
    # Something to know about Radiative heatflux. If you actually have ablation it's going to severely change the effective absorptivity and emissivity because typically the outgoing gas blocks the radiation. 
    # Just be aware that even if Radiation heatflux is present it can still be incredibly difficult to predict the effects in the presence of ablation and you'll be lucky to get +/- 20% of reality
    
    Constant_Heatflux_left = Fixed_Heatflux_left*q_in_left # Just a simple heatflux in/out of the boundary. This is where my calculator started a long time ago
    Constant_Heatflux_right = Fixed_Heatflux_right*q_in_right
    
    # Transient Boundary Conditions from the excel file
    # Down here we handle all the transient boundary conditions but only if a transient boundary condition is switched on. Otherwise this is ignored which speeds up the code
    if (Transient_Conduction_left+Transient_Convection_left+Transient_Radiation_to_Boundary_left+Transient_Radiation_from_Boundary_left+Transient_Fixed_Heatflux_left > 0): #Test if any transient boundary conditions are set to True (1)
        
        # If any of the transient boundary conditions are turned on, all of this will end up getting calculated. So make sure the transient boundary conditions is filled out. Otherwise I might have to break up the if statements and that's a lot of work
        # Everything below here is exactly the same as the steady state, but with the words "Transient_" in front of it. Amazing how that works
        Conduction_Heatflux_left += Transient_Conduction_left*(Transient_T_wall_left(sim_time)-Temperature_array[0])*(K_array[0]/x_step)
        Convection_Heatflux_left += Transient_Convection_left*(Transient_h_left(sim_time)*(Transient_T_conv_left(sim_time)-Temperature_array[0]))
        Radiation_Heatflux_to_Boundary_left += Transient_Radiation_to_Boundary_left*((Transient_absorptivity_left(sim_time)*Transient_View_Factor_1_left(sim_time)*steffan_boltzman_constant*Transient_T_rad_1_left(sim_time)**4)+(Transient_absorptivity_left(sim_time)*Transient_View_Factor_2_left(sim_time)*steffan_boltzman_constant*Transient_T_rad_2_left(sim_time)**4))
        Radiation_Heatflux_from_Boundary_left += Transient_Radiation_from_Boundary_left*(-Transient_emissivity_left(sim_time)*steffan_boltzman_constant*Temperature_array[0]**4)
        Constant_Heatflux_left += Transient_Fixed_Heatflux_left*Transient_q_in_left(sim_time)
   
    if (Transient_Conduction_right+Transient_Convection_right+Transient_Radiation_to_Boundary_right+Transient_Radiation_from_Boundary_right+Transient_Fixed_Heatflux_right > 0):
    
        Conduction_Heatflux_right += Transient_Conduction_right*(Transient_T_wall_right(sim_time)-Temperature_array[-1])*(K_array[-1]/x_step)
        Convection_Heatflux_right += Transient_Convection_right*(Transient_h_right(sim_time)*(Transient_T_conv_right(sim_time)-Temperature_array[-1]))
        Radiation_Heatflux_to_Boundary_right += Transient_Radiation_to_Boundary_right*((Transient_absorptivity_right(sim_time)*Transient_View_Factor_1_right(sim_time)*steffan_boltzman_constant*Transient_T_rad_1_right(sim_time)**4)+(Transient_absorptivity_right(sim_time)*Transient_View_Factor_2_right(sim_time)*steffan_boltzman_constant*Transient_T_rad_2_right(sim_time)**4))
        Radiation_Heatflux_from_Boundary_right += Transient_Radiation_from_Boundary_right*(-Transient_emissivity_right(sim_time)*steffan_boltzman_constant*Temperature_array[-1]**4)
        Constant_Heatflux_right += Transient_Fixed_Heatflux_right*Transient_q_in_right(sim_time)
    
    # Total all the heat fluxes together for the left and right boundary conditions
    # Notice that the time-step hasn't come into the equation yet. So this is just the per-second heatflux
    Left_Heatflux = Conduction_Heatflux_left+Convection_Heatflux_left+Radiation_Heatflux_to_Boundary_left+Radiation_Heatflux_from_Boundary_left+Constant_Heatflux_left
    Right_Heatflux = Conduction_Heatflux_right+Convection_Heatflux_right+Radiation_Heatflux_to_Boundary_right+Radiation_Heatflux_from_Boundary_right+Constant_Heatflux_right
    
    return(Left_Heatflux, Right_Heatflux) # Return the left and right heatflux

# A simple function to convert temperature from Kelvin, the best unit, to fahrenheit, objectively the worst temperature unit. It's 2021 and this is still necessary to properly communicate information. I cry    
def K_to_F(Temperature_array):
    return ((1.8*(Temperature_array-273.15))+32)

#%%

# =============================================================================
# Initialize all the necessary Arrays
# These are all just empty arrays/dataframes which will hold information during the calculation process
# =============================================================================

Temperature_array = np.full(len(cell_materials),T_initial) # This numpy array contains the actual temperature values for each cell during the simulation process
# Note that it's initialized at T_initial. In the future I should make this available to be overwritten by some initial temperature function
Heatflux_array = np.zeros(len(cell_materials)) # This stores all the heatfluxes during each calculation step (the q's if you will).
# If you want to know the heat in/out of one of the sides, just look at Heatflux_array[0] or Heatflux_array[-1] (Left and Right BC)

# These arrays store the material properties for each cell. They will get recalculated at some set time interval to keep from gumming up the solver
K_array = np.ones_like(Temperature_array)
rho_array = np.ones_like(Temperature_array)
Cp_array = np.ones_like(Temperature_array)
alpha_array = K_array/(rho_array*Cp_array) # thermal diffusivity for you cool kids

# percent_ablated stores how far along the right most cell is from maximum ablation. Ablation is tracked as a reduction in density based on values in the materials database
# Once percent_ablated goes to 100%, a copy of the Temperature_array is made without the last cell. This represents the new materials, sans the last cell which is dead
percent_ablated = 0.0

# This pandas dataframe took a while to figure out. Basically it has an index of the x location for each cell and column names for the time.
# The time is initially set to zero and a copy of the Temperature_array is stored
# Temperature_data is used for plotting but not part of the calculations so feel free to mess with it
Temperature_data = pd.DataFrame(index=np.linspace(0,sum(Thickness),int(sum(cell_count))), columns=[0], copy=True) 
# Create a dataframe named "Temperature_data" with an index (left column for you non-panda's folks)
# That index starts at 0 and goes to the total thickness of the original material sample. I guess in theory you could maybe add material? Another day
# Columns = 0 just means that the first column of stored data has a time-stamp of 0s
Temperature_data.iloc[:,0].iloc[0:Temperature_array.size]=Temperature_array
# Fill in that initial t = 0s with the data initialized into Temperature_array

#%%

# =============================================================================
# In this section we setup some basic parameters for the simulation. In the future I want to have a dedicated sheet for solve parameters that comes out with output.xlsx 
# That should make it easier to track what the simulation was for and understanding the results
# =============================================================================

sim_time = 0.0 # Resetting the simulation time back to 0. If you wanted to simulate a rocket there's technically no reason you couldn't have this negative
interval = plotting_interval + 1 # Interval is poorly named. It actually tracks the time since the last time data was output. So it's related to plotting interval, thus the name
material_interval = Material_Properties_Calculation_Interval+1.0 # How often material properties will be recalculated. By adding the +1 it will trigger the
# material property calculator at the start of the simulation
iteration = 0 # Tracks how many iterations the simulation takes. Typical simulations can be 100,000 ~ 1,000,000 iterations and runs about 200,000 iterations/minute
t = time.time() # Store the current actual physical time for estimating how much longer is left in the simulation

#%%

# =============================================================================
# This section is the actual calculation loop. Here we will take the current temperatures and heatfluxes and apply a small time_step to move energy around propotional to the heatflux * time_step
# There are a few logic loops that mostly track time for updating material properties and exporting data for plotting
# The ablation logic loop is the most complicated and I'll talk about it when I get to that section
# =============================================================================

while sim_time < time_final: # Continue iterations until the current simulation time is greater than the stated final time
    iteration += 1 # at the start of each iteration, increase the iteration count by 1
    
    # This logic loop checks if the time since the last time material properties were calculated is greater than the designated interval
    if (material_interval > Material_Properties_Calculation_Interval): # check if the material_interval is greater than the Material_Properties_Calculation_Interval
        material_interval = 0.0 # Once this loop has been triggered, reset the material_interval back to zero
        for i in range(0,len(Temperature_array)): # Go cell-by-cell recalculating the material properties. I really need to switch this from linear to doing the whole array at once
            K_array[i]     = K(cell_materials[i],Temperature_array[i],P_initial) # Pass the cell material, temperature, and P_initial (this should switch to pressure instead of P_initial. Might need a loop for Pressure(time))
            rho_array[i]   = rho(cell_materials[i],Temperature_array[i],P_initial) # Same as above
            Cp_array[i]    = Cp(cell_materials[i],Temperature_array[i],P_initial) # Same as above
        alpha_array = K_array/(rho_array*Cp_array) # Calculate the thermal diffusivity from the new material properties

    Conductive_Heatflux = Temperature_Gradient(Temperature_array)*K_array/x_step # Calculate the heatflux from conduction by calculating the temperature difference with "Temperature_Gradient" (not a gradient) and then scaling by K/dx to get heatflux
    Boundary_Heat_Flux = Boundary_Fluxes(Temperature_array, sim_time) # Calculate the boundary heat flux from the temperature array and the simulation time

    Heatflux_array = Conductive_Heatflux/(rho_array*Cp_array*x_step) # This is once again miss-labeled but it calculates the actual temperature rise per second for each cell from the heatflux we calculated earlier
    Heatflux_array[0] += (Boundary_Heat_Flux[0])/(rho_array[0]*Cp_array[0]*x_step) # This adds the effects of the boundary conditions on the temperature rise at the boundaries
    Heatflux_array[-1] += (Boundary_Heat_Flux[1])/(rho_array[-1]*Cp_array[-1]*x_step) # same as above

    time_step = (0.5*x_step**2)/np.max(alpha_array) # Calculate the appropriate timestep for a forward Euler solver that maintains stability. This is all based on the previously mentioned github document

    Temperature_array += time_step*Heatflux_array # Take the calculated dT/dt and multiply by the calculated time_step to get the change in temperature for each cell for this time step

# =============================================================================
#     The following section controls the logic for ablation. First ablation has to be turned on, thus the Ablation == True
#     Second, the material which makes up the outer right boundary has to actually be able to ablate. That's something set in the materials database
#     Finally, that cell also needs to be above it's ablation temperature for ablation to begin, If all three conditions are satisfied, then we calculate
#     ablation for the outermost cell layer
# =============================================================================
    
    if (Ablation == True) and (Materials[cell_materials[-1]]["Ablates"]==True) and (Temperature_array[-1]>Materials[cell_materials[-1]]["T_ab"]): # Check to see if ablation is on, if the last cell on the right can ablate, and if its temperature is high enough to ablate
        
# =============================================================================
#         This section actually tracks the ablation. What is happening here is, however much hotter than the ablation temperature the cell becomes, an equivalent mass of the material is removed
#         such that the temperature difference above the ablation temperature times the heat of ablation times the mass removal would leave the remaining material at the ablation temperature.
#         i.e. any heat above the ablation temperature just removes mass such that the heat of ablation would cancel out that extra heat
# =============================================================================
        
        # percent_ablated is cumulative starting from 0.00 so each iteration maybe 0.1% of the outer layer ablates and percent_ablated keeps track of that change
        percent_ablated += ((Temperature_array[-1]-Materials[cell_materials[-1]]["T_ab"])*Cp_array[-1]*rho_array[-1]/Materials[cell_materials[-1]]["h_ab"])/rho_array[-1]
        
        # Update the last cells density using the outer cells material, temperature, and pressure to calculate the baseline density and then multiply by the percentage ablated to get the new current density
        # I update the material properties of the final cell every iteration just because its the cell most likely to undergo rapid temperature changes that might not get caught by the default material properties update loop
        rho_array[-1] = rho(cell_materials[-1],Temperature_array[-1],P_initial)*(1.0-percent_ablated)
        K_array[-1] = K(cell_materials[-1],Temperature_array[-1],P_initial)
        Cp_array[-1] = Cp(cell_materials[-1],Temperature_array[-1],P_initial)
        
        # Update the outer layer Temperature to just be the ablation temperature
        Temperature_array[-1] = Materials[cell_materials[-1]]["T_ab"]
        
        # This section actually destroys the outermost cell if the density of the material gets low enough
        # It compared the current percent ablated against the maximum allowable percent ablated and, if the current percent_ablated is above the material limit for percent ablated, destroys that cell
        if (percent_ablated > Materials[cell_materials[-1]]["percent_available"]): # check if ablation is more than the limit
            percent_ablated = 0.0 # reset percent ablated 
            Temperature_array = Temperature_array[:-1] # Set temperature array to be a copy of itself but excluding the outer most cell. This effectively destroys that cell since everything is driven by temperature array
            K_array = K_array[:-1] # recalculate K array for the new temperature_array size
            rho_array = rho_array[:-1] # same as above
            Cp_array = Cp_array[:-1] # same as above
            alpha_array = alpha_array[:-1] # same as above
    
    interval += time_step # Add the current time_step to interval for tracking when data should be plotted and saved
    material_interval += time_step # Same as above but for recalculating material properties
    sim_time += time_step # Same as above but for total simulation time
    
    if interval > plotting_interval: # If the interval is larger than the plotting interval, save the data out
        Temperature_data[str(np.around(sim_time))] = np.nan # Create a pandas dataframe of nans (this triggers when the simulaiton states because I set interval = plotting_interval + 1 earlier)
        Temperature_data.iloc[:,-1].iloc[:Temperature_array.size]=Temperature_array.tolist() # Store the current Temperature_array into the Temperature_data dataframe
        interval = 0.0 # Reset the interval to zero
        Remaining_iterations = (time_final-sim_time)/time_step # Estimate the remaining iterations by calculating the remaining simulation time and dividing by the current time_steps (not super accurate)
        Remaining_time = Remaining_iterations*((time.time()-t)/iteration) # Estimate the remaining time from the remaining iteraitons multiplied by the average iteration time since starting (First time estimate is garbage but subsequent estimates are very accurate)
        print("sim_time = %3.4f (s), time_step = %3.6f (s), remaining_iterations = %6.0f, remaining_time = %4.1f minutes"%(sim_time,time_step, Remaining_iterations, Remaining_time/(60))) # print out some data to let the user know the calculator is running and
        # generally just to understand if this simulation will take forever


#%%
# =============================================================================
# At this point all the temperature data has been written to a pandas database named "Temperature_data" so now we plot the information
# =============================================================================

# Convert temperature to heretic units
Temperature_data_F = K_to_F(Temperature_data) # convert K to °F

fig1 = Temperature_data_F.plot(figsize = [11,8],legend=False); # Plot the Temperature_data on a figure sized 11" by 8" without a legend
plt.ylabel("Temperature (F)") # the y axis is temperature in °F
plt.xlabel("thickness (mm)") # the x axis is the location within the material stackup in mm
if (Temperature_data.shape[1]<20.0): # if there is more than 20 time series just don't plot the legend otherwise the figure gets messed up
    labels = ["t = " + s for s in fig1.get_legend_handles_labels()[1]] # get the current legend labels from fig1 ([0] is handles which is each line and [1] is the label of said line)
    # for whatever reason we have to use "string comprehension" because fig1.get_legend_handles_labels is a method and you cant just work with them directly so this code morphs them into a list of strings (I guess)
    fig1, plt.legend(labels = labels, bbox_to_anchor=(1.0, 1.01), loc='upper left') # Create an invisible box whos upper left corner starts at x = 1 plot width, y = 1.01 plot width
    # Put the lengend in this invisible box, anchored to the upper left corner

fig1, plt.tight_layout() # Add a tight layout to the plot because space is at a premium

for j in range(0,Thickness.size): # This code calculates how many materials you have specified by figuring out how many discreet entries there are in the "Thickness" array
    fig1, plt.axvline(np.sum(Thickness[0:j+1]), color = "black", lw=1) # Once per entry in the thickness array draw a black, vertical line at the x-position corresponding to the last cell of that particular material
    # This just makes it a lot easier to visualize where each material is in the stackup. Maybe I should fill in the background different ot something as well
    
# These are some handy plotting functions for maximuum temperatures and material thicknesses at various points. I turn them on/off with comments like a pleb
fig1, plt.axvline(Temperature_data_F.iloc[:,-1].count()*x_step*1000, color = "black", linestyle = '--', lw = 1) # Add a dashed vertical line (axvline) at the x position of the last cell in the simulation. Makes it easier to see the final point of ablation
#fig1, plt.axhline(Temperature_data_F.iloc[0,: ].max()              , color = "red"  , linestyle = '--', lw = 1) # Add a red, dashed horizontal line at the maximum temperature point of the left most cell. Good for checking on peak temperatures that made it through the materials
#fig1, plt.axhline(Temperature_data_F.max().max()                   , color = "red"  , linestyle = '--', lw = 1) # Add a red, dashed horizontal line at the location of maximum temperature in the entire simulation. Good for checking peak temperatures during transient runs
#fig1, plt.axhline(Temperature_data_F.iloc[-1,: ].max()              , color = "red"  , linestyle = '--', lw = 1) # Add a red, dashed horizontal line at the peak temperature of the right boundary condition. Doesn't work with ablation because most of those temps are nan
#fig1, plt.yticks(list(plt.yticks()[0])+[Temperature_data_F.iloc[-1,: ].max()]) # This takes the list of y-axis labels and adds a label for the peak temperature reached by the right boundary condition cell. Also doesn't work with ablation because most of the cells are nan
fig1, plt.yticks(list(plt.yticks()[0])+[Temperature_data_F.iloc[0,:].max()]+[Temperature_data_F.max().max()]) # Takes the y-axis labels and adds a label for the maximum temperature of the left boundary and the maximum overall simulation temperature
fig1, plt.xticks(list(plt.xticks()[0])+[Temperature_data.iloc[:,0].last_valid_index()]+[Temperature_data.iloc[:,-1].last_valid_index()]) # Adds a label on the x-axis for the last valid temperature data on the last recorded timestep. This should be equal to the ablation depth
fig1, plt.xlim([0,np.sum(Thickness)]) # Set the x-axis limits to be 0 and the maximum material thickness
#fig1, plt.ylim([np.min(np.min(np.floor(Temperature_data_F/100)))*100,np.max(np.max(np.ceil(Temperature_data_F/100)))*100]) # I tried to add some fancy code for setting the y-axis limit to nearest 100 for low and high temperature but just the default y-axis is better

#%%
Temperature_data.to_excel("output.xlsx", sheet_name = "output", na_rep = '') # Output the data as an xlsx file named "output", on a sheet named "output", throwing away any nan values (makes it easier to plot since the temp series will be different lengths if there is ablation)
Temperature_data.to_csv("output.csv",sep=',') # Output the data as a .csv with each value seperated by commas