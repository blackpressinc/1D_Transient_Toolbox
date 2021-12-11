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
time_final = 181.0 #[s]
plotting_interval = 20.0 #[s] Data will ***saved*** be plotted at this interval. Note that this controls when data is saved
T_initial = 274.0 #[K] All cells in the simulation will be initialized with this value. Once the Temperature_array has been created, you can add a custom temperature function or something
# but I have not yet coded that
P_initial = 1.0 #[atm] #pressure needs to be re-worked so that it can be made a function of time. This just acts as a place holder
steffan_boltzman_constant = 5.67*10**-8 #used for radiation calculations
# =============================================================================
# This solver supports transient material properties where available but calculating them each iteration makes the solver incredibly slow
# So you can set a time to recalculate the material properties every X seconds
# =============================================================================
Material_Properties_Calculation_Interval = 1.0 #[s]

# =============================================================================
# This is the actual material stackup you are simulating. Materials are assigned left to right based on the thicknesses. 
# Material list is just strings that have to match with the some name in the first row of the Material_Properties.xlsx
# Thickness is the total thickness of each material in millimeters (I will never use english to define a property god help me)
# =============================================================================
Material_list = ["T792","A12","P50"]
Thickness = np.array([3.15,0.25,10.0]) #[mm] Array defining the thickness of the individual layers

# =============================================================================
# These parameters actually setup the simulation mesh. Whichever is the tightest requirement will be used
# You can make a layer as thin as you want but timestep gets smaller as a function of the smallest cell thickness squared dt ≈ 1/(dx^2)
# =============================================================================
minimum_cell_count_per_layer = 4 #Thinnest Layer will be used to calculate dx for the mesh
maximum_cell_size = 0.05 #[mm] Cells will not be made larger than this value

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
T_wall_left = 270 #[K]
h_left = 20 #[W/m²·K]
T_conv_left = 266.483 #[K]
emissivity_left = 0.03
absorptivity_left = 0.03
View_Factor_1_left = 0.5
View_Factor_2_left = 1 - View_Factor_1_left
T_rad_1_left = 300 #[K]
T_rad_2_left = 300 #[K]
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
Transient_Convection_left = False
Transient_Radiation_to_Boundary_left = False
Transient_Radiation_from_Boundary_left = False
Transient_Fixed_Heatflux_left = False

# =============================================================================
# If you are loading in transient boundary conditions, you'll need to point them to an excel file somewhere in the directory you are running this script from
# It shouldn't run if all of the Transient_ booleans are set to false
# =============================================================================
#Left Transient Boundary Conditions
if (Transient_Conduction_left+Transient_Convection_left+Transient_Radiation_to_Boundary_left+Transient_Radiation_from_Boundary_left+Transient_Fixed_Heatflux_left > 0): #Test if any transient boundary conditions are set to True (1)
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
T_wall_right = 310 #[K]
h_right = 25 #[W/m²·K]
T_conv_right = 310 #[K]
emissivity_right = 0.9
absorptivity_right = 0.9
View_Factor_1_right = 0.44
View_Factor_2_right = 1 - View_Factor_1_right
T_rad_1_right = 2200 #[K]
T_rad_2_right = 300 #[K]
q_in_right = 600000 #[W/m²]

Conduction_right = False
Convection_right = False
Radiation_to_Boundary_right = False
Radiation_from_Boundary_right = False
Fixed_Heatflux_right = False

Transient_Conduction_right = False
Transient_Convection_right = True
Transient_Radiation_to_Boundary_right = True
Transient_Radiation_from_Boundary_right = True
Transient_Fixed_Heatflux_right = False

if (Transient_Conduction_right+Transient_Convection_right+Transient_Radiation_to_Boundary_right+Transient_Radiation_from_Boundary_right+Transient_Fixed_Heatflux_right > 0):
#Right Transient Boundary Conditions
    Right_Transient_BC = pd.read_excel('Example_Transient_Boundary_Conditions.xlsx',index_col=0,skiprows=1,header=1)

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

Materials = pd.read_excel('Material_Properties.xlsx',index_col=(0),usecols=("C:W"),skiprows=(0)) #Note that you can use "usecols" to just select a subset of 
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
"""
This section creates functions to quickly calculate the K, rho, and Cp using the materials name
Temperature, and pressure. 

I need to add in reference pressure and temperature for materials in the future. Otherwise T and P are a bit abstract
"""

def K(name, Temperature, Pressure):
    index = 0
    K_as_a_function_of_Temperature = 0
    for i in range (0,5):
        K_as_a_function_of_Temperature += globals()[name][i+index]*Temperature**i
        
    Pressure_Modification = 0
    for i in range (0,5):
        Pressure_Modification += globals()[name][i+index+15]*Pressure**i
    
    return(K_as_a_function_of_Temperature*Pressure_Modification)
    
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

def Temperature_Gradient(Temperature_array):
    Thermal_Gradient = np.zeros_like(Temperature_array)
    Thermal_Gradient[1:-1] = Temperature_array[:-2]+Temperature_array[2:]-2*Temperature_array[1:-1]
    Thermal_Gradient[0] = Temperature_array[1]-Temperature_array[0]
    Thermal_Gradient[-1] = Temperature_array[-2]-Temperature_array[-1]
    
    return(Thermal_Gradient)
    
def Boundary_Fluxes(Temperature_array, sim_time):
    
    # Static Boundary Fluxes from the python file
    Conduction_Heatflux_left = Conduction_left*(T_wall_left-Temperature_array[0])*(K_array[0]/x_step)
    Conduction_Heatflux_right = Conduction_right*(T_wall_right-Temperature_array[-1])*(K_array[-1]/x_step)
    
    Convection_Heatflux_left = Convection_left*(h_left*(T_conv_left-Temperature_array[0]))
    Convection_Heatflux_right = Convection_right*(h_right*(T_conv_right-Temperature_array[-1]))
    
    Radiation_Heatflux_to_Boundary_left = Radiation_to_Boundary_left*((absorptivity_left*View_Factor_1_left*steffan_boltzman_constant*T_rad_1_left**4)+(absorptivity_left*View_Factor_2_left*steffan_boltzman_constant*T_rad_2_left**4))
    Radiation_Heatflux_to_Boundary_right = Radiation_to_Boundary_right*((absorptivity_right*View_Factor_1_right*steffan_boltzman_constant*T_rad_1_right**4)+(absorptivity_right*View_Factor_2_right*steffan_boltzman_constant*T_rad_2_right**4))
    
    Radiation_Heatflux_from_Boundary_left = Radiation_from_Boundary_left*(-emissivity_left*steffan_boltzman_constant*Temperature_array[0]**4)
    Radiation_Heatflux_from_Boundary_right = Radiation_from_Boundary_right*(-emissivity_right*steffan_boltzman_constant*Temperature_array[-1]**4)
    
    Constant_Heatflux_left = Fixed_Heatflux_left*q_in_left
    Constant_Heatflux_right = Fixed_Heatflux_right*q_in_right
    
    #Transient Boundary Conditions from the excel file
    if (Transient_Conduction_left+Transient_Convection_left+Transient_Radiation_to_Boundary_left+Transient_Radiation_from_Boundary_left+Transient_Fixed_Heatflux_left > 0): #Test if any transient boundary conditions are set to True (1)

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
    
    Left_Heatflux = Conduction_Heatflux_left+Convection_Heatflux_left+Radiation_Heatflux_to_Boundary_left+Radiation_Heatflux_from_Boundary_left+Constant_Heatflux_left
    Right_Heatflux = Conduction_Heatflux_right+Convection_Heatflux_right+Radiation_Heatflux_to_Boundary_right+Radiation_Heatflux_from_Boundary_right+Constant_Heatflux_right
    
    return(Left_Heatflux, Right_Heatflux)
    
def K_to_F(Temperature_array):
    return ((1.8*(Temperature_array-273.15))+32)

#%%
"""
Initialize all the necessary Arrays
"""
Temperature_array = np.full(len(cell_materials),T_initial)
Temperature_history = np.copy(Temperature_array)
Heatflux_array = np.zeros(len(cell_materials))
K_array = np.ones_like(Temperature_array)
rho_array = np.ones_like(Temperature_array)
Cp_array = np.ones_like(Temperature_array)
alpha_array = K_array/(rho_array*Cp_array)
percent_ablated = 0.0

Temperature_data = pd.DataFrame(index=np.linspace(0,sum(Thickness),sum(cell_count)), columns=[0])
Temperature_data.iloc[:,0].iloc[0:Temperature_array.size]=Temperature_array

#%%
"""
Initializing of the solution
"""
sim_time = 0.0
interval = plotting_interval+1.0
material_interval = Material_Properties_Calculation_Interval+1.0
iteration = 0
data_storage_index = 0
phi = 0.5
iteration_time = 0.0
data = []

for i in range(0,len(cell_materials)):
    K_array[i]     = K(cell_materials[i],Temperature_array[i],P_initial)
    rho_array[i]   = rho(cell_materials[i],Temperature_array[i],P_initial)
    Cp_array[i]    = Cp(cell_materials[i],Temperature_array[i],P_initial)
    
alpha_array = K_array/(rho_array*Cp_array)

t = time.time()


#%%
"""
Looping solution begins here
"""

while sim_time < time_final:
    if (material_interval > Material_Properties_Calculation_Interval):
        material_interval = 0.0
        for i in range(0,len(Temperature_array)):
            K_array[i]     = K(cell_materials[i],Temperature_array[i],P_initial)
            rho_array[i]   = rho(cell_materials[i],Temperature_array[i],P_initial)
            Cp_array[i]    = Cp(cell_materials[i],Temperature_array[i],P_initial)
        alpha_array = K_array/(rho_array*Cp_array)
        data.append(Temperature_array[0])
        
    
    iteration += 1
    """
    Initialize the Material properties arrays
    """
    
    """
    Calculate the conductive heat fluxes
    """
    
    Conductive_Heatflux = Temperature_Gradient(Temperature_array)*K_array/x_step
    Boundary_Heat_Flux = Boundary_Fluxes(Temperature_array, sim_time)
        
    Heatflux_array = Conductive_Heatflux/(rho_array*Cp_array*x_step)
    Heatflux_array[0] += (Boundary_Heat_Flux[0])/(rho_array[0]*Cp_array[0]*x_step)
    Heatflux_array[-1] += (Boundary_Heat_Flux[1])/(rho_array[-1]*Cp_array[-1]*x_step)
    
    time_step = (0.5*x_step**2)/np.max(alpha_array)
    
    Temperature_array += time_step*Heatflux_array  
    
    if (Ablation == True) and (Materials[cell_materials[-1]]["Ablates"]==True) and (Temperature_array[-1]>Materials[cell_materials[-1]]["T_ab"]):
        percent_ablated += ((Temperature_array[-1]-Materials[cell_materials[-1]]["T_ab"])*Cp_array[-1]*rho_array[-1]/Materials[cell_materials[-1]]["h_ab"])/rho_array[-1]
        Temperature_array[-1] = Materials[cell_materials[-1]]["T_ab"]
        
        if (percent_ablated > Materials[cell_materials[-1]]["percent_available"]):
            percent_ablated = 0.0
            Temperature_array = Temperature_array[:-1]
            K_array = K_array[:-1]
            rho_array = rho_array[:-1]
            Cp_array = Cp_array[:-1]
            alpha_array = alpha_array[:-1]
    
    interval += time_step
    material_interval += time_step
    sim_time = sim_time + time_step
    
    if interval > plotting_interval:
        Temperature_data[str(np.around(sim_time))] = np.nan
        Temperature_data.iloc[:,-1].iloc[:Temperature_array.size]=Temperature_array.tolist()
        data_storage_index += 1
        interval = 0.0
        Remaining_iterations = (time_final-sim_time)/time_step
        Remaining_time = Remaining_iterations*((time.time()-t)/iteration)
        print("sim_time = %3.4f (s), time_step = %3.6f (s), remaining_iterations = %6.0f, remaining_time = %4.1f minutes"%(sim_time,time_step, Remaining_iterations, Remaining_time/(60)))


#%%
Temperature_data_F = K_to_F(Temperature_data)

fig1 = Temperature_data_F.plot(figsize = [11,8],legend=False);
plt.ylabel("Temperature (F)")
plt.xlabel("thickness (mm)")
if (Temperature_data.shape[1]<20.0):
    fig1, plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
fig1, plt.tight_layout()
for j in range(0,Thickness.size):
    fig1, plt.axvline(np.sum(Thickness[0:j+1]), color = "black", lw=1)
fig1, plt.axvline(Temperature_data_F.iloc[:,-1].count()*x_step*1000, color = "black", linestyle = '--', lw = 1)
fig1, plt.axhline(Temperature_data_F.iloc[0,: ].max()              , color = "red"  , linestyle = '--', lw = 1)
fig1, plt.yticks(list(plt.yticks()[0])+[Temperature_data_F.iloc[0,:].max()])
fig1, plt.xticks(list(plt.xticks()[0])+[Temperature_data.iloc[:,0].last_valid_index()]+[Temperature_data.iloc[:,-1].last_valid_index()])
fig1, plt.xlim([0,np.sum(Thickness)])
fig1, plt.ylim([0,1200])

#%%
Temperature_data.to_excel("output.xlsx", sheet_name = "output", na_rep = '')