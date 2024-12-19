#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 10:48:13 2024

@author: hopeangotti
"""
'''
Projected Lionfish Population Dynamics in the Mediterranean as it Relates 
to Sea Surface Temperature Data 

BACKGROUND: 
    
    Lionfish are known as the ideal invaders due to their lack of natural 
    predators outside of their home range, their generalist diets, their 
    tolerance to a wide range of depths, temperatures, and salinity ranges, and 
    high reproductive rates. This species has already caused ecological 
    decimation to the Gulf of Mexico and the Caribbean since the 1980s but now 
    in more recent years, they have been sighted in the Mediterranean and with 
    increasing sea surface temperatures, better suited to their preference, 
    are projected to grow more rapidly. 


PURPOSE: 
    
    This model is simulating the spread and population dynamics of lionfish 
    species invasive to a specific geographic area (the Mediterranean) using a 
    a two-dimensional diffusion equation to model their spatial distribution 
    and population dynamics over time. It uses Sea Surface Temperature (SST)
    data from the National Oceanic and Atmospheric Administration (NOAA) to 
    adjust the population growth rate of the lionfish based on the environmental
    conditions. The model also incorporates a logistic growth model to account 
    for the carrying capacity of the species.


MASTERKEY:
    
    Section 1: 
        file_path -- path to the NetCDF file containing Sea Surface Temperature 
        (SST) data. 
        sst_dataset -- this is the dataset that has been loaded from NOAA’s 
        NetCDF file on Reconstructed SST records
        temps -- the extracted SST variable from the dataset which represents 
        sea surface temperatures for each spatial grid 
        
    Section 2: 
        nx -- represents the number of grid points in the x-direction (longitude)
        ny -- represents the number of grid points in the y-direction (latitude)
        dx -- represents the grid spacing in the x-direction (longitude) in 
        kilometers 
        dy -- represents the grid spacing in the y-direction (latitude) in 
        kilometers 
        x -- represents the 1D array of x-coordinates from 0 to nx*dx with 
        spacing dx 
        y -- represents the 1D array of y-coordinates from 0 to ny*dy with 
        spacing dy
        X -- represents the 2D meshgrid array of x-coordinates for the grid 
        Y -- represents the 2D meshgrid array of y-coordinates for the grid
        D_rate -- the constant diffusion rate of the lionfish population in kilometers 
        per year
        dt -- the time step for the simulation in years (how much time between 
        each update of populations)
        carrying_cap -- the carrying capacity of the lionfish population 
        (maximum population density of lionfish per one square kilometer)
        sst_specific -- a 1D array which stores the subset of SST data for the 
        specific geographic location in the Mediterranean (adjusted latitude 
        point = 28, adjusted longitude point = 110)
        C_lionfish -- the 2D array that represents the lionfish population 
        density at each grid point
        c_lionfish -- the flattened 1D version of C_lionfish which can be used for 
        necessary computations
        
    Section 3: 
        sx -- the stability constant in the x-direction for the von Neumann 
        stability condition
        sy -- the stability constant in the y-direction for the von Neumann 
        stability condition 
        
    Section 4: 
        A -- the square matrix used for the diffusion equation, representing the 
        diffusion coefficients for the model at each individual grid point 
        ik -- the index for the 2D array 
        i -- common loop variable representing x-direction 
        k -- common loop variable representing y-direction
        
    Section 6: 
        totaltime -- the total time for the simulation in years 
        time -- the current time for the simulation 
        newc_lionfish -- the new population values after multiplying the diffusion A matrix 
        growthrate -- the growth rate of the lionfish population depending on SST
        c_lionfish -- the flattened array of lionfish population 


INPUTS:
    
    1. Sea Surface Temperature Data (SST) -- from NOAA's NetCDF file ('NOAA 
    Extended Reconstructed SST V5') it is a global monthly SST analysis from 
    1854 to the present 
        - sst_dataset, sst_specific
    
    2. Grid parameters 
        - nx, ny, dx, dy, X, Y
    
    3. Diffusion and growth parameters 
        - D_rate, dt, carrying_cap, growthrate
    
    4. Initial population distribution of lionfish 
        - C_lionfish, c_lionfish
        
    5. Simulation parameters 
        - totaltime, time

OUTPUTS:
    
    1. Lionfish population distribution over time (initially and after simulated
       time)
    
    2. 2D colormap plots showing spatial distribution of lionfish population 
    
    3. 3D surface plot showing spatial distribution of lionfish population 
    
    4. Stability of model 

'''

#   SECTION 1
#---------------
# Import all packages needed
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


# Import xarray to handle NetCDF data file 
import xarray as xr


# Load the NetCDF data file that contains the Sea Surface Temperature data (SST)
file_path = 'sst.mnmean.nc copy' # File path to the SST data
sst_dataset = xr.open_dataset(file_path) # Open dataset using xarray


temps = sst_dataset.sst # Extract the 'sst' variable from the dataset and call it 'temps'
print(temps) # Display the sst data 


#   SECTION 2
#---------------
###-- SET INITIAL CONDITIONS --###
nx = 100 # Define number of gridpoints in the x-direction (longitude)
ny = 100 # Define number of gridpoints in the y-direction (latitude)


# Define the distance between the gridpoints in kilometers (spatial resolution of grid)
dx = 0.3 # Grid spacing in x-direction in kilometers 
dy = 0.3 # Grid spacing in y-direction in kilometers 


# Create 1D arrays of x-coordinates and y-coordinates based on grid-size and grid-spacing
x = np.arange(0, nx*dx, dx) # Array of x-coordinates from 0 to nx*dx in kilometers (longitude)
y = np.arange(0, ny*dy, dy) # Array of x-coordinates from 0 to nx*dx in kilometers (longitude)


# Create a meshgrid of x and y to represent a 2D spatial grid 
X, Y = np.meshgrid(x, y, indexing = 'ij') # create a 2-D coordinate system for plotting 


# Define the constant diffusion rate at which lionfish populations spread in kilometers per year 
D_rate = 0.5 # km/year 


# Define timestep for simulation in years 
dt = 0.009 # years 


# Define constant carrying capacity of lionfish (the maximum popultion density of lionfish per squre kilometer)
carrying_cap = 50000 # the carrying capacity of lionfish is approximately  50,000 lionfish per square km


# Extract a specific portion of the SST data for a region of interest (latitude by longitude)
# the format for calling the data at these specific points is temps[time, longitude, latitude]
# longitude boundaries of Mediterranean converted to matrix points = 90:124
# latitude boundaries of Mediterranean converted to matrix points = 23:28
sst_specific = temps[:,28,110]  # using points 28 latitude by 110 longitude 


# Define initial concentrations of lionfish in the 2D array of zeros 
C_lionfish = np.zeros((nx, ny)) # Starting with an empty grid of no lionfish 
C_lionfish[int(nx/2), int(ny/2)] = 20 # Introducing an intial population of 20 lionfish at the center 
c_lionfish = C_lionfish.flatten() # Flatten the 2D array into a 1D array by rows for further computation 


#   SECTION 3
#---------------
###-- STABILITY CHECK --###
# Check the stability of the problem using the von Neumann stability formula for the diffusion equation 
sx = dt * D_rate / dx**2 # von Neumann stability constant in x-direction 
sy = dt * D_rate / dy**2 # von Neumann stability constant in y-direction


# Test the stability for the computations 
import sys # Use to exit the program if the conditions are unstable 
if sx>0.5:
    print('x is unstable') # Indicate instability in x-direction 
    sys.exit() # Exit simulation if x is unstable 
elif sy>0.5:
    print('y is unstable') # Indicate instability in y-direction 
    sys.exit() # Exit simulation if y is unstable 


#   SECTION 4
#---------------
###-- CREATE THE A MATRIX --###
A = np.zeros((nx*ny, nx*ny)) # Create square A matrix initialized to zero 


# Populate the A matrix based on the diffusion equation 
for i in range(nx):
    for k in range (ny): 
        ik = i*ny + k # Calculates index for 2D array 
        
        ## ----- BOUNDARY CONDITIONS ----- ##
        # Set boundary points where ther is no diffusion (edges of grid)
        if i == 0:
            A[ik, ik] = 1 # no change in scenario/ diffusion at boundary points
        elif i == (nx-1): 
            A[ik, ik] = 1 
        elif k == 0: 
            A[ik, ik] = 1
        elif k == (ny-1):
            A[ik, ik] = 1
        else: 
            ##-----MATRIX COEFFICIENT-----## 
            A[ik,ik] = 1-2*sx - 2*sy # Central coefficient for the current point
            A[ik, (i+1)*ny + k] = sx # Diffusion in x-direction (right)
            A[ik, (i-1)*ny + k] = sx # Diffusion in x-direction (left)
            A[ik, i*ny + k + 1] = sy # Diffusion in y-direction (top)
            A[ik, i*ny + k - 1] = sy # Diffusion in y-direction (bottom)
             

#   SECTION 5
#---------------        
###-- PLOTTING INITIAL CONDITIONS --###
## Plot the intial lionfish population distribution 
# Create a 3D Surface Plot 
fig = plt.figure() # Create new figure
ax = fig.add_subplot(111, projection = '3d') # Add subplot 3D
ax.plot_surface(X, Y, C_lionfish) # Plot intial distribution 


# Create a 2D Colormesh Grid 
fig1, ax1 = plt.subplots(1,1) # Create a subplot 2D
im = ax1.pcolormesh(X, Y, C_lionfish) # Plot the initial population on the grid 
ax1.set_title('Lionfish Presence after Initial Introduction') # Set the title 
fig1.colorbar(im, ax=ax1, label = 'Concentration of Lionfish (fish/km)') # Add a colorbar for scale 


#   SECTION 6
#---------------
###-- TIME LOOP --###
# Create a time loop to plot the simulation modeling lionfish population dynamics over time 
totaltime = 10 # Set the total time for the simulation in years 
time = 0 # Initialize the current time 

while time <= totaltime:
    newc_lionfish = np.dot(A, c_lionfish) # Multiply dot product of the A matrix and the current population 
    c_lionfish[:] = newc_lionfish # Update the population with new values 
    
    # Determine the growth rate of the lionfish population based on the Sea Surface Temperature of the area
    # 23.33 degrees Celsius is the minimum prime water temperature for lionfish to thrive in
    if sst_specific[i] >= 23.33: # If SST is greater than or equal to 23.33 degrees celsius, the growth rate will increase 
        growthrate = 2.0 # Population will experience a higher growth rate 
    else: # If SST does not exceed 23.33 degrees celsius, lionfish population will experience lower growth rate 
        growthrate = 1.1 # Population will experience lower growth rate 
    
    # Apply the logistic growth model to the population of lionfish: 
    ## growth * population * (1-pop/carrying capacity)
    c_lionfish += (growthrate * c_lionfish * (1-c_lionfish/carrying_cap))*dt
    
    # Use timestep to increment time
    time +=dt


#   SECTION 7
#---------------
###-- PLOT THE LIONFISH POPULATION DYNAMICS AFTER TIME --###
# Reshape the population vector back to 2D grid for plotting purposes 
C_lionfish = c_lionfish.reshape(X.shape) # Reshape 1D vector into a 2D array 

# Plot the distribution of lionfish after the simulation 
# Create a 3D Surface Plot 
fig = plt.figure() # Create new figure
ax = fig.add_subplot(111, projection = '3d')  # Add subplot 3D
ax.plot_surface(X, Y, C_lionfish) # Plot final distribution 

# Create a 2D Colormesh Grid 
fig1, ax1 = plt.subplots(1,1) # Create a subplot 2D
im = ax1.pcolormesh(X, Y, C_lionfish) # Plot the initial population on the grid 
ax1.set_title('Lionfish Presence 10 years after Initial Introduction') # Set the title 
fig1.colorbar(im, ax=ax1, label = 'Concentration of Lionfish (fish/km)') # Add a colorbar for scale 



