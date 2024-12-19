# FinalProject_LionfishConcentrations
This model simulates the spread and population dynamics of lionfish in the Mediterranean Sea using Sea Surface Temperature (SST) data and a two-dimensional diffusion equation. 

OVERVIEW: 
This model simulates the spread and population dynamics of lionfish in the Mediterranean Sea using Sea Surface Temperature (SST) data and a two-dimensional diffusion equation. The simulation takes into account environmental conditions (temperature) and some biological characteristics of the species such as their high reproductive rates and high tolerance for a wide range of environmental conditions. 
So far lionfish as an invasive species have caused a significant amount of ecological damage to the Gulf of Mexico and the Caribbean. This model aims to explore their potential future distribution in the Mediterranean Region particularly as sea surface temperatures increase due to climate change. 
The code models lionfish diffusion and population growth over time, using SST data from the National Oceanic and Atmospheric Administration’s (NOAA) Extended Reconstructed SST V5 dataset and applies a logistic growth model to account for the species’ carrying capacity. 

BACKGROUND: 
Lionfish originate in the Indo-Pacific however aquarium releases have allowed them to become the ideal invasive species because of their natural lack of predators outside their home range, generalist diets, tolerance to a wide range of environmental conditions (depth, temperature, salinity), and high reproductive rates. Since the mid 1980s lionfish have sprawled across the Gulf of Mexico and the Caribbean. More recent sightings in the Mediterranean and the continuous rise in sea surface temperatures suggest there may be a rapid increase in this region. 

The purpose of this code is to model the spatial distribution and growth rates of lionfish over time as affected by Sea Surface Temperatures. In this model we are assuming that lionfish are introduced initially at a fixed value at one central location in the Mediterranean, the population is spreading only according to this diffusion pattern, and the population growth rates are affected  mostly by Sea Surface Temperatures. 

The structure of the code is split up into 7 sections following: section 1: loading the sst data, section 2: setting the initial conditions, section 3: checking diffusion problem stability, section 4: computing the diffusion equation in the A matrix, section 5: plotting the initial concentrations of lionfish, section 6: creating a time loop for the simulation, and section 7: plotting post-simulation concentrations of lionfish. 

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

Language: Python 
Libraries: numpy, matplotlib, xarray, and  mpl_toolkits.mplot3d
