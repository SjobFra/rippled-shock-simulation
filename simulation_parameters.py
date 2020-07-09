"""
This file contains the parameters used in the rippled shock simulation.
"""


#Name of the file where the results are saved. 
#Must be .hdf5 format. (Should rework so user cant edit format) 
filename = "test.hdf5"

#Inside the -hdf5 file data can be saved in groups. 
#Changing just the groupname when running multiple simulations 
#is preferred so user does not have to contend with multiple files.
groupname = "test1"

#Number of data writes to file
#if 1 then only the starting and ending steps are saved
#if more than 1 then the simulation time is divided to that many parts and
#results are written at these times.
savesteps = 1

#Choose how many cores to use in the simulation. 
#If incompatible number is given, then uses all but one core
sim_cores = 1

#Number of particles to be simulated
particle_number = 10

#Maximum simulation time
t_max = 50.0  

#Gas compression ratio
r = 8.0

#mean free path
meanfreepath = 1.0

#Upstream plasma velocity
v_up = 2.0   

#The mean obliquity of the shock
obliquity_mean = 0.0

#The chance of an ripple (from 0.0 with no chance of ripple to 1.0 where the ripple is always present)
ripple_chance = 1.0

#The amplitude of the ripple in the obliquity
#At the beginning of the simulation a test is made to run only acceptable ripples
#(So that the shock obliquity does not ever exceed 90 degrees)
ripple_amplitude = 0.1

#The wavelength of the ripple in the obliquity
ripple_wavelength = 1.0

#The ripple velocity in the x-plane
ripplevx = 1.0

#The ripple velocity in the y-plane
ripplevy = 0.0

#Particle position at injection
#The shock is located at infitesimally close to zero on the positive side
x = 0.0

#Particle speed at injection
v = 6.0