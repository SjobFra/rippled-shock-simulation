#File containing the simulation parameters for the single rippled shock simulation

#Simulation file and group name inside file

filename = "errortest.hdf5"
groupname = "test3"

#Simulation parameters

unitless = True

if unitless:
    
    #Number of particles to be simulated
    particle_number = 1

    #Maximum simulation time
    t_max = 500.0  
    
    #Number of data writes to file
    savesteps = 1

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
    ripple_amplitude = 0.1
    
    #The wavelength of the ripple in the obliquity
    ripple_wavelength = 1.0
    
    #The ripple velocity in the x-plane
    ripplevx = 1.0
    
    #The ripple velocity in the y-plane
    ripplevy = 0.0

    #Particle position at injection
    x = 0.0
    
    #Particle speed at injection
    v = 6.0

else:
  
    #Number of particles to be simulated
    particle_number = 10000     

    #Maximum simulation time
    t_max = 300.0 #s

    #Number of data writes to file
    savesteps = 1

    #Gas compression ratio
    r = 8.0

    #Magnetic compression ratio
    r_mag = 4.0 
    
    #mean free path of the particles in kilometers
    meanfreepath = 1.0 #km

    #Upstream plasma velocity
    v_up = 2.0 #km/s
    
    #The mean obliquity of the shock
    obliquity_mean = 0.0
    
    #The chance of an ripple 
    ripple_chance = 0.0
    
    #The amplitude of the ripple in the obliquity
    ripple_amplitude = 0.0
    
    #The wavelength of the ripple in the obliquity
    ripple_wavelength = 1.0
    
    #Particle position at injection
    x = 0.0

    #Particle speed at injection
    v = 8.0 #km/s


    #Changing to an unitless system for the simulation

    t_max = t_max*(v_up/2)/meanfreepath
    