#
#Single shock simulation main loop
#
import random as rnd #random() generates a random float uniformly in the semi-open range [0.0, 1.0)
import numpy as np
import h5py
import time
import warnings
from multiprocessing import Pool
from simulation_functions import changeFrame
from simulation_functions import scatterParticle
from simulation_functions import moveParticle
from simulation_functions import calculateCompressionRatios
import simulation_parameters


def simulateParticle(particle, x, v, t, mu, t_max, count, totalparticles, v_up, obliquity_mean, ripple_amplitude, ripple_wavelength, ripple_chance, ripplevx, ripplevy, r, mfp):
    
    #Seed random with each particle since multiprocessing breaks higher seedings
    rnd.seed((particle+(totalparticles*(count-1))))
    #Testing with taking compression ratio calculation random numbers from numpy random
    np.random.seed(int((particle+(totalparticles*(count-1)))))
    
    #Calculate random variables needed in the compression ratio calculations
    randx = np.random.uniform()
    randy = np.random.uniform()
    
    r, r_mag, deltau = calculateCompressionRatios(randx, randy, obliquity_mean, ripple_amplitude, ripple_wavelength, ripplevx, ripplevy, r)
    
    #Follow particle until maximum simulation time is reached
    while t < t_max:
        
        #Check whether the particle is in the upstream or downstream region
        if x <= 0:
            u = v_up
        else:
            u = ((v_up*r_mag) / r)
        
        #Shock frame to plasma frame
        v, mu = changeFrame(v, mu, u)
        
        #Calculate new timestep according to particle speed
        dt = 0.1 * mfp / v
        
        #Calculate random variables needed in scattering so numba can be used
        thetarn = rnd.random()
        phirn = rnd.random()
        
        #Scatter in plasma frame
        mu = scatterParticle(dt, v,  mu, mfp, thetarn, phirn)
        
        #Plasma frame to shock frame
        v, mu = changeFrame(v, mu, (-u))
        
        #Calculate the random number needed in checking if the ripple is there
        ripple_rand = np.random.uniform()
        
        #Calculate threshold value for mu
        mu_0 = (-x)/(v*dt)
        
        #Check if the particle is in the upstream
        if x < 0:
            
            #Also if the pitch angle is above threshold then shock crossing happens
            if mu > mu_0:
                
                #Check if there is a ripple
                if ripple_rand <= ripple_chance:
                    
                    #Start loop to recalculate deltau if the ripple is escaping the particle
                    while True:
                    
                        #Calculate random variables needed in the compression ratio calculations
                        randx = np.random.uniform()
                        randy = np.random.uniform()
                        
                        #Calculate new magnetic compression ratio for downstream
                        r, r_mag, deltau = calculateCompressionRatios(randx, randy, obliquity_mean, ripple_amplitude, ripple_wavelength, ripplevx, ripplevy, r)
                        
                        if v*mu + deltau > 0:
                            break
                    
                    ###THIS VERSION JUST REPLACES INCORRECT RIPPLES WITH 0.0
                    ##Calculate random variables needed in the compression ratio calculations
                    #randx = np.random.uniform()
                    #randy = np.random.uniform()
                    #
                    ##Calculate new magnetic compression ratio for downstream
                    #r, r_mag, deltau = calculateCompressionRatios(randx, randy, obliquity_mean, ripple_amplitude, ripple_wavelength, ripplevx, ripplevy, r)
                    #
                    #if v*mu + deltau < 0:
                    #    deltau = 0.0
                
                            
                else:
                    
                    #Calculate random variables needed in the compression ratio calculations
                    randx = np.random.uniform()
                    randy = np.random.uniform()
                    
                    #If no ripple is present then calculate the compression ratios from the mean obliquity
                    r, r_mag, deltau = calculateCompressionRatios(randx, randy, obliquity_mean, 0.0, ripple_wavelength, ripplevx, ripplevy, r)
            
        #calculate new position
        x, mu, r_mag = moveParticle(x, v, mu, dt, deltau, r_mag)
        #add timestep to particle time
        t += dt
        
        
    final_info = []
    #Add final values of pitch angle, position and velocity of each particle    
    final_info = [particle, x, v, t, mu]
    return final_info
        

list_of_results = []
def writeResult(result):
    #Append result to a list of lists
    list_of_results.append(result)
        
#    #Also print status each time 1000 particles are simulated
#    if (len(list_of_results)) % 100 == 0:
#            print("\rNumber of particles simulated so far: " + str((len(list_of_results))), end='\r', flush=True)


def parallelSim(particle_number, t_max, v_up, starting_x, starting_v):
    #Seed random in starting parameter calculations
    np.random.seed(0)
    
    #Build the parameter matrix from arrays
    #Start with particle number, time and mu arrays
    particle = np.array(range(0, particle_number) )

    #Starting time is picked randomly from a unifrom distribution from 0.0 to t_max
    t_start = np.random.uniform(0.0, t_max, particle_number) 

    #Cosine of pitch angle is between [-1,1]
    mu = 2 * np.random.uniform(0.0, 1.0, particle_number) - 1.0

    #Start filling x and v with first allocating right size arrays
    x = np.zeros(particle_number)
    v = np.zeros(particle_number)
    
    #Give each particle a starting position and a starting velocity
    x.fill(starting_x)
    v.fill(starting_v)
    
    v, mu = changeFrame(v, mu, -1.0 * group.attrs['v_up'])

    #Finally stack the parameter arrays to form the parameter matix
    parameters = np.column_stack((particle, x, v, t_start, mu))

    #Write starting parameters to file  
    data = group.create_dataset("T=0", data=parameters)
    
#    ###For testing purposes
#    wantedparticle = 0
#    print(data[wantedparticle])
#    attrib = np.append(data[wantedparticle], [group.attrs['t_max'], 1, group.attrs['particles_total'], group.attrs['v_up'], group.attrs['obliquity_mean'], group.attrs['ripple_amplitude'],  group.attrs['ripple_wavelength'],  group.attrs['ripple_chance'], group.attrs['ripple_velocity_x'], group.attrs['ripple_velocity_y'], group.attrs['r'], group.attrs['meanfreepath']])
#    final = simulateParticle(attrib[0], attrib[1], attrib[2], attrib[3], attrib[4], attrib[5], attrib[6], attrib[7], attrib[8], attrib[9], attrib[10], attrib[11], attrib[12], attrib[13], attrib[14], attrib[15], attrib[16])
#    print(final)
    
    count = 1
    while count <= group.attrs['savesteps']:
        
        
        t_step = group.attrs['t_max']/group.attrs['savesteps']*count
        #open a pool of workers (without parameters uses os.cpu_count amount of workers)
        pool = Pool()
        
        #simulate all particles
        for param in data:
            
            attrib = param
            attrib = np.append(attrib, [t_step, count, group.attrs['particles_total'], group.attrs['v_up'], group.attrs['obliquity_mean'], group.attrs['ripple_amplitude'],  group.attrs['ripple_wavelength'],  group.attrs['ripple_chance'], group.attrs['ripple_velocity_x'], group.attrs['ripple_velocity_y'], group.attrs['r'], group.attrs['meanfreepath']])
            #Use pool of workers to distribute each particle simulation to free workers
            #use function simulatePartice with attributes that are give
            #Give the result of the run to function writeResult
            pool.apply_async(simulateParticle, attrib, callback = writeResult)
    
        
        #Close workerpool to free cores            
        pool.close()
        
        #wait for final runs to finish before continuing
        pool.join()
        
        #Take the list of results and format it to a numpy array
        particles = np.array(list_of_results)
        list_of_results[:] = []
        
        #Write the result to the .hdf5 file
        data = group.create_dataset("T=%d" %t_step, data=particles)
        count += 1
   

#__name__ == '__main__' needed by the multiprocessing library to not start infinately many programs
if __name__ == '__main__':
    #Simulation parameters
    r = simulation_parameters.r
    t_max = simulation_parameters.t_max
    savesteps = simulation_parameters.savesteps
    v_up = simulation_parameters.v_up
    obliquity_mean = simulation_parameters.obliquity_mean
    ripple_amplitude = simulation_parameters.ripple_amplitude
    ripple_wavelength = simulation_parameters.ripple_wavelength
    ripple_chance = simulation_parameters.ripple_chance
    ripplevx = - v_up * np.tan(obliquity_mean)
    ripplevy = simulation_parameters.ripplevy
    particle_number = simulation_parameters.particle_number 
    filename = simulation_parameters.filename
    groupname = simulation_parameters.groupname
    starting_x = simulation_parameters.x
    starting_v = simulation_parameters.v
    meanfreepath = simulation_parameters.meanfreepath

    
    try:
        
        testx = 0.0
        testy = np.pi / 2.0
        
        #Calculate the ripple strength
        a = 2.0 * np.pi * ripple_amplitude / ripple_wavelength
        
        #Cosines and Sines are needed in the obliquity calculation
        cosx = np.cos(testx)
        cosy = np.cos(testy)
        sinx = np.sin(testx)
        siny = np.sin(testy)
        
        #Calculate the semi-random obliquity of the shock
        numerator = np.cos(obliquity_mean) - a * np.sin(obliquity_mean) * cosx * siny
        denominator = np.sqrt( a*a*(cosx*cosx*siny*siny + cosy*cosy*sinx*sinx) + 1.0 )
        
        obliquity_test = np.arccos( numerator/denominator )
        obliquity_test = np.degrees(obliquity_test)
        
        if obliquity_test >= 90.0:
            raise Exception("Mean obliquity and ripple too big. Chance of {} angle when 90 is maximum. Change one or both!".format(obliquity_test))
        elif obliquity_test > 80.0:
            warnings.warn("Mean olbiquity and ripple large. Chance of over 80 degree obliquity. Program will continue to run.")
            
        
        #Open .hdf5 binary file to save starting parameters and create group for run results
        file = h5py.File(filename, "a")
    
        #Simulation runs can be saved to the same file using a different group each time
        group = file.create_group(groupname)
        group.attrs['r'] = r
        group.attrs['meanfreepath'] = meanfreepath
        group.attrs['t_max'] = t_max
        group.attrs['savesteps'] = savesteps
        group.attrs['v_up'] = v_up
        group.attrs['obliquity_mean'] = obliquity_mean
        group.attrs['ripple_amplitude'] = ripple_amplitude
        group.attrs['ripple_wavelength'] = ripple_wavelength
        group.attrs['particles_total'] = particle_number
        group.attrs['ripple_chance'] = ripple_chance
        group.attrs['ripple_velocity_x'] = ripplevx
        group.attrs['ripple_velocity_y'] = ripplevy
        
        #startt is used to time the simulation time
        startt = time.time()
        
        #Start the parallel simulation
        parallelSim(particle_number, t_max, v_up, starting_x, starting_v) 
        
        group.attrs['runtime'] = time.time() - startt

        #Close .hdf5 file
        file.close()
        
    except:           
        #Close .hdf5 file even if an error occures but before raising the error
        file.close()
        
        #Raise the error
        raise
        
    #Finally print runtime of the simulation
    print("\nRuntime: " + str(time.time() - startt) + "s")