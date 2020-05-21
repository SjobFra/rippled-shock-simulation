#
#Functions called in the main simulation program
#
import numpy as np
from numba import jit


@jit(nopython=True)
def changeFrame(v, mu, u):
    
    #Caluclate new velocity
    v_prime = np.sqrt((v*mu - u)**2 + (v**2)*(1 - mu**2))
    
    #Calculate new pitch angle
    mu_prime = (v * mu - u) / v_prime
    
    return v_prime, mu_prime


@jit(nopython=True)
def scatterParticle(dt, v, mu, mfp, thetaranvalue, phiranvalue):
    
    #calculate scattering angles
    theta = np.sqrt( -2.0 * dt *  (v / mfp) * np.log( 1.0 - thetaranvalue ) )
    phi = 2.0 * np.pi * phiranvalue

    #Calculate new pitch angle
    mu = mu * np.cos(theta) + np.sqrt(1.0 - (mu**2)) * np.sin(theta) * np.cos(phi)
    
    return mu

@jit(nopython=True)
def calculateRippleX(randx, randy, obliquity_mean, ripple_amplitude):
    
    lastx = randx
    
    #Should be a/2*pi, but with ripple_wavelength = 1 then a/2*pi = ripple_amplitude
    constants = ripple_amplitude * np.tan(obliquity_mean) * np.sin(2 * np.pi * randy)
    
    x = randx + constants * np.sin(2 * np.pi * lastx)
    
    #Keep iterating calculating x until it is the same as the last one to 6 decimals
    while (np.abs(x - lastx) > 1e-6):
        
        lastx = x
        x = randx + constants * np.sin(2 * np.pi * lastx)
    
    return x


@jit(nopython=True)
def calculateCompressionRatios(x, y, obliquity_mean, ripple_amplitude, ripple_wavelength, ripplevx, ripplevy, r):
    
    x = calculateRippleX(x, y, obliquity_mean, ripple_amplitude)
    
    #Calculate the ripple strength
    a = 2.0 * np.pi * ripple_amplitude / ripple_wavelength
    
    #Cosines and Sines are needed in the obliquity calculation
    cosx = np.cos(2 * np.pi * x)
    cosy = np.cos(2 * np.pi * y)
    sinx = np.sin(2 * np.pi * x)
    siny = np.sin(2 * np.pi * y)
    
    #Calculate the semi-random obliquity of the shock
    numerator = np.cos(obliquity_mean) - a * np.sin(obliquity_mean) * cosx * siny
    denominator = np.sqrt( a*a*(cosx*cosx*siny*siny + cosy*cosy*sinx*sinx) + 1.0 )
    
    obliquity = np.arccos( numerator/denominator )
    
    #Tangent of the obliquity is needed in the magnetic compression calculation
    obtan = np.tan(obliquity)
    
    #Calculate the new magnetic compression ratio
    r_mag = np.sqrt(1.0 + r*r*obtan*obtan) / np.sqrt(1.0 + obtan*obtan)
    
    #Calculate the change in velocity
    numr = a * (ripplevx * cosx * siny + ripplevy * sinx * cosy)
    denr = a * np.sin(obliquity_mean) * cosx * siny - np.cos(obliquity_mean)
    
    deltau = numr/denr
    
    return r, r_mag, deltau

@jit(nopython=True)
def moveParticle(x, v, mu, dt, deltau, r_mag):
    #Calculate threshold value for mu
    mu_0 = (-x)/(v*dt)
    
    #Check if the particle is in the upstream
    if x < 0:
        
        #Also if the pitch angle is above threshold then shock crossing happens
        if mu > mu_0:
            
            #Check if there is a ripple
            if deltau != 0.0:                    
                
                #Change from average shock frame to rippled shock frame
                v, mu = changeFrame(v, mu, -deltau)
            
            #Calculate time spent on this timestep before crossing
            dt_before = (-x)/(v*mu)
            
            #Particle transimtted through the shock
            if mu > np.sqrt( 1 - (1/r_mag)):
                mu = np.sqrt(1 - (1 - mu**2)*r_mag)
                x_new = v*mu*(dt - dt_before)
                
            #Paricle is reflected at the shock
            else:
                mu = -mu
                x_new = v*mu*(dt - dt_before)
                
            #Check if there was a ripple
            if deltau != 0.0:
                
                #Change from rippled shock frame to average shock frame
                v, mu = changeFrame(v, mu, deltau)
        
        #If no shock crossing has happened then just move the particle
        else:
            x_new = x + v*mu*dt
            
    #If particle is in the downstream
    else:
        
        #Again if pitch angle is above threshold then shock crossing happens
        #Going back upstream there is no possibility of reflection
        if mu < mu_0:
            dt_before = -x/(v*mu)
            mu = (-1)*np.sqrt(1 - (1 - mu**2)*(1/r_mag))
            x_new = v*mu*(dt - dt_before)
            
        #If no shock crossing has happened then just move the particle
        else:
            x_new = x + v*mu*dt
    
    #Return new position and pitch angle
    return x_new, mu, r_mag