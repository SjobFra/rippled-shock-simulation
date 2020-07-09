"""
Functions called by the rippled_shock.py simulation.
"""

import numpy as np
from numba import jit


@jit(nopython=True)
def changeFrame(v, mu, u):
    """
    Changes the frame of reference for the particle from unprimed to primed frame.

    Parameters
    ----------
    v : float
        Velocity of the particle in the starting frame.
    mu : float
        Pitch angle cosine of the particle in the starting frame.
    u : float
        The velocity difference between the two frames.

    Returns
    -------
    v_prime : float
        The velocity of the particle in the second frame.
    mu_prime : float
        The pitch angle cosine of the particle in the second frame.

    """
    
    #Caluclate new velocity
    v_prime = np.sqrt((v*mu - u)**2 + (v**2)*(1 - mu**2))
    
    #Calculate new pitch angle
    mu_prime = (v * mu - u) / v_prime
    
    return v_prime, mu_prime


@jit(nopython=True)
def scatterParticle(dt, v, mu, mfp, thetaranvalue, phiranvalue):
    """
    Scatters the particle by calculating a new pitch angle cosine for the particle.

    Parameters
    ----------
    dt : float
        The length of the current timestep taken by the particle.
    v : float
        The velocity of the particle.
    mu : float
        Current pitch angle cosine of the particle.
    mfp : float
        Mean free path of the particle.
    thetaranvalue : float
        A random float between (0, 1]. Numba changes the function to machine code and so random values cannot be taken inside this function.
    phiranvalue : float
        A random float between (0, 1]. Numba changes the function to machine code and so random values cannot be taken inside this function.

    Returns
    -------
    mu : float
        New pitch angle cosine for the particle.

    """
    
    #calculate scattering angles
    theta = np.sqrt( -2.0 * dt *  (v / mfp) * np.log( 1.0 - thetaranvalue ) )
    phi = 2.0 * np.pi * phiranvalue

    #Calculate new pitch angle
    mu = mu * np.cos(theta) + np.sqrt(1.0 - (mu**2)) * np.sin(theta) * np.cos(phi)
    
    return mu

@jit(nopython=True)
def calculateRippleX(randx, randy, obliquity_mean, ripple_amplitude):
    """
    Calculates the x-coordinate of the rippled crossing position iteratively.

    Parameters
    ----------
    randx : float
        A random float between (0, 1]. Numba changes the function to machine code and so random values cannot be taken inside this function.
    randy : float
        A random float between (0, 1]. Numba changes the function to machine code and so random values cannot be taken inside this function.
    obliquity_mean : float
        Angle between the smooth shock normal and the magnetic field.
    ripple_amplitude : float
        Ripple amplitude.

    Returns
    -------
    x : float
        The x-coordinate of the rippled shock crossing point.

    """
    
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
    """
    Calculates the magnetic compression ratio and the shock surface movement speed along the particle trajectory at a randomly chosen point at the ripple surface.

    Parameters
    ----------
    x : float
        A random float between (0, 1]. Numba changes the function to machine code and so random values cannot be taken inside this function.
    y : float
        A random float between (0, 1]. Numba changes the function to machine code and so random values cannot be taken inside this function.
    obliquity_mean : float
        Angle between the smooth shock normal and the magnetic field.
    ripple_amplitude : float
        Ripple amplitude.
    ripple_wavelength : float
        Ripple wavelength.
    ripplevx : float
        Ripple surface movement speed in the x-direction.
    ripplevy : float
        Ripple surface movement speed in the y-direction.
    r : float
        Gas compression ratio of the shock.

    Returns
    -------
    r : float
        Gas compression ratio of the shock.
    r_mag : float
        Magnetic compression ratio of the shock.
    deltau : float
        Movement speed of the rippled shock surface along the particle trajectory.

    """
    
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
    """
    Moves the particle. If the particle can cross the shock during the movement checks if the particle is transmitted through or mirrored at the shock.

    Parameters
    ----------
    x : float
        The current position of the particle.
    v : float
        The current velocity of the particle.
    mu : float
        The current pitch angle cosine of the particle.
    dt : float
        The timestep of the simulation.
    deltau : float
        The movement speed of the rippled shock surface along the particle trajectory.
    r_mag : float
        The magnetic compression ratio of the shock.

    Returns
    -------
    x_new : float
        The position of the particle after movement.
    mu : float
        The pitch angle cosine of the particle after movement.
    r_mag : float
        The magnetic compression ratio of the shock.

    """
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