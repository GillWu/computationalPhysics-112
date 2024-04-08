import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from .particles import Particles 
from numba import jit, njit, prange, set_num_threads

"""
The N-Body Simulator class is responsible for simulating the motion of N bodies



"""

class NBodySimulator:

    def __init__(self, particles: Particles):

        # set up the properties of the particles
        self.particles = particles
        self.time      = particles.time
        # call defult settings in setup()
        self.setup()
        return

    def setup(self, G=1,
                    rsoft=0.01,
                    method="RK4",
                    io_freq=10,
                    io_header="nbody",
                    io_screen=True,
                    visualization=False):
        """
        Customize the simulation enviroments.

        :param G: the graivtational constant
        :param rsoft: float, a soften length
        :param meothd: string, the numerical scheme
                       support "Euler", "RK2", and "RK4"

        :param io_freq: int, the frequency to outupt data.
                        io_freq <=0 for no output. 
        :param io_header: the output header
        :param io_screen: print message on screen or not.
        :param visualization: on the fly visualization or not. 
        """
        self.G = G
        self.rsoft = rsoft
        self.method = method
        self.io_freq = io_freq
        self.io_header = io_header
        self.io_screen = io_screen
        self.visualization = visualization
        
        # set up the advance_particles method
        if method.lower() == "Euler":
            self._advance_particles = self._advance_particles_Euler
        elif method.lower() == "RK2":
            self._advance_particles = self._advance_particles_RK2
        elif method.lower() == "RK4":
            self._advance_particles = self._advance_particles_RK4
        else:
            print("Please enter the method (Euler, RK2, RK4)")
            return

        return

    def evolve(self, dt:float, tmax:float):
        """
        Start to evolve the system

        :param dt: float, the time step
        :param tmax: float, the total time to evolve
        
        """
        time = self.particles.time
        tsteps = np.ceil(tmax/dt)

        for n in range(tsteps):

            # make sure the last step will arrive at tmax
            if (time+dt) > tmax:
                dt = tmax - time

            # complicate physics model
            particles = self._advance_particles(dt, particles)

            # Output
            # output is not ready in particles.py
            #if (n % self.io_freq == 0):
            #    if self.io_screen:
            #        print("n=",n , "Time: ", time, " dt: ", dt)
            #    fn = self.io_header+"_"+str(n).zfill(6)+".dat"
            #    fn = io_folder+"/"+fn
            #    particles.output(fn)

            #update the time
            time += dt

        print("Simulation is done!")
        return

    def _calculate_acceleration(self, nparticles, masses, positions):
        """
        Calculate the acceleration of the particles
        """
        accelerations = np.zeros_like(positions)
        

        




        return accelerations
        
    def _advance_particles_Euler(self, dt, particles):
        nparticles = particles.nparticles
        mass = particles.masses
        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles, mass, pos)

        # do the Euler update
        pos = pos + vel*dt
        vel = vel + acc*dt
        acc = self._calculate_acceleration(nparticles, mass, pos)

        particles = particles.set_particles(mass, pos, vel, acc)

        return particles

    def _advance_particles_RK2(self, dt, particles):

        # TODO





        return particles

    def _advance_particles_RK4(self, dt, particles):
        
        #TODO








        return particles



if __name__ == "__main__":
    
    pass