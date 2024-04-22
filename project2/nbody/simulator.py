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
        


        return

    def evolve(self, dt:float, tmax:float):
        """
        Start to evolve the system

        :param dt: float, the time step
        :param tmax: float, the total time to evolve
        
        """
        time = self.particles.time
        tsteps = int(np.ceil(tmax/dt))

        # set up the advance_particles method
        method = self.method
        if method.lower() == "euler":
            self._advance_particles = self._advance_particles_Euler
        elif method.lower() == "rk2":
            self._advance_particles = self._advance_particles_RK2
        elif method.lower() == "rk4":
            self._advance_particles = self._advance_particles_RK4
        else:
            print("Please enter the method (Euler, RK2, RK4)")
            return
        
        # initialize the output
        io_folder = "data_" + self.io_header
        Path(io_folder).mkdir(parents=True, exist_ok=True)

        # start main loop
        for n in range(tsteps):

            # make sure the last step will arrive at tmax
            if (time+dt) > tmax:
                dt = tmax - time
         
            # Output
            if (n % self.io_freq == 0):
                if self.io_screen:
                    print("n=",n , "Time: ", time, " dt: ", dt)
                fn = self.io_header+"_"+str(n).zfill(6)+".dat"
                fn = io_folder + "/"+fn
                self.particles.output(fn)

            #update the time
            time += dt

        print("Simulation is done!")
        return

    def _calculate_acceleration(self, nparticles, masses, positions):
        """
        Calculate the acceleration of the particles
        """
        accelerations = np.zeros_like(positions)
        # insert acceleration from outside function
        acceleration = acc_gravity(nparticles,masses,positions,accelerations,self.G,self.rsoft)

        return acceleration
        
    def _advance_particles_Euler(self, dt, particles):
        nparticles = particles.nparticles
        mass = particles.masses
        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles, mass, pos)

        # do the Euler update
        y0 = np.array([pos, vel])
        yderv = np.array([vel, acc])
        ynext = y0 + yderv * dt
        pos = ynext[0]
        vel = ynext[1]
        acc = self._calculate_acceleration(nparticles, mass, pos)

        # update the particles
        particles = particles.set_particles(mass, pos, vel, acc)

        return particles

    def _advance_particles_RK2(self, dt, particles):
        nparticles = particles.nparticles
        mass = particles.masses
        # y0
        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles, mass, pos)
        
        # do the RK2 update
        # y1
        pos1 = pos + vel*dt
        vel1 = vel + acc*dt
        acc1 = self._calculate_acceleration(nparticles, mass, pos1)
        # y2
        pos2 = pos1 + vel1*dt
        vel2 = vel1 + acc1*dt
        # RK2 sum up
        pos = 0.5*(pos + pos2)
        vel = 0.5*(vel + vel2)
        acc = self._calculate_acceleration(nparticles, mass, pos)
        # update the particles
        particles = particles.set_particles(mass, pos, vel, acc)

        return particles

    def _advance_particles_RK4(self, dt, particles):
        
        nparticles = particles.nparticles
        mass = particles.masses
        dt2 = dt/2
        # do the RK4 update
        # y0
        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles, mass, pos)
        # y1
        pos1 = pos + vel*dt2
        vel1 = vel + acc*dt2
        acc1 = self._calculate_acceleration(nparticles, mass, pos1)
        # y2
        pos2 = pos + vel1*dt2
        vel2 = vel + acc1*dt2
        acc2 = self._calculate_acceleration(nparticles, mass, pos2)
        # y3
        pos3 = pos + vel2*dt
        vel3 = vel + acc2*dt
        acc3 = self._calculate_acceleration(nparticles, mass, pos3)
        # RK4 sum up
        pos = pos + (vel + 2*vel1 + 2*vel2 + vel3)*dt/6
        vel = vel + (acc + 2*acc1 + 2*acc2 + acc3)*dt/6
        acc = self._calculate_acceleration(nparticles, mass, pos)
        # update the particles
        particles = particles.set_particles(mass, pos, vel, acc)

        return particles

@njit(parallel=True)
def acc_gravity(nparticles, mass, pos, acc, G, rsoft):
    """
    Calculate the acceleration due to gravity.
    """
    for i in prange(nparticles):
        for j in prange(nparticles):
            if i != j:
                rij = pos[i,:] - pos[j,:]
                r = np.sqrt(np.sum(rij**2) + rsoft**2)
                runit = rij/r
                force = -(G * mass[i,0] * mass[j,0] /r**2 )*runit
                acc[i,:] += force[:]/mass[i,0]
                acc[j,:] -= force[:]/mass[j,0]

if __name__ == "__main__":
    
    pass