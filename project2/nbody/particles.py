import numpy as np
import matplotlib.pyplot as plt


class Particles:
    """
    Particle class to store particle properties
    """
    def __init__(self, N:int = 100):
        """
        :param N: number of particles.
        """
        self.nparticles     = N
        self._masses        = np.ones((N,1))
        self._positions     = np.zeros((N,3))
        self._velocities    = np.zeros((N,3))
        self._accelerations = np.zeros((N,3))
        self._tags          = np.arange(N)
        self.time = 0
        return


    @property
    def tags(self):
        return self._tags
    
    @tags.setter
    def tags(self,some_name):
        if len(some_name) != self.nparticles:
            print("Name is too long or too short!!")
            raise ValueError
        self._tags = some_name
    
    @property
    def masses(self):
        return self._masses
    
    @masses.setter
    def masses(self,some_masses):
        size = some_masses.shape
        if size[1] != 1:
            print("Mass should be a N*1 array.")
            raise ValueError
        self._tags = some_masses
        return
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self,some_x):
        if len(some_x) != self.nparticles or some_x.shape[1] != 3:
            print("positions is too long or too short!!")
            raise ValueError
        self._positions = some_x

    @property
    def velocities(self):
        return self._velocities
    
    @velocities.setter
    def velocities(self,some_v):
        if len(some_v) != self.nparticles or some_v.shape[1] != 3:
            print("velocities is too long or too short!!")
            raise ValueError
        self._velocities = some_v
        return

    @property
    def accelerations(self):
        return self._accelerations
    
    @accelerations.setter
    def accelerations(self,some_a):
        if len(some_a) != self.nparticles or some_a.shape[1] != 3:
            print("accelerations is too long or too short!!")
            raise ValueError
        self._accelerations = some_a

    def add_particles(self, masses, positions, velocities, accelerations):
        M      = self._masses
        X      = self._positions
        V      = self._velocities
        A      = self._accelerations

        number = len(masses)
        self.nparticles += number
        self._masses = np.concatenate((M,masses),axis=0)
        self._positions = np.concatenate((X , positions),axis=0)
        self._velocities = np.concatenate((V, velocities),axis=0)
        self._accelerations = np.concatenate((A, accelerations),axis=0)
        return
    
    def set_particles(self, masses, positions, velocities, accelerations):
        self._masses = masses
        self._positions = positions
        self._velocities = velocities
        self._accelerations = accelerations
        return
    
    def output(self, filename):
        """
        Output particle properties to a file

        :param filename: output file name
        """
        masses = self.masses
        pos = self.positions
        vel = self.velocities
        acc = self.accelerations
        tags = self.tags
        time = self.time

        header = """
                ----------------------------------------------------
                Data from a 3D direct N-body simulation. 

                rows are i-particle; 
                coumns are :mass, tag, x ,y, z, vx, vy, vz, ax, ay, az

                NTHU, Computational Physics 

                ----------------------------------------------------
                """
        header += "Time = {}".format(time)
        np.savetxt(filename,(tags[:],masses[:,0],pos[:,0],pos[:,1],pos[:,2],
                            vel[:,0],vel[:,1],vel[:,2],
                            acc[:,0],acc[:,1],acc[:,2]),header=header)
        return


    
    def draw(self, dim=2):
        """
        Draw particles in 2,3D space
        """
        fig = plt.figure()

        if dim == 2:
            ax = fig.add_subplot(111)
            ax.scatter(self.positions[:,0], self.positions[:,1])
            ax.set_xlabel('X [code unit]')
            ax.set_ylabel('Y [code unit]')
            
        elif dim == 3:
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.positions[:,0], self.positions[:,1], self.positions[:,2])
            ax.set_xlabel('X [code unit]')
            ax.set_ylabel('Y [code unit]')
            ax.set_zlabel('Z [code unit]')
        else:
            print("Invalid dimension!")
            return

if __name__ =='__main__':
    time          = 0    # the starting  time
    num_particles = 100  # number of particles
    masses        = np.ones((num_particles,1))
    positions     = np.zeros((num_particles,3)) # 3 directions
    velocities    = np.zeros((num_particles,3))
    accelerations = np.zeros((num_particles,3))
    tags          = np.linspace(1,num_particles,num_particles)
    particles = Particles(N=num_particles)
    particles.output(filename='data.txt')