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
    
    def output(self,filename):
        """
        This function is not finished yet.
        It will writes down the properties of particles in a txt file.
        """
        # 將所有數組轉換為字符串
        tags = np.array(self.tags, dtype=str).reshape(-1, 1)
        masses = np.array(self.masses, dtype=str).reshape(-1, 1)
        pos = np.array(self.positions, dtype=str).reshape(-1, 1)
        vel = np.array(self.velocities, dtype=str).reshape(-1, 1)
        acc = np.array(self.accelerations, dtype=str).reshape(-1, 1)

        header = "# time,tag,mass,x,y,z,vx,vy,vz,ax,ay,az\n"
        header+= "# s,,kg,m,m,m,m/s,m/s,m/s,m/s^2,m/s^2,m/s^2\n"
        # 儲存數組
        np.savetxt(filename, np.hstack((np.ones((self.nparticles,1))*self.time, 
                                             tags, masses, pos, vel, acc)), 
                                             fmt='%s', header=header)
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
    num_particles = 100
    masses        = np.ones((num_particles,1))
    print(len(masses))