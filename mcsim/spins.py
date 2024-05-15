import numbers
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import mcsim




class Spins:
    """Field of spins on a two-dimensional lattice.

    Each spin is a three-dimensional vector s = (sx, sy, sz). Underlying data
    stucture (``self.array``) is a numpy array (``np.ndarray``) with shape
    ``(nx, ny, 3)``, where ``nx`` and ``ny`` are the number of spins in the x
    and y directions, respectively, and 3 to hold all three vector components of
    the spin.

    Parameters
    ----------
    n: Iterable

        Dimensions of a two-dimensional lattice ``n = (nx, ny)``, where ``nx``
        and ``ny`` are the number of atoms in x and y directions, respectively.
        Values of ``nx`` and ``ny`` must be positive integers.

    value: Iterable

        The value ``(sx, sy, sz)`` that is used to initialise all spins in the
        lattice. All elements of ``value`` must be real numbers. Defaults to
        ``(0, 0, 1)``.

    """

    def __init__(self, n, value=(0, 0, 1)):
        # Checks on input parameters.
        if len(n) != 2:
            raise ValueError(f"Length of iterable n must be 2, not {len(n)=}.")
        if any(i <= 0 or not isinstance(i, int) for i in n):
            raise ValueError("Elements of n must be positive integers.")

        if len(value) != 3:
            raise ValueError(f"Length of iterable value must be 3, not {len(value)=}.")
        if any(not isinstance(i, numbers.Real) for i in value):
            raise ValueError("Elements of value must be real numbers.")

        self.n = n
        self.array = np.empty((*self.n, 3), dtype=np.float64)
        self.array[..., :] = value

        if not np.isclose(value[0] ** 2 + value[1] ** 2 + value[2] ** 2, 1):
            # we ensure all spins' magnitudes are normalised to 1.
            self.normalise()

    @property
    def mean(self):

        """ Calculates the mean of each spin component (Sx, Sy, Sz), from all the
        spins in the 2D lattice. 
        
        It gives us an idea of the overall orientation of the spins in the lattice. 
        This aggregate measure, will aid us in the inspection of the emerged magnetisation state 
        at the end of the simulation

        Returns
        -------
        numpy.ndarray

            A 3 element numpy array containing the mean value for each spin component across
            all the spins in the lattice: ``(mean_Sx, mean_Sy, mean_Sz)``

        Unoptimised Approach 
        --------------------
        Computes the mean method using ``for`` loops. It may not be optimally efficient for
        larger lattices

            counter=np.zeros((3)) #(0,0,0)
            for i in range(self.n[0]):
                for j in range(self.n[1]):
                    counter=counter+self.array[i,j,:]
            m_c = 1/(self.n[0]*self.n[1])*counter
            return m_c
        """

        return np.mean(self.array, (0,1)) #Reference [2]

        

    def __abs__(self):
        """"Calculates the norm (magintude) of each spin in the 2D lattice

        As we are using the simplified and modified Heisenberg's spin model, one of the
        approximations made is that all the magnitudes of the spins must remain constant.
        This method will help us when computing the norm. This method calculates the 
        magintude of each spin, facilitating the normalisation process


        Returns
        -------
        numpy.ndarray
            A 3D array with a shape of (nx,ny,1) which each `(i,j,0)` element represents the magnitude
            of the spin at the position (i,j) in the lattice. 

        Unoptimised way 
        ---------------
        norms=np.zeros((self.n[0], self.n[1],1)) #(lattice)
        for i in range(self.n[0]):
            for j in range(self.n[1]):
                spins_values=self.array[i,j]
                norms[i,j,0]=(spins_values[0]**2+spins_values[1]**2+spins_values[2]**2)**(1/2)
        return norms
        """

        return np.sqrt(np.sum(self.array**2, axis=2, keepdims=True)) # Reference [3] for keepdims




    def normalise(self):
        """Normalise the magnitude of all spins to 1.

        This method modifies the `self.array` in place. After calling this method, all spins in the
        lattice will have a magnitude of 1.

        Returns
        -------
        It doesn't return any value, it just modifies the array
        """

        self.array = self.array / abs(self)  # This computation will be failing until you implement __abs__.




    def randomise(self):
        """Initialise the lattice with random spins.

        Components of each spin are between -1 and 1: -1 <= si <= 1, and all
        spins are normalised to 1.

        """
        self.array = 2 * np.random.random((*self.n, 3)) - 1
        self.normalise()





    def plot(self):
        """Plots the 3 components of the spins in a 2D lattice

        As we are in 2D dimensions, the Sz component values of the spins
        will be indicated with a colorbar.

        If run after applying driver method (Monte Carlo algorithm), the plot
        will show the configuration of the spins after reaching the equilibrium state

        If run before, the plot will show the random configuration of spins (initial state)


         Notes:
            - The colorbar represents the magnitude of the Sz component.
            - Arrows represent the direction and magnitude of the spins in the (nx,ny) plane.

        Returns
        -------
        None
            The function directly visualizes the spin configuration using matplotlib 
            but doesn't return any value.

        """
        nx, ny = np.meshgrid(range(self.n[0]), range(self.n[1]))
        Sx_final = self.array[..., 0]
        Sy_final = self.array[..., 1]
        Sz_final = self.array[..., 2]
        quiver=plt.quiver(nx, ny, Sx_final, Sy_final, Sz_final, pivot='middle', angles='xy', scale_units='xy', scale=1, cmap='plasma')
        plt.title("2D Lattice spins after appplying Metropolis algorithm")
        plt.colorbar(quiver, ax=plt.gca(), label="$S_z$ Component")
        plt.xlabel('nx')
        plt.ylabel('ny')
        plt.tight_layout()
        plt.show()
