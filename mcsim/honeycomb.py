import numbers
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.patches as patches
import numpy as np


class HoneycombSpins:
    def __init__(self, n, value=(0,0,1)):
        """Exactly same as 2D lattice
        But when need to take into account that there are 2 atoms (A, B)
        in each unit cell, so we need to distinguish two spins.
        """
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
        self.array = np.empty((self.n[0], self.n[1], 2, 3), dtype=np.float64)
        self.array[..., :] = value

        """
        So here self.array[i,j,0,:] = Spin of the atom A in (i,j)
                self.array[i,j,1,:] = Spin of the atom B in (i,j)

        """

    def plot_honeycomb(self):
        """Plots the 3 components of the spins on a honeycomb lattice"""

        fig, ax = plt.subplots()
        a = 1 #Lattice constant
        delta_y = a * np.sqrt(3)/2

        nx_A, ny_A = np.meshgrid(range(self.n[0]), range(self.n[1]))
        nx_B, ny_B = nx_A + 0.5, ny_A + delta_y  #HoneyComb basis 



        Sx_A, Sy_A, Sz_A = self.array[..., 0, 0], self.array[..., 0, 1], self.array[..., 0, 2]
        Sx_B, Sy_B, Sz_B = self.array[..., 1, 0], self.array[..., 1, 1], self.array[..., 1, 2]

    
        quiver_A = plt.quiver(nx_A, ny_A, Sx_A, Sy_A, Sz_A, pivot='middle', angles='xy', scale_units='xy', scale=1, cmap='plasma', alpha=0.6)
        quiver_B = plt.quiver(nx_B, ny_B, Sx_B, Sy_B, Sz_B, pivot='middle', angles='xy', scale_units='xy', scale=1, cmap='plasma', alpha=0.6)

        plt.title("Honeycomb lattice spins")
        plt.colorbar(quiver_A, ax=plt.gca(), label="$S_z$ Component")  # Colorbar based on A sublattice (but it's the same for both)
        plt.xlabel('nx')
        plt.ylabel('ny')
        plt.tight_layout()
        plt.show()

    def randomise(self):
        """Randomise the spins on the lattice and normalize them."""
        self.array = 2 * np.random.random((*self.n, 2, 3)) - 1
        norms = np.sqrt(np.sum(self.array ** 2, axis=-1, keepdims=True))
        self.array /= norms

    def mean(self):
        """Calculate the mean spin value (Same as in Spins class)"""
        return np.mean(self.array, axis=(0, 1, 2))

    def magnetization(self):
        """Calculate the net magnetization of the lattice (Same as in Spins class)."""
        return np.sum(self.array, axis=(0, 1, 2))