import numpy as np

def random_spin(s0, alpha=0.1):
    """Generate a new random spin based on the original one.

    Parameters
    ----------
    s0: np.ndarray

        The original spin that needs to be changed.

    alpha: float

        Larger alpha, larger the modification of the spin. Defaults to 0.1.

    Returns
    -------
    np.ndarray

        New updated spin, normalised to 1.

    """
    delta_s = (2 * np.random.random(3) - 1) * alpha
    s1 = s0 + delta_s
    return s1 / np.linalg.norm(s1)


class Driver:
    """Driver class.

    Driver class does not take any input parameters at initialisation.

    """

    def __init__(self):
        pass

    def drive(self, system, n, alpha=0.1):
        """Executes the Metropolis Monte Carlo algorithm.
        The purpose is to simulate the system's behaviour and bring it closer to the equilibrium by 
        probabilistically accepting or rejecting proposed spin updates based on their impact on 
        the system's energy.
    
        For each iteration, a random spin in the system is chosen, and its value is slightly
        perturbed. The system's energy is computed before and after the perturbation. The 
        perturbation is accepted if it lowers the system's energy otherwise, it is rejected.
        (We are not implementing the probabilistic acceptation with Boltzmann factor, as T=0)

        Parameters
        ----------
        system 
            An instance of the system whose spins are to be updated. This object has a 'energy' 
            method to compute its energy and an 's' attribute representing its spins (object)
        
         
            Number of iterations the Metropolis Monte Carlo algorithm should be run for (int)

        alpha 
            A parameter passed to the 'random_spin' function which determines the magnitude
            of the spin perturbation. Larger 'alpha' results in larger spin perturbations.
            Defaults to 0.1 (float, optional)

        Returns
        -------

        The method updates the spins in-place and does not return any value. The implementation 
        provided accepts changes that lower the energy but does not include probabilistic acceptance 
        based on the Boltzmann factor.
        """

        for _ in range(n):
            E_0=system.energy() 
            i = np.random.randint(0, system.s.n[0])
            j = np.random.randint(0, system.s.n[1])
            s0 = system.s.array[i,j].copy() #We create the random spin
            s_0_changed = random_spin(s0, alpha)
            system.s.array[i,j] = s_0_changed #We change the value of that random spin
            E_1 = system.energy()
            delta_E = E_1 - E_0
            if delta_E < 0:
                E_0 = E_1 #We accept the change
            else:
                system.s.array[i,j]=s0

