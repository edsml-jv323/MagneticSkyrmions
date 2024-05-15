import numpy as np


class System:
    """System object with the spin configuration and necessary parameters.

    Parameters
    ----------
    s: mcsim.Spins

        Two-dimensional spin field.

    B: Iterable

        External magnetic field, length 3.

    K: numbers.Real

        Uniaxial anisotropy constant.

    u: Iterable(float)

        Uniaxial anisotropy axis, length 3. If ``u`` is not normalised to 1, it
        will be normalised before the calculation of uniaxial anisotropy energy.

    J: numbers.Real

        Exchange energy constant.

    D: numbers.Real

        Dzyaloshinskii-Moriya energy constant.

    """

    def __init__(self, s, B, K, u, J, D):
        self.s = s
        self.J = J
        self.D = D
        self.B = B
        self.K = K
        self.u = u

    def energy(self):
        """Total energy of the system.

        The total energy of the system is computed as the sum of all individual
        energy terms.

        Returns
        -------
        float

            Total energy of the system.

        """
        return self.zeeman() + self.anisotropy() + self.exchange() + self.dmi()
    


    def zeeman(self):
        """" Calculates the Zeeman energy of the entire 2D lattice.

        The Zeeman energy represents the interaction of the magnetic moments
        in the lattice with an external magnetic field. It is given by the following
        expression:

        E_Zeeman = - Σ s_i,j · B

        where:
            Σ : represents the sum over all spins in the lattice
            s_i,j : represents the spin at the (i,j) position in the lattice
            B : is the external magentic field (vector)


        Attributes Used
        ---------------
        self.s.array
            3D array representing the spins in the lattice (numpy.ndarray)
        self.B
            A vector representing the external magnetic field applied (numpy.ndarray)

        Returns
        -------
        float
            The total Zeeman energy of the system
   
        Unoptimised approach
        --------------------
        The nested for loops are much less efficient than the vectorised approach
        using numpy

        E_Zeeman=0
        for i in range(self.s.n[0]): #Reference [1] Composition relationship
            for j in range(self.s.n[1]):
                spins_values=self.s.array[i,j]
                E_Zeeman=E_Zeeman-(spins_values[0]*self.B[0]+spins_values[1]*self.B[1]+spins_values[2]*self.B[2])
        return E_Zeeman
        """

        E_Zeeman = -np.sum(self.s.array * self.B)
        return E_Zeeman




    def anisotropy(self):
        """"Calculates the anisotropy energy of the system.

        The magentic anisotropy energy represents the tendency of spins
        to favor magentisation along specific directions (anysotropy axes).

        It is given by:

            e_a = -K * Σ (s_i,j · u)^2

        where:
            Σ : represents the sum over all the spins in the lattice
            s_i,j : represents the spin at the (i,j) position in the lattice
            u : is the preferred direction for magnetisation (anisotropy axis) (unit vector)
            K : is the anisotropy constant

        Attributes used
        ---------------
        self.s.array
            3D array representing the spins in the lattice (numpy.ndarray)
        self.u 
            Vector representing the preferred direction for magnetisation (numpy.ndarray)
        self.k
            Anisotropy constant (float)
        

        Returns
        -------
        float
            The total magnetic anisotropy energy of the system

        Unoptimised way
        ---------------
        norm_u=np.sqrt(self.u[0]**2+self.u[1]**2+self.u[2]**2) #We normalise
        if norm_u != 1:
            self.u=(self.u[0]/norm_u, self.u[1]/norm_u, self.u[2]/norm_u) 
        
        e_a=0
        for i in range(self.s.n[0]): #Reference [1] Composition relationship
            for j in range(self.s.n[1]):
                spins_values=self.s.array[i,j]
                e_a=e_a-self.K*(spins_values[0]*self.u[0]+spins_values[1]*self.u[1]+spins_values[2]*self.u[2])**(2)
        return e_a
        """

        if not np.isclose(np.linalg.norm(self.u), 1):  #Reference [4]
            self.u = self.u / np.linalg.norm(self.u)

        e_a = -self.K * np.sum((np.sum(self.s.array * self.u, axis=2))**2)
        return e_a


    def exchange(self):
        """"Calculates the exchange energy of the 2D lattice.

        The exchange energy arises from interactions between nearest
        neighbours spins. For a given pair of neigbouring spins, the exchange 
        energy is:

            E_ex = -J * (s_i,j · s_k,l)

        where: 
            s_i,j and s_k,l are first neighboring spins in the lattice
            J is the exchange coupling constant
            The dot product indicates we are performing the scalar product between
            two spins

        The total exchange energy is obtained by summing up the exchange energies of all 
        nearest-neighbor pairs in the lattice. In this implementation, the method loops
        over all spin pairs to compute the sum. In further versions, the exchange energy 
        calculation will be optimised.

         Attributes used
         ---------------
         self.s.array
            3D array representing the spins in the lattice  (numpy.ndarray)
         self.J 
            The exchange coupling constant (float)
           

         Returns
         -------
         float
            The total exchange energy of the system

        """

    #We just apply distributive property
        E_ex = 0
        for i in range(self.s.n[0]):
            for j in range(self.s.n[1]-1):
                first_part = self.s.array[i,j][0]*self.s.array[i, j+1][0]+self.s.array[i,j][1]*self.s.array[i, j+1][1]+self.s.array[i,j][2]*self.s.array[i, j+1][2]
                E_ex = E_ex - first_part * self.J

        for j in range(self.s.n[1]):
            for i in range(self.s.n[0]-1):
                second_part = self.s.array[i,j][0]*self.s.array[i+1, j][0]+self.s.array[i,j][1]*self.s.array[i+1, j][1]+self.s.array[i,j][2]*self.s.array[i+1, j][2]
                E_ex = E_ex - second_part * self.J

        return E_ex
        

    
    def dmi(self):
        """" Calculates the Dzyaloshinskii-Moriya interaction (DMI) energy
        for the 2D lattice.
        
        The DMI energy represents the interaction energy due to the antisymmetric exchange 
        between neighboring spins. It arises from the combination of spin-orbit coupling 
        and broken inversion symmetry in magnetic structures.

        For the mathematical expression check Reference [5]

        Attributes Used
        ---------------
        self.s.array 
            A 3D array representing the spins in the lattice (numpy.ndarray)
        self.D 
                The DMI constant (float)

        Returns
        -------
        float
            The total DMI energy of the system

        Unoptimised way
        ---------------
        res1 = 0
        res2 = 0
        r_ij_x=np.array([1,0,0]) #Vector from s[i,j] to s[i, j+1]
        for i in range(self.s.n[0]):
            for j in range(self.s.n[1]-1):
                cross_product_1=-np.dot(r_ij_x, np.cross(self.s.array[i,j], self.s.array[i,j+1]))
                res1 = res1 + cross_product_1 

        r_ij_y=np.array([0,1,0]) #Vector from s[i,j] to s[i+1, j]
        for j in range(self.s.n[1]):
            for i in range(self.s.n[0]-1):
                cross_product_2=-np.dot(r_ij_y, np.cross(self.s.array[i,j], self.s.array[i+1,j]))
                res2 = res2 + cross_product_2 


        return -self.D*(res1 + res2)   
        """
        
        r_ij_x = np.array([1, 0, 0]) 
        r_ij_y = np.array([0, 1, 0])  

        first_cross_product = np.cross(self.s.array[:, :-1, :], self.s.array[:, 1:, :])
        second_cross_product = np.cross(self.s.array[:-1, :, :], self.s.array[1:, :, :])

        result1 = -np.sum(r_ij_x * first_cross_product)
        result2 = -np.sum(r_ij_y * second_cross_product)

        return -self.D * (result1 + result2)
      
