import numpy as np

import mcsim

rtol = 0.01
atol = 0.01


class TestIndividualEnergies:
    def test_zeeman(self):
        n = (5, 5)
        s = mcsim.Spins(n=n)
        s.randomise()

        B = (1, 0, 0)
        K = 0
        u = (0, 1, 0)
        J = 0
        D = 0

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        driver = mcsim.Driver()
        driver.drive(system, n=10_000)

        assert np.allclose(system.s.mean, (1, 0, 0), rtol=rtol, atol=atol)
        assert np.allclose(abs(system.s), 1, rtol=rtol, atol=atol)

    def test_anisotropy(self):
        n = (5, 5)
        s = mcsim.Spins(n=n)
        s.randomise()

        B = (0, 0, 0)
        K = 1
        u = (0, 0, 1)
        J = 0
        D = 0

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        driver = mcsim.Driver()
        driver.drive(system, n=10_000)

        assert np.allclose(abs(system.s.array[..., -1]), 1, rtol=rtol, atol=atol)
        assert np.allclose(abs(system.s), 1, rtol=rtol, atol=atol)

    def test_exchange(self):
        n = (5, 5)
        s = mcsim.Spins(n=n)
        s.randomise()

        B = (0, 0, 0)
        K = 0
        u = (0, 0, 1)
        J = 1
        D = 0

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        driver = mcsim.Driver()
        driver.drive(system, n=100_000)

        assert np.allclose(system.s.array, system.s.mean, rtol=rtol, atol=atol)
        assert np.allclose(abs(system.s), 1, rtol=rtol, atol=0.1)
