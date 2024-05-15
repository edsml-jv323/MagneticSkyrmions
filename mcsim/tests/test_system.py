import numbers

import numpy as np

import mcsim


class TestInitialisation:
    def test_init(self):
        n = (100, 101)
        value = (0, 0, 1)
        s = mcsim.Spins(n=n, value=value)

        B = (0, 0, 1)
        K = 0.01
        u = (0, 1, 0)
        J = 0.6
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert np.allclose(system.s.mean, (0, 0, 1))
        assert np.isclose(system.J, J)
        assert np.isclose(system.D, D)
        assert np.isclose(system.K, K)
        assert np.allclose(system.u, u)
        assert np.allclose(system.B, B)


class TestZeeman:
    def test_zeeman_zero(self):
        n = (100, 101)
        s = mcsim.Spins(n=n)
        s.randomise()

        B = (0, 0, 0)
        K = 0.01
        u = (0, 1, 0)
        J = 0.6
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert isinstance(system.zeeman(), numbers.Real)
        assert np.isclose(system.zeeman(), 0)

    def test_zeeman_uniform(self):
        n = (10, 10)
        value = (0, 0, 1)
        s = mcsim.Spins(n=n, value=value)

        B = (0, 0, 0.5)
        K = 0.01
        u = (0, 1, 0)
        J = 0.6
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert np.isclose(system.zeeman(), -50)

    def test_zeeman_non_uniform(self):
        n = (10, 10)
        s = mcsim.Spins(n=n)
        s.array[0:5:, ...] = (0, 0, 1)
        s.array[5::, ...] = (0, 0, -1)

        B = (0, 0, 1)
        K = 0.01
        u = (0, 1, 0)
        J = 0.6
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert np.isclose(system.zeeman(), 0)


class TestAnisotropy:
    def test_anisotropy_zero(self):
        n = (50, 51)
        s = mcsim.Spins(n=n)
        s.randomise()

        B = (0, 0, 1)
        K = 0
        u = (0, 1, 0)
        J = 0.6
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert isinstance(system.zeeman(), numbers.Real)
        assert np.isclose(system.anisotropy(), 0)

    def test_anisotropy_uniform(self):
        n = (10, 10)
        value = (0, 0, 1)
        s = mcsim.Spins(n=n, value=value)

        B = (0, 0, 5)
        K = 1
        u = (0, 1, 0)
        J = 0.6
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert np.isclose(system.anisotropy(), 0)

    def test_anisotropy_non_uniform(self):
        n = (10, 10)
        s = mcsim.Spins(n=n)
        s.array[0:5:, ...] = (0, 0, 1)
        s.array[5::, ...] = (0, 0, -1)

        B = (0, 0, -2)
        K = 5
        u = (0, 0, 2)  # This value should be normalised before computing the energy.
        J = 0.6
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert np.isclose(system.anisotropy(), -500)


class TestExchange:
    def test_exchange_zero(self):
        n = (10, 10)
        s = mcsim.Spins(n=n)
        s.randomise()

        B = (0, 0, 1)
        K = 1
        u = (0, 1, 0)
        J = 0.0
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert isinstance(system.exchange(), numbers.Real)
        assert np.isclose(system.exchange(), 0)

    def test_exchange_uniform(self):
        n = (5, 6)
        value = (0, 0, 1)
        s = mcsim.Spins(n=n, value=value)

        B = (0, 0, 1)
        K = 1
        u = (0, 1, 0)
        J = 1
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert np.isclose(system.exchange(), -49)

    def test_exchange_non_uniform(self):
        n = (5, 6)
        s = mcsim.Spins(n=n)
        s.array[:, :3:, :] = (0, 1, 0)
        s.array[:, 3::, :] = (1, 0, 0)

        B = (0, 0, 1)
        K = 1
        u = (0, 1, 0)
        J = 1
        D = 0.7

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert np.isclose(system.exchange(), -44)


class TestDMI:
    def test_dmi_zero(self):
        n = (10, 10)
        s = mcsim.Spins(n=n)
        s.randomise()

        B = (0, 0, 1)
        K = 1
        u = (0, 1, 0)
        J = 1.5
        D = 0

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert isinstance(system.dmi(), numbers.Real)
        assert np.isclose(system.dmi(), 0)

    def test_dmi_uniform(self):
        n = (6, 5)
        s = mcsim.Spins(n=n)

        B = (0, 0, 1)
        K = 1
        u = (0, 1, 0)
        J = 1.5
        D = 1

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert isinstance(system.dmi(), numbers.Real)
        assert np.isclose(system.dmi(), 0)

    def test_dmi_non_uniform(self):
        n = (5, 6)
        s = mcsim.Spins(n=n)
        s.array[:, :3:, :] = (0, 1, 0)
        s.array[:, 3::, :] = (0, 0, 1)

        B = (0, 0, 1)
        K = 1
        u = (0, 1, 0)
        J = 1
        D = 1

        system = mcsim.System(s=s, B=B, K=K, u=u, J=J, D=D)

        assert np.isclose(system.dmi(), 5)
