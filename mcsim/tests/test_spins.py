import numpy as np
import pytest
import mcsim


class TestInitialisation:
    def test_init(self):
        n = (12, 15)
        value = (1, 0, 0)
        s = mcsim.Spins(n=n, value=value)

        assert s.n == n

        assert isinstance(s.array, np.ndarray)
        assert s.array.shape == (*n, 3)
        assert np.allclose(s.array, value)

    def test_init_wrong_n(self):
        n = (12, 0)
        value = (1, 0, 0)
        with pytest.raises(ValueError):
            s = mcsim.Spins(n=n, value=value)

    def test_init_wrong_value(self):
        n = (12, 0)
        value = (5, 6)
        with pytest.raises(ValueError):
            s = mcsim.Spins(n=n, value=value)


class TestRandomise:
    def test_randomise_component_values(self):
        n = (25, 25)
        s = mcsim.Spins(n=n)
        s.randomise()

        for i in range(s.array.shape[-1]):
            # Are all values of components between -1 and 1?
            assert np.greater_equal(s.array[..., i], -1).all()
            assert np.less_equal(s.array[..., i], 1).all()

    def test_randomise_mean(self):
        n = (100, 100)
        s = mcsim.Spins(n=n)
        s.randomise()

        # Is the mean somewhere around zero?
        assert all(-0.25 < s.mean[i] < 0.25 for i in range(3))

    def test_randomise_normalised(self):
        n = (15, 15)
        s = mcsim.Spins(n=n)
        s.randomise()

        assert np.allclose(abs(s), 1)


class TestMean:
    def test_mean_type(self):
        n = (22, 41)
        value = (0, 0, 1)
        s = mcsim.Spins(n=n, value=value)

        assert isinstance(s.mean, np.ndarray)
        assert len(s.mean) == 3

    def test_mean_uniform(self):
        n = (20, 20)
        value = (1, 1, 1)
        s = mcsim.Spins(n=n, value=value)

        assert np.allclose(s.mean, np.sqrt(3) / 3)

    def test_mean_non_uniform(self):
        n = (50, 50)
        s = mcsim.Spins(n=n)

        # Assign a non-uniform value to array.
        s.array = np.ones((*n, 3))
        s.array[:25:, :, :] = (0, 0, 1)
        s.array[25::, :, :] = (0, 0, -1)

        assert np.allclose(s.mean, (0, 0, 0))


class TestNormalise:
    def test_normalise_type(self):
        n = (100, 101)
        s = mcsim.Spins(n=n)

        assert abs(s).shape == (*n, 1)

    def test_normalise_rational(self):
        n = (22, 41)
        value = (0, 5, 0)
        s = mcsim.Spins(n=n, value=value)

        assert np.allclose(s.array, (0, 1, 0))
        assert np.allclose(abs(s), 1)

    def test_normalise_irrational(self):
        n = (15, 11)
        value = (1, 1, 1)
        s = mcsim.Spins(n=n, value=value)

        assert np.allclose(s.array, np.sqrt(3) / 3)
        assert np.allclose(abs(s), 1)


class TestPlot:
    def test_plot(self):
        # In this test, we are only ensuring we can run the plot method.
        n = (5, 6)
        s = mcsim.Spins(n=n)
        s.randomise()

        s.plot()  # There is no assert statement here.
