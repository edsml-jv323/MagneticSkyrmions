from setuptools import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name="mcsim",  # Package name must be mcsim.
    version = "1.0.0",
    packages = ["mcsim"],
    install_requires=required,
    description = "A Python package for finding a magnetic skyrmion in a 2D lattice spins using the Metropolis algorithm",
    author = "Jorge Veiras Yanes",
    author_email = "jv323@ic.ac.uk"
)


