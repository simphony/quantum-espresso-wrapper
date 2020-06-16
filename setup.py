from setuptools import setup, find_packages

#Read description
with open('README.md', 'r') as readme:
    README_TEXT = readme.read()

#Main setup configuration class
setup(
    name = 'quantum-espresso',
    version = '0.1',
    author = 'Materials Informatics team, Fraunhofer IWM',
    description = 'Simulation wrapper for Quantum Espresso/SimPhoNy',
    long_description = README_TEXT,
    packages = find_packages(),
    test_suite = 'tests',
    
)