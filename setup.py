from setuptools import setup, find_packages
import sys, os.path


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'trappedionsqsim',
    'description': 'Quantum and Classical Simulations for Trapped Ions',
    'author': 'Omid Khosravani, Daniel Murphy, Lauren Lopez',
    'url': 'https://github.com/trxw/trappedionsqsim',
    'version': '0.1',
    'install_requires': ['qutip', 'numpy', 'matplotlib', 'scipy'],
    'packages': ['trappedionsqsim', 'trappedionsqsim.utils']
}

setup(**config)
