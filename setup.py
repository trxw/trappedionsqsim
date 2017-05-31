try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Quantum and Classical Simulations for Trapped Ions',
    'author': 'Omid Khosravani',
    'url': 'https://github.com/trxw/trappedionsqsim',
    'version': '0.1',
    'install_requires': ['qutip'],
    'packages': ['trappedionsqsim',],
    'scripts': ['bin/trappedionsqsim'],
    'name': 'trappedionsqsim'
}

setup(**config)
