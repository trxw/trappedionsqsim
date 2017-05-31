try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'trappedionsqsim',
    'author': 'Omid Khosravani',
    'url': 'https://github.com/trxw/trappedionsqsim',
    'version': '0.1',
    'install_requires': ['qutip'],
    'packages': ['trappedionsqsim',],
    'name': 'trappedionsqsim'
}

setup(**config)
