try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='OxidizeME',
      version='0.01',
      description='ME model with ROS stress damage and repair',
      author='Laurence Yang',
      author_email='',
      url='https://github.com/SBRG/oxidizeme',
      packages=['oxidizeme'],
      )
