
from setuptools import setup, find_packages
#import "/home/yosef/Desktop/Organizing DataSearch Code"
import os

# parse dependencies and their links from requirements.txt files
install_requires = [line.rstrip() for line in open('requirements.txt')]
#tests_require = [line.rstrip() for line in open('tests/requirements.txt')]

# install package
setup(
    name='mypackage',
    packages = ['mypackage'],
    version = '0.1',
    description="finds most relevant kinetic data",
    url="https://github.com/KarrLab/Kinetic-Datanator",
    download_url='https://github.com/KarrLab/Kinetic-Datanator/tarball/0.1',
    author="Yosef Roth",
    author_email="yosefdroth@gmail.com",
    license="MIT",
    keywords=['kinetic data', 'biology', 'reactions'],
   # packages=find_packages(exclude=['tests', 'tests.*']),
    #install_requires=install_requires,
    #tests_require=tests_require,
    classifiers=[]
)