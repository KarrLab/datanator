
from setuptools import setup, find_packages
#import "/home/yosef/Desktop/Organizing DataSearch Code"
import os

# parse dependencies and their links from requirements.txt files
#install_requires = [line.rstrip() for line in open('requirements.txt')]
#tests_require = [line.rstrip() for line in open('tests/requirements.txt')]

# install package

setup(
    name='KineticDatanator',
    packages = ['KineticDatanator'],
    version = '0.1.18',
    description="finds most relevant kinetic data",
    url="https://github.com/KarrLab/Kinetic-Datanator",
    download_url='https://github.com/KarrLab/Kinetic-Datanator/tarball/0.1',
    author="Yosef Roth",
    author_email="yosefdroth@gmail.com",
    license="MIT",
    keywords=['kinetic data', 'biology', 'reactions'],
    package_data = {},
    include_package_data=True,
   # packages=find_packages(exclude=['tests', 'tests.*']),
    install_requires=['requests==2.2.1', 'urllib3==1.7.1', 'openpyxl==2.3.5', 'ete3==3.0.0b35', 'argparse', 'jxmlease==1.0.1', 'numpy==1.8.2', 'cement==2.10.2'],
    #tests_require=tests_require,
    classifiers=[]
)
