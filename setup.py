
from setuptools import setup, find_packages
#import "/home/yosef/Desktop/Organizing DataSearch Code"
import os

# parse dependencies and their links from requirements.txt files
install_requires = [line.rstrip() for line in open('requirements.txt')]
#tests_require = [line.rstrip() for line in open('tests/requirements.txt')]

# install package
setup(
    name="kineticDatanator2",
    version="1.0.0",
    description="finds most relevant kinetic data",
    url="https://github.com/KarrLab/Kinetic-Datanator",
    #download_url='https://github.com/KarrLab/nose2unitth/tarball/{}'.format(nose2unitth.__version__),
    author="Yosef Roth",
    author_email="yosefdroth@gmail.com",
    license="MIT",
    keywords='kinetic data biology reactions',
    packages=find_packages(exclude=['tests', 'tests.*']),
    install_requires=install_requires,
    #tests_require=tests_require,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
    ],
    entry_points={
        'console_scripts': [
            'kineticDatanator2 = __main__:main',
        ],
    },
)