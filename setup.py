from setuptools import setup, find_packages
import os
import pip
import re

# get long description
if os.path.isfile('README.rst'):
    with open('README.rst', 'r') as file:
        long_description = file.read()
else:
    long_description = ''

# get version
with open('kinetic_datanator/VERSION', 'r') as file:
    version = file.read().strip()

# parse dependencies and their links from requirements.txt files
install_requires = []
tests_require = []
dependency_links = []

for line in open('requirements.txt'):
    pkg_src = line.rstrip()
    match = re.match('^.+#egg=(.*?)$', pkg_src)
    if match:
        pkg_id = match.group(1)
        dependency_links.append(pkg_src)
    else:
        pkg_id = pkg_src
    install_requires.append(pkg_id)
    pip.main(['install', line])

for line in open('tests/requirements.txt'):
    pkg_src = line.rstrip()
    match = re.match('^.+#egg=(.*?)$', pkg_src)
    if match:
        pkg_id = match.group(1)
        dependency_links.append(pkg_src)
    else:
        pkg_id = pkg_src
    tests_require.append(pkg_id)

dependency_links = list(set(dependency_links))

# install package
setup(
    name='kinetic_datanator',
    version=version,
    description='Finds relevant kinetic data for biochemical models',
    long_description=long_description,

    url='https://github.com/KarrLab/kinetic_datanator',
    download_url='https://github.com/KarrLab/kinetic_datanator',
    license='MIT',

    author='Yosef Roth',
    author_email='yosefdroth@gmail.com',

    keywords=['kinetic data', 'systems biology', 'computational biology', ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],

    packages=find_packages(exclude=['tests', 'tests.*']),
    package_data={
        'kinetic_datanator': [
            'VERSION',
            'config/schema.cfg',
            'config/default.cfg',
            'data_source/*.txt',
            'data/*.txt',
            'data/*.xlsx',
        ],
    },
    entry_points={
        'console_scripts': [
            'kinetic_datanator = kinetic_datanator.__main__:main',
        ],
    },

    install_requires=install_requires,
    tests_require=tests_require,
    dependency_links=dependency_links,
)

# download NCBI taxonomy database
from kinetic_datanator.util import taxonomy_util
taxonomy_util.setup_database()
