from setuptools import setup, find_packages
import kinetic_datanator

# parse dependencies and their links from requirements.txt files
install_requires = [line.rstrip() for line in open('requirements.txt')]
tests_require = [line.rstrip() for line in open('tests/requirements.txt')]
dependency_links = []

# install package
setup(
    name='kinetic_datanator',
    version = kinetic_datanator.__version__,
    description='Finds relevant kinetic data for biochemical models',

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
            'data/*.txt',
            'data/*.xlsx',
        ],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'kinetic_datanator = kinetic_datanator.__main__:main',
        ],
    },

    install_requires=install_requires,
    tests_require=tests_require,
    dependency_links=dependency_links,
)
