import re
import setuptools
import subprocess
import sys
try:
    result = subprocess.run(
        [sys.executable, "-m", "pip", "show", "pkg_utils"],
        check=True, capture_output=True)
    match = re.search(r'\nVersion: (.*?)\n', result.stdout.decode(), re.DOTALL)
    assert match and tuple(match.group(1).split('.')) >= ('0', '0', '5')
except (subprocess.CalledProcessError, AssertionError):
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "-U", "pkg_utils"],
        check=True)
import os
import pkg_utils

name = 'datanator'
dirname = os.path.dirname(__file__)

# get package metadata
md = pkg_utils.get_package_metadata(dirname, name)

# install package
setuptools.setup(
    name=name,
    version=md.version,
    description='Finds relevant kinetic data for biochemical models',
    long_description=md.long_description,

    url='https://github.com/KarrLab/' + name,
    download_url='https://github.com/KarrLab/' + name,
    license='MIT',

    author='Karr Lab',
    author_email='members@karrlab.org',

    keywords=['kinetic data', 'systems biology', 'computational biology', ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],

    packages=setuptools.find_packages(exclude=['tests', 'tests.*']),
    package_data={
        name: [
            'config/core.schema.cfg',
            'config/core.default.cfg',
            'data_source/*.txt',
            'data/*.txt',
            'data/*.xlsx',
        ],
    },
    entry_points={
        'console_scripts': [
            'datanator = datanator.__main__:main',
        ],
    },

    install_requires=md.install_requires,
    extras_require=md.extras_require,
    tests_require=md.tests_require,
    dependency_links=md.dependency_links,
)
