import setuptools
try:
    import pkg_utils
except ImportError:
    import pip._internal
    pip._internal.main(['install', '--process-dependency-links', 'git+https://github.com/KarrLab/pkg_utils.git#egg=pkg_utils'])
    import pkg_utils
import os

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
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],

    packages=setuptools.find_packages(exclude=['tests', 'tests.*']),
    package_data={
        name: [
            'VERSION',
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
