Installation
============

Install dependencies
--------------------
* Open Babel
* Python
* Pip

Install Open Babel
^^^^^^^^^^^^^^^^^^

Install Openbabel from `http://openbabel.org/wiki/Category:Installation <http://openbabel.org/wiki/Category:Installation>`_::

    apt-get install \
        libcairo2-dev \
        libeigen3-dev \
        libxml2-dev \
        swig \
        zlib1g-dev

    cd /tmp
    wget https://sourceforge.net/projects/openbabel/files/openbabel/2.4.1/openbabel-2.4.1.tar.gz/download -O /tmp/openbabel-2.4.1.tar.gz
    tar -xvvf /tmp/openbabel-2.4.1.tar.gz
    cd openbabel-2.4.1
    mkdir build
    cd build
    cmake ..
    make
    make install
    ldconfig

Latest revision from GitHub
---------------------------
Run the following command to install the latest version from GitHub::

    pip install git+git://github.com/KarrLab/kinetic_datanator.git#egg=kinetic_datanator

Latest release From PyPI
---------------------------
Run the following command to install the latest release from PyPI::

    pip install kinetic_datanator
