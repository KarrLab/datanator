Installation
============

*Install dependencies*

Install Openbabel from <http://openbabel.org/wiki/Category:Installation>::

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

*Install kinetic_datanator*::

    $ pip install kinetic_datanator
