Installation
============
The following instructions describe how to install ``datanator`` onto Ubuntu Linux 16.04.
Datanator only supports Python 3.

Install dependencies
--------------------
First, please install the following dependencies:

* `Git <https://git-scm.com>`_
* `Open Babel <http://openbabel.org>`_
* `Python 3 <https://www.python.org>`_
* `Pip <https://pip.pypa.io>`_

The following shell commands can be used to install these dependencies onto Ubuntu Linux 16.04::

    apt-get install \
        git \
        libcairo2-dev \
        libeigen3-dev \
        libxml2-dev \
        python \
        python-pip \
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


Install ``datanator``
-----------------------------
Second, please run the following shell commands to clone and install ``datanator`` from GitHub::

    git clone git@github.com:KarrLab/datanator.git
    pip3 install -e datanator

Because ``datanator`` is under active development, we recommend regularly pulling the latest revision of ``datanator`` from GitHub.
