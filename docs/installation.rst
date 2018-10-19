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

Run ``datanator``
-----------------------------
The API for datanator can be run with a test and production server.

In order to run the test server, run the following command::

    python3 manage.py runserver

NOTE: You will need to have the correct configuration in the datanator/__init__.py
file. Configurations can be found in datanator/config.py include:

* LocalDevelopmentConfig - Local server for database
* CircleTestingConfig - CircleCI server for database
* BuildConfig - Docker Compose/UCONN HPC server for database
* ProductionConfig - AWS RDS server for database (PRIVATE) 

In order to run the production server, run the following command::

    gunicorn -w 4 -b localhost:5000 --timeout 120 manage:app
        
This command will create a gunicorn production server with 4 workers at the localhost:5000 address with a timeout of 2 min 




Contact `Saahith <mailto:saahith116@gmail.com>`_ for any questions regarding installation and running the server
