Installation
============
The following instructions describe how to install ``datanator`` onto a Debian-based Linux OS
using the docker image `wc_env <https://hub.docker.com/r/karrlab/wc_env>`

Install dependencies
--------------------
First, please install the following dependencies:

* `Docker <https://docs.docker.com/get-docker/>`_
* `Docker-compose <https://docs.docker.com/compose/install/>`_

The following shell commands can be used to install these dependencies onto Ubuntu Linux 16.04::

    apt-get update
    
    apt-get install \
        apt-transport-https \
        ca-certificates \
        curl \
        gnupg-agent \
        software-properties-common

    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

    sudo add-apt-repository \
        "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
        $(lsb_release -cs) \
        stable"

    sudo apt-get update

    sudo apt-get install docker-ce docker-ce-cli containerd.io

    sudo curl -L "https://github.com/docker/compose/releases/download/1.26.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose

    sudo chmod +x /usr/local/bin/docker-compose

    sudo ln -s /usr/local/bin/docker-compose /usr/bin/docker-compose

Install ``datanator``
-----------------------------
Second, please run the following shell commands to clone and install ``datanator`` from GitHub::

    mkdir karr_lab
    mkdir ~/.wc
    cd ./karr_lab
    git clone git@github.com:KarrLab/pkg_utils.git
    git clone git@github.com:KarrLab/wc_utils.git
    git clone git@github.com:KarrLab/karr_lab_aws_manager.git
    git clone git@github.com:KarrLab/datanator_query_python.git
    git clone git@github.com:KarrLab/datanator.git
    cd ./datanator
    nano docker-compose.yml # change ``zl`` on line 13 to the proper username. Save and exit by pressing ``Ctrl + X`` followed by ``Y``
    docker-compose up -d


Run ``datanator``
-----------------------------
One needs to find the docker container ID in order to use Datanator package::

    docker ps
    docker exec -it <container_id> bash
    cd karr_lab/datanator

All python scripts in ``datanator`` dicrectory can be run with python3, for example::

    python3 datanator/data_source/corum_nosql.py

Running the command above will parse `Corum <http://mips.helmholtz-muenchen.de/corum/>` and
store the parsed data in KarrLab's MongoDB.


Contact `Yang <mailto:zhouyang.lian@familian.life>`_ for any questions regarding installation and running the package.
