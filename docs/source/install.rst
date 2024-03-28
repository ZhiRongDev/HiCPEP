Installation
=============

All the programs were tested in Ubuntu 22.04.4 LTS, HiCPAP requires ``python3``, ``pip`` and ``libcurl4-openssl-dev`` installed on your system. 

For example (Paste these commands in Bash or Zsh):

.. code:: bash

    sudo apt-get update
    sudo apt-get install -y libcurl4-openssl-dev
    sudo apt-get install -y python3
    sudo apt-get install -y pip
    sudo apt-get install -y git 
    git clone git@github.com:ZhiRongDev/HiCPAP.git
    cd HiCPAP
    python3 -m pip install -e .

If you have already installed the requirements, just paste these commands:

.. code:: bash

    git clone git@github.com:ZhiRongDev/HiCPAP.git
    cd HiCPAP
    python3 -m pip install -e .