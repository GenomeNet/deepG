#!/bin/bash

# pip
sudo apt-get install python3-setuptools
sudo easy_install3 pip

# venv
sudo apt-get install python3-venv
pip install --user virtualenv

# hdf5
cd ..
wget "$HDF5_RELEASE_URL/hdf5-${HDF5_VERSION%.*}/hdf5-$HDF5_VERSION/src/hdf5-$HDF5_VERSION.tar.gz"
tar -xzf "hdf5-$HDF5_VERSION.tar.gz"
cd "hdf5-$HDF5_VERSION"
./configure --prefix=/usr/local
sudo make install
cd ..
