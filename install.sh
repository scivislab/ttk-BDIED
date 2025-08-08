#!/bin/bash

apt-get install -y cmake-qt-gui libboost-system-dev libpython3.10-dev libxt-dev libxcursor-dev libopengl-dev
apt-get install -y qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools 
apt-get install -y python3-sklearn python3-seaborn
apt-get install -y libsqlite3-dev 
apt-get install -y gawk
apt-get install -y git p7zip-full wget

git clone https://github.com/topology-tool-kit/ttk-paraview.git
cd ttk-paraview
git checkout 5.10.1
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DCMAKE_INSTALL_PREFIX=../install ..
make -j4
make -j4 install

PV_PREFIX=`pwd`/../install
export PATH=$PATH:$PV_PREFIX/bin
export LD_LIBRARY_PATH=$PV_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PV_PREFIX/lib/python3.10/site-packages

cd ../../ttk

mkdir build && cd build
paraviewPath=`pwd`/../../ttk-paraview/install/lib/cmake/paraview-5.10
cmake -DCMAKE_INSTALL_PREFIX=../install -DParaView_DIR=$paraviewPath ..
make -j4
make -j4 install

TTK_PREFIX=`pwd`/../install
export PV_PLUGIN_PATH=$TTK_PREFIX/bin/plugins/TopologyToolKit
export LD_LIBRARY_PATH=$TTK_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:$TTK_PREFIX/lib/python3.10/site-packages

cd ../..
