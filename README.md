# Branch Decomposition-Independent Edit Distances for Merge Trees

This repository contains the source code and examples for the following research papers:

"Branch Decomposition-Independent Edit Distances for Merge Trees"  
Florian Wetzels, Heike Leitte and Christoph Garth.
Computer Graphics Forum, 2022.
<!---[Link to paper]()-->

"A Deformation-based Edit Distance for Merge Trees"  
Florian Wetzels and Christoph Garth.
TopoInVis 2022.
<!---[Link to paper]()-->

"Merge Tree Geodesics and Barycenters with Path Mappings"  
Florian Wetzels, Mathieu Pont, Julien Tierny and Christoph Garth.
IEEE Transactions on Visualization and Computer Graphics, 2024.
<!---[Link to paper]()-->

"Accelerating Computation of Stable Merge Tree Edit Distances using Parameterized Heuristics"  
Florian Wetzels, Heike Leitte and Christoph Garth  
IEEE Transactions on Visualization and Computer Graphics, 2025.  
<!---[Link to paper]()-->

The sections below contain detailed instructions on how to compile and run the code on a vanilla Ubuntu 22.04 to reproduce the specific figures from each paper.
For simple usage, you can clone this repository and run the installation script, the data script and the python code (after setting the environment variables correctly) to render the figure:

```bash
./install.sh
./get_data_vortexstreet.sh

PV_PREFIX=ttk-paraview/install
export PATH=$PATH:$PV_PREFIX/bin
export LD_LIBRARY_PATH=$PV_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PV_PREFIX/lib/python3.10/site-packages

TTK_PREFIX=ttk/install
export PV_PLUGIN_PATH=$TTK_PREFIX/bin/plugins/TopologyToolKit
export LD_LIBRARY_PATH=$TTK_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:$TTK_PREFIX/lib/python3.10/site-packages

python3 scripts/compute_matrices_vortexstreet.py
```

## Usage

The directory `ttk` contains a complete ttk source with the methods from all four papers.
The installation script `install.sh` installs all dependencies, downloads the paraview source code and compiles paraview as well as ttk.
The data script `scripts/get_data_vortexstreet.sh` downloads the [2D vortex street dataset](https://www.csc.kth.se/~weinkauf/notes/cylinder2d.html) and converts it to vtk files.
The python script `scripts/get_data_vortexstreet.sh` loads the vtk dataset, computes merge trees and the distance matrix with the path mapping distance for look-ahead values of 0-3.

This produces four pdf files containing the distance matrices from the paper teaser, as well as four pdf files containing the t-SNE embeddings from the same figure.
All outputs are located in the root directory of the repository.
For a look-ahead value k (0-3), the two pdf are named `dm_la<k>.pdf` and `embedding_la<k>.pdf`.

Other datasets can be downloaded through the script `scripts/get_data_synthetic.sh`.
Alternatively, they can be manually downloaded from the data publications listed in the respective section below.
For these datasets, distance matrices can be computed using the scripts `scripts/compute_matrices_four_clusters.py` and `scripts/compute_matrices_four_clusters.py`.
A merge tree barycenter can be computed using `scripts/compute_barycenters_four_clusters.py`.

## Installation Note

Tested on Ubuntu 22.04.3 LTS.

### Install the dependencies

```bash
apt-get install -y cmake-qt-gui libboost-system-dev libpython3.10-dev libxt-dev libxcursor-dev libopengl-dev
apt-get install -y qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools 
apt-get install -y python3-sklearn python3-seaborn
apt-get install -y libsqlite3-dev 
apt-get install -y gawk
apt-get install -y git p7zip-full wget
```

### Install Paraview

First, go in the root of this repository and run the following commands:
(replace the `4` in `make -j4` by the number of available cores on your system)

```bash
git clone https://github.com/topology-tool-kit/ttk-paraview.git
cd ttk-paraview
git checkout 5.10.1
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DCMAKE_INSTALL_PREFIX=../install ..
make -j4
make -j4 install
```

Some warnings are expected when using the `make` command, they should not cause any problems.

Stay in the build directory and set the environment variables:
(replace `3.10` in `python3.10` by your version of python)

```bash
PV_PREFIX=`pwd`/../install
export PATH=$PATH:$PV_PREFIX/bin
export LD_LIBRARY_PATH=$PV_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PV_PREFIX/lib/python3.10/site-packages
```

### Install TTK

Go in the `ttk-lookahead` directory then run the following commands:
(replace the `4` in `make -j4` by the number of available cores on your system)

```bash
mkdir build && cd build
paraviewPath=`pwd`/../../ttk-paraview/install/lib/cmake/paraview-5.10
cmake -DCMAKE_INSTALL_PREFIX=../install -DParaView_DIR=$paraviewPath ..
make -j4
make -j4 install
```

Stay in the build directory and set the environment variables:
(replace `3.10` in `python3.10` by your version of python)

```bash
TTK_PREFIX=`pwd`/../install
export PV_PLUGIN_PATH=$TTK_PREFIX/bin/plugins/TopologyToolKit
export LD_LIBRARY_PATH=$TTK_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:$TTK_PREFIX/lib/python3.10/site-packages
```

## Datasets
"Vertical Instability Example — Four Clusters"  
Florian Wetzels, Heike Leitte, Christoph Garth.  
Zenodo, 2025.  
DOI: [10.5281/zenodo.16756130](https://doi.org/10.5281/zenodo.16756130)

"Vertical Instability Example — Outlier"  
Florian Wetzels, Heike Leitte, Christoph Garth.  
Zenodo, 2025.  
DOI: [10.5281/zenodo.16755706](https://doi.org/10.5281/zenodo.16755706)
