# Wigner-Seitz-Py
Compact Wigner-Seitz script for LAMMPS dump files. Can handle several atomic species if a host lattice species is defined. 

## Getting Started

git clone https://github.com/mb4512/Wigner-Seitz-Py.git

### Prerequisites

This project was written and tested only with Python 3.9.1. Likely any Python 3+ version will work.

Required Python libraries:
```
numpy
scipy
sklearn
```

### Installing

It is highly recommended to run this on a [virtual environment](https://docs.python.org/3/tutorial/venv.html) for python. All required libraries are easily installed with 
```
pip3 install numpy
pip3 install scipy
pip3 install sklearn
```

## Running the code 

The script only supports serial mode. The LAMMPS dump file is expected to contain the columns `ITEM: ATOMS id type x y z`. The script was tested only for orthogonal boxes, but in principle also supports non-orthogonal boxes.

The script is run
```
python3 ws-analysis.py REFERENCE.DUMP DISTORTED.DUMP -host ID -export EXPORT.DUMP
```

where `REFERENCE.DUMP` is the path to a LAMMPS dump file of the pristine reference lattice used for Wigner-Seitz analysis, `DISTORTED.DUMP` is the dump file containing distorted or damaged atomic structure. The optional `-host ID` flag tells the script which atomic ID to use as the host lattice (default: first occurring species in `REFERENCE.DUMP`) and the optional `-export EXPORT.DUMP` flag tells the script to what path the Wigner-Seitz output coordinates are stored at (default: `ws-output.dump`)

With the sample files gives in the `test` subdirectory, the script can be executed as
```
python3 ws-analysis.py test/W.dump test/W-H.856856.dump -host 1 -export test/export.dump
```

The output file `test/export.dump` contains all atomic coordinates of the host lattice with an occupancy other than 1, as well as all atomic coordinates of the other atomic species (mapped onto the reference lattice) with an occupancy other than 0.

## To be implemented

Passing `LAMMPS` data files, more flexibility with importing dump files. Automatic determination of the reference structure based on simulated scattering intensity. Currently only 3D periodic boundaries are supported.

## Authors

* **Max Boleininger**, [UKAEA](http://www.ccfe.ac.uk/) 
