# initialise
import os, sys, glob
import numpy as np

from scipy.spatial import cKDTree 

from lib.affine import AffineTransform
from lib.dumpio import ReadDump
from lib.dumpio import write_dump

def announce(string):
    print ("\n=====================================================")
    print (string)
    print ("=====================================================\n")
    return 0

def main():
    announce(" Wigner-Seitz analysis for non-monoatomic structures\n Max Boleininger (2021), mboleininger@gmail.com")

    print ("Command line arguments:", sys.argv)
    print ()

    path0 = sys.argv[1] # path to reference structure
    path  = sys.argv[2] # asd

    host = None
    epath = None
    wsmode = None
    for _i,_arg in enumerate(sys.argv):

        # atom type id of host lattice
        if _arg == "-host":
            host = int(sys.argv[_i+1])

        # WS output file export path 
        if _arg == "-export":
            epath = sys.argv[_i+1]

        # WS output file export path 
        if _arg == "-wsmode":
            wsmode = sys.argv[_i+1]


    # import reference cell
    print ("Importing reference cell data %s... " % path0, end='')
    data0 = ReadDump(path0) 
    print ("imported %d atoms." % data0.natoms)

    ref_species = np.array(np.unique(data0.data[:,0]), dtype=np.int64)
    print ("Reference file contains the atomic species:", ref_species)
    
    if not host:
        host = ref_species[0]
        print ("Defaulting to using atom type %d as host lattice." % host)
    print ()

    if not epath:
        epath = "ws-output.dump"
        print ("Defaulting to using %s as output path." % epath)
    print ()

    if not wsmode:
        wsmode = "default"
        print ("Defaulting to using %s Wigner-Seitz mode for implanted species." % wsmode)
    print ()



    # import distorted cell
    print ("Importing distorted cell data %s... " % path, end='')
    datap = ReadDump(path)
    print ("imported %d atoms." % datap.natoms)

    species = np.array(np.unique(datap.data[:,0]), dtype=np.int64)
    print ("Distorted file contains the atomic species:", species)
    print ()

    print ("Affine mapping of distorted cell to reference cell... ", end='')
    # set up affine transformations
    atrans = AffineTransform(data0.celldim, datap.celldim)

    # transform distorted coordinates to reference cell
    datap.data[:,1:] = atrans.go2to1_all(datap.data[:,1:])

    # wrap atoms back into reference cell if outside
    # todo: implement fast version of pbcwrap1 
    datap.data[:,1:] = np.r_[[atrans.pbcwrap1(_x) for _x in datap.data[:,1:]]]
    print ("done.")

    # get cell vectors of reference cell and create 1st neb. periodic images
    print ("Creating periodic images of reference cell... ", end='')
    cvecs = atrans.get_cell_vectors(data0.celldim)

    # define unique atomic ids of periodically repeated atoms
    repeated_ids = np.r_[[i for i in range(data0.natoms)]]
    dcopy = np.copy(data0.data)
    
    _ids = np.copy(repeated_ids)  
    for ix in [-1,0,1]:
        for iy in [-1,0,1]:
            for iz in [-1,0,1]:
                if ix == iy == iz == 0:
                    continue
                _dcopy = np.copy(data0.data)
                _dcopy[:,1:] += ix*cvecs[0] + iy*cvecs[1] + iz*cvecs[2]
                dcopy = np.r_[dcopy, _dcopy]
                repeated_ids = np.r_[repeated_ids, _ids]
    print ("done.")

    # build KDTree of periodically repeated reference structure for nearest neighbour search
    print ("Building k-tree of reference cell atoms... ", end='')
    ktree0 = cKDTree(dcopy[:,1:])
    print ("done.")

    announce ("Wigner-Seitz analysis")

    # first, find defect content of host lattice type
    amax = 8.0
    if host:
        _hostxyz = datap.data[datap.data[:,0] == host, 1:]
    else:
        _hostxyz = datap.data[:,1:]

    neb_indices = ktree0.query(_hostxyz, k=1, distance_upper_bound=amax)[1] 
    neb_indices = repeated_ids[neb_indices]

    # find number of vacant sites
    uvals, ucounts = np.unique(neb_indices, return_counts=True)

    vac_indices = np.setdiff1d(repeated_ids, uvals)
    nvac = len(vac_indices)

    # find number of interstitial sites 
    int_indices = uvals[ucounts > 1]
    nint = np.sum(ucounts[ucounts>1]-1)

    print("Host species %d,      n_vac n_int = %d %d" % (host, nvac, nint))

    # arrays storing W-S analysis results for exporting
    xyz_array   = np.r_[data0.data[int_indices,1:], data0.data[vac_indices,1:]]
    types_array = np.r_[[host]*(len(int_indices)+len(vac_indices))]
    occ_array   = np.r_[ucounts[ucounts>1], [0]*len(vac_indices)]

    # next, find defect content of other atomic species
    if host and (wsmode == "default"):
        for _sp in species:
            if _sp == host:
                continue

            _restxyz = datap.data[datap.data[:,0] == _sp, 1:]
            neb_indices = ktree0.query(_restxyz, k=1, distance_upper_bound=amax)[1]
            neb_indices = repeated_ids[neb_indices]

            # find number of vacant sites
            uvals, ucounts = np.unique(neb_indices, return_counts=True)

            vac_indices = np.setdiff1d(repeated_ids, uvals)
            nvac = len(vac_indices)

            # since the base occupancy of the implanted species is 0, only look at occupied sites
            occ_indices = uvals[ucounts > 0] 
            nocc = np.sum(ucounts[ucounts>0])

            print("Implanted species %d, n_vac n_occ = %d %d" % (_sp, nvac, nocc))

            # add W-S analysis of implanted species for exporting
            xyz_array   = np.r_[xyz_array, data0.data[occ_indices,1:]]
            types_array = np.r_[types_array, [_sp]*len(occ_indices)]
            occ_array   = np.r_[occ_array, ucounts[ucounts>0]]

    if host and (wsmode == "host"):
        for _sp in species:
            if _sp == host:
                continue

            _restxyz = datap.data[datap.data[:,0] == _sp, 1:] 

            print("Implanted species %d, n_atoms = %d" % (_sp, len(_restxyz)))
            
            # add implanted species for exporting
            xyz_array   = np.r_[xyz_array, _restxyz] 
            types_array = np.r_[types_array, [_sp]*len(_restxyz)]
            occ_array   = np.r_[occ_array,  [1]*len(_restxyz)]

    types_array = np.array(types_array, dtype=np.int64)
    occ_array = np.array(occ_array, dtype=np.int64)

    print ()
    print ("Exporting W-S output to file %s... " % epath, end="")
    write_dump(data0.celldim, types_array, occ_array, xyz_array, epath, 0)
    print ("done.")
    print ()


if __name__=="__main__":
    main()

