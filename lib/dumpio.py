import numpy as np

class ReadDump:
    '''Class for importing, processing, and storing dump LAMMPS files'''
    def __init__(self, path, collist=None,  nlines=None):
        '''Import dump file and save content as attributes.'''
        reffile = []
        # possibly faster: enter into dictionary as file is read line-by-line
        with open(path) as file:
            for i, line in enumerate(file):
                reffile.append(line.rstrip())
                if nlines:
                    if i > nlines:
                        break

        # enter input file data
        ref = {}
        self.step = int(reffile[1])
        self.natoms = int(reffile[3]) 
        self.celldim = np.r_[[_.split() for _ in reffile[5:8]]].astype(np.float)
        self.data = np.r_[[_.split() for _ in reffile[9:]]].astype(np.float)

        # if orthorhombic cell is imported, add zero-valued xy, xz, yz values 
        if self.celldim.shape == (3,2):
            self.celldim = np.c_[self.celldim, [0.,0.,0.]]

        self.ids = self.data[:,0].astype(int)

        # sort data according to atomic id
        iorder = np.argsort(self.ids)
        self.ids = self.ids[iorder]
        self.data = self.data[iorder]
        
        if collist:
            self.data = self.data[:,collist]
        else:
            self.data = self.data[:,1:]


def write_dump(cell, types, occ, xyz, path, frame):

    # exctract cell data
    N = len(xyz) 
    xlo_bound, xhi_bound, xy = cell[0]
    ylo_bound, yhi_bound, xz = cell[1]
    zlo_bound, zhi_bound, yz = cell[2]

    # write file header
    wfile = open(path, 'w')

    wfile.write("ITEM: TIMESTEP\n")
    wfile.write("%d\n" % frame)
    wfile.write("ITEM: NUMBER OF ATOMS\n")
    wfile.write("%d\n" % N)

    # simulation box info specific for orthogonal box with PBC
    wfile.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
    wfile.write("%f %f %f\n" % (xlo_bound, xhi_bound, xy))
    wfile.write("%f %f %f\n" % (ylo_bound, yhi_bound, xz))
    wfile.write("%f %f %f\n" % (zlo_bound, zhi_bound, yz))

    wfile.write("ITEM: ATOMS id type x y z occ\n")

    for _i,_xyz in enumerate(xyz):
        wfile.write("%5d %3d %14.8f %14.8f %14.8f %3d\n" % (_i+1, types[_i], _xyz[0], _xyz[1], _xyz[2], occ[_i])) 

    wfile.close()

    return 0

