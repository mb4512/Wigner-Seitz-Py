import numpy as np

class AffineTransform:
    '''Define the affine transformation between two sytems.'''
    
    def __init__(self, celldim1, celldim2):
        # initialise affine transformation matrix and vector
        
        # get offset vectors
        self.r1 = self.get_origin(celldim1)
        self.r2 = self.get_origin(celldim2)
        
        # get basis matrices
        self.c1mat = self.get_cell_vectors(celldim1)
        self.c2mat = self.get_cell_vectors(celldim2)

        # get inverse matrices for lattice vector transformation
        self.c1mati = np.linalg.inv(self.c1mat)
        self.c2mati = np.linalg.inv(self.c2mat)

        # build transformation matrices
        self.amatrix  = np.matmul(self.c2mat.T, np.linalg.inv(self.c1mat.T))
        self.amatrixi = np.matmul(self.c1mat.T, np.linalg.inv(self.c2mat.T))
    
    
    def get_cell_dimensions(self, celldim):
        xlo_bound, xhi_bound, xy = celldim[0]
        ylo_bound, yhi_bound, xz = celldim[1]
        zlo_bound, zhi_bound, yz = celldim[2]

        xlo = xlo_bound - min(0., xy, xz, xy+xz)
        xhi = xhi_bound - max(0., xy, xz, xy+xz)
        ylo = ylo_bound - min(0., yz)
        yhi = yhi_bound - max(0., yz)
        zlo = zlo_bound
        zhi = zhi_bound

        return (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz)

    
    def get_origin(self, celldim):
        (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz) = self.get_cell_dimensions(celldim)
        return np.r_[xlo, ylo, zlo]

    
    def get_cell_vectors(self, celldim):
        (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz) = self.get_cell_dimensions(celldim)

        c1 = np.r_[xhi-xlo, 0., 0.]
        c2 = np.r_[xy, yhi-ylo, 0.]
        c3 = np.r_[xz, yz, zhi-zlo]

        return np.c_[[c1,c2,c3]]
    
    def go1to2(self, xyz):
        '''Transform vector xyz from reference frame 1 to reference frame 2.
        
        Only for testing single transformations! It is slow.'''
        return  self.amatrix.dot(xyz - self.r1) + self.r2
    
    def go2to1(self, xyz):
        '''Transform vector xyz from reference frame 2 to reference frame 1.
        
        Only for testing single transformations! It is slow.'''
        return self.amatrixi.dot(xyz - self.r2) + self.r1
 
    def go1to2_all(self, xyz):
        '''Transform all vectors xyz from reference frame 1 to reference frame 2.
        
        This is reasonably efficient.'''
        return np.r_[[self.amatrix[0][0]*(xyz[:,0] - self.r1[0]) +
                      self.amatrix[0][1]*(xyz[:,1] - self.r1[0]) +
                      self.amatrix[0][2]*(xyz[:,2] - self.r1[0]) + self.r2[0]],
                     [self.amatrix[1][0]*(xyz[:,0] - self.r1[1]) +
                      self.amatrix[1][1]*(xyz[:,1] - self.r1[1]) +
                      self.amatrix[1][2]*(xyz[:,2] - self.r1[1]) + self.r2[1]],
                     [self.amatrix[2][0]*(xyz[:,0] - self.r1[2]) +
                      self.amatrix[2][1]*(xyz[:,1] - self.r1[2]) +
                      self.amatrix[2][2]*(xyz[:,2] - self.r1[2]) + self.r2[2]]].T
    
    def go2to1_all(self, xyz):
        '''Transform all vectors xyz from reference frame 2 to reference frame 1.
        
        This is reasonably efficient.'''
        return np.r_[[self.amatrixi[0][0]*(xyz[:,0] - self.r2[0]) +
                      self.amatrixi[0][1]*(xyz[:,1] - self.r2[0]) +
                      self.amatrixi[0][2]*(xyz[:,2] - self.r2[0]) + self.r1[0]],
                     [self.amatrixi[1][0]*(xyz[:,0] - self.r2[1]) +
                      self.amatrixi[1][1]*(xyz[:,1] - self.r2[1]) +
                      self.amatrixi[1][2]*(xyz[:,2] - self.r2[1]) + self.r1[1]],
                     [self.amatrixi[2][0]*(xyz[:,0] - self.r2[2]) +
                      self.amatrixi[2][1]*(xyz[:,1] - self.r2[2]) +
                      self.amatrixi[2][2]*(xyz[:,2] - self.r2[2]) + self.r1[2]]].T

    def pbcwrap1(self, xyz):
        '''Check if a vector falls outside the box, and if so, wrap it back inside.'''

        fcoords = np.matmul(self.c1mati, xyz - self.r1)
        gcoords = fcoords - np.floor(fcoords)

        return self.r1 + np.matmul(self.c1mat, gcoords)

    def pbcdistance(self, xyz1, xyz2):
        '''Compute the distance between two vectors, taking into account pbc.
        Only for testing single transformations! It is slow.'''

        fcoords1 = np.matmul(self.c1mati, xyz1 - self.r1)
        fcoords2 = np.matmul(self.c1mati, xyz2 - self.r1)
        
        df = fcoords2 - fcoords1
        dg = df - np.sign(df)*(np.abs(df) > .5).astype(int)

        return np.linalg.norm(np.matmul(self.c1mat, dg))


    def pbcdistance_all(self, xyz, xyz_all):
        '''Compute the distance between one vectors and a list of vectors, taking into account pbc.'''

        fcoords     = np.matmul(self.c1mati, xyz - self.r1)
        
        xyzr_all = xyz_all - self.r1
        fcoords_all = np.r_[[self.c1mati[0][0]*(xyzr_all[:,0]) +
                             self.c1mati[0][1]*(xyzr_all[:,1]) +
                             self.c1mati[0][2]*(xyzr_all[:,2])], 
                            [self.c1mati[1][0]*(xyzr_all[:,0]) +
                             self.c1mati[1][1]*(xyzr_all[:,1]) +
                             self.c1mati[1][2]*(xyzr_all[:,2])],
                            [self.c1mati[2][0]*(xyzr_all[:,0]) +
                             self.c1mati[2][1]*(xyzr_all[:,1]) +
                             self.c1mati[2][2]*(xyzr_all[:,2])]].T
         
        df = fcoords_all - fcoords
        dg = df - np.sign(df)*(np.abs(df) > .5).astype(int)

        dxyz_all = np.r_[[self.c1mat[0][0]*dg[:,0] +
                          self.c1mat[0][1]*dg[:,1] +
                          self.c1mat[0][2]*dg[:,2]], 
                         [self.c1mat[1][0]*dg[:,0] +
                          self.c1mat[1][1]*dg[:,1] +
                          self.c1mat[1][2]*dg[:,2]],
                         [self.c1mat[2][0]*dg[:,0] +
                          self.c1mat[2][1]*dg[:,1] +
                          self.c1mat[2][2]*dg[:,2]]]
        
        dxyz_norm = np.sqrt(dxyz_all[0]*dxyz_all[0] + dxyz_all[1]*dxyz_all[1] + dxyz_all[2]*dxyz_all[2])
        return dxyz_norm 

