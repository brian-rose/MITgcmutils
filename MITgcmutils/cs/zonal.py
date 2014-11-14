#  alternative to zonalAvg.py:
#  construct an array of weights of size (6*nc,nc,nr,ny)
#   for each latitude in the zonal average, 
#  and for each level,
#    the array is (6*nc,nc)
#   and just needs to multiplied by the field to be averaged, and summed.
#  i.e. we construct a kernel for the zonal average
#  only need to do this once for each grid, then zonal averages are trivial
#  for any field on that grid.

import numpy as np

def zonalKernel(csgrid,ny):
    '''Construct a kernel for zonal averaging.'''
    ar = csgrid['rA']  # area of grid cells
    hc = csgrid['HFacC']
    yc = csgrid['YC']
    
    nc = np.size(yc,axis=0) # size of cube face
    nr = np.size(hc,axis=2)  # number of vertical levels
    
    # Define latitude axis (with regular spacing) to average to.
    dy = 180./ny;   ylat = (np.arange(ny)+0.5)*dy - 90.
    
    # Compute zonal line weights.
    nMax = 6*nc*3           # Max number of cells near zonal line
    #ijzav = np.zeros((ny,nMax), dtype=int)  # Index of cells near band line (within dy)
    izav = np.zeros((ny,nMax), dtype=int)  # first index of cells near band line (within dy)
    jzav = np.zeros((ny,nMax), dtype=int)  # second index of cells near band line (within dy)
    alzav = np.zeros((ny,nMax))  # Area weights for cells near band line
    npzav = np.zeros(ny)     # Number of cells close to band line
    #for ij = 1:6*nc*nc                # Loop over all cells
    #for ij in range(npoints):
    for i in range(6*nc):
        for j in range(nc):
            delta = (yc[i,j]-ylat[0])/dy  # Latitude difference from bottom
            jz = 1 + np.floor(delta)              # Number of zonal lines away (round up)
            delta = delta % 1               # Fracional distance from upper line
            if jz < 1:
                n = npzav[jz] + 1   # Contribution from cell NORTH of line
                #ijzav[jz,n-1] = ij
                izav[jz,n-1] = i
                jzav[jz,n-1] = j
                alzav[jz,n-1] = 1     # Weight = 1 because all area accounted for
                npzav[jz] = n 
            elif jz >= ny:
                n = npzav[jz-1] + 1     # Contribution from cell NORTH of line
                #ijzav[jz-1,n-1] = ij
                izav[jz,n-1] = i
                jzav[jz,n-1] = j
                alzav[jz-1,n-1] = 1       # Weight = 1 because all area accounted for
                npzav[jz-1] = n
            else:
                n = npzav[jz-1] + 1;     # Contribution from cell SOUTH of line
                #ijzav[jz-1,n-1] = ij
                izav[jz-1,n-1] = i
                jzav[jz-1,n-1] = j
                alzav[jz-1,n-1] = 1 - delta
                npzav[jz-1] = n
                n = npzav[jz] + 1   # Contribution from cell NORTH of line
                #ijzav[jz,n-1] = ij
                izav[jz,n-1] = i
                jzav[jz,n-1] = j
                alzav[jz,n-1] = delta
                npzav[jz] = n
    NbMx = np.max(npzav)  # Reduce size
    #ijzav = ijzav[:,:NbMx]
    izav = izav[:,:NbMx]
    jzav = jzav[:,:NbMx]
    alzav = alzav[:,:NbMx]
    
    # Area associated with zonal lines:
    areazon = np.zeros(ny)
    for j in range(ny):
        #ijLoc = ijzav[j,:npzav[j]]
        iLoc = izav[j,:npzav[j]]
        jLoc = jzav[j,:npzav[j]]
        vvLoc = alzav[j,:npzav[j]]  #';
        areazon[j] = np.sum(vvLoc * ar[iLoc,jLoc])
    
    # in matlab code, mskzon was calculate in the basin block
    #  it's an area mask per basin and per level. with nBas = 0 it's identical to areazon
    #  but only for a single flat bottomed basin!
    #  will need to do more work here if we want to take proper basin averages
    #mskLoc = areazon
    
    kernel = np.zeros(6*nc,nc,ny,nr)
    for k in range(nr):
        var = ar * hc[:,:,k]
        for j in range(ny):
            iLoc = izav[j,:npzav[j]]
            jLoc = jzav[j,:npzav[j]]
            vvLoc = alzav[j,:npzav[j]]
            kernel[iLoc,jLoc,j,k] = vvLoc * var[iLoc,jLoc] / areazon[j]
    
    ## Compute zonal average:
    #for it in range(nt):
    #    fld1t = fld2[:,:,it]
    #    #for nb = 1:1+abs(nBas),
    #        #if nb == 1, vv1 = ar;
    #        #else        vv1 = ar.*mskB(:,nb-1); end
    #    for k in range(nr):
    #            #if nBas < 0, mskLoc=mskzon(:,nb);
    #            #else         mskLoc=mskzon(:,k,nb); end
    #        var = ar2 * hc2[:,k]
    #        var = var * fld1t[:,k]
    #        for j in range(ny):
    #            if mskLoc[j] > frcZmn*areazon[j]:
    #                ijLoc=ijzav[j,:npzav[j]]
    #                vvLoc=alzav[j,:npzav[j]]
    #                fldzon[j,k,it] = np.sum(vvLoc*var[ijLoc])/mskLoc[j]
    #            else:
    #                fldzon[j,k,it] = np.nan
    
    return kernel, ylat, areazon