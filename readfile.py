import numpy as np

def readfile(pre, n, p, library_scale, theta, option, z_scale, end):
    filename = pre+str(int(n))+'_'+str(int(p))+'_'+str(float(library_scale))+'_'+str(float(theta))+'_'+option+'_'+str(float(z_scale))+end
    print(filename)
    return np.load(filename)