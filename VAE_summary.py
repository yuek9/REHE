import shelve
import numpy as np
import VAE_function
import pandas as pd

n=300
p=100
theta = 8.0
option = 'logisticnormal'
option = 'dirichlet'
z_scale= 15


import glob, os
res_file = []
os.chdir('/home/students/yuek/Desktop/VAE/VAEres')
for file in glob.glob("*"):
    res_file.append(file)

res_out = pd.DataFrame([], columns = ['method','n','p','distribution','theta','z_scale','time','frob','l1','shannon', 'simpson', 'KL', 'mean_sparsity'])



for filename in res_file:
    # filename = 'shelve'+str(n)+'_'+str(p)+'_'+str(theta)+'_'+option+'_'+str(z_scale)+'.out'
    print(filename)
    my_shelf = shelve.open(filename, 'r')
    #for key in my_shelf:
    #    globals()[key]=my_shelf[key]
    true_comp = my_shelf['true_comp']
    true_count = my_shelf['true_count']
    VAE_comp = my_shelf['VAE_comp']
    naive_comp = my_shelf['naive_comp']
    nreps = my_shelf['nreps']
    VAE_get_l1 = []
    VAE_get_frob = []
    VAE_get_shannon = []
    VAE_get_simpson = []
    VAE_get_KL = []
    sparsity = []
    VAE_time = my_shelf['VAE_time']
    naive_get_l1 = []
    naive_get_frob = []
    naive_get_shannon = []
    naive_get_simpson = []
    naive_get_KL = []
    np.mean(VAE_time)
    for i in range(nreps):
        sparsity.append(np.mean(true_count[i][0]==0))
        naive_alr = VAE_function.alr_trans(naive_comp[i][0])
        VAE_alr = VAE_function.alr_trans(VAE_comp[i][0])
        VAE_get_frob.append(VAE_function.eval_frob(VAE_comp[i][0], true_comp[i][0]))
        VAE_get_l1.append(VAE_function.eval_l1(VAE_comp[i][0], true_comp[i][0]))
        VAE_get_shannon.append(VAE_function.eval_shannon(VAE_comp[i][0], true_comp[i][0]))
        VAE_get_simpson.append(VAE_function.eval_simpson(VAE_comp[i][0], true_comp[i][0]))
        VAE_get_KL.append(VAE_function.eval_KL(VAE_comp[i][0], true_comp[i][0]))
        naive_get_frob.append(VAE_function.eval_frob(naive_comp[i][0], true_comp[i][0]))
        naive_get_l1.append(VAE_function.eval_l1(naive_comp[i][0], true_comp[i][0]))
        naive_get_shannon.append(VAE_function.eval_shannon(naive_comp[i][0], true_comp[i][0]))
        naive_get_simpson.append(VAE_function.eval_simpson(naive_comp[i][0], true_comp[i][0]))
        naive_get_KL.append(VAE_function.eval_KL(naive_comp[i][0], true_comp[i][0]))

    if 'z_scale' in my_shelf.keys():
        z_scale = my_shelf['z_scale']
    else:
        z_scale=1.2

    res_out = pd.concat([res_out, pd.DataFrame([['VAE', my_shelf['n'], my_shelf['p'], my_shelf['option'], my_shelf['theta'],z_scale ,np.mean(VAE_time),np.mean(VAE_get_frob), np.mean(VAE_get_l1), np.mean(VAE_get_shannon), np.mean(VAE_get_simpson), np.mean(VAE_get_KL), np.mean(sparsity)],['naive', my_shelf['n'], my_shelf['p'], my_shelf['option'], my_shelf['theta'], z_scale, 0,np.mean(naive_get_frob), np.mean(naive_get_l1), np.mean(naive_get_shannon), np.mean(naive_get_simpson), np.mean(naive_get_KL), np.mean(sparsity)] ], columns = ['method','n','p','distribution','theta','z_scale','time','frob','l1','shannon', 'simpson', 'KL', 'mean_sparsity'])], axis=0)

    my_shelf.close()

os.chdir('/home/students/yuek/Desktop/VAE')
filename='VAE_summary.out'
my_shelf = shelve.open(filename,'n') # 'n' for new
for key in dir():
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()

np.savetxt(r'VAE_res_out.txt', res_out.values, fmt='%s')