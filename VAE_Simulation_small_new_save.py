from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from keras.layers import Lambda, Input, Dense, advanced_activations, BatchNormalization, Dropout
from keras.models import Model
from keras.utils import plot_model
from keras import backend as K
from keras import optimizers as op
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
import time
from keras.callbacks import LearningRateScheduler
#from rpy2.robjects.packages import importr
#import rpy2.robjects as ro
#import rpy2.robjects.numpy2ri
import sys
import VAE_function

w_sampling = VAE_function.w_sampling
x_sampling = VAE_function.x_sampling
normalization_transform = VAE_function.normalization_transform
logistic_transform = VAE_function.logistic_transform
add_epsilon = VAE_function.add_epsilon
lr_scheduler = VAE_function.lr_scheduler
activation_wrapper = VAE_function.activation_wrapper
generate_zi_NN = VAE_function.generate_zi_NN
generate_xi = VAE_function.generate_xi
generate_yi = VAE_function.generate_yi
eval_mse = VAE_function.eval_mse
eval_frob = VAE_function.eval_frob
eval_l1 = VAE_function.eval_l1
y_z_loss = VAE_function.y_z_loss
data_generation = VAE_function.data_generation
naive_compute = VAE_function.naive_compute
alr_trans = VAE_function.alr_trans
generate_zi_Cao = VAE_function.generate_zi_Cao
comp_shannon = VAE_function.comp_shannon
comp_simpson = VAE_function.comp_simpson 
eval_KL = VAE_function.eval_KL
eval_simpson = VAE_function.eval_simpson
eval_shannon = VAE_function.eval_shannon 

#--------------
# build a simple neural network
#--------------

def build_network_and_train(x_train, x_test, intermediate_dim = 32):
    original_dim = x_train.shape[1] # input dimension
    input_shape = (original_dim, )
    intermediate_dim = int(intermediate_dim)
    # intermediate_dim = 32 # neural network layer
    batch_size = 32
    latent_dim = 5
    epochs = 1000

    ## VAE model = encoder + decoder; first build encoder model
    inputs = Input(shape=input_shape, name='encoder_input')

    # in SAVER-X, inputs are normalized with this extra step.
    # inputs_2 = Lambda(normalization_transform, output_shape=input_shape, name='normalization')(inputs)

    inputs_drop = Dropout(0.2)(inputs)

    ## only have one main layer: counts data -> latent layer (of dimension 10) -> latent mean and variance for wi
    ## it is possible that activation function is important in convergence

    wrapped_relu = activation_wrapper(tf.keras.activations.relu)

    latent_layer = Dense(latent_dim, activation='relu')(inputs_drop)
    latent_layer = BatchNormalization(center=True, scale=False)(latent_layer)

    wi_mean = Dense(latent_dim, activation=wrapped_relu(alpha=0.1, threshold=-1, max_value=5))(latent_layer)
    wi_log_var = Dense(latent_dim, activation=wrapped_relu(alpha=0.1, threshold=0, max_value=5))(latent_layer)

    ## use reparameterization trick to get wi
    wi_layer = Lambda(w_sampling, output_shape=(latent_dim,), name='wi')([wi_mean, wi_log_var])

    ## instantiate encoder model
    encoder = Model(inputs, [wi_mean, wi_log_var, wi_layer], name='encoder')
    encoder.summary() # encoder has 2570 parameters
    # plot_model(encoder, to_file='vae_mlp_encoder.png', show_shapes=True)

    ## VAE model = encoder + decoder; next build decoder model
    ## latent wi -> one intermediate layer -> output
    latent_inputs = Input(shape=(latent_dim,), name='wi_sampling')
    latent_inputs_drop = Dropout(0.2)(latent_inputs)
    layer1 = Dense(intermediate_dim, activation='relu')(latent_inputs_drop)
    layer1 = BatchNormalization(center=True, scale=False)(layer1)
    layer1 = Dropout(0.2)(layer1)
    ## option 1: use z_i as the neural net output, assume xi is from Dirichlet(z_i), then loss function P(Y|Z) is closed form, and xi of interest will be posterior mode of x|y
    # output a posterior sample of zi based on input yi; computed with one layer of neural network
    # notice that we must have positive zi in the Dirichlet distribution; will add machine epsilon to zero values
    zii = Dense(original_dim, activation='relu', kernel_initializer='random_normal')(layer1)
    zi_layer = Lambda(add_epsilon, output_shape=(original_dim,))(zii)
    # instantiate decoder model
    decoder = Model(latent_inputs, zi_layer, name='decoder')
    decoder.summary()
    # plot_model(decoder, to_file='vae_mlp_decoder.png', show_shapes=True)

    # instantiate VAE model
    outputs = decoder(encoder(inputs)[2])

    # construct VAE and train the model
    vae = Model(inputs, outputs, name='vae_mlp')



    vae.add_loss(y_z_loss(inputs, outputs))
    vae.compile(optimizer=op.RMSprop(lr = 0.001))
    vae.summary()
    # plot_model(vae,to_file='vae_mlp.png',show_shapes=True)

    callbacks = [LearningRateScheduler(lr_scheduler, verbose=1)] # this is the learning rate automatic reduction schedule
    callbacks = [tf.keras.callbacks.ReduceLROnPlateau(
        monitor='val_loss', factor=0.5, patience=20, verbose=1, mode='auto',
        min_delta=0.0001, cooldown=0, min_lr=0),
        tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=100)] # this is a scheduler based on validation loss; reduce learning rate if validation loss does not decrease in 20 epoches

    history = vae.fit(x=x_train,
                      callbacks = callbacks,
                      epochs=1000,
                      batch_size=batch_size,
                      validation_data=(x_test, None))



    # ## take a look at the trace of training loss


    # history.history.keys()
    # track = plt.figure()
    # f1 = track.add_subplot(1, 2, 1)
    # f2 = track.add_subplot(1, 2, 2)
    # f1.plot(history.history['loss'])
    # f2.plot(history.history['val_loss'])
    # f1.set_title('loss')
    # f2.set_title('val_loss')
    # f2.ticklabel_format(useOffset=False)
    # track.show()

    ## compute outputs (zi) based on input yi
    impute_train = vae.predict(x_train)
    # impute_train[0,1:5]
    # impute_train[1,1:5] # but it is weird that the predicted output is almost the same for all inputs, and values dont make sense


    # compute posterior mean of xi given zi and yi (assuming Dirichlet distribution for xi)
    impute_xi = (impute_train+x_train) / np.sum(impute_train+x_train, 1)[:,None]

    return {'impute_xi': impute_xi, 'impute_zi':impute_train, 'history_train':history.history['loss'], 'history_test': history.history['val_loss']}



##--------------------------
## start simulation
##--------------------------

#-----
# generate data, with sample size n and gene size p; count data follows:
# 1. multinomial logistic normal distribution;
# 2. Dirichlet multinomial distribution;
# both hyperparameter level zi based on Cao's low rank X_star
#-----

# input as: VAE_Simulation_small.py n p library_scale theta option z_scale
# use input parameter values from input_arg.txt

print('input for: VAE_Simulation_small.py n p library_scale theta option z_scale')
tmp = pd.read_csv('input_arg.txt', sep="\t", comment='#')
sys.argv = int(sys.argv[1])-1
print(tmp.iloc[sys.argv])

n=int(tmp.iloc[sys.argv,0]) # total number of cells
p=int(tmp.iloc[sys.argv,1])  # total number of gene entries-1
library_scale = float(tmp.iloc[sys.argv,2]) # relative scale of the total counts for each sample
theta = float(tmp.iloc[sys.argv,3])
option = tmp.iloc[sys.argv,4]
z_scale = float(tmp.iloc[sys.argv,5])


nreps=2

true_comp = []
true_count = []
VAE_comp = []
naive_comp = []

VAE_time = []


for i in range(nreps):
    data = data_generation(n=n,p=p,library_scale=library_scale,seed=2020+i, option=option, theta=theta, z_scale=z_scale)
    
    true_comp.append([data['xi']])
    true_count.append([data['yi']])

    tic = time.clock()
    VAE_res = build_network_and_train(x_train = data['yi'], x_test = data['yi'], intermediate_dim = 32)
    toc = time.clock()
    VAE_comp.append([VAE_res['impute_xi']])
    VAE_time.append(toc - tic)

    # rpy2.robjects.numpy2ri.activate()
    # ro.r("""
    #       filepath  = 'E:/Dropbox/Jing_Ali/other_code/composition-estimate-master/composition-estimate-master/'
    #       setwd(filepath)
    #       source(paste0(filepath, "proximal_gradient.R"))
    #       source(paste0(filepath,"tune_proximal_gradient.R"))
    #       source(paste0(filepath,"naive_methods.R"))
    # """)

    # autoTuneProxGradient = ro.r('''autoTuneProxGradient''')
    # tic = time.clock()
    # R_res = autoTuneProxGradient(data['yi'], 5)
    # Cao_res  = rpy2.robjects.numpy2ri.ri2py(R_res.rx('X_hat')[0])
    # toc = time.clock()
    # Cao_comp.append([Cao_res])
    # Cao_time.append(toc-tic)


    naive_xi = naive_compute(data['yi'])
    naive_comp.append([naive_xi])

    eval_frob(VAE_res['impute_xi'], data['xi'])
    # eval_frob(Cao_res, data['xi'])
    eval_frob(naive_xi, data['xi'])


filename='VAEres/'+str(n)+'_'+str(p)+'_'+str(library_scale)+'_'+str(theta)+'_'+option+'_'+str(z_scale)
print(filename)
np.save(filename+'_true_comp.npy', true_comp)
np.save(filename+'_true_count.npy', true_count)
np.save(filename+'_VAE_comp.npy', VAE_comp)
np.save(filename+'_naive_comp.npy', naive_comp)
np.save(filename+'_VAE_time.npy', VAE_time)