from keras.layers import Lambda, Input, Dense, advanced_activations, BatchNormalization, Dropout
from keras.models import Model
from keras.utils import plot_model
from keras import backend as K
from keras import optimizers as op
import tensorflow as tf
import numpy as np
import pandas as pd
import time
from keras.callbacks import LearningRateScheduler


##-------------------------
## Definition of functions
##-------------------------

# reparametrization trick for generating a normal random vector
# Here is for generating the posterior latent normal vector wi, which will be passed to a neural net to generate zi
def w_sampling(args):
    w_mean, w_log_var = args
    batch = K.shape(w_mean)[0]
    dim = K.int_shape(w_mean)[1]
    # by default, random_normal has mean = 0 and std = 1.0
    epsilon = K.random_normal(shape=(batch, dim))
    return w_mean + K.exp(0.5 * w_log_var) * epsilon

# this defines generation of logistic normal vector xi of dimension p+1, with a mean vector, and a single scalar value for the common variance (the covariance matrix is diagonal)
# later will input zi as the mean vector, and set common variance theta2 as a parameter to estimate
def x_sampling(args):
    x_mean, x_log_var = args
    batch = K.shape(x_mean)[0]
    dim = K.int_shape(x_mean)[1]
    # by default, random_normal has mean = 0 and std = 1.0
    epsilon = K.random_normal(shape=(batch, dim))
    x_mean + K.exp(0.5 * x_log_var) * epsilon
    return 0

# this transformation of the input counts data was used in one reference paper for VAE; it might help
def normalization_transform(x):
    tmp = K.sum(x, 1)
    return K.log(x/K.reshape(tmp,(-1, 1)) * 10000 + 1)

# this is logistic transformation without additional weight parameters
def logistic_transform(z):
    tmp = K.concatenate([z, [K.constant(0)]], 0)
    x = tmp/(K.sum(K.exp(tmp)))
    return x

# this try to make zi output positive
def add_epsilon(z):
    z = K.abs(z) + K.constant(np.finfo(float).eps)
    return z


# This is a sample of a learning rate scheduler from online reference; decrease learning rate per 10 epoch
def lr_scheduler(epoch, lr):
    decay_rate = 0.85
    decay_step = 10
    if epoch % decay_step == 0 and epoch:
        return lr * pow(decay_rate, np.floor(epoch / decay_step))
    return lr

# This helps to change activation function paramters in the layers
class activation_wrapper(object):
    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        def _func(x):
            return self.func(x, *args, **kwargs)
        return _func

# this is generation of zi from a randomly defined neural network
def generate_zi_NN(wi, input_node, layer_1_node, layer_2_node, output_node, last_activation):
    n = wi.shape[0]
    weight_layer_1 = np.random.normal(1, 0.1, (input_node+1)*layer_1_node).reshape((input_node+1,layer_1_node))
    bias_layer_1 = np.random.gamma(1, 1, 1).reshape((1,1))
    layer_1_y = np.dot(np.concatenate([wi, np.repeat(bias_layer_1,n).reshape((n,1))], 1), weight_layer_1)
    layer_1_res = tf.keras.activations.relu(layer_1_y, alpha=0.1)

    weight_layer_2 = np.random.normal(-1, 0.1, (layer_1_node+1)*layer_2_node).reshape((layer_1_node+1,layer_2_node))
    bias_layer_2 = np.random.normal(0,2,1).reshape((1,1))
    layer_2_y = np.dot(np.concatenate([layer_1_res, np.repeat(bias_layer_2,n).reshape((n,1))], 1), weight_layer_2)
    layer_2_res = tf.keras.activations.relu(layer_2_y, alpha=0.1)

    weight_layer_output = np.random.normal(-1, 0.1, (layer_2_node+1)*output_node).reshape((layer_2_node+1,output_node))
    bias_layer_output = np.random.normal(-2,2,1).reshape((1,1))
    layer_output_y =  np.dot(np.concatenate([layer_2_res, np.repeat(bias_layer_output,n).reshape((n,1))], 1), weight_layer_output)
    # scale the output
    layer_output_y = layer_output_y / np.sqrt(np.var(layer_output_y))
    layer_output_res = last_activation(layer_output_y)
    return K.eval(layer_output_res)

# this follows idea of Cao's paper for genrating zi
# r is controlling low rank of the 'composition'
def generate_zi_Cao(n, p, r=10, option='logisticnormal', z_scale = 2):
    X_star = np.zeros((n, p))
    if np.sum(X_star<=0)>0:
        U = np.abs(np.random.normal(size=n*r).reshape((n, r)))
        V1 = np.random.uniform(low=0, high=1, size=p*r).reshape((p,r))
        index = V1<0.55
        V1[index] = 1
        V1[~index] = 0
        np.fill_diagonal(V1, 1)
        V2 = np.random.normal(loc=0, scale=10**(-3), size=p*r).reshape((p, r)) # changed scale value
        V = V1+V2
        Z = np.dot(U, V.transpose())   
        # should not have zero entries in the Z or X_star, but if continue sampling until positive it is stuck. 
        # here need modification  
        X_star = Z / np.sum(Z, 1)[:, None] 
    if option == 'logisticnormal':
        # use alr transformation
        zi = alr_trans(X_star)
        return zi
    if option == 'dirichlet':
        # multiply a scaling parameter 
        zi = X_star * z_scale
        return zi
    print('zi generation error')


# this is generation of xi vector from 1. Logistic normal; 2. dirichlet
# we want to have all xi being positive
def generate_xi(n, p, zi, option, theta=1):
    xi = 0
    if (option == 'logisticnormal'):
        while (np.sum(xi <= 0) > 0):
            logi_normal = np.random.normal(0, 1, size=n * (p - 1)).reshape((n, (p - 1))) * np.sqrt(theta) + zi
            tmp = np.concatenate([logi_normal, np.repeat(0, n).reshape((n, 1))], 1)
            xi = np.exp(tmp) / np.sum(np.exp(tmp), 1)[:, None]
        return xi
    if (option == 'dirichlet'):
        while (np.sum(xi <= 0) > 0):
            dirichlet_zi = zi
            xi = np.vstack(pd.Series(range(n)).apply(lambda i: np.random.dirichlet(alpha=dirichlet_zi[i, :])))
        return xi
    print('xi generation error')

# this is generation of counts from given total count vector and probability vector
def generate_yi(n, counts, xi):
    x_all = pd.Series(range(n)).apply(lambda x: np.random.multinomial(counts[x], xi[x,]))
    x_all = np.vstack(x_all.values)
    return x_all

# evaluation metrics for outputs
def eval_mse(a, b):
    return 'need repeated simulation'

def eval_frob(a, b):
    return (np.sqrt(np.sum((a - b) ** 2)))

def eval_l1(a, b):
    return np.mean(np.sum(np.abs(a-b), 1))

def comp_shannon(a):
    tmp = np.array(np.sum(-a*np.log(a), 1))
    return tmp

def comp_simpson(a):
    return np.array(np.sum(a**2, 1))

def eval_KL(a, b):
    return np.sum(a*(np.log(a)-np.log(b)))

def eval_shannon(a, b):
    return np.sum((comp_shannon(a) - comp_shannon(b))**2)

def eval_simpson(a, b):
    return np.sum((comp_simpson(a) - comp_simpson(b))**2)

# define VAE loss based on output
# with zi output, use log likelihood of the counts
def y_z_loss(y_actual, y_pred):
    # y_actual = y_actual + K.constant(np.finfo(float).eps)
    # y_pred = y_pred + K.constant(np.finfo(float).eps)
    m = K.sum(y_actual)
    part1 = tf.math.lgamma(m + 1) - K.sum(tf.math.lgamma(y_actual + 1))
    part2 = K.sum(tf.math.lgamma(y_pred)) - tf.math.lgamma(K.sum(y_pred))
    part3 = tf.math.lgamma(K.sum(y_pred)+K.sum(y_actual)) - K.sum(tf.math.lgamma(y_pred + y_actual))
    log_y_z = part1 - part2 - part3
    return -log_y_z / K.sqrt(K.constant(n))  # try no regularization with KL distance for the latent layer; scale the negative log likelihood

def data_generation(n=300,p=100,library_scale=5,seed=2020, option='logisticnormal', theta=8, z_scale=2, r=10):
    # this library_scale will be the gamma parameter as in Cao's paper
    np.random.seed(int(seed))
    n = int(n)
    p = int(p)
    library_scale = int(library_scale)

    # total counts for each cell (mi); follow Cao's paper setting
    counts = pd.Series(range(n)).apply(lambda x: np.random.poisson(np.random.uniform(low=library_scale*p, high=10*library_scale*p, size=1), size=1))
    counts = np.vstack(counts)

    # # generate hyperparameter wi from normal(0,Id)
    # d=10
    # wi = np.random.normal(0, 1, d*n).reshape((n, d))

    # # generate a simple two layer neural net with randomly generated weight and bias parameters
    # zi = generate_zi_NN(wi, input_node=d, layer_1_node=64, layer_2_node=64, output_node=p, last_activation = tf.keras.activations.tanh)

    # We follow Cao's paper for generation of an underlying composition matrix X*, and from that we generate zi:

    zi = generate_zi_Cao(n, p, r=r, option=option, z_scale = z_scale)

    # generate xi from zi;
    # logistic normal will have diagonal covariance with var=theta; dimension p+1
    # Dirichlet will have zi positive; dimension p

    xi = generate_xi(n, p, zi, option=option, theta=theta)



    ## then generate count data from multinomial
    yi = generate_yi(n, counts, xi)
    print('count sparsity: ',  np.mean(yi==0)) # sparsity of yi

    ## we do not need test set here, only for imputation purpose
    x_train = yi
    prob_vec_train = xi
    x_test = yi
    prob_vec_test = xi

    # np.savetxt("x_all.csv", yi, delimiter=',', fmt="%s")

    return {'yi':yi, 'xi':xi, 'zi':zi}

def naive_compute(yi, add=0.5):
    yi = np.array(yi).astype(np.float)
    yi[yi==0] = add
    return (yi)/np.sum(yi, 1)[:,None] # naive imputed prob vecs

def alr_trans(impute_xi):
    xi_len = impute_xi.shape[1]
    return np.log(impute_xi[:,0:(xi_len-1)]) - np.log(impute_xi[:,xi_len-1])[:,None]


n=300 # total number of cells
p=100  # total number of gene entries-1
library_scale = 5 # relative scale of the total counts for each sample
theta = 8
option = 'logisticnormal'

i=1

# data = data_generation(n=n,p=p,library_scale=library_scale,seed=2020+i, option=option, theta=theta, z_scale=1.2)

# print(data)
