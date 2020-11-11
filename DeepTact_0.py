#!/usr/bin/env python
#keras version: keras-1.2.0

import sys
import os, re
import time
import random
import itertools
import datetime
import string
import numpy as np
import hickle as hkl
from sklearn import metrics

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Convolution2D, MaxPooling2D, Flatten
from tensorflow.keras.layers import LSTM, Bidirectional
from tensorflow.keras.layers import Dense, Activation

from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import BatchNormalization
# from tensorflow.keras.layers import Reshape, Merge, Permute
from tensorflow.keras.layers import Reshape, Concatenate, Permute
from tensorflow.keras import optimizers
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
import tensorflow.keras.backend as K
from tensorflow.keras.models import load_model
# from tensorflow.keras.engine.topology import Layer, InputSpec
from tensorflow.keras.layers import Layer
# from tensorflow.keras import initializations

import graphtools as gt
"""
DeepTACT.py

Training DeepTACT for P-P/P-E interactions

@author: abhiag
"""

######################## GPU Settings #########################
# gpu_use = str(0)#raw_input('Use gpu: ')
# gpu_cnmem = str(0.9)#raw_input('CNMeM: ')
# os.environ['THEANO_FLAGS'] = "warn.round=False,device=cuda"+gpu_use+",lib.cnmem="+gpu_cnmem
# os.environ['THEANO_FLAGS'] = "warn.round=False,device=gpu"+gpu_use+",lib.cnmem="+gpu_cnmem


########################### Input #############################
if len(sys.argv)<3:
    print('[USAGE] python DeepTACT.py cell interaction_type num_DNase_experiments')
    print('For example, python DeepTACT.py demo P-E 3')
    sys.exit()
CELL = sys.argv[1]
TYPE = sys.argv[2]
NUM_REP = int(sys.argv[3])
if TYPE == 'P-P':
    filename1 = 'promoter1'
    filename2 = 'promoter2'
    RESIZED_LEN = 1000 #promoter
elif TYPE == 'P-E':
    filename1 = 'enhancer'
    filename2 = 'promoter'
    RESIZED_LEN = 2000 #enhancer
else:
    print('[USAGE] python DeepTACT.py cell interaction_type num_DNase_experiments')
    print('For example, python DeepTACT.py demo P-E 3')
    sys.exit()


######################## Initialization #######################
NUM_SEQ = 4
# NUM_ENSEMBL = int(sys.argv[4])

########################### Training ##########################
# Attention GRU network
class AttLayer(Layer):
    def __init__(self, **kwargs):
        # self.init = initializations.get('normal')
        #self.input_spec = [InputSpec(ndim=3)]
        super(AttLayer, self).__init__(**kwargs)

    def build(self, input_shape):
        assert len(input_shape)==3
        self.W = self.add_weight(shape=(input_shape[-1],), initializer="random_normal", trainable=True)
        super(AttLayer, self).build(input_shape)  # be sure you call this somewhere!

    def call(self, x):
        M = K.tanh(x)
        alpha = K.dot(M,K.expand_dims(self.W, axis=-1))
        alpha = K.squeeze(alpha,axis=-1)#.dimshuffle(0,2,1)

        ai = K.exp(alpha)
        weights = ai/K.expand_dims(K.sum(ai, axis=1),axis=-1)
        weighted_input = x*K.expand_dims(weights,axis=-1)
        return K.tanh(K.sum(weighted_input,axis=1))

    def get_output_shape_for(self, input_shape):
        return (input_shape[0], input_shape[-1])



def model_def():
    inp_region1_seq = Input(shape=(1, NUM_SEQ, RESIZED_LEN),name='enh_seq')
    inp_region2_seq = Input(shape=(1, NUM_SEQ, 1000),name='pr_seq')
    inp_region1_expr = Input(shape=(1, NUM_REP, RESIZED_LEN),name='enh_dnase')
    inp_region2_expr = Input(shape=(1, NUM_REP, 1000),name='pr_dnase')
    drop_rate = 0.5 
    conv_enhancer_seq = Sequential()
    conv_enhancer_seq.add(Convolution2D(1024, (NUM_SEQ, 40), activation = 'relu', padding = 'valid', data_format = 'channels_first', input_shape = (1, NUM_SEQ, RESIZED_LEN)))
    conv_enhancer_seq.add(MaxPooling2D(pool_size = (1, 20), padding = 'valid', data_format = 'channels_first'))
    conv_enhancer_seq.add(Reshape((1024, (RESIZED_LEN-40+1)//20)))
    out_enh_seq = conv_enhancer_seq(inp_region1_seq)

    conv_promoter_seq = Sequential()
    conv_promoter_seq.add(Convolution2D(1024, (NUM_SEQ, 40), activation = 'relu', padding = 'valid', data_format = 'channels_first', input_shape = (1, NUM_SEQ, 1000)))
    conv_promoter_seq.add(MaxPooling2D(pool_size = (1, 20), padding = 'valid', data_format = 'channels_first'))
    conv_promoter_seq.add(Reshape((1024, 48)))
    out_pr_seq = conv_promoter_seq(inp_region2_seq)

    merged_seq = Concatenate()([out_enh_seq, out_pr_seq])

    conv_enhancer_DNase = Sequential()
    conv_enhancer_DNase.add(Convolution2D(1024, (NUM_REP, 40), activation = 'relu', padding = 'valid', data_format = 'channels_first', input_shape = (1, NUM_REP, RESIZED_LEN)))
    conv_enhancer_DNase.add(MaxPooling2D(pool_size = (1, 20), padding = 'valid', data_format = 'channels_first'))
    conv_enhancer_DNase.add(Reshape((1024, (RESIZED_LEN-40+1)//20)))
    out_enh_DNase = conv_enhancer_DNase(inp_region1_expr)

    conv_promoter_DNase = Sequential()
    conv_promoter_DNase.add(Convolution2D(1024, (NUM_REP, 40), activation = 'relu', padding = 'valid', data_format = 'channels_first', input_shape = (1, NUM_REP, 1000)))
    conv_promoter_DNase.add(MaxPooling2D(pool_size = (1, 20), padding = 'valid', data_format = 'channels_first'))
    conv_promoter_DNase.add(Reshape((1024, 48)))
    out_pr_DNase = conv_promoter_DNase(inp_region2_expr)

    merged_DNase = Concatenate()([out_enh_DNase, out_pr_DNase])
    #   
    merged = Concatenate(axis=-2)([merged_seq, merged_DNase]) 

    merged = Permute((2, 1))(merged)
    merged = BatchNormalization()(merged)
    merged = Dropout(drop_rate)(merged)
    merged = Bidirectional(LSTM(100, return_sequences=True), merge_mode="concat")(merged)
    merged = AttLayer()(merged)
    merged = BatchNormalization()(merged)
    merged = Dropout(drop_rate)(merged)
    merged = Dense(925)(merged)
    merged = BatchNormalization()(merged)
    merged = Activation('relu')(merged)
    merged = Dropout(drop_rate)(merged)
    merged = Dense(1, activation = 'sigmoid')(merged)

    model = Model(inputs=[inp_region1_seq, inp_region2_seq, inp_region1_expr, inp_region2_expr], outputs=merged)
    return model

def f1(y_true, y_pred):
    cast = lambda x:tf.keras.backend.cast(x,dtype='float64')
    TP = K.sum(cast(K.equal(y_true, 1) & K.equal(K.round(y_pred), 1)))
    FP = K.sum(cast(K.equal(y_true, 0) & K.equal(K.round(y_pred), 1)))
    FN = K.sum(cast(K.equal(y_true, 1) & K.equal(K.round(y_pred), 0)))
    TN = K.sum(cast(K.equal(y_true, 0) & K.equal(K.round(y_pred), 0)))
    P = TP / (TP + FP + K.epsilon())
    R = TP / (TP + FN + K.epsilon())
    F1 = 2 * P * R / (P + R + K.epsilon())
    return F1

def true_labels(path,index_file,lim_data=None):
    # path = f"{CELL}/{TYPE}/data"
    # index_file = f"{CELL}/{TYPE}/bagData/train_{t}.hkl"
    labels = []
    for i in itertools.islice(hkl.load(index_file),lim_data):
        fpath = f"{path}/{i}.npz"
        labels.append(np.load(fpath)['label'])
    return np.array(labels)

def data_gen(path,index_file,lim_data=None, check_num_rep = True):
    # path = f"{CELL}/{TYPE}/data"
    # index_file = f"{CELL}/{TYPE}/bagData/train_{t}.hkl"
    for i in itertools.islice(hkl.load(index_file),lim_data):
        fpath = f"{path}/{i}.npz"
        obj = np.load(fpath)
        input_dict = {}
        input_dict['enh_seq'] = tf.convert_to_tensor(obj['enh_seq'])
        input_dict['pr_seq'] = tf.convert_to_tensor(obj['pr_seq'])
        if check_num_rep:
            assert obj['enh_dnase'].shape[1] == NUM_REP
            assert obj['pr_dnase'].shape[1] == NUM_REP
        input_dict['enh_dnase'] = tf.convert_to_tensor(obj['enh_dnase'])[:,:NUM_REP,:]
        input_dict['pr_dnase'] = tf.convert_to_tensor(obj['pr_dnase'])[:,:NUM_REP,:]
        label = tf.convert_to_tensor(obj['label'])
        # enh_seq.set_shape([1,NUM_SEQ,RESIZED_LEN])
        # pr_seq.set_shape([1,NUM_SEQ,1000])
        # enh_dnase.set_shape([1,NUM_REP,RESIZED_LEN])
        # pr_dnase.set_shape([1,NUM_REP,1000])
        # label.set_shape([None])
        yield input_dict, label


def split_train_val_bootstrap(test=0.2):
    labels = np.load(f"{CELL}/{TYPE}/{filename1}_Seq.npz")['label']
    NUM = labels.shape[0]
    testsize = int(NUM*test)
    trainsize = NUM - testsize
    random_perm = np.random.permutation(range(NUM))

    test_index = random_perm[0:testsize]
    train_index = random_perm[testsize:]

    hkl.dump(train_index, CELL+'/'+TYPE+'/train.hkl')
    hkl.dump(test_index, CELL+'/'+TYPE+'/test.hkl')
    return 0

def bagging():
    traingen_callable = lambda:data_gen(f"{CELL}/{TYPE}/data", f"{LOG_DIR}/train.hkl",lim_data=None, check_num_rep = False)
    valgen_callable = lambda:data_gen(f"{CELL}/{TYPE}/data", f"{LOG_DIR}/val.hkl",lim_data=None, check_num_rep = False)
    output_types = {'enh_seq':tf.float64, 'pr_seq':tf.float64, 'enh_dnase':tf.float64, 'pr_dnase':tf.float64}
    output_shapes = {'enh_seq':[1,NUM_SEQ,RESIZED_LEN], 'pr_seq':[1,NUM_SEQ,1000], 'enh_dnase':[1,NUM_REP,RESIZED_LEN], 'pr_dnase':[1,NUM_REP,1000]}
    train_set = tf.data.Dataset.from_generator(traingen_callable, output_types=(output_types, tf.int64), output_shapes = (output_shapes,[]))
    shuffle_buffer_size = 256#len(hkl.load(f"{LOG_DIR}/train.hkl"))
    train_set = train_set.shuffle(shuffle_buffer_size, reshuffle_each_iteration=True).batch(BATCH_SIZE)
    val_set = tf.data.Dataset.from_generator(valgen_callable, output_types=(output_types, tf.int64), output_shapes = (output_shapes,[]))
    val_set = val_set.shuffle(shuffle_buffer_size, reshuffle_each_iteration=True).batch(BATCH_SIZE)

    model = model_def()
    print('compiling...')
    model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers.Adam(lr = 0.00001),
                  metrics = ['acc', f1])
    filename = CELL+'/'+TYPE+'/models_'+time_append+'/best_model.h5'
    modelCheckpoint = ModelCheckpoint(filename, monitor = 'val_acc', save_best_only = True, mode = 'max')

    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=LOG_DIR, histogram_freq=1)

    print('fitting...')
    model.fit(train_set, validation_data=val_set, epochs = 40, callbacks = [modelCheckpoint,tensorboard_callback])

def train():
    ### CHECK GPU usage
    print(tf.config.list_physical_devices('GPU'))
    print(tf.test.gpu_device_name())

    VAL_FRAC = 0.2
    time.sleep(np.random.uniform(np.random.uniform(high=3.0))) #avoid file lock for .hkl files
    train_index = hkl.load(CELL+'/'+TYPE+'/train.hkl')
    random_perm = np.random.permutation(range(len(train_index)))
    valsize = int(len(train_index)*VAL_FRAC)
    val_index = random_perm[0:valsize]
    train_index = random_perm[valsize:]

    hkl.dump(train_index, LOG_DIR+'/train.hkl')
    hkl.dump(val_index, LOG_DIR+'/val.hkl')

    bagging()


########################### Evaluation ##########################
def bag_pred(label, bag_pred, bag_score):
    vote_pred = np.zeros(bag_pred.shape[1])
    vote_score = np.zeros(bag_score.shape[1])
    for i in range(bag_pred.shape[1]):
        vote_pred[i] = stats.mode(bag_pred[:,i]).mode
        vote_score[i]= np.mean(bag_score[:,i])
    f1 = metrics.f1_score(label, vote_pred)
    auprc = metrics.average_precision_score(label, vote_score)
    return f1, auprc

def evaluate(eval_cell, bootstrap_time):
    eval_cell_path = '/'.join(CELL.split('/')[:-1])+'/'+eval_cell + '/' + TYPE

    NUM_ENSEMBL=len(bootstrap_time)

    output_types = {'enh_seq':tf.float64, 'pr_seq':tf.float64, 'enh_dnase':tf.float64, 'pr_dnase':tf.float64}
    output_shapes = {'enh_seq':[1,NUM_SEQ,RESIZED_LEN], 'pr_seq':[1,NUM_SEQ,1000], 'enh_dnase':[1,NUM_REP,RESIZED_LEN], 'pr_dnase':[1,NUM_REP,1000]}


    model = model_def()
    model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers.Adam(lr = 0.00001),
                  metrics = ['acc', f1])

    test_labels = true_labels(f"{eval_cell_path}/data", f"{eval_cell_path}/test.hkl",lim_data=None)

    avg_score = np.zeros((test_labels.shape[0],1))

    for time_append in bootstrap_time[:NUM_ENSEMBL]:
        model.load_weights(CELL+'/'+TYPE+f"/models_{time_append}/best_model.h5")

        testgen_callable = lambda:data_gen(f"{eval_cell_path}/data", f"{eval_cell_path}/test.hkl",lim_data=None, check_num_rep = False)

        test_set = tf.data.Dataset.from_generator(testgen_callable, output_types=(output_types, tf.int64), output_shapes = (output_shapes,[]))
        test_set = test_set.batch(BATCH_SIZE)
        
        score = model.predict(test_set)
        avg_score = avg_score + score

    avg_score = avg_score/NUM_ENSEMBL    
    test_pred = (avg_score>0.5).astype(int).reshape(-1)
    outpath = CELL+'/'+TYPE+f'/bootstrap_out_{eval_cell}.png'
    gt.plot_roc(outpath, f"train = {CELL}; test = {eval_cell}", test_labels, test_pred)

    return outpath
# def evaluate():
#     region1_seq, region2_seq, region1_expr, region2_expr, label = load_test_data()#the same format as training data
#     bag_pred = np.zeros((NUM_ENSEMBL,label.shape[0]))
#     bag_score = np.zeros((NUM_ENSEMBL,label.shape[0]))
#     for t in range(NUM_ENSEMBL):
#         model.load_weights('./models_'+time_append+'/best_model_'+str(t)+'.h5')
#         score = model.predict([region1_seq, region2_seq, region1_expr, region2_expr], batch_size = 100)
#         bag_pred[t,:] = (score > 0.5).astype(int).reshape(-1)
#         bag_score[t,:] = score.reshape(-1)
#     f1, auprc = bag_pred(label, bag_pred, bag_score)
#     return f1, auprc


############################ MAIN ###############################
BATCH_SIZE=32
time_append = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")+'-'+''.join(random.choices(string.ascii_uppercase + string.digits,k=10))
os.makedirs(CELL+'/'+TYPE+'/models_'+time_append, exist_ok=True)
LOG_DIR = CELL+'/'+TYPE+"/logs_" + time_append
os.makedirs(LOG_DIR, exist_ok=True)
train()

# split_train_val_bootstrap()


# bootstrap_time = ['20201030-102007', '20201030-062033', '20201030-052652', '20201029-140022', '20201029-100734', 
#     '20201029-080806', '20201029-052739', '20201029-045235', '20201029-003927', '20201028-050134', '20201028-025051', 
#     '20201028-021157', '20201028-001550', '20201027-194725', '20201027-143432', '20201026-133133', '20201026-132933', 
#     '20201026-133436', '20201026-131102', '20201026-131111', '20201023-150244', '20201023-150253', '20201023-143322', 
#     '20201023-142612', '20201023-142204', '20201023-142053']
# EVAL_CELL = sys.argv[4]
# evaluate(EVAL_CELL, bootstrap_time)
