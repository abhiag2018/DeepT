#!/usr/bin/env python
#keras version: keras-1.2.0

import sys
import os, re
import random
import itertools
import datetime
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
NUM_ENSEMBL = int(sys.argv[4])

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

def data_gen(path,index_file,lim_data=None):
    # path = f"{CELL}/{TYPE}/data"
    # index_file = f"{CELL}/{TYPE}/bagData/train_{t}.hkl"
    for i in itertools.islice(hkl.load(index_file),lim_data):
        fpath = f"{path}/{i}.npz"
        obj = np.load(fpath)
        input_dict = {}
        input_dict['enh_seq'] = tf.convert_to_tensor(obj['enh_seq'])
        input_dict['pr_seq'] = tf.convert_to_tensor(obj['pr_seq'])
        input_dict['enh_dnase'] = tf.convert_to_tensor(obj['enh_dnase'])
        input_dict['pr_dnase'] = tf.convert_to_tensor(obj['pr_dnase'])
        label = tf.convert_to_tensor(obj['label'])
        # enh_seq.set_shape([1,NUM_SEQ,RESIZED_LEN])
        # pr_seq.set_shape([1,NUM_SEQ,1000])
        # enh_dnase.set_shape([1,NUM_REP,RESIZED_LEN])
        # pr_dnase.set_shape([1,NUM_REP,1000])
        # label.set_shape([None])
        yield input_dict, label


def split_train_val_bootstrap(val=0.2):
    labels = np.load(f"{CELL}/{TYPE}/{filename1}_Seq.npz")['label']
    NUM = labels.shape[0]
    valsize = int(NUM*val)
    trainsize = NUM - valsize
    val_index = np.random.choice(range(NUM),valsize,replace=False)
    train_index = [item for item in range(NUM) if item not in val_index]
    for t in range(NUM_ENSEMBL):
        indices = np.random.choice(train_index,trainsize,replace=True)
        hkl.dump(indices, CELL+'/'+TYPE+'/bagData/train_'+str(t)+'.hkl')
    hkl.dump(val_index, CELL+'/'+TYPE+'/bagData/val.hkl')
    return 0

def bagging(t):
    traingen_callable = lambda:data_gen(f"{CELL}/{TYPE}/data",f"{CELL}/{TYPE}/bagData/train_{t}.hkl",lim_data=1000)
    valgen_callable = lambda:data_gen(f"{CELL}/{TYPE}/data",f"{CELL}/{TYPE}/bagData/val.hkl",lim_data=200)
    output_types = {'enh_seq':tf.float64, 'pr_seq':tf.float64, 'enh_dnase':tf.float64, 'pr_dnase':tf.float64}
    output_shapes = {'enh_seq':[1,NUM_SEQ,RESIZED_LEN], 'pr_seq':[1,NUM_SEQ,1000], 'enh_dnase':[1,NUM_REP,RESIZED_LEN], 'pr_dnase':[1,NUM_REP,1000]}
    train_set = tf.data.Dataset.from_generator(traingen_callable, output_types=(output_types, tf.int64), output_shapes = (output_shapes,[]))
    train_set = train_set.shuffle(50000).batch(32)
    val_set = tf.data.Dataset.from_generator(valgen_callable, output_types=(output_types, tf.int64), output_shapes = (output_shapes,[]))
    val_set = val_set.shuffle(50000).batch(32)

    # ## load data: sequence
    # region1 = np.load(CELL+'/'+TYPE+'/bagData/'+filename1+'_Seq_'+str(t)+'.npz')
    # region2 = np.load(CELL+'/'+TYPE+'/bagData/'+filename2+'_Seq_'+str(t)+'.npz')
    # label = region1['label']
    # region1_seq = region1['sequence']
    # region2_seq = region2['sequence']
    
    # ## load data: DNase
    # region1 = np.load(CELL+'/'+TYPE+'/bagData/'+filename1+'_DNase_'+str(t)+'.npz')
    # region2 = np.load(CELL+'/'+TYPE+'/bagData/'+filename2+'_DNase_'+str(t)+'.npz')
    # region1_expr = region1['expr']
    # region2_expr = region2['expr']

    model = model_def()
    print('compiling...')
    model.compile(loss = 'binary_crossentropy',
                  optimizer = optimizers.Adam(lr = 0.00001),
                  metrics = ['acc', f1])
    filename = CELL+'/'+TYPE+'/models_'+time_append+'/best_model_' + str(t) + '.h5'
    modelCheckpoint = ModelCheckpoint(filename, monitor = 'val_acc', save_best_only = True, mode = 'max')

    log_dir = CELL+'/'+TYPE+"/logs_" + time_append
    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)

    print('fitting...')
    model.fit(train_set, validation_data=val_set, epochs = 40, callbacks = [modelCheckpoint,tensorboard_callback])

def train():
    ### CHECK GPU usage
    print(tf.test.gpu_device_name())
    # for t in range(NUM_ENSEMBL):
    t=int(sys.argv[5])
    print(t)
    bagging(t)


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

def evaluate():
    region1_seq, region2_seq, region1_expr, region2_expr, label = load_test_data()#the same format as training data
    bag_pred = np.zeros((NUM_ENSEMBL,label.shape[0]))
    bag_score = np.zeros((NUM_ENSEMBL,label.shape[0]))
    for t in range(NUM_ENSEMBL):
        model.load_weights('./models_'+time_append+'/best_model_'+str(t)+'.h5')
        score = model.predict([region1_seq, region2_seq, region1_expr, region2_expr], batch_size = 100)
        bag_pred[t,:] = (score > 0.5).astype(int).reshape(-1)
        bag_score[t,:] = score.reshape(-1)
    f1, auprc = bag_pred(label, bag_pred, bag_score)
    return f1, auprc


############################ MAIN ###############################
time_append = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
os.system('mkdir -p '+CELL+'/'+TYPE+'/models_'+time_append)
os.system('mkdir -p '+CELL+'/'+TYPE+"/logs_" + time_append)
train()


