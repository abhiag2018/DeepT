#!/usr/bin/env python
#keras version: keras-1.2.0

import sys
import os, re
import random
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
# from tensorflow.keras import initializations



########################### Input #############################
if len(sys.argv)<3:
    print( '[USAGE] python DeepTACT.py cell interaction_type num_DNase_experiments')
    print( 'For example, python DeepTACT.py demo P-E 3')
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
    print( '[USAGE] python DeepTACT.py cell interaction_type num_DNase_experiments')
    print( 'For example, python DeepTACT.py demo P-E 3')
    sys.exit()


######################## Initialization #######################
NUM_SEQ = 4
NUM_ENSEMBL = int(sys.argv[4])

########################### Training ##########################

def model_def():
    inp_region1_seq = Input(shape=(1, NUM_SEQ, RESIZED_LEN))
    inp_region2_seq = Input(shape=(1, NUM_SEQ, 1000))
    inp_region1_expr = Input(shape=(1, NUM_REP, RESIZED_LEN))
    inp_region2_expr = Input(shape=(1, NUM_REP, 1000))
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
    # merged.add(Merge([merged_seq, merged_DNase], mode = 'concat', concat_axis = -2))
    merged = Concatenate(axis=-2)([merged_seq, merged_DNase]) 
    ################# FIX THIS ?

    merged = Permute((2, 1))(merged)
    merged = BatchNormalization()(merged)
    merged = Dropout(drop_rate)(merged)

    # merged = Bidirectional(LSTM(100, return_sequences = True), merge_mode = 'concat')(merged)
    # merged.add(AttLayer())
    merged_1 = Bidirectional(LSTM(100, return_sequences=True), merge_mode="concat")(merged)
    merged_2 = Bidirectional(LSTM(100, return_sequences=True), merge_mode="concat")(merged)
    merged = tf.keras.layers.Attention()([merged_1,merged_2])

    # merged = Permute((2, 1))(merged)
    merged = Dense(1)(merged)
    merged = tf.keras.layers.Flatten()(merged)
    # merged = K.sum(merged,axis=-2)

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

# # load data: sequence
t=0
# CELL ='/projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.1/TrainingData/testCell_py3'
# TYPE='P-E'
# filename1='enhancer'
# filename2='promoter'
region1 = np.load(CELL+'/'+TYPE+'/bagData/'+filename1+'_Seq_'+str(t)+'.npz')
region2 = np.load(CELL+'/'+TYPE+'/bagData/'+filename2+'_Seq_'+str(t)+'.npz')
label = region1['label']
region1_seq = region1['sequence']
region2_seq = region2['sequence']

## load data: DNase
region1 = np.load(CELL+'/'+TYPE+'/bagData/'+filename1+'_DNase_'+str(t)+'.npz')
region2 = np.load(CELL+'/'+TYPE+'/bagData/'+filename2+'_DNase_'+str(t)+'.npz')
region1_expr = region1['expr']
region2_expr = region2['expr']

print(tf.test.gpu_device_name())

model = model_def()
print('compiling...')
model.compile(loss = 'binary_crossentropy',
              optimizer = optimizers.Adam(lr = 0.00001),
              metrics = ['acc', f1])
filename = CELL+'/'+TYPE+'/models/best_model_' + str(t) + '.h5'
modelCheckpoint = ModelCheckpoint(filename, monitor = 'val_acc', save_best_only = True, mode = 'max')
print('fitting...')
model.fit([region1_seq, region2_seq, region1_expr, region2_expr], label, epochs = 40, batch_size = 100, validation_split = 0.1, callbacks = [modelCheckpoint])


