'''
Descripttion: 
version: 
Author: sueRimn
Date: 2022-09-16 15:09:19
LastEditors: sueRimn
LastEditTime: 2022-09-16 16:09:41
'''
from keras.models import load_model
from tensorflow.keras.utils import plot_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np

input_sequence = 'GGGACCAGAT'
# Replace this with your custom sequence

context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 2))
models = [load_model(resource_filename('spliceai', x)) for x in paths]
# plot_model(models[1], to_file='/home/kfwang/002.png', show_shapes=True)
x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
y = np.mean([models[m].predict(x) for m in range(1)], axis=0)

print(y)
acceptor_prob = y[0, :, 1]
donor_prob = y[0, :, 2]


print(acceptor_prob)
print(donor_prob)