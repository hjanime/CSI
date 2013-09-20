'''
Created on 17 Sep, 2013

@author: caofan
'''

import numpy as np


class VDPGM:
    def __init__(self):
        pass
    
def mk_log_likelihood(data, hp_posterior, hp_prior, opts):
    D,N = data.shape()