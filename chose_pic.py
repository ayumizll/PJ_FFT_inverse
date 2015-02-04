# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 12:48:35 2015

@author: Lili_ZHENG
"""
#==============================================================================
# function used to detect spectrum signal pics
#==============================================================================
import numpy as np

def chose_pic(spectre,phase,nbr_pic):
    n = np.size(spectre)
    out_spectre = np.zeros(n)
    out_phase = np.zeros(n)
    loc = np.where((spectre[1:-1]>spectre[2:]) & (spectre[1:-1]>spectre[:-2])) 
    out_spectre[loc]=spectre[loc]
    out_phase[loc]=phase[loc]
    return [loc,out_spectre,out_phase]