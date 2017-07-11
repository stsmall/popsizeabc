'''
Module for GenomeWide Association Studies analysis
'''
import numpy as np


missing=np.int16(-1)

def _complete_cases(x):
    return x!=missing

complete_cases=np.vectorize(_complete_cases)

def _is_na(x):
    return x==missing

is_na=np.vectorize(_is_na)

