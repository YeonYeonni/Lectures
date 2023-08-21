#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# PringMatrix.py
'''
    선형 연립방정식의 출력
'''

import numpy as np

def printEqs(A, b):
    n = b.size
    for i in range(n):
        for j in range(n):
            print("{0:10.3e} ".format(A[i][j]), end = "  ")
        print("|   {0:10.3e}".format(b[i]))
        
def printMat(A):
    n, m = A.shape
    for i in range(n):
        for j in range(m):
            print("{0:10.3e} ".format(A[i][j]), end = "")
        print("")

