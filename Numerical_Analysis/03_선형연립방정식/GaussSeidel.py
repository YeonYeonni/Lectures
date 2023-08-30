#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# GaussSeidel.py
# Gauss-Seidel 반복법 모듈
'''
    x, it = GaussSeidel(A, b, tol = 1.0e-7)
        Gauss-Seidel 법으로 [A]{x} = {b} 풀기
'''
import math
import numpy as np

def GaussSeidel(A, b, tol = 1.0e-7):
    rows, cols = A.shape
    x = np.zeros((rows), dtype=float)
    xo = np.zeros((rows), dtype=float)
    for k in range(1, 501):
        for i in range(rows):
            x[i] = b[i]
            for j in range(i):
                x[i] -= A[i][j] * x[j]
            for j in range(i+1, cols):
                x[i] -= A[i][j] * xo[j]
            x[i] = x[i] / A[i][i]
            
        # 수렴성 검토
        for i in range(rows):
            if math.fabs((x[i] - xo[i]) > tol):
                break
        # 수렴함
        else:
            return x, k
        
        xo = x.copy()
    else:
        print('Gauss-Seidel법은 수렴하지 않음')

