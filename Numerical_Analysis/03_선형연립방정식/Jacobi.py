#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jacobi.py
# Jacobi 반복법 모듈
'''
    x, it = Jacobi(A, b, tol = 1.0e-7)
        Jacobi 법으로 [A]{x} = {b} 풀기
'''
import math
import numpy as np

def Jacobi(A, b, tol = 1.0e-7):
    rows, cols = A.shape
    x = b.copy()
    xo = b.copy()
    for k in range(1, 501):
        for i in range(rows):
            x[i] = b[i]
            for j in range(cols):
                if (i == j) : continue
                x[i] -= A[i][j] * xo[j]
            x[i] /= A[i][i]
                
        # 수렴성 검토
        for i in range(rows):
            if math.fabs((x[i] - xo[i]) > tol):
                break
        # 수렴함
        else:
            return x, k
        
        xo = x.copy()
    else:
        print('Jacobi법은 수렴하지 않음')

