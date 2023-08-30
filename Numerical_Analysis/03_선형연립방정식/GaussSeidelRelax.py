#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# GaussSeidelRelax.py
# 이완이 있는 Gauss-Seidel 반복법 모듈
'''
    x, it = GaussSeidelRelax(A, b, omega = 1.0, tol = 1.0e-7)
        이완이 있는 Gauss-Seidel법으로 [A]{x} = {b} 풀기
'''
import math
import numpy as np

def GaussSeidelRelax(A, b, omega = 1.0, tol = 1.0e-7):
    rows, cols = A.shape
    x = np.zeros((rows), dtype=float)
    xo = np.zeros((rows), dtype=float)
    for k in range(1, 501):
        for i in range(rows):
            dx = b[i]
            for j in range(i):
                dx -= A[i][j] * x[j]
            for j in range(i, cols):
                dx -= A[i][j] * xo[j]
            x[i] = xo[i] + omega * dx / A[i][i]
            
        # 수렴성 검토
        for i in range(rows):
            # 수렴 안함
            if math.fabs((x[i] - xo[i]) > tol):
                break
        # 수렴함
        else:
            return x, k
        
        xo = x.copy()
        
    else:
        print('Gauss-Seidel법은 수렴하지 않음')

