#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# LudCholesky.py
# Cholesky법에 의한 LU 분해와 연립방정식 풀이
'''
    분해 과정 중 대각요소가 음수이면 프로그램을 종료
    
    L = LudCholesky(A)
        Cholesky 분해: [L][L]T = [A]
    
    x = solveCholesky(L, b)
        Cholesky 분해에 의한 연립방정식 풀이
'''
import math
import sys
import numpy as np

# Cholesky법에 의한 LU 분해
def decompCholesky(A):
    n = len(A)
    for k in range(n):
        try:
            A[k, k] = math.sqrt(A[k, k] - np.dot(A[k, 0:k], A[k, 0:k]))
        except ValueError:
            print("행렬이 정부호가 아님")
            sys.exit()
        for i in range(k+1, n):
            A[i, k] = (A[i, k] - np.dot(A[i, 0:k], A[k, 0:k])) / A[k, k]
        for k in range(1, n):
            A[0:k, k] = 0.0
        return A
    
def solveCholesky(L, b):
    n = len(b)
    # [L]{y} = {b}의 풀이
    for k in range(n):
        b[k] = (b[k] - np.dot(L[k, 0:k], b[0:k])) / L[k, k]
    # [L_t]{x} = {y}의 풀이
    for k in range(n-1, -1, -1):
        b[k] = (b[k] - np.dot(L[k+1:n, k], b[k+1:n]))/L[k, k]
    return b

