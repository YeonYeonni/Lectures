#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Lud3.py 모듈
'''
    x = Lud3A(A, b, prt)
        삼대각행렬을 입력하여 풀이
        [A] : 삼대각행렬
        {b} : 우변벡터
        {x} : 삼대각 선형연립방정식의 풀이
    
    x = Lud3V(c, d, e, b, prt):
        삼대각행렬을 3개의 벡터로 입력하여 풀이
        
    c, d, e = decompLud3(c, d, e, prt)
        삼대각행렬 [A]의 LU 분해
        여기서 {c}, {d}, {e}는 [A]의 대각선
        출력 {c}, {d}, {e}는 분해행렬의 대각선
    
    x = solveLud3(c, d, e, b)
        [A]{x} {b}의 해
        여기서 {c}, {d}, {e}는 LUdecomp3()의 반환벡터
    
    c, d, e = Tri2Vecs(A) : 삼대각행렬을 3개의 벡터로 분리
    
    A = Vecs2Tri(c, d, e) : 3개의 벡터를 삼대각행렬로 합병
'''
import numpy as np
from PrintMatrix import *

# 삼대각행렬을 입력하여 방정식 풀이
def Lud3A(A, b, prt=False):
    c, d, e = Tri2Vecs(A)
    c, d, e = decompLud3(c, d, e, prt)
    b = solveLud3(c, d, e, b)
    return b

# 삼대각행렬을 벡터로 입력하여 방정식 풀이
def Lud3V(c, d, e, b, prt=False):
    c, d, e = decompLud3(c, d, e, prt)
    b = solveLud3(c, d, e, prt)
    return b

# 삼대각 벡터의 LU 분해
def decompLud3(c, d, e, prt=False):
    n = len(d)
    for k in range(1, n):
        lam = c[k] / d[k-1]
        d[k] = d[k] - lam * e[k-1]
        c[k] = lam
        
    # 분해결과 출력
    if (prt):
        print('LU 분해된 삼대각행렬')
        printMat(Vecs2Tri(c, d, e))
        
    return c, d, e

# LU 분해된 벡터의 풀이
def solveLud3(c, d, e, b):
    n = len(d)
    for k in range(1, n):
        b[k] = b[k] - c[k] * b[k-1]
    b[n-1] = b[n-1] / d[n-1]
    for k in range(n-2, -1, -1):
        b[k] = (b[k] - e[k] * b[k+1]) / d[k]
    return b

# 삼대각행렬을 3개의 벡터로 분리
def Tri2Vecs(A):
    n, _ = A.shape
    c = np.zeros(n)
    d = np.zeros(n)
    e = np.zeros(n)
    
    d[0] = A[0][0]
    e[0] = A[0][1]
    
    for i in range(1, n-1):
        c[i] = A[i][i-1]
        d[i] = A[i][i]
        e[i] = A[i][i+1]
    c[n-1] = A[n-1][n-2]
    d[n-1] = A[n-1][n-1]
    return c, d, e

# 3개의 벡터를 삼대각행렬로 결합
def Vecs2Tri(c, d, e):
    n = d.size
    A = np.zeros((n, n))
    
    A[0][0] = d[0]
    A[0][1] = e[0]
    
    for i in range(1, n-1):
        A[i][i-1] = c[i]
        A[i][i] = d[i]
        A[i][i+1] = e[i]
    A[n-1][n-2] = c[n-1]
    A[n-1][n-1] = d[n-1]
    return A

