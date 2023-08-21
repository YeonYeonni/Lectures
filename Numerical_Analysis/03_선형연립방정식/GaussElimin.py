#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# GaussElimin()
'''
x = GaussElimin(A, b, prt)
    방정식 [A]{b} = {x}를 Gauss 소거법으로 푼다.
    A : 계수행렬
    b : 우변벡터
    prt : 계산의 중간 출력
'''

import numpy as np
from PrintMatrix import *

def GaussElimin(A, b, prt = False):
    n = len(b)
    
    # 전방 소거 단계
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            if A[i, k] != 0.0:
                lam = A[i, k] / A[k, k]
                A[i, k + 1 : n] -= lam * A[k, k + 1 : n]
                b[i] -= lam * b[k]
    
    # 계산중간 출력
    if (prt):
        print("전방 소거 후")
        printEqs(A, b)
    
    # 후방대입 단계
    for k in range(n - 1, -1, -1):
        b[k] = (b[k] - np.dot(A[k, k + 1 : n], b[k + 1 : n])) / A[k, k]
        
    return b

