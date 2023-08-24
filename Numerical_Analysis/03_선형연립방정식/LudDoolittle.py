#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# LudDoolittle.py
# Doolittle법에 의한 LU 분해와 연립방정식 풀이
'''
    분해와 풀이 단계를 포함
    분해 단계 : 행렬[L / U]반환
    풀이 단계 : 
        전방대입 : b의 요소들은 y로 치환
        후방대입 : y는 해 x로 덮어쓰기
    
    A = decompDoolittle(A)
        LU 분해: [L][U] = [A]
    
    x = solveDoolittle(A, b)
        해석 단계: [L][U]{x} = {b}
'''
import numpy as np

# LU 분해
def decompDoolittle(A):
    n = len(A)
    for k in range(0, n-1):
        for i in range(k+1, n):
            if A[i, k] != 0.0:
                lam = A[i, k] / A[k, k]
                A[i, k+1:n] = A[i, k+1:n] - lam * A[k, k+1:n]
                A[i, k] = lam
    return A

# 후방대입법
def solveDoolittle(A, b):
    n = len(A)
    for k in range(1, n):
        b[k] = b[k] - np.dot(A[k, 0:k], b[0:k])
    b[n-1] = b[n-1] / A[n-1, n-1]
    
    for k in range(n-2, -1, -1):
        b[k] = (b[k] - np.dot(A[k, k+1:n], b[k+1:n])) / A[k, k]
    return b

