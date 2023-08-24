#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# LudSym5.py 모듈
'''
    대칭오대각행렬 풀이
    decompLudSym5() : 대칭오대각행렬 A를 A = [f / e / d / e / f] 형태로 분해
    원벡터 d, e, f는 변경. 분해된 행렬의 벡터로 대치
    분해 뒤에 연립방정식 Ax = b의 해는 solveLudSym5()에 의해 얻을 수 있음
    전방대입 동안, 원벡터 b는 y로 치환
    후방대입 동안, y는 x로 덮어 씌움
    solveLudSym5()를 나갈 때 b는 해벡터를 반환
    
    x = LudFSym5A(A, b, prt)
        대칭오대각행렬을 입력하여 풀이
        [A]: 오대각행렬
        {b}: 우변벡터
        {x}: 오대각 선형연립방정식의 풀이
        
    d,e,f = decompLudSym5(d,e,f)
        대칭오대각행렬 [A]의 LU 분해
        여기서 {f}, {e}, {d}는 [A]의 대각선
        출력에서 {d}, {e}, {f}는 분해행렬의 대각선

    x = solveLudSym5(d,e,f,b)
        [A]{x} = {b}의 해
        여기서 {d}, {e}, {f}는 decompLud5()에서 반환된 벡터
'''

import numpy as np

def decompLudSym5(d, e, f):
    n = len(d)
    for k in range(n-2):
        lam = e[k] / d[k]
        d[k+1] = d[k+1] - lam*e[k]
        e[k+1] = e[k+1] - lam*f[k]
        e[k] = lam
        lam = f[k] / d[k]
        d[k+2] = d[k+2] - lam*f[k]
        f[k] = lam
    lam = e[n-2] / d[n-2]
    d[n-1] = d[n-1] - lam*e[n-2]
    e[n-2] = lam
    return d, e, f

def solveLudSym5(d, e, f, b):
    n = len(d)
    b[1] = b[1] - e[0]*b[0]
    for k in range(2,n):
        b[k] = b[k] - e[k-1]*b[k-1] - f[k-2]*b[k-2]
        
    b[n-1] = b[n-1] / d[n-1]
    b[n-2] = b[n-2] / d[n-2] - e[n-2]*b[n-1]
    for k in range(n-3, -1, -1):
        b[k] = b[k] / d[k] - e[k]*b[k+1] - f[k]*b[k+2]
    return b

# 오대각행렬을 3개의 벡터로 분리
def Penta2Vecs(A):
    n, _ = A.shape
    d = np.zeros(n)
    e = np.zeros(n)
    f = np.zeros(n)
    
    for i in range(1, n-2):
        d[i] = A[i][i]
        e[i] = A[i][i+1]
        f[i] = A[i][i+2]
    d[n-2] = A[n-2][n-2]
    e[n-2] = A[n-2][n-1]
    d[n-1] = A[n-1][n-1]
    return d, e, f

# 3개의 벡터를 오대각행렬로 결합
def Vecs2Penta(d, e, f):
    n = d.size
    A = np.zeros((n,n))
    
    # 상삼각행렬 만들기
    for i in range(n-2):
        A[i][i] = d[i]
        A[i][i+1] = e[i]
        A[i][i+1] = f[i]
    A[n-2][n-2] = d[n-2]
    A[n-2][n-1] = e[n-2]
    A[n-1][n-1] = d[n-1]
    
    # 상삼각 -> 대칭
    for i in range(1, n):
        for j in range(i):
            A[i][j] = A[j][i]
    return A

