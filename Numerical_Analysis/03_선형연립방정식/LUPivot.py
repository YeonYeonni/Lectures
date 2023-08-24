#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## LUPivot 모듈
'''
    LUdecomp() : 분해 단계 동안 행교환의 기록을 seq 배열에 유지
        초기 seq은 [0, 1, 2, ...]
        두 행이 교환되면 대응하는 요소를 seq에서도 교환
        seq는 원래의 행들이 재배열된 순서
        이는 풀이 단계인 LUsolve()에 전달
    LUsolve() : 전방대입과 후방대입에 앞서 같은 순서로 상수벡터의 요소를 재배열
        마지막에 행피봇팅을 이용해 LU 분해를 한 뒤, 역행렬을 구하는 루틴 추가
        
    a, seq = LUdecompPiv(a, tol=1.0e-9)
        행피봇팅을 이용해 행렬 [a]의 LU 분해
        반환된 행렬 [a]는 상삼각에 [U], 하삼각에 [L]을 담음
        [L][U]는 [a]의 행단위 순열
        이 순열은 벡터 {seq}에 저장
    x = LUsolve(a, b, seq)
        [L][U]{x} = {b}의 풀이
        여기서 행렬 [a]와 벡터 {seq}는 LUdecompPiv()의 반환값

'''
import sys
import numpy as np
import SwapRC

def LUdecompPiv(a, tol=1.0e-9):
    n = len(a)
    seq = np.array(range(n))
    
    # 크기 인수 설정
    s = np.zeros((n))
    for i in range(n):
        s[i] = max(abs(a[i, :]))
        
    for k in range(0, n-1):
        # 필요시 행교환
        p = np.argmax(np.abs(a[k:n, k]) / s[k:n]) + k
        if abs(a[p, k]) < tol:
            print('특이행렬임')
            sys.exit()
        if p != k:
            SwapRC.swapRows(s, k, p)
            SwapRC.swapRows(a, k, p)
            SwapRC.swapRows(seq, k, p)
            
        # 소거
        for i in range(k+1, n):
            if a[i, k] != 0.0:
                lam = a[i, k] / a[k, k]
                a[i, k+1:n] = a[i, k+1:n] - lam*a[k, k+1:n]
                a[i, k] = lam
    return a, seq

def LUsolve(a, b, seq):
    n = len(a)
    
    # 상수벡터를 재정렬하여 {x}에 저장
    x = b.copy()
    for i in range(n):
        x[i] = b[seq[i]]
    
    # 풀이
    for k in range(1, n):
        x[k] = x[k] - np.dot(a[k, 0:k], x[0:k])
    x[n-1] = x[n-1] / a[n-1, n-1]
    for k in range(n-2, -1, -1):
        x[k] = (x[k] - np.dot(a[k, k+1:n], x[k+1:n])) / a[k, k]
    return k

