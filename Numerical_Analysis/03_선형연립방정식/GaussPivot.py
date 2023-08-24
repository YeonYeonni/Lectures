#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## GaussPivot
'''
    GaussPivot() : 행피봇팅과 Gauss 소거법 수행
    
    x = GaussPivot(a, b, tol=1.0e-12).
    선형연립방정식 [a]{x} = {b}을
    행피봇팅을 하고 Gauss 소거법으로 푼다.
'''
import numpy as np
import SwapRC

def GaussPivot(a, b, tol=1.0e-12):
    n = len(b)
    
    # 크기 척도 설정
    s = np.zeros(n)
    for i in range(n):
        s[i] = max(np.abs(a[i, :]))
    
    for k in range(0, n-1):
        # 필요시 행교환
        p = np.argmax(np.abs(a[k:n, k]) / s[k:n]) + k
        if abs(a[p, k]) < tol:
            print('특이행렬임')
            sys.exit()
        if p != k:
            SwapRC.swapRows(b, k, p)
            SwapRC.swapRows(s, k, p)
            SwapRC.swapRows(a, k, p)
        
        # 소거
        for i in range(k+1, n):
            if a[i, k] != 0.0:
                lam = a[i, k] / a[k, k]
                a[i, k+1:n] = a[i, k+1:n] - lam*a[k, k+1:n]
                b[i] = b[i] - lam*b[k]
    
    if abs(a[n-1, n-1]) < tol:
        print('특이행렬임')
        sys.exit()
    
    # 후방대입
    b[n-1] = b[n-1] / a[n-1, n-1]
    for k in range(n-2, -1, -1):
        b[k] = (b[k] - np.dot(a[k, k+1:n], b[k+1:n])) / a[k, k]
    return b

