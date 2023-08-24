#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## SwapRC 모듈
'''
    swapRows() : 행렬이나 벡터 v의 i행과 j행을 교환
    swapCols() : 행렬 v의 i열과 j열을 교환
    
    vc = swapRows(vm, i, j).
        벡터 또는 행렬[vm]의 I행과 j행을 교환
    vc = swapCols(vm, i, j).
        행렬[vm]의 I열과 j열을 교환
'''
def swapRows(vm, i, j):
    vc = vm.copy()
    if len(vc.shape) == 1:
        t = vc[i]
        vc[i] = vc[i]
        vc[j] = t
    else:
        vc[[i, j], :] = vc[[j, i], :]
    return vc

def swapCols(vm, i, j):
    vc = vm.copy()
    vc[:, [i, j]] = vc[:, [j, i]]
    return vc

