#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
from itertools import accumulate
from cvxopt import spmatrix, matrix

def conelp_form(file):
    file = open(file)
    mdim=int(file.readline())# the number of matrix
    nblocks=int(file.readline())# the number of blocks in every matrix
    block_size=[abs(int(n)) for n in file.readline().split(' ') if n!='\n'] # the size of each block
    obj=np.array([float(n) for n in file.readline().split(' ') if n!='\n'])# objective
    lst=[0]+list(accumulate([x for x in block_size[0:-1]])) #index of blocks
    G_col=mdim
    G_row=(lst[-1]+block_size[-1])**2
    matrix_size=lst[-1]+block_size[-1]
    h=spmatrix([],[],[], (G_row,1))
    G=spmatrix([],[],[], (G_row, G_col))

    line=file.readline().split(' ') #read the matrix entries
    while (line[0]=='0'):
        add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
        h[(add+int(line[2]))*matrix_size+add+int(line[3])]=float(line[4])
        line=file.readline().split(' ') 

    add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
    G[(add+int(line[2]))*matrix_size+add+int(line[3]),int(line[0])-1]=float(line[4])

    for line in file:
        line=line.split(' ')
        add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
        G[(add+int(line[2]))*matrix_size+add+int(line[3]),int(line[0])-1]=float(line[4])
    G=-G
    h=-matrix(h)
    c=matrix(obj)
    dims = {'l': 0, 'q': [], 's': [matrix_size]}
    return c,G,h,dims



