from scipy.sparse import csc_matrix
import scs
import math
import numpy as np
from itertools import accumulate
from cvxopt import spmatrix, matrix
a=math.sqrt(2)
def G_half(n):
    lst=[]
    for i in range(1,n):
        lst=lst+list(range(i*n,i+n*i))
    return lst
def G_rem(n):
    lst=[]
    for i in range(n):
        lst=lst+list(range(i*n+i,n-i+n*i+i))
    return lst

#convert the sdpa form into the form that can be used by cvxopt linear cone solver
def conelp_form(file):
    file = open(file)
    mdim=int(file.readline())# the number of matrix
    nblocks=int(file.readline())# the number of blocks in every matrix
    block_size=[abs(int(n)) for n in file.readline().split(' ') if n!='\n' and n!=''] # the size of each block
    obj=np.array([float(n) for n in file.readline().split(' ') if n!='\n' and n!=''])# objective
    lst=[0]+list(accumulate([x for x in block_size[0:-1]])) #index of blocks
    G_col=mdim
    G_row=(lst[-1]+block_size[-1])**2
    matrix_size=lst[-1]+block_size[-1]
    h=spmatrix([],[],[], (G_row,1))
    G=spmatrix([],[],[], (G_row, G_col))

    line=file.readline().split(' ') #read the matrix entries
    line =[ele for ele in line if ele!=''] 
    while (line[0]=='0'):
        add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
        h[(add+int(line[2]))*matrix_size+add+int(line[3])]=float(line[4])
        line=file.readline().split(' ') 
        line =[ele for ele in line if ele!=''] 

    add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
    G[(add+int(line[2]))*matrix_size+add+int(line[3]),int(line[0])-1]=float(line[4])

    for line in file:
        line=line.split(' ')
        line =[ele for ele in line if ele!=''] 
        add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
        G[(add+int(line[2]))*matrix_size+add+int(line[3]),int(line[0])-1]=float(line[4])
    G=-G
    h=-matrix(h)
    c=matrix(obj)
    dims = {'l': 0, 'q': [], 's': [matrix_size]}
    return c,G,h,dims


#convert the sdpa form into the form that can be used by scs (splitting conic solver)
def scs_form(file): 
    file=open(file)
    mdim=int(file.readline())# the number of matrix
    nblocks=int(file.readline())# the number of blocks in every matrix
    block_size=[abs(int(n)) for n in file.readline().split(' ')  if n!='\n' and n!=''] # the size of each block
    obj=np.array([float(n) for n in file.readline().split(' ') if n!='\n' and n!=''])# objective
    lst=[0]+list(accumulate([x for x in block_size[0:-1]])) #index of blocks
    G_col=mdim
    G_row=(lst[-1]+block_size[-1])**2
    matrix_size=lst[-1]+block_size[-1]
    h=np.zeros(G_row,dtype=np.double)
    G=csc_matrix((G_row, G_col), dtype=np.double)
    line=file.readline().split(' ') #read the matrix entries
    line =[ele for ele in line if ele!=''] 
    while (line[0]=='0'):
        add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
        x=add+int(line[2])
        y=add+int(line[3])
        if x!=y:
            h[x*matrix_size+y]=a*float(line[4]) 
        else:
            h[x*matrix_size+y]=float(line[4]) 
        line=file.readline().split(' ') 
        line =[ele for ele in line if ele!=''] 

    add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
    G[(add+int(line[2]))*matrix_size+add+int(line[3]),int(line[0])-1]=float(line[4])

    for line in file:
        line=line.split(' ')
        line =[ele for ele in line if ele!=''] 
        add=lst[int(line[1])-1]-1 #add this to column and row to  get the position in the matrix 
        x=add+int(line[2])
        y=add+int(line[3])
        if x!=y:
            G[x*matrix_size+y,int(line[0])-1]=a*float(line[4])
        else:
            G[x*matrix_size+y,int(line[0])-1]=float(line[4])

    cone= { 's':[matrix_size]}
    rem=G_rem(matrix_size)
    data={'A':-G[rem,:],'b':-h[rem],'c':obj}  
    return data, cone

