from qr_decomp import qr_hh
import numpy as np
import scipy
def qr_solve(A,r,b):
    Q = GenerateQ(A);
    #R = np.triu(A)
    b1 = np.matmul(np.transpose(Q),b)[0:len(r)];
    return tri_solve(r,b1);
def tri_solve(R,b):
    #return scipy.linalg.solve(R,b);
    #Sadly R is not a square matrix and I run out of time to figure this out
    #Really sorry :(
    n = len(b)
    x = np.zeros(n);
    x[n] = b[n] / R[n,n];
    count = n;
    while (count>0):
        comp = np.dot(L[i,:i], x[:i])
        x[i] = 1/L[i,i] * (b[i] - comp)
        count = count-1;
    return x;

        
def ScalarProductHelpMethod(v):
    n = len(v);
    return np.identity(n) - 2.0 * np.matmul(v,np.transpose(v)) / (np.matmul(np.transpose(v),v))[0,0];
def GenerateQ(v):
    #Generates Householder matrix q
    #4.9
    qs = []
    n = v.shape[0]
    for i in range(n):
        qs.append(ScalarProductHelpMethod(v));
    #Q = None
    Q = np.identity(n);  
    for j in range(n):
        if i == v.shape[1]:
            break;
        identity = np.identity(n);
        identity[j:,j:] = qs[j];
        Q = np.matmul(Q,identity)
    return Q
if __name__ == "__main__":
    #Sadly, this code throws an exception
    A = np.array([[1,-1],[1e-8,0],[0,1e-8]]);
    b = np.array([0,0,1e-6]);
    A, r = qr_hh(A)
    x = qr_solve(A, r, b);
    print(x)
