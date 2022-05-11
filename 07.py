import numpy as np
from numpy import ndarray
from scipy import linalg
from qr_decomp import *
import matplotlib.pyplot as plt
def create_q(M):
    n = M.shape[0]
    Q = np.eye(n)
    for i in range(n):
        if i != M.shape[1]:
            one = np.eye(n)
            q = construct_q(M[i:n, i])
            one[i:, i:] = q
        else:
            break
        Q = Q @ one
    return Q
def construct_q(v):
    h = v.copy()
    h = h.reshape((len(v), 1))
    return np.eye(len(v)) - 2 * h @ h.T / (h.T @ h)[0, 0]
def QR_eigsolve(A,eps,kmax):
    a_k = []
    a_k.append(A)
    for i in range(kmax):
        M, r = qr_hh(A)
        Q = create_q(M)
        a_h=Q.T@a_k[i]@Q
        a_k.append(a_h)
        
        if i > 2:
            diags = np.diag(np.diag(a_h))
            diagsum = np.cumsum(diags)[-1]
            offdiagsum = np.cumsum(a_h - diags)[-1]
            if offdiagsum/diagsum < eps or i > kmax:
                break
    return a_k

if __name__ == "__main__":
    n=16
    A=np.random.rand(n,n); A=np.tril(A)@np.tril(A).T
    ev = linalg.eig(A)[0]
    
    alist = QR_eigsolve(A, 0.005, 50)
    for i in alist:
        Am=np.min(abs(i))
        AM=np.max(abs(i))
        plt.imshow(i,cmap='ocean_r',vmin=Am,vmax=AM)
    
    plt.colorbar()
    
    print("Relative error of the eigenvalues:")
    relError = (linalg.norm(alist[-1].diagonal() - ev))/linalg.norm(ev)
    print(relError)