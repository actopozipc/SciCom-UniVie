from scipy import linalg
import numpy as np
import matplotlib.pyplot as plt
def dir_veciter(A,eps):
    #random unit vector
    v_1 = np.random.rand(A.shape[0])
    v_1 = v_1 / linalg.norm(v_1)
    #array of v
    V = []
    V.append(v_1)
    #eigenvalues
    lam = []
    k = 0
    while True:
        v_dach = A@V[k]
        lambda_k = v_dach@V[k]
        lam.append(lambda_k)
        v_k1 = v_dach / linalg.norm(v_dach)
        V.append(v_k1)
        if k>1:
            if abs((lambda_k-lam[k-1])/lambda_k)<eps:
                return lam,V
        k = k+1
def inv_veciter(A,lest,eps):
    #random unit vector
    v_1 = np.random.rand(A.shape[0])
    v_1 = v_1 / linalg.norm(v_1)
    #array of v
    V = []
    V.append(v_1)
    #eigenvalues
    lam = []
    k = 0
    while True:
        v_dach = np.linalg.inv(A-(lest*np.identity(A.shape[0]))) @ V[k]
        lambda_k = lest + 1/(v_dach@V[k])
        lam.append(lambda_k)
        v_k1 = v_dach / linalg.norm(v_dach)
        V.append(v_k1)
        if k>1:
            if abs((lambda_k-lam[k-1])/lambda_k)<eps:
                return lam,V
        k = k+1
    return lam,V

def PlotEigenvalueDifference(maxv, eigenvalues):
        deltas = []
        i = []
        count = 0;
        for k in eigenvalues:
            deltas.append(k-maxv)
            count = count+1
            i.append(count)
        return i,deltas

if __name__ == "__main__":
    eps = 1e-5 #epsilon
    fig, axs = plt.subplots(2, 2)
    A = np.array([[2, 0, 0.2], [0, -2, 1], [0.2, 1, -2]]) #A
    lam, V = dir_veciter(A, eps) 
    richtig = linalg.eig(A) #richtige Eigenwerte
    foundEigenvalues = lam
    biggestFoundEigenvalue = lam[-1]
    print("Gefundener größter Eigenwert für A:", biggestFoundEigenvalue)
    print("Gefunde Eigenwerte:",foundEigenvalues)#gefundene Eigenwerte von linalg
    x,y = PlotEigenvalueDifference(biggestFoundEigenvalue, foundEigenvalues)
    axs[0,0].plot(x,y)
    lam, V = inv_veciter(A, -4, eps) #inverse methode
    foundEigenvalues = lam
    biggestFoundEigenvalue = lam[-1]
    x,y = PlotEigenvalueDifference(biggestFoundEigenvalue, foundEigenvalues)
    axs[0,1].plot(x,y)
    print("Größter gefundener Eigenwert mit Inverser Methode:",biggestFoundEigenvalue) #letzter eigenwert von der inversen methode
    B = np.array([[2, 0, 0.2], [0, 2.02, 1], [0.2, 1, -2]]) #B
    lam, V = dir_veciter(B, eps)
    richtig = linalg.eig(B) #richtige eigenwerte von linalg
    foundEigenvalues = lam
    biggestFoundEigenvalue = lam[-1]
    x,y = PlotEigenvalueDifference(biggestFoundEigenvalue, foundEigenvalues)
    axs[1,0].plot(x,y)
    print("Größter gefundener Eigenwert für B:",biggestFoundEigenvalue)
    print("Gefundene Eigenwerte:", foundEigenvalues)
    lam, V = inv_veciter(B, 4, eps)
    print("Größter gefundener Eigenwert mit Inverser Methode:",biggestFoundEigenvalue)
    foundEigenvalues = lam
    biggestFoundEigenvalue = lam[-1]
    x,y = PlotEigenvalueDifference(biggestFoundEigenvalue, foundEigenvalues)
    axs[1,1].plot(x,y)
    fig.tight_layout()
    plt.show()


