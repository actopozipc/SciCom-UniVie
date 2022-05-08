import numpy as np
from numpy import ndarray
from scipy import linalg
import matplotlib.pyplot as plt




def pol_fit(x, y, n):
    M = np.vander(x, N=n, increasing=True);
    a = linalg.solve(np.dot(np.transpose(M), M), np.dot(np.transpose(M), y));
    return a;


def pol(c, x):
    y = 0;
    for i in range(c.size):
        y += c[:][i] * x ** i;
    return y;


if __name__ == "__main__":
    #NOTE: cd into the directory before using python file
    Z = np.genfromtxt('ex05.dat');
    x = Z[:, 0];
    y = Z[:, 1];
    plt.plot(Z[:, 0], Z[:, 1]);
    x1 = np.linspace(0, 5);
    for i in [2, 3, 6, 9, 20]:
        plt.plot(x1, pol(pol_fit(x, y, i), x1));
    plt.legend();
    plt.ylim((3, -3));
    plt.show();
    for i in [2, 3, 6, 9, 20]:
        print("Predictet value of y for x=5 for the different fits:",{pol(pol_fit(x, y, i), 5)})


