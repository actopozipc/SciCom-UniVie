import numpy as np
import scipy as sp
from numpy import ndarray
from scipy import linalg

if __name__ == "__main__":
    a = [];
    b = [];
    for i in range(1,7):
        a.append(float(i/10));
        b.append(float(i));
    a = np.array(a);
    b = np.array(b);
    temp = [0,0,0,0.001,0,0];
    b_s = b + np.array(temp);
    V = np.vander(a);
    sol1 = sp.linalg.solve(V, b);
    sol2 = sp.linalg.solve(V, b_s);
    #errors
    errorSolutions = sp.linalg.norm(sol2 - sol1, ord=2) / sp.linalg.norm(sol1, ord=2);
    errorVectors = sp.linalg.norm(b_s - b, ord=2) / sp.linalg.norm(b, ord=2);
    print("Relative error of the solutions:", errorSolutions);
    print("Relative error of the vectors:", errorVectors);
    #lower bound
    lower_bound = errorSolutions / errorVectors;
    print("Lower bound:",lower_bound);
    realCondition = sp.linalg.norm(V, ord=2) * sp.linalg.norm(np.linalg.inv(V), ord=2);
    print("The difference is", realCondition-lower_bound);