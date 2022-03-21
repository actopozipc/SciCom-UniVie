
from math import prod
import numpy as np
from scipy import linalg
import timeit
import matplotlib.pyplot as plt
import lu_decomp
'''
1.a
Lx = b
forward substitution method
params: 
    L: np.array
    b: np.array

returns:
    x: np.array
'''
def tri_solve(L,b):
    #x1 = b1
    #x_i = b_i - sum(L_iy * y_j)
    n = len(b)
    x = np.zeros(n)
    x[0] = b[0]/L[0,0];
    comp = 0;
    for i in range(1,n):
        
        # for k in range(0,i):
        #     index = L[i,k]
        #     preSolution = x[k]
        comp = np.dot(L[i,:i], x[:i])
        x[i] = 1/L[i,i] * (b[i] - comp)
    return x;
'''
1.b
Measure and visualize time of my function and linalg.solve for different values
params:
    None
Returns:
    None
'''
def measure_time():
    #Do it for 5,10,20,50,100,200,400 n
    ns = [5,10,20,50,100,200,400]
    average_trisolve_times = [] #my times
    average_linalgsolve_times = [] #linalg times
    for n in ns:
        trisolve_times = []
        linalgsolve_times = []
        for i in range(10):

            L,b = generate_random_lower_triangle_matrix(n)
            t1 = timeit.default_timer()
            tri_solve(L, b)
            trisolve_times.append(timeit.default_timer()-t1)
            t2 = timeit.default_timer()
            linalg.solve(L, b, lower=True, check_finite=False)
            linalgsolve_times.append(timeit.default_timer() - t2)
        #Save time in array
        average_trisolve_times.append(sum(trisolve_times) / len(trisolve_times))
        average_linalgsolve_times.append(sum(linalgsolve_times) / len(linalgsolve_times))
    
    #Plot graph
    #linalg
    x = ns
    y = average_linalgsolve_times
    fig, ax = plt.subplots()
    ax.plot(x, y, ".-", linewidth=2.0, color = 'blue', label='linalg.solve' )
    #tri_solve
    y = average_trisolve_times
    ax.plot(x, y, ".-", linewidth=2.0, color = 'red', label='tri_solve')
    #legende
    plt.legend()
    #log scale
    plt.yscale("log")
    plt.xscale("log")

    plt.show()
    #Think: have you implemented tri_solve efficiently?
    #if I wanted efficiency, I wouldn't take python
    #Beside that, its hard to compete with a library written in a language that is compiled and runs bare metal
    #and scipy is not unlikely multithreaded
    #I still assume my code is okay, since the runtime increases approximately linearly
    pass
'''
2
Solves system of equations Ax=b by gaussian elimination
params:
    a: n*n matrix, np.array
    b: n-d vector, np.array
    pv: pivoting is used or not, bool
returns:
    x: n-d vector, np.array
'''
def lu_solve(A,b,pv):
    #1. LU decomposition of the matrix A.
    M, z = lu_decomp.lu(A,pv)
    #2. Determination of y in Ly=b by forward substitution
    #3 Determination of x in Rx = y by backward substitution
    print(M)
    print()
    ##return x;
'''
Helper method for 1.
Generates a random lower triangle matrix with n dimensions
params:
    n: integer
returns:
    tuple: 
        triangle matrix with n dimensions: np.array
        vector with n elements: np.array
'''
def generate_random_lower_triangle_matrix(n):
    L = np.tril(np.random.randint(1,10,size=(n,n)))
    #print(L)
    b = np.random.randint(1,10,size=(n))
    #print(b)
    return [L,b]

k = 8.6173324 * 10**(-5);
def fermi_dirac(T, u, e):

    return 1/(math.exp((e-u)/k*T)+1)
    pass
if __name__ == "__main__":
    #2
    A = np.matrix([[3,4,2],[1,5,7],[2,3,1]])#np.matrix([[10*10**-17,1],[2,1]])
    b = [1,3,3]
    lu_solve(A, b, False)
    #1.a
    L = np.matrix([[3,0,0],[1,5,0],[2,3,1]])
    b = [3,2,1];
    x = tri_solve(L,b) #(1.0,0.2,-1.6)
    y = linalg.solve(L, b) #(1.0,0.2,-1.6)
    #Generate random matrices to test
    L, b = generate_random_lower_triangle_matrix(3)
    x = tri_solve(L,b)
    solution = linalg.solve(L,b)
    if (x == solution).all(): #if result from my function and linalg.solve are the same
        print("Works! :)")
    #1.b
    measure_time()

    #3
    