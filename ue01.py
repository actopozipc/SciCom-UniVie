
from math import prod
import numpy as np
from scipy import linalg
import timeit
import matplotlib.pyplot as plt
import lu_decomp
import scipy
import math
'''
Ruft alle Methoden für das erste Beispiel auf
'''

def ErsteAufgabeAufruf():
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
def ZweiteAufgabeAufruf():
    #2
    A = np.matrix([[10*10**-17,1],[2,1]])
    b = [1,3]
    lu_solve(A, b, True)
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
    #     (a1, a2, a3)                (a1, a2, a3)         (1 , 0 , 0 )
    # A = (a4, a5, a6) , dann ist R = (0 , a5, a6) und L = (a4, 1 , 0 )
    #      (a7, a8, a9)                (0 , 0 , a9)         (a7, a8, 1 )
    #1. LU decomposition of the matrix A.
    #1. LU (often LR in German) decomposition of the matrix A.
    M, z = lu_decomp.lu(A,pv)
    R = []
    linecount = 0
    #die i-te Zeile von R ist die M Zeile mit den letzen i Elementen 0
    for i in M: #für jede Zeile in M
        line = [] #array für R-Zeile
        arr = np.squeeze(np.asarray(i)) #über np.matrix iterieren gibt wieder eine matrix, deswegen konvertierung zu array
        for l in range(linecount): #fülle restliche Zeile mit 0
            line.append(0)
        diff = len(arr)-linecount #diff = wie viel alte werte
        for j in range(diff): #von 0 bis diff element restliche elemnete hinzufügen
            line.append(arr[j+linecount])
        R.append(line[:]) #einzelne zeile zu R matrix hinzufügen
        line.clear() #platz für die nächste zeile machen
        linecount = linecount+1 #nächste zeile
    
    L = []
    linecount = 0
    #die i-te Zeile von L ist i Einträge von der i-ten M-Zeile, dann eine 1, dann eine 0
    for i in M: #für jede Zeile in M
        line = [] #array für L-Zeile
        arr = np.squeeze(np.asarray(i)) #über np.matrix iterieren gibt wieder eine matrix, deswegen konvertierung zu array
        for l in range(linecount): #fülle anfang der zeile mit werten
            line.append(arr[l])
        diff = len(arr)-linecount #differenz ist, wie viele restliche werte mit 1 und 0 befüllt werden
        for j in range(diff): #fülle rest der zeile mit 1 und dann nur 0
            if j == 0:
                line.append(1)
            else:
                line.append(0)
        L.append(line[:]) #einzelne Zeile zu L Matrix hinzufügen
        line.clear() #platz für neue zeile machen
        linecount = linecount+1 #nächste zeile
    
    #2. Determination of y in Ly = b by forward substitution      
    y = linalg.solve(L , b, lower=True, check_finite=False)
        
    P,Lower,Upper = scipy.linalg.lu(A)
    #2. Determination of y in Ly=b by forward substitution
    #3 Determination of x in Rx = y by backward substitution
    x =  linalg.solve(R , y)
    print(M)
    print(np.array(R).reshape(M.shape)) 
    print(np.array(L).reshape(M.shape))    
    # print(L)
    # print(U)
    print()
    ##return x;
def DritteAufgabeAufruf():
 #für u=-1 eV and e = -3 to 0 und T = 10,300,1000
    PlotProbabilityForThreeTemps();
    #     Now calculate explicitly p(e) for energies e = 0.0, −0.2, . . . , −2.8, −3.0 eV,
    # T = 300 K, and the same value of μ. 
    T = 300     #for T = 300K 
    start = 0
    es = []
    u = -1
    while start >=-3.0:
        es.append(start)
        start = start-0.2
    print()
'''
Fermi Dirac function
'''
k = 8.6173324 * 10**(-5);
def fermi_dirac(T, u, e):
    x = (e-u)/(k*T)
    return 1/(np.exp(x)+1)
'''
3.a
Draws fermi-dirac function for 3 given Temp values and a given energy and given u
params:
    None
Return:
    None
'''
def PlotProbabilityForThreeTemps():
    #Plot the probability p(e) = 1 − f (e) that a state is unoccupied 
    Ts = [10,300,1000] #for T = 10, 300 and 1000 K
    u = -1 #and u = −1 eV,
    #or energies e between -3 and 0 eV.
    es =  []
    start = -3.0
    while start <=0:
        es.append(start)
        start = start+0.1
    p = [] #array for probabilities
    print("energy    p=1-f       T") #print table header
    for i in Ts: #for each of the three temps
        for j in es: #for each of the energy
            p1 = 1 - fermi_dirac(i, u, j) #p = 1-f(e)
            p.append(p1) #add p to the array for probabilities
            print(j, "      ", np.round(p1,decimals=2), "       ", i) #print result with tab between so it fits the table, also rounded for the table

        plt.xlim([-3.0, 0]) #set axis range for x from -3 to 0 
        plt.plot(es,p) #plot energies with probabilities
        p.clear() #clear probability array so it can get filled with the next probabilities for the next temp
    plt.show()
if __name__ == "__main__":

    #ErsteAufgabeAufruf();
    ZweiteAufgabeAufruf();
    #DritteAufgabeAufruf();
    
    

   
    