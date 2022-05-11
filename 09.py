coeff = [10,-0.98,0.34, -0.027, 0.001, -0.000016]

def phi1(t):
    y = t
    for i in range(len(coeff)):
        y += coeff[i] * t ** i
    return y
def fprime(t):
    y = 0
    for i in range(len(coeff)):
        y += coeff[i]*i*t**(i-1)
    return y
def f(t):
    y = 0
    for i in range(len(coeff)):
        y += coeff[i]*t**i
    return y
def phi2(t):
    return t - f(t)/fprime(t)
def g(t):
    return f(t+f(t))/f(t)-1
def phi3(t):
    return t-f(t) / g(t)
def fixiter(phi,x0,kmax, epsilon):
    x = []
    counter = 0
    while True:
       x.append(x0) 
       x0 = phi(x0)
       if counter>3:
            diff = (x[counter]-x[counter-1]) / x[counter-1]
            if abs(diff)<epsilon:
                return x
            elif counter == kmax:
                break;
       counter = counter + 1;
    return x;

    pass
if __name__ == "__main__":
    x0 = 27
    kmax = 5
    epsilon = 0.005
    res1 = fixiter(phi1, x0, kmax, 0.005)
    res2 = fixiter(phi2, x0, kmax, epsilon)
    res3 = fixiter(phi3, x0, kmax, epsilon)
    print("Iter.    Phi1        Phi2          Phi3")
    for i in range (5):
        print(i, "      ", res1[i], "       ", res2[i], "       ", res3[i])