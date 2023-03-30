import numpy as np
def vetor(A):
    if A['n'] == 'i':
        A = [1,0,0]
    elif A['n'] == 'j':
        A = [0,1,0]
    else:
        A = [0,0,1]
    return A
def fator_paralelo(X, Y, L):
    x = X/L
    y = Y/L
    Fij = (2/(np.pi*x*y))*(np.log((((1 + x**2)*(1 + y**2))/(1 + x**2 + y**2))**(1/2))
                           + x*(1 + y**2)**(1/2)*np.arctan(x/(1+y**2)**(1/2))
                           + y*(1 + x**2)**(1/2)*np.arctan(y/(1+x**2)**(1/2))
                           - x*np.arctan(x) - y*np.arctan(y))
    return Fij

def fator_perpendicular(X, Y, Z):
    h = Z/X
    w = Y/X

    Fij = 1/(np.pi*w) * (w*np.arctan(1/w) + h*np.arctan(1/h) - (h**2 + w**2)**(1/2)*np.arctan(1/(h**2 + w**2)**(1/2))
                      + 0.25*np.log((((1 + w**2)*(1 + h**2))/(1 + w**2 + h**2))*((w**2*(1 + w**2 + h**2))/((1 + w**2)*(w**2 + h**2)))**(w**2) *
                                    ((h**2*(1 + h**2 + w**2))/((1 + h**2)*(h**2 + w**2)))**(h**2)))
    return Fij

def fator_paralelo_dict(A, B, X, Y):
    #X = A['L1']
    #Y = A['L2']
    a = vetor(A)
    b = vetor(B)
    v1 = np.array([A['x'] * a[0], A['y'] * a[1], A['z'] * a[2]])
    v2 = np.array([B['x'] * b[0], B['y'] * b[1], B['z'] * b[2]])
    L = np.linalg.norm(v1-v2)
    x = X/L
    y = Y/L
    Fij = (2/(np.pi*x*y))*(np.log((((1 + x**2)*(1 + y**2))/(1 + x**2 + y**2))**(1/2))
                           + x*(1 + y**2)**(1/2)*np.arctan(x/(1+y**2)**(1/2))
                           + y*(1 + x**2)**(1/2)*np.arctan(y/(1+x**2)**(1/2))
                           - x*np.arctan(x) - y*np.arctan(y))
    return Fij

def fator_perpendicular_dict(A, B, X, Y, Z):
    #X = A['L2']
    #Y = A['L1']
    #Z = B['L1']
    h = Z/X
    w = Y/X
    Fij = 1/(np.pi*w) * (w*np.arctan(1/w) + h*np.arctan(1/h) - (h**2 + w**2)**(1/2)*np.arctan(1/(h**2 + w**2)**(1/2))
                      + 0.25*np.log((((1 + w**2)*(1 + h**2))/(1 + w**2 + h**2))*((w**2*(1 + w**2 + h**2))/((1 + w**2)*(w**2 + h**2)))**(w**2) *
                                    ((h**2*(1 + h**2 + w**2))/((1 + h**2)*(h**2 + w**2)))**(h**2)))
    return Fij

def fator_paralelo_dict(A, B, X, Y):
    #X = A['L1']
    #Y = A['L2']
    a = vetor(A)
    b = vetor(B)
    v1 = np.array([A['x'] * a[0], A['y'] * a[1], A['z'] * a[2]])
    v2 = np.array([B['x'] * b[0], B['y'] * b[1], B['z'] * b[2]])
    L = np.linalg.norm(v1-v2)
    x = X/L
    y = Y/L
    Fij = (2/(np.pi*x*y))*(np.log((((1 + x**2)*(1 + y**2))/(1 + x**2 + y**2))**(1/2))
                           + x*(1 + y**2)**(1/2)*np.arctan(x/(1+y**2)**(1/2))
                           + y*(1 + x**2)**(1/2)*np.arctan(y/(1+x**2)**(1/2))
                           - x*np.arctan(x) - y*np.arctan(y))
    return Fij

def fator_perpendicular_dict(A, B, X, Y, Z):
    #X = A['L2']
    #Y = A['L1']
    #Z = B['L1']
    h = Z/X
    w = Y/X
    Fij = 1/(np.pi*w) * (w*np.arctan(1/w) + h*np.arctan(1/h) - (h**2 + w**2)**(1/2)*np.arctan(1/(h**2 + w**2)**(1/2))
                      + 0.25*np.log((((1 + w**2)*(1 + h**2))/(1 + w**2 + h**2))*((w**2*(1 + w**2 + h**2))/((1 + w**2)*(w**2 + h**2)))**(w**2) *
                                    ((h**2*(1 + h**2 + w**2))/((1 + h**2)*(h**2 + w**2)))**(h**2)))
    return Fij

if __name__ == '__main__':
    print(fator_perpendicular(1,1,1))
    print(fator_paralelo(1,1,1))
    print(4*fator_perpendicular(1,1,1) + fator_paralelo(1,1,1))