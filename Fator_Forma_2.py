#Fator_Forma_2.py>
import numpy as np


def fator_forma_montecarlo(altitude):
    import random
    import math

    Quant_Teste=1000000
    d= altitude*1000 #Altitude
    Raio_Terra=6371000 #Raio da terra

    Counter = 0
    for i in range(0, Quant_Teste):

        R_Teta = random.uniform(0, 1)  # Randomico para teta
        R_Fi = random.uniform(-1, 1)  # Randomico para fi

        Teta = math.asin(math.sqrt(R_Teta))  # Calculo de teta
        Fi = math.fabs(2 * math.pi * R_Fi)

        Teta_max = math.asin(Raio_Terra / d)

        if (Teta <= Teta_max and Teta >= 0):
            Counter += 1

    FS_T = Counter / Quant_Teste  # Fator de forma do satelite para a terra

    return FS_T

#Fator de Forma Sat - Terra - Plano perpendicular do satelite
import random
import math

def fator_forma_perpendicular(altitude):
    Quant_Teste=1000000
    d= altitude #Altitude
    Raio_Terra= 6371.000 #Raio da terra
    Orbita= d - Raio_Terra

    Counter=0
    for i in range(0,Quant_Teste):

        R_Teta=random.uniform(0,1) # Randomico para teta
        R_Fi=random.uniform(0,1) # Randomico para fi

        Teta=math.asin(math.sqrt(R_Teta)) # Calculo de teta
        Fi=2*math.pi*R_Fi

    #    print(math.degrees(Fi))

        Alpha_max=math.asin(Raio_Terra/d)
        Psi_max=math.pi-math.pi/2-Alpha_max
        h_max=Raio_Terra*math.sin(Psi_max)
        Dist_Z=Raio_Terra*math.cos(Psi_max)
        Alpha=math.pi/2-Teta

        d1=d-Dist_Z
        X1=d1*math.tan(Fi)
        S1=d1/math.cos(Fi)
        Y1=S1*math.tan(Alpha)

        Raio_Circulo=math.sqrt(Y1**2+X1**2)

        if (Raio_Circulo<=h_max):
            if (Fi<=math.pi and Fi>=0):
                Counter+=1

    FS_T=Counter/Quant_Teste # Fator de forma do satelite para a terra
    return FS_T
def fator_forma_classico(vetor_posicao, angulo):

    import numpy as np

    raio_terra = 6371
    k = raio_terra/(vetor_posicao)
    d = angulo # angulo da face do satelite

    if 0<=d<=np.arccos(k):
        FS = k**2*np.cos(d)

    if np.arccos(k) < d < np.pi - np.arccos(k):
        FS = k**2*np.cos(d) + (1/np.pi)*(np.pi/2 - ((1 - k**2)*(k**2 - np.cos(d)**2))**(0.5) - np.arcsin((1 - k**2)**0.5/(np.sin(d)))
         - k**2*np.cos(d)*np.arccos((1/k**2 - 1)**0.5*np.tan(d)**-1))
    if np.pi - np.arccos(k) <= d <= np.pi:
        FS = 0

    return FS

if __name__ == '__main__':
    print(f'fator de forma frontal: {fator_forma_montecarlo(7000.000)}')
    print(f'fator de forma modelo analitico: {fator_forma_classico(7000.000, np.pi/2)}')
    print(f'fator de forma perpendicular: {fator_forma_perpendicular(7000.00)}')