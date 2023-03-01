"""
    Universidade Federal de Santa Catarina
    Laboratory of Applications and Research in Space - LARS
    Orbital Mechanics Division

    Título do Algoritmo = Codigo principal de propagacao orbital e analise termica
    Autor= Rodrigo S. Cardozo
    Versão= 0.0.1
    Data= 24/10/2022

"""
import matplotlib.pyplot as plt

from Propagador_Orbital import propagador_orbital
from Calor_Incidente import calor_incidente
from Periodo_Orbital import periodo_orbital
from datetime import datetime
import numpy as np
import pandas as pd
import Plots
import os, sys

def resource_path(relative_path):

    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

df = pd.read_csv(resource_path("Dados\dados_entrada.csv"), sep='=', engine='python', on_bad_lines='skip')
SMA = float(df.iloc[0, 0])
a = float(df.iloc[1, 0])  # ecentricidade da orbita
if a < 0.002:
    ecc = 0.002
else:
    ecc = a
Raan = float(df.iloc[2, 0])  # ascencao direita do nodo ascendente
arg_per = (float(df.iloc[3, 0]))  # argumento do perigeu
true_anomaly = (float(df.iloc[4, 0]))  # anomalia verdadeira
b = (float(df.iloc[5, 0]))  # inclinacao
if b < 1.0:
    inc = 1.0
else:
    inc = b
mu = float(df.iloc[6, 0])  # constante gravitacional da terra
J2 = float(df.iloc[7, 0])  # zona harmonica 2
Raio_terra = float(df.iloc[8, 0])  # raio da terra
num_orbita = float(df.iloc[9, 0])  # numero de obitas
rp = SMA * (1 - ecc)
T_orbita = periodo_orbital(SMA)

PSIP = float(df.iloc[11, 0])
TETAP = float(df.iloc[12, 0])
PHIP = (2*np.pi)/T_orbita #
psi0 = float(df.iloc[14, 0])
teta0 = float(df.iloc[15, 0])
phi0 = float(df.iloc[16, 0])
input_string = df.iloc[18, 0]
data = datetime.strptime(input_string, " %m/%d/%Y %H:%M:%S")
delt = float(df.iloc[19, 0])
n = int(df.iloc[20, 0])
massa = float(df.iloc[22, 0])  # massa do cubesat
largura = float(df.iloc[23, 0])  # comprimento do sat
comprimento = float(df.iloc[24, 0])  # largura do sat
altura = float(df.iloc[25, 0])  # altura do sat

# Intensidade radiante do sol e terra e valores de emissividade

Is = float(df.iloc[27, 0])    # radiacao solar
Ir = float(df.iloc[28, 0])    # radiacao IR Terra
e = float(df.iloc[29, 0])     # emissividade Terra
ai = float(df.iloc[30, 0])    # absortividade do satelite
gama = float(df.iloc[31, 0])  # refletividade da Terra

Propagacao_orbital = propagador_orbital(data, SMA, ecc, Raan, arg_per, true_anomaly, inc, num_orbita, delt, psi0, teta0,
                                        phi0, PSIP, TETAP, PHIP, massa, largura, comprimento, altura)


calor_total = calor_incidente(Propagacao_orbital, Is, Ir, e, ai, gama, num_orbita)
size = SMA*1.1
#plot_animacao = Plots.plot_animacao_orbita(Propagacao_orbital, size, num_orbita)
#plot_groundtrack3d = Plots.plot_groundtrack_3D(Propagacao_orbital)
#plot_groundtrack2d = Plots.plot_groundtrack_2D(Propagacao_orbital)
plot_calor_sol = Plots.calor_solar(calor_total)
plot_calor_Terra = Plots.calor_IR_Terra(calor_total)
plot_calor_albedo = Plots.calor_albedo(calor_total)
plot_calor_total = Plots.calor_total(calor_total)


import os

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


# Example
createFolder('./data/')
# Creates a folder in the current directory called data