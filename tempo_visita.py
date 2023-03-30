"""
Determinacao do tempo de visita a uma ground station

analise da ground station de Floripa

coordenadas floripa: lat(-27.6), long(-48.5).




"""

import pandas as pd
import os, sys

def resource_path(relative_path):

    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)
df = pd.read_csv(resource_path("data\posicao_tempo.csv"), sep='=', engine='python', on_bad_lines='skip')

import numpy as np

""" Refazendo o exmplo do livro """
H = 1000.0 # altitude do satelite
R_E = 6371.00 # Raio da Terra
a = H + R_E
i = 50.0 # inclinacao
e = np.radians(5.0) # elevacao minima
lat_gs = np.radians(33.5) # latitude ground station
long_gs = np.radians(247.6)   # longitude ground station
lat_pole = np.radians(40.0)
long_pole = np.radians(230)
import Periodo_Orbital as PO

# raio angular da Terra
rho = np.arcsin(R_E/(R_E + H))
print(f'Raio angular da Terra: {np.degrees(rho)} deg')

# Periodo orbital
P = PO.periodo_orbital(a)
print(f'Periodo Orbital: {P/60} min')

# Max angulo nadir
n_max = np.arcsin(np.sin(rho) * np.cos(e))
print(f'Maximo angulo nadir: {np.degrees(n_max)} deg')

# Maximo angulo central da Terra
lamb_max = 90 - np.degrees(e) - np.degrees(n_max)
print(f'Maximo angulo central: {(lamb_max)} deg')

# Distancia maxima ao alvo
D_max = R_E*(np.sin(np.radians(lamb_max))/np.sin(n_max))
print(f'Distancia maxima ao alvo: {D_max} km')

# Minimo angulo central
lamb_min = np.pi/2 - np.arcsin(np.sin((lat_pole))*np.sin((lat_gs)) +
                     np.cos((lat_pole))*np.cos((lat_gs))*np.cos((long_pole - long_gs)))
print(f'Minimo angulo central: {np.degrees(lamb_min)} deg')

# Minima angulo nadir eta_min
n_min = np.arctan((np.sin(rho) * np.sin(lamb_min)) / (1 - np.sin(rho) * np.cos(lamb_min)))
print(f'Minimo angulo nadir: {np.degrees(n_min)} deg')

# Máxima elevação epsilon
e_max = 90 - np.degrees(lamb_min) - np.degrees(n_min)
print(f'Máximo angulo de elevação: {e_max} deg')

# Distância mínima D_min
D_min = R_E*(np.sin(lamb_min) / np.sin(n_min))
print(f'Distância mínima D: {D_min} km')

# Máxima taxa angular theta_p_max
theta_p_max = (360 * (R_E + H))/(P*D_min)
print(f'Máxima taxa angular theta_p_max: {theta_p_max*60} deg/min')

# faixa azimutal
Delta_phi = 2 * np.arccos((np.tan(lamb_min) / np.tan(lamb_max)))
print(f'Faixa azimutal: {np.degrees(Delta_phi)} deg')

# Azimute da passagem central phi_center
phi_center = np.arccos(((np.sin(lamb_min) * np.sin(lat_gs)) - np.sin(lat_pole)) / (np.cos(lamb_min) * np.cos(lat_gs)))
print(f'Azimute do centro da passagem: {np.degrees(phi_center)} deg')

# Tempo de visada T
T = (P / 180) * np.arccos(np.cos(lamb_max) / np.cos(lamb_min))
print(f'Tempo de visada: {T} min')









