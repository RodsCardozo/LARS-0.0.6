"""
Determinação do contato, tempo e posição do satélite na visada de um alvo ou ground station.
A ideia é calcular a partir dos dados obtidos do dataframe gerado a partir da integração do movimento do
satélite, advindo da função "propagador_orbital_2.py"



"""
import pandas
import pandas as pd
import os, sys

def resource_path(relative_path):

    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)
df = pd.read_csv(resource_path("data\posicao_tempo2.csv"), sep=',', engine='python', on_bad_lines='skip')

def tempo_visada(dataframe, coordenadas_alvo: list, angulo_elevacao: float, inclinacao: float):
    """
    :param dataframe: pandas.Dataframe
    :param coordenadas_alvo: list(latitude, longitude)
    :param angulo_elevacao: float(elevation angle minimun)
    :param inclinacao: float(Orbit inclination)
    :return: pandas.Dataframe with all communication time between a target and a satellite
    """
    import numpy as np
    import pandas as pd
    import Periodo_Orbital as PO
    from tqdm import tqdm
    R_E = 6371.00 # Raio da Terra em km
    e = np.radians(angulo_elevacao)
    i = np.radians(inclinacao)
    lat_gs = np.radians(coordenadas_alvo[0]) # latitude ground station
    #print(lat_gs)
    long_gs = np.radians(coordenadas_alvo[1])   # longitude ground station
    #print(long_gs)
    # raio angular da Terra
    df = dataframe
    lat_pole = np.pi / 2 - i
    posicoes = pd.DataFrame()
    for i in tqdm(range(0, len(df)), colour='#00a6ff'):
        H = np.linalg.norm([df.iloc[i, 0], df.iloc[i, 1], df.iloc[i, 2]])
        #print(H)
        rho = np.arcsin(R_E / (H))
        #print(f'Raio angular da Terra: {np.degrees(rho)} deg')
        # Periodo orbital
        P = PO.periodo_orbital(H)
        #print(f'Periodo Orbital: {P / 60} min')
        long_pole = df.iloc[i, df.columns.get_loc('Raan')] - np.pi / 2
        # Max angulo nadir
        n_max = np.arcsin(np.sin(rho) * np.cos(e))
        #print(n_max)
        #print(f'Maximo angulo nadir: {np.degrees(n_max)} deg')
        posicao_gs = [R_E * np.cos(lat_gs) * np.cos(long_gs), R_E * np.cos(lat_gs) * np.sin(long_gs), R_E * np.sin(lat_gs)]
        #print(posicao_gs)
        lamb = np.dot(posicao_gs, [df.iloc[i, 0],
                                   df.iloc[i, 1],
                                   df.iloc[i, 2]]) / (np.linalg.norm([df.iloc[i, 0], df.iloc[i, 1], df.iloc[i, 2]])
                                                      * np.linalg.norm(posicao_gs))
        posicoes = []
        vetor_posicoes = []
        # Maximo angulo central da Terra
        lamb_max = np.pi/2 - (e) - (n_max)
        #print(lamb_max)
        #print(f'Maximo angulo central: {np.degrees(lamb_max)} deg')
        if np.linalg.norm(np.arccos(lamb)) < (lamb_max):

            nova_lina = {'X_ECI': df.iloc[i, 0],
                         'Y_ECI': df.iloc[i, 1],
                         'Z_ECI': df.iloc[i, 2],
                         'longitude': np.degrees(np.arctan2(df.iloc[i, 1], df.iloc[i,0])) ,
                         'latitude': np.degrees(np.arcsin(df.iloc[i, 2] / H)),
                         'end': [0]}
            df_nova_linha = pd.DataFrame(nova_lina)
            posicoes = pd.concat([posicoes, df_nova_linha], ignore_index=True)


        '''# Distancia maxima ao alvo
        D_max = R_E*(np.sin(np.radians(lamb_max))/np.sin(n_max))
        #print(f'Distancia maxima ao alvo: {D_max} km')

        # Minimo angulo central
        lamb_min = np.pi/2 - np.arcsin(np.sin((lat_pole))*np.sin((lat_gs)) +
                             np.cos((lat_pole))*np.cos((lat_gs))*np.cos((long_pole - long_gs)))
        #print(f'Minimo angulo central: {np.degrees(lamb_min)} deg')

        # Minima angulo nadir eta_min
        n_min = np.arctan((np.sin(rho) * np.sin(lamb_min)) / (1 - np.sin(rho) * np.cos(lamb_min)))
        #print(f'Minimo angulo nadir: {np.degrees(n_min)} deg')

        # Máxima elevação epsilon
        e_max = 90 - np.degrees(lamb_min) - np.degrees(n_min)
        #print(f'Máximo angulo de elevação: {e_max} deg')

        # Distância mínima D_min
        D_min = R_E*(np.sin(lamb_min) / np.sin(n_min))
        #print(f'Distância mínima D: {D_min} km')

        # Máxima taxa angular theta_p_max
        theta_p_max = (360 * (R_E + H))/(P*D_min)
        #print(f'Máxima taxa angular theta_p_max: {theta_p_max*60} deg/min')

        # faixa azimutal
        Delta_phi = 2 * np.arccos((np.tan(lamb_min) / np.tan(lamb_max)))
        #print(f'Faixa azimutal: {np.degrees(Delta_phi)} deg')

        # Azimute da passagem central phi_center
        phi_center = np.arccos(((np.sin(lamb_min) * np.sin(lat_gs)) - np.sin(lat_pole)) / (np.cos(lamb_min) * np.cos(lat_gs)))
        #print(f'Azimute do centro da passagem: {np.degrees(phi_center)} deg')

        # Tempo de visada T
        T = (P / 180) * np.arccos(np.cos(lamb_max) / np.cos(lamb_min))
        #print(f'Tempo de visada: {T} min')'''

    return posicoes

from Plots import plot_groundtrack_3D as plt3d

a = tempo_visada(dataframe=df, coordenadas_alvo=[0.0, 0.0], angulo_elevacao=10.0, inclinacao= 52.0)

import os.path
a.to_csv(os.path.join('./dados_comunicacao/', 'posicoes_comunicacao.csv'), sep=',', index = False, columns=a.columns[0:])
import plotly.express as px
'''fig = px.line_3d(a, x='X_ECI', y= 'Y_ECI', z='Z_ECI')
fig.show()'''
print(a)
print(df)
plt3d(a)