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
df = pd.read_csv(resource_path("data\posicao_tempo.csv"), sep='=', engine='python', on_bad_lines='skip')


def tempo_visada(dataframe, coordenadas_alvo: list, angulo_elevacao: float):
    """
    :param dataframe: pandas.Dataframe
    :param coordenadas_alvo: list(latitude, longitude)
    :param angulo_elevacao: float(elevation angle minimun)
    :return: pandas.Dataframe with all communication time between a target and a satellite
    """
    import numpy as np
    import pandas as pd
    import Periodo_Orbital as PO
    R_E = 6371.00 # Raio da Terra em km
    e = angulo_elevacao
    lat_gs = np.radians(coordenadas_alvo[0]) # latitude ground station
    long_gs = np.radians(coordenadas_alvo[1])   # longitude ground station
    # raio angular da Terra
    df = dataframe

    for i in range(0, len(df)):
        H = np.linalg.norm((df.iloc[i,0], df.iloc[i,1], df.iloc[i,2]))
        rho = np.arcsin(R_E / (R_E + H))
        print(f'Raio angular da Terra: {np.degrees(rho)} deg')
        # Periodo orbital
        P = PO.periodo_orbital(a)
        print(f'Periodo Orbital: {P / 60} min')

        # Max angulo nadir
        n_max = np.arcsin(np.sin(rho) * np.cos(e))
        print(f'Maximo angulo nadir: {np.degrees(n_max)} deg')
        posicao_gs = [R_E * np.cos(lat_gs) * np.cos(long_gs), R_E * np.cos(lat_gs) * np.sin(long_gs), R_E * np.sin(lat_gs)]
        lamb = np.dot(posicao_gs, ([df.iloc[i,0], df.iloc[i,1], df.iloc[i,2]]))
        # Maximo angulo central da Terra
        lamb_max = np.pi/2 - (e) - (n_max)
        print(f'Maximo angulo central: {(lamb_max)} deg')
        if lamb < (lamb):
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
    return 'batata'
a = tempo_visada()