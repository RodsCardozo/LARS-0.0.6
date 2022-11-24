
"""
    Universidade Federal de Santa Catarina
    Laboratory of Applications and Research in Space - LARS
    Orbital Mechanics Division

    Título do Algoritmo= Dados de entrada para simulação orbital e radiação incidente
    Autor= Rodrigo S. Cardozo
    Versão= 0.0.1
    Data= 24/10/2022

"""
def calor_incidente(posicao_orientacao, Radiacao_solar, Radiacao_Terra, emissividade_Terra, absortividade_satelite,
                    refletividade_terra, numero_divisoes):

    """
    & Posicao_orientacao = Dataframe com a orientacao do cubesat e a sua posicao
    & Radiacao_solar = Intensidade da irradiacao solar
    & Radiacao_Terra = Intensidade da IR da Terra
    & emissividade_Terra = valor da emissidade médio da Terra
    & absortividade_satelite = valor média da absortividade de cada face do satelite
    & refletividade_Terra = valor da refletividade médio da Terra
    & numero_divisoes = divisao da terra em n elementos de area
    """
    print("Calculando calor")
    import Fator_Forma_2
    import matplotlib.pyplot as plt
    import icosphere
    from Fator_Forma_2 import fator_forma_classico as FS
    import numpy as np
    import pandas as pd
    def cos(angulo):
        import numpy as np
        return np.cos(angulo)
    def sin(angulo):
        import numpy as np
        return np.sin(angulo)
    ''' Dados iniciais da orbita '''

    # propagador orbital

    prop_orb = posicao_orientacao  #propagador_orbital(data, SMA, ecc, Raan, arg_per, true_anomaly, inc, 1, delt, psi0, teta0, phi0, PSIP, TETAP,
                                   #PHIP, massa, largura, comprimento, altura, omegap)

    # Intensidade radiante do sol e terra e valores de emissividade

    Is = Radiacao_solar          # radiacao solar
    Ir = Radiacao_Terra          # radiacao IR Terra
    e = emissividade_Terra       # emissividade Terra
    ai = absortividade_satelite  # absortividade do satelite
    gama = refletividade_terra   # refletividade da Terra

    Vs = np.array([1, 0, 0]) # vetor solar

    # divisao da terra em elementos de area utilizando um icosaedro
    Raio_terra = 6371
    nu = numero_divisoes
    vertices, faces = icosphere.icosphere(nu)
    center = []
    for i in range(0, len(faces), 1):

        a = faces[i][0]
        b = faces[i][1]
        c = faces[i][2]
        A = vertices[a]*(Raio_terra)
        B = vertices[b]*(Raio_terra)
        C = vertices[c]*(Raio_terra)
        x = float(A[0] + B[0] + C[0])
        y = float(A[1] + B[1] + C[1])
        z = float(A[2] + B[2] + C[2])
        center.append([x / 3, y / 3, z / 3])

    As = (4*np.pi*(Raio_terra)**2)/len(center) # area de cada elemento

    vet_terra = pd.DataFrame(center, columns=['Terra_X', 'Terra_Y', 'Terra_Z'])
    Ai = [a * c,
          b * c,
          a * c,
          b * c,
          a * b,
          a * b]
    Ni = [[1, 0, 0],
          [0, 1, 0],
          [-1, 0, 0],
          [0, -1, 0],
          [0, 0, -1],
          [0, 0, 1]]
    df1 = pd.DataFrame(Ni, columns=['x', 'y', 'z'])
    Posicao_orientacao = pd.concat([prop_orb, df1], axis=1)

    # determinacao da orientacao de cada face

    names = [['N1_X', 'N1_Y', 'N1_Z'],
             ['N2_X', 'N2_Y', 'N2_Z'],
             ['N3_X', 'N3_Y', 'N3_Z'],
             ['N4_X', 'N4_Y', 'N4_Z'],
             ['N5_X', 'N5_Y', 'N5_Z'],
             ['N6_X', 'N6_Y', 'N6_Z']]
    R = []
    for j in range(0, len(Ni), 1):
        for i in range(0, len(Posicao_orientacao), 1):
            A = Ni[j]

            psi = float(Posicao_orientacao.iloc[i, 0])
            teta = float(Posicao_orientacao.iloc[i, 1])
            phi = float(Posicao_orientacao.iloc[i, 2])

            vetor = A


            Q_rot = np.array([[np.cos(psi) * np.cos(phi) - np.sin(psi) * np.sin(phi) * np.cos(teta),
                               np.cos(psi) * np.sin(phi) + np.sin(psi) * np.cos(teta) * np.cos(phi),
                               np.sin(psi) * np.sin(teta)],
                              [-np.sin(psi) * np.cos(phi) - np.cos(psi) * np.sin(phi) * np.cos(teta),
                               -np.sin(psi) * np.sin(phi) + np.cos(psi) * np.cos(teta) * np.cos(phi),
                               np.cos(psi) * np.sin(teta)],
                              [np.sin(teta) * np.sin(phi), -np.sin(teta) * np.cos(phi), np.cos(teta)]])

            Q = np.dot(np.transpose(Q_rot), vetor)
            R1 = Q[0]
            R2 = Q[1]
            R3 = Q[2]
            R.append([np.array(R1), np.array(R2), np.array(R3)])

        df2 = pd.DataFrame(R, columns=names[j])
        R = []
        Posicao_orientacao = pd.concat([Posicao_orientacao, df2], axis=1)

    Posicao_orientacao = pd.concat([Posicao_orientacao, vet_terra], axis=1)
    Posicao_orientacao['r'] = np.sqrt(Posicao_orientacao['X_ECI']**2 + Posicao_orientacao['Y_ECI']**2 + Posicao_orientacao['Z_ECI']**2)
    Posicao_orientacao['fim'] = 1
    import os.path
    Posicao_orientacao.to_csv(os.path.join('./data/', 'Posicao_orientacao.csv'), sep=',')

    vetor_terra = []
    for i in range(0, len(vet_terra), 1):
        vetor_terra.append(np.array([(vet_terra.iloc[i, 0]), (vet_terra.iloc[i, 1]), (vet_terra.iloc[i, 2])]))

    vetor_posicao = []
    for i in range(0, len(prop_orb), 1):
        vetor_posicao.append(np.array([np.array(Posicao_orientacao.iloc[i, 3]), np.array(Posicao_orientacao.iloc[i, 4]),
                      np.array(Posicao_orientacao.iloc[i, 5])]))

    '''Inicio do calculo de radiacao'''
    print('Calculando radiacao solar')
    Qs1 = []
    Qs2 = []
    Qs3 = []
    Qs4 = []
    Qs5 = []
    Qs6 = []

    for i in range(0, len(vetor_posicao), 1):

        PSI = np.arccos(np.dot(vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), Vs / np.linalg.norm(Vs)))
        QSI = np.arcsin(Raio_terra / np.linalg.norm((vetor_posicao[i])))

        if PSI + QSI < np.pi:

            A1 = np.array([Posicao_orientacao.iloc[i, 13], Posicao_orientacao.iloc[i, 14], Posicao_orientacao.iloc[i, 15]])
            A2 = np.array([Posicao_orientacao.iloc[i, 16], Posicao_orientacao.iloc[i, 17], Posicao_orientacao.iloc[i, 18]])
            A3 = np.array([Posicao_orientacao.iloc[i, 19], Posicao_orientacao.iloc[i, 20], Posicao_orientacao.iloc[i, 21]])
            A4 = np.array([Posicao_orientacao.iloc[i, 22], Posicao_orientacao.iloc[i, 23], Posicao_orientacao.iloc[i, 24]])
            A5 = np.array([Posicao_orientacao.iloc[i, 25], Posicao_orientacao.iloc[i, 26], Posicao_orientacao.iloc[i, 27]])
            A6 = np.array([Posicao_orientacao.iloc[i, 28], Posicao_orientacao.iloc[i, 29], Posicao_orientacao.iloc[i, 30]])

            k1 = np.dot(A1/np.linalg.norm(A1), Vs)
            k2 = np.dot(A2/np.linalg.norm(A2), Vs)
            k3 = np.dot(A3/np.linalg.norm(A3), Vs)
            k4 = np.dot(A4/np.linalg.norm(A4), Vs)
            k5 = np.dot(A5/np.linalg.norm(A5), Vs)
            k6 = np.dot(A6/np.linalg.norm(A6), Vs)

            if k1 > 0:
                qs1 = ai * Is * k1
                Qs1.append(qs1)
            else:
                Qs1.append(0)
            if k2 > 0:
                qs2 = ai * Is * k2
                Qs2.append(qs2)
            else:
                Qs2.append(0)

            if k3 > 0:
                qs3 = ai * Is * k3
                Qs3.append(qs3)
            else:
                Qs3.append(0)
            if k4 > 0:
                qs4 = ai * Is * k4
                Qs4.append(qs4)
            else:
                Qs4.append(0)
            if k5 > 0:
                qs5 = ai * Is * k5
                Qs5.append(qs5)
            else:
                Qs5.append(0)
            if k6 > 0:
                qs6 = ai * Is * k6
                Qs6.append(qs6)
            else:
                Qs6.append(0)

        else:
            Qs1.append(0)
            Qs2.append(0)
            Qs3.append(0)
            Qs4.append(0)
            Qs5.append(0)
            Qs6.append(0)

    '''Radiacao de albedo incidente'''

    print('Calculando radiacao de albedo')

    Halb1 = 0
    Halb2 = 0
    Halb3 = 0
    Halb4 = 0
    Halb5 = 0
    Halb6 = 0

    Qalb1 = []
    Qalb2 = []
    Qalb3 = []
    Qalb4 = []
    Qalb5 = []
    Qalb6 = []

    from tqdm import tqdm

    for i in tqdm(range(0, len(vetor_posicao), 1), colour='green'):

        A1 = np.array([Posicao_orientacao.iloc[i, 13], Posicao_orientacao.iloc[i, 14], Posicao_orientacao.iloc[i, 15]])
        A2 = np.array([Posicao_orientacao.iloc[i, 16], Posicao_orientacao.iloc[i, 17], Posicao_orientacao.iloc[i, 18]])
        A3 = np.array([Posicao_orientacao.iloc[i, 19], Posicao_orientacao.iloc[i, 20], Posicao_orientacao.iloc[i, 21]])
        A4 = np.array([Posicao_orientacao.iloc[i, 22], Posicao_orientacao.iloc[i, 23], Posicao_orientacao.iloc[i, 24]])
        A5 = np.array([Posicao_orientacao.iloc[i, 25], Posicao_orientacao.iloc[i, 26], Posicao_orientacao.iloc[i, 27]])
        A6 = np.array([Posicao_orientacao.iloc[i, 28], Posicao_orientacao.iloc[i, 29], Posicao_orientacao.iloc[i, 30]])
        phi = (np.dot(Vs, vetor_posicao[i] / np.linalg.norm(vetor_posicao[i])))

        if np.pi > phi >= 0:

            d1 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A1))
            d2 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A2))
            d3 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A3))
            d4 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A4))
            d5 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A5))
            d6 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A6))
            VP = (np.linalg.norm(-vetor_posicao[i]))

            FS1 = FS(VP, d1)
            FS2 = FS(VP, d2)
            FS3 = FS(VP, d3)
            FS4 = FS(VP, d4)
            FS5 = FS(VP, d5)
            FS6 = FS(VP, d6)

            Qalb1.append(ai * gama * Is * FS1 * phi) #ai * gama * Is *

            Qalb2.append(ai * gama * Is * FS2 * phi)

            Qalb3.append(ai * gama * Is * FS3 * phi)

            Qalb4.append(ai * gama * Is * FS4 * phi)

            Qalb5.append(ai * gama * Is * FS5 * phi)

            Qalb6.append(ai * gama * Is * FS6 * phi)
        else:

            Qalb1.append(0.0)

            Qalb2.append(0.0)

            Qalb3.append(0.0)

            Qalb4.append(0.0)

            Qalb5.append(0.0)

            Qalb6.append(0.0)

    ''' Radiacao da terra '''
    print('Calculando radiacao da terra')

    Qrad1 = []
    Qrad2 = []
    Qrad3 = []
    Qrad4 = []
    Qrad5 = []
    Qrad6 = []

    for i in tqdm(range(0, len(vetor_posicao), 1), colour='cyan'):
        A1 = np.array([Posicao_orientacao.iloc[i, 13], Posicao_orientacao.iloc[i, 14], Posicao_orientacao.iloc[i, 15]])
        A2 = np.array([Posicao_orientacao.iloc[i, 16], Posicao_orientacao.iloc[i, 17], Posicao_orientacao.iloc[i, 18]])
        A3 = np.array([Posicao_orientacao.iloc[i, 19], Posicao_orientacao.iloc[i, 20], Posicao_orientacao.iloc[i, 21]])
        A4 = np.array([Posicao_orientacao.iloc[i, 22], Posicao_orientacao.iloc[i, 23], Posicao_orientacao.iloc[i, 24]])
        A5 = np.array([Posicao_orientacao.iloc[i, 25], Posicao_orientacao.iloc[i, 26], Posicao_orientacao.iloc[i, 27]])
        A6 = np.array([Posicao_orientacao.iloc[i, 28], Posicao_orientacao.iloc[i, 29], Posicao_orientacao.iloc[i, 30]])

        d1 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A1))
        d2 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A2))
        d3 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A3))
        d4 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A4))
        d5 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A5))
        d6 = np.arccos(np.dot(-vetor_posicao[i] / np.linalg.norm(vetor_posicao[i]), A6))
        VP = np.linalg.norm(vetor_posicao[i])

        FS1 = FS(VP, d1)
        FS2 = FS(VP, d2)
        FS3 = FS(VP, d3)
        FS4 = FS(VP, d4)
        FS5 = FS(VP, d5)
        FS6 = FS(VP, d6)

        Qrad1.append(e * Ir * (FS1)) #e * Ir *

        Qrad2.append(e * Ir * (FS2))

        Qrad3.append(e * Ir * (FS3))

        Qrad4.append(e * Ir * (FS4))

        Qrad5.append(e * Ir * (FS5))

        Qrad6.append(e * Ir * (FS6))

    rad_sol = []
    for i in range(0, len(Qs1), 1):
        rad_sol.append([Qs1[i], Qs2[i], Qs3[i], Qs4[i], Qs5[i], Qs6[i]])

    Q_sol = pd.DataFrame(rad_sol, columns=['Solar 1', 'Solar 2', 'Solar 3', 'Solar 4', 'Solar 5', 'Solar 6'])

    rad_alb = []
    for i in range(0, len(Qalb1), 1):
        rad_alb.append([Qalb1[i], Qalb2[i], Qalb3[i], Qalb4[i], Qalb5[i], Qalb6[i]])
    Q_alb = pd.DataFrame(rad_alb, columns=['Albedo 1', 'Albedo 2', 'Albedo 3', 'Albedo 4', 'Albedo 5', 'Albedo 6'])

    rad_terra = []
    for i in range(0, len(Qrad1), 1):
        rad_terra.append([Qrad1[i], Qrad2[i], Qrad3[i], Qrad4[i], Qrad5[i], Qrad6[i]])
    Q_terra = pd.DataFrame(rad_terra, columns=['IR Terra 1', 'IR Terra 2', 'IR Terra 3', 'IR Terra 4', 'IR Terra 5', 'IR Terra 6'])

    QT = pd.concat([Q_sol, Q_alb], axis=1)
    QT = pd.concat([QT, Q_terra], axis=1)
    QT['Total 1'] = QT['Solar 1'] + QT['Albedo 1'] + QT['IR Terra 1']
    QT['Total 2'] = QT['Solar 2'] + QT['Albedo 2'] + QT['IR Terra 2']
    QT['Total 3'] = QT['Solar 3'] + QT['Albedo 3'] + QT['IR Terra 3']
    QT['Total 4'] = QT['Solar 4'] + QT['Albedo 4'] + QT['IR Terra 4']
    QT['Total 5'] = QT['Solar 5'] + QT['Albedo 5'] + QT['IR Terra 5']
    QT['Total 6'] = QT['Solar 6'] + QT['Albedo 6'] + QT['IR Terra 6']

    import os.path
    QT.to_csv(os.path.join('./data/', 'Calor_Incidente.csv'), sep=',')
    print('Fim')

    return QT