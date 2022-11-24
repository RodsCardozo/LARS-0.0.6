# -*- coding: utf-8 -*-
# propagador_orbital_mk4.py>
"""
    Universidade Federal de Santa Catarina
    Laboratory of Applications and Research in Space - LARS
    Orbital Mechanics Division

    Título do Algoritmo = Codigo principal de propagacao orbital e analise termica
    Autor= Rodrigo S. Cardozo
    Versão= 0.0.1
    Data= 24/10/2022

"""
def propagador_orbital(data, semi_eixo, excentricidade, Raan, argumento_perigeu, anomalia_verdadeira,
        inclinacao, num_orbitas, delt, psi, teta, phi, psip, tetap, phip, massa, largura, comprimento, altura):

    """
    & Data = inicio da simulacao
    & Semi_eixo = altitude no periapse da orbita
    & excentricidade = e
    & Raan= Angulo da posicao do nodo ascendente
    & argumento_perigeu = Angulo da orientacao da linha dos apses
    & anomalia_verdadeira = algulo do vetor posicao e a linha dos apses com origem no foco
    & inclinacao = inclinacao da orbita
    & num_orbitas = numero de orbitas a serem simuladas
    & psi = primeiro angulo de Euler
    & teta = segundo angulo de Euler
    & psi = terceiro angulo de Euler
    & psip = velocidade angular do primeiro angulo de Euler
    & tetap = velocidade angular do segundo angulo de Euler
    & phip = velocidade angular do terceiro angulo de Euler
    & massa = massa do cubesat ex: 3 kg
    & largura = largura do cubsat ex: 0.1m
    & comprimento = comprimento do cubesat ex: 0.1m
    & altura = altura do cubesat ex: 0.2

    """
    print("Propagando o movimento")
    import numpy as np
    import pandas as pd
    from scipy.integrate import odeint
    from datetime import datetime
    from datetime import timedelta
    from nrlmsise00 import msise_model
    from Periodo_Orbital import periodo_orbital
    def propagador(q, t, Rho, velocidade, massa, largura, comprimento, altura, CD, posicao, Area_transversal):  # funcao para integrar
        import numpy as np
        mu = 398600
        J2 = 1.08263e-3
        R_terra = 6371
        rho = Rho
        r = posicao
        m = float(massa)  # massa do cubesat
        a = float(largura)  # comprimento do sat
        b = float(comprimento)  # largura do sat
        c = float(altura)  # altura do sat
        Ix3 = (m / 12) * (b ** 2 + c ** 2)  # momento de inercia na direcao x
        Iy3 = (m / 12) * (a ** 2 + c ** 2)  # momento de inercia na direcao y
        Iz3 = (m / 12) * (a ** 2 + b ** 2)  # momento de inercia na direcao z
        h, ecc, anomalia_verdadeira, raan, inc, arg_per, q0, q1, q2, q3, wx3, wy3, wz3 = q

        dMdt = [r*((-1/(2*r))*h*rho*velocidade*((CD*Area_transversal)/m) - 1.5 * ((J2*mu*R_terra**2)/r**4) * np.sin(inc)**2*np.sin(2*(arg_per + anomalia_verdadeira))),

                (h/mu)*np.sin(anomalia_verdadeira)*((-1/(2*h))*mu*ecc*rho*velocidade*((CD*Area_transversal)/m)*np.sin(anomalia_verdadeira)
                - 1.5*((J2*mu*R_terra**2)/r**4)*(1 - 3*np.sin(inc)**2*np.sin(arg_per + anomalia_verdadeira)**2))
                + (((-1/(2*r))*h*rho*velocidade*((CD*Area_transversal)/m) - 1.5*((J2*mu*R_terra**2)/r**4)*np.sin(inc)**2*np.sin(2*(arg_per
                + anomalia_verdadeira)))/(mu*h))*((h**2 + mu*r)*np.cos(anomalia_verdadeira) + mu*ecc*r),

                (h/r**2 + ((h**2*np.cos(anomalia_verdadeira))/(mu*ecc*h))*((-1/(2*h))*mu*ecc*rho*velocidade*((CD*Area_transversal)/m)*np.sin(anomalia_verdadeira)
                - 1.5*((J2*mu*R_terra**2)/r**4)*(1 - 3*np.sin(inc)**2*np.sin(arg_per + anomalia_verdadeira)**2))
                - (r + h**2/mu)*(np.sin(anomalia_verdadeira)/(ecc*h))*((-1/(2*r))*h*rho*velocidade*((CD*Area_transversal)/m)
                - (1.5 * (J2*mu*R_terra**2)/r**4) * np.sin(inc)**2 * np.sin(2*(arg_per + anomalia_verdadeira)))),

                (r/(h*np.sin(inc)))*np.sin(arg_per + anomalia_verdadeira)*(- 1.5*((J2*mu*R_terra**2)/r**4)*np.sin(2*inc)*np.sin(arg_per + anomalia_verdadeira)),

                (r / (h)) * np.cos(arg_per + anomalia_verdadeira) * (- 1.5 * ((J2 * mu * R_terra ** 2) / r ** 4) * np.sin(2 * inc) * np.sin(arg_per + anomalia_verdadeira)),

                (-1/(ecc*h))*((h**2/mu)*np.cos(anomalia_verdadeira)*((-1/(2*h))*mu*ecc*rho*velocidade*((CD*Area_transversal)/m)*np.sin(anomalia_verdadeira)
                - 1.5*((J2*mu*R_terra**2)/r**4)*(1 - 3*np.sin(inc)**2*np.sin(arg_per + anomalia_verdadeira)**2))
                - (r + h**2/mu)*np.sin(anomalia_verdadeira)*((-1/(2*r))*h*rho*velocidade*((CD*Area_transversal)/m)
                - 1.5 * ((J2*mu*R_terra**2)/r**4) * np.sin(inc)**2 * np.sin(2*(arg_per + anomalia_verdadeira))))
                - ((r*np.sin(arg_per + anomalia_verdadeira))/(h*np.tan(inc)))*(- 1.5 * (J2*mu*R_terra**2)/r**4 * np.sin(2*inc) * np.sin(arg_per + anomalia_verdadeira)),

                0.5 * (- q1 * wx3 - q2 * wy3 - q3 * wz3),
                0.5 * (+ q0 * wx3 + q3 * wy3 - q2 * wz3),
                0.5 * (- q3 * wx3 + q0 * wy3 + q1 * wz3),
                0.5 * ( q2 * wx3 - q1 * wy3 + q0 * wz3),
                ((Iy3 - Iz3) / Ix3) * wy3 * wz3,
                ((Iz3 - Ix3) / Iy3) * wz3 * wx3,
                ((Ix3 - Iy3) / Iz3) * wx3 * wy3]
        return dMdt

    def rho(data, altitude, latitude, longitude):
        densidade = msise_model(data, altitude, latitude, longitude, 150, 150, 4, lst=16)
        rho = densidade[0][5] * 1000
        return rho

    # condicoes iniciais

    SMA = float(semi_eixo)  # semi eixo maior
    ecc0 = float(excentricidade)  # ecentricidade da orbita
    Raan0 = np.radians(float(Raan))  # ascencao direita do nodo ascendente
    arg_per0 = np.radians(float(argumento_perigeu))  # argumento do perigeu
    true_anomaly0 = np.radians(float(anomalia_verdadeira))  # anomalia verdadeira
    inc0 = np.radians(float(inclinacao))  # inclinacao
    rp0 = SMA*(1-ecc0)
    T_orb = periodo_orbital(SMA)
    mu = 398600
    J2 = 1.08263e-3
    R_terra = 6371
    h0 = np.sqrt(rp0*mu*(1 - ecc0))
    psi = float(np.radians(psi))  # angulo inicial de PSI
    teta = float(np.radians(teta))  # angulo inicial de TETA
    phi = float(np.radians(phi))  # angulo inicial de PHI
    psip = float(psip)  # velocidade angular do angulo PSI
    tetap = float(tetap)  # velocidade angular do angulo TETA
    phip = float(phip)  # velocidade angular do angulo PHI
    wx3_i = float(-psip * np.sin(teta) * np.cos(phi) + tetap * np.sin(phi))  # velocidade angular do corpo em x
    wy3_i = float(psip * np.sin(teta) * np.sin(phi) + tetap * np.cos(phi))  # velocidade angular do corpo em y
    wz3_i = float(psip * np.cos(teta) + phip)  # velocidade angular do corpo em z

    q0 = float((np.cos(psi / 2) * np.cos(teta / 2) * np.cos(phi / 2) - np.sin(psi / 2) * np.cos(teta / 2) * np.sin(
              phi / 2)))  # quaternion q0
    q1 = float((np.cos(psi / 2) * np.sin(teta / 2) * np.sin(phi / 2) - np.sin(psi / 2) * np.sin(teta / 2) * np.cos(
                phi / 2)))  # quaternion q1
    q2 = float((np.cos(psi / 2) * np.sin(teta / 2) * np.cos(phi / 2) + np.sin(psi / 2) * np.sin(teta / 2) * np.sin(
                phi / 2)))  # quaternion q2
    q3 = float((np.sin(psi / 2) * np.cos(teta / 2) * np.cos(phi / 2) + np.cos(psi / 2) * np.cos(teta / 2) * np.sin(
                phi / 2)))  # quaternion q3

    x_rot = np.cos(argumento_perigeu) * np.cos(Raan) - np.cos(inclinacao) * np.sin(argumento_perigeu) * np.sin(Raan)
    y_rot = np.cos(argumento_perigeu) * np.sin(Raan) + np.cos(inclinacao) * np.sin(argumento_perigeu) * np.cos(Raan)
    z_rot = np.sin(inclinacao) * np.sin(argumento_perigeu)

    Posi_ini = int(rp0)*[x_rot, y_rot, z_rot]

    lamb_e = (np.arctan2(Posi_ini[1], Posi_ini[0]))

    lat = []
    long = []
    # comeco da integracao

    DELTAT = delt
    mu = 398600
    J2 = 1.08263e-3
    R_terra = 6371
    Time_step = delt
    passo = 10000
    ini_date = data
    n = num_orbitas
    T = T_orb*n
    t = np.linspace(0, Time_step, passo)
    data = []
    solution = []
    time_simu = []
    cont = 0
    while cont < T:
        qi = [h0, ecc0, true_anomaly0, Raan0, inc0, arg_per0, q0, q1, q2, q3, wx3_i, wy3_i, wz3_i]
        altitude = rp0 - R_terra
        xp = (h0 ** 2 / mu) * (1 / (1 + ecc0 * np.cos(true_anomaly0))) * np.cos(true_anomaly0)
        yp = (h0 ** 2 / mu) * (1 / (1 + ecc0 * np.cos(true_anomaly0))) * np.sin(true_anomaly0)
        zp = 0
        X_ECI = ((np.cos(Raan0) * np.cos(arg_per0) - np.sin(Raan0) * np.sin(arg_per0) * np.cos(inc0)) * xp
                 + (-np.cos(Raan0) * np.sin(arg_per0) - np.sin(Raan0) * np.cos(inc0) * np.cos(arg_per0)) * yp
                 + np.sin(Raan0) * np.sin(inc0) * zp)

        Y_ECI = ((np.sin(Raan0) * np.cos(arg_per0) + np.cos(Raan0) * np.cos(inc0) * np.sin(arg_per0)) * xp
                 + (-np.sin(Raan0)*np.sin(arg_per0) + np.cos(Raan0)*np.cos(inc0)*np.cos(arg_per0)) * yp
                 - np.cos(Raan0)*np.sin(inc0) * zp)

        Z_ECI = (np.sin(inc0) * np.sin(arg_per0) * xp
                 + np.sin(inc0) * np.cos(arg_per0) * yp
                 + np.cos(inc0)*zp)
        posicao = np.linalg.norm(np.array([X_ECI, Y_ECI, Z_ECI]))

        RAAN = lamb_e - ((2*np.pi)/(24*3600 + 56*60 + 4))*DELTAT

        lamb_e = RAAN

        X_ECEF = ((np.cos(RAAN) * np.cos(arg_per0) - np.sin(RAAN) * np.sin(arg_per0) * np.cos(inc0)) * xp
                 + (-np.cos(RAAN) * np.sin(arg_per0) - np.sin(RAAN) * np.cos(inc0) * np.cos(arg_per0)) * yp
                 + np.sin(RAAN) * np.sin(inc0) * zp)

        Y_ECEF = ((np.sin(RAAN) * np.cos(arg_per0) + np.cos(RAAN) * np.cos(inc0) * np.sin(arg_per0)) * xp
                 + (-np.sin(RAAN) * np.sin(arg_per0) + np.cos(RAAN) * np.cos(inc0) * np.cos(arg_per0)) * yp
                 - np.cos(RAAN) * np.sin(inc0) * zp)

        Z_ECEF = (np.sin(inc0) * np.sin(arg_per0) * xp
                 + np.sin(inc0) * np.cos(arg_per0) * yp
                 + np.cos(inc0) * zp)

        r = np.array([X_ECEF, Y_ECEF, Z_ECEF])

        latitude = np.degrees((np.arctan2(r[2], np.sqrt(r[0] ** 2 + r[1] ** 2))))
        longitude = np.degrees((np.arctan2(r[1], r[0])))

        lat.append(latitude)
        long.append(longitude)

        Rho = rho(ini_date, altitude, latitude, longitude)
        velocidade = (mu/h0)*np.sqrt(np.sin(true_anomaly0)**2 + (ecc0 + np.cos(true_anomaly0))**2)
        massa = massa
        CD = 2.2
        Area_transversal = 0.1*0.1
        largura = largura
        comprimento = comprimento
        altura = altura
        sol = odeint(propagador, qi, t, args=(Rho, velocidade, massa, largura, comprimento, altura, CD, posicao, Area_transversal))
        solution.append(sol[passo - 1])
        h0 = sol[passo-1][0]
        ecc0 = sol[passo-1][1]
        true_anomaly0 = sol[passo-1][2]
        Raan0 = sol[passo-1][3]
        inc0 = sol[passo-1][4]
        arg_per0 = sol[passo-1][5]
        posicao = (h0**2/mu)*(1/(1-ecc0*np.cos(true_anomaly0)))
        q0 = sol[passo-1][6]
        q1 = sol[passo-1][7]
        q2 = sol[passo-1][8]
        q3 = sol[passo-1][9]
        wx3_i = sol[passo-1][10]
        wy3_i = sol[passo-1][11]
        wz3_i = sol[passo-1][12]
        cont = cont + Time_step
        time_simu.append(cont)
        final_date = timedelta(seconds=Time_step)
        ini_date = ini_date + final_date
        data.append(ini_date)

    solucao = pd.DataFrame(solution, columns=['h', 'ecc', 'anomalia_verdadeira', 'raan', 'inc', 'arg_per', 'q0', 'q1', 'q2', 'q3', 'wx3', 'wy3', 'wz3'])
    solucao['X_perifocal'] = (solucao['h']**2/mu)*(1/(1 + solucao['ecc']*np.cos(solucao['anomalia_verdadeira'])))*np.cos(solucao['anomalia_verdadeira'])
    solucao['Y_perifocal'] = (solucao['h']**2/mu)*(1/(1 + solucao['ecc']*np.cos(solucao['anomalia_verdadeira'])))*np.sin(solucao['anomalia_verdadeira'])
    solucao['Z_perifocal'] = 0
    solucao['distancia'] = np.sqrt(solucao['X_perifocal']**2 + solucao['Y_perifocal']**2)

    df = pd.DataFrame()
    phi = []
    for i in range(0, len(solution), 1):
        a = (np.arctan2(2*(solucao.iloc[i, 8]*solucao.iloc[i, 9] - solucao.iloc[i, 6]*solucao.iloc[i, 7]), 2*(solucao.iloc[i, 7]*solucao.iloc[i, 9] + solucao.iloc[i, 6]*solucao.iloc[i, 8])))
        if np.linalg.norm(a) < 0.000000001:
            phi.append(0.0)
        elif np.linalg.norm(np.linalg.norm(a) - np.pi) < 0.000001:
            phi.append(0.0)
        else:
            phi.append(a)

    dfphi = pd.DataFrame(np.unwrap(phi), columns=['PHI'])
    df = pd.concat([df, dfphi], axis=1)

    teta = []
    for i in range(0, len(solution), 1):
        a = np.arccos(2*(solucao.iloc[i, 6]**2 + solucao.iloc[i, 9]**2) - 1)

        if np.linalg.norm(a) < 0.000000001:
            teta.append(0.0)
        elif a == -np.pi:
            teta.append(0.0)
        else:
            teta.append(a)
    dfteta = pd.DataFrame(teta, columns=['TETA'])
    df = pd.concat([df, dfteta], axis=1)

    psi = []
    for i in range(0, len(solution), 1):
        a = (np.arctan2(2*(solucao.iloc[i, 8]*solucao.iloc[i, 9] + solucao.iloc[i, 6]*solucao.iloc[i, 7]), -2*(solucao.iloc[i, 7]*solucao.iloc[i, 9] - solucao.iloc[i, 6]*solucao.iloc[i, 8])))
        if np.linalg.norm(a) < 0.000000001:
            psi.append(0.0)
        elif np.linalg.norm(np.linalg.norm(a) - np.pi) < 0.000001:
            psi.append(0.0)
        else:
            psi.append(a)

    dfpsi = pd.DataFrame(np.unwrap(psi), columns=['PSI'])
    df = pd.concat([df, dfpsi], axis=1)

    '''df = pd.DataFrame()
    df['PHI'] = np.unwrap(np.arctan2(2*(solucao['q2']*solucao['q3'] - solucao['q0']*solucao['q1']), 2*(solucao['q1']*solucao['q3'] + solucao['q0']*solucao['q2'])))
    df['TETA'] = np.arccos(2*(solucao['q0']**2 + solucao['q3']**2) - 1)
    df['PSI'] = np.unwrap(np.arctan2(-2*(solucao['q2']*solucao['q3'] + solucao['q0']*solucao['q1']), 2*(solucao['q1']*solucao['q3'] - solucao['q0']*solucao['q2']) ))'''

    df['X_ECI'] = ((np.cos(solucao['raan'])*np.cos(solucao['arg_per']) - np.sin(solucao['raan'])*np.sin(solucao['arg_per'])*np.cos(solucao['inc']))*solucao['X_perifocal']

                        + (-np.cos(solucao['raan'])*np.sin(solucao['arg_per']) - np.sin(solucao['raan'])*np.cos(solucao['inc'])*np.cos(solucao['arg_per']))*solucao['Y_perifocal']

                        + np.sin(solucao['raan'])*np.sin(solucao['inc'])*solucao['Z_perifocal'])

    df['Y_ECI'] = ((np.sin(solucao['raan'])*np.cos(solucao['arg_per']) + np.cos(solucao['raan'])*np.cos(solucao['inc'])*np.sin(solucao['arg_per']))*solucao['X_perifocal']

                        + (-np.sin(solucao['raan'])*np.sin(solucao['arg_per']) + np.cos(solucao['raan'])*np.cos(solucao['inc'])*np.cos(solucao['arg_per']))*solucao['Y_perifocal']

                        - np.cos(solucao['raan'])*np.sin(solucao['inc'])*solucao['Z_perifocal'])

    df['Z_ECI'] = (np.sin(solucao['inc'])*np.sin(solucao['arg_per'])*solucao['X_perifocal']
                        + np.sin(solucao['inc'])*np.cos(solucao['arg_per'])*solucao['Y_perifocal']
                        + np.cos(solucao['inc'])*solucao['Z_perifocal'])

    df1 = pd.DataFrame(lat, columns=['latitude'])
    df = pd.concat([df, df1], axis=1)

    df2 = pd.DataFrame(long, columns=['longitude'])
    df = pd.concat([df, df2], axis=1)

    df3 = pd.DataFrame(data, columns=['Data'])
    df = pd.concat([df, df3], axis=1)

    df4 = pd.DataFrame(time_simu, columns=['Tempo'])
    df = pd.concat([df, df4], axis=1)

    import os.path
    df.to_csv(os.path.join('./data/', 'dados_ECI.csv'), sep=',')
    import os.path
    solucao.to_csv(os.path.join('./data/', 'solver.csv'), sep=',')

    return df