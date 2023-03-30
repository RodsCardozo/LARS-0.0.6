"""primeira versáo de um modelo lumped para cubesats

Teste 1 - Criando um modelo de caixa com radiação interna

comparação com o artigo de Garzón et al.


"""
"""
Primeira etapa: criação dos nós

"""
import fator_forma_interno as ff
import numpy as np
import pandas as pd

lump = pd.DataFrame()
A1 = {'nome': 1,
      'x': 5.0,
      'y': 0.0,
      'z': 0.0,
      'Lx': 0.2,
      'Ly': 10,
      'Lz': 30,
      'n': 'i',
      'emissividade': 0.99,
      'densidade': 1000.0,
      'condutividade térmica': 200
      }
lump = pd.DataFrame(A1, index=[0])
A2 = {'nome': 2,
      'x': 0.0,
      'y': 5.0,
      'z': 0.0,
      'Lx': 10.0,
      'Ly': 0.2,
      'Lz': 30,
      'n': 'j',
      'emissividade': 0.99,
      'densidade': 1000.0,
      'condutividade térmica': 200
      }
lump.loc[len(lump)] = A2
#print(lump)
A3 = {'nome': 3,
      'x': -5.0,
      'y': 0.0,
      'z': 0.0,
      'Lx': 0.2,
      'Ly': 10,
      'Lz': 30,
      'n': 'i',
      'emissividade': 0.99,
      'densidade': 1000.0,
      'condutividade térmica': 200
      }
lump.loc[len(lump)] = A3
#print(lump)
A4 = {'nome': 4,
      'x': 0.0,
      'y': -5.0,
      'z': 0.0,
      'Lx': 10.0,
      'Ly': 0.2,
      'Lz': 30.0,
      'n': 'j',
      'emissividade': 0.99,
      'densidade': 1000.0,
      'condutividade térmica': 200
      }
lump.loc[len(lump)] = A4
#print(lump)
A5 = {'nome': 5,
      'x': 0.0,
      'y': 0.0,
      'z': -15.0,
      'Lx': 10.0,
      'Ly': 10.0,
      'Lz': 0.2,
      'n': 'k',
      'emissividade': 0.99,
      'densidade': 1000.0,
      'condutividade térmica': 200
      }
lump.loc[len(lump)] = A5
#print(lump)
A6 = {'nome': 6,
      'x': 0.0,
      'y': 0.0,
      'z': 15.0,
      'Lx': 10.0,
      'Ly': 10.0,
      'Lz': 0.2,
      'n': 'k',
      'emissividade': 0.99,
      'densidade': 1000.0,
      'condutividade térmica': 200
      }
lump.loc[len(lump)] = A6
print(lump)

"""
Segunda etapa: junção das interfaces. Designar quem está tocando quem e quem está vendo quem. 
"""

""" Fator de forma interno para cavidades quadradas """

def fator_forma(A, B, tipo_fator_forma, **kwargs):
      import fator_forma_interno as ff
      tipo = str(tipo_fator_forma)
      if tipo.lower() == 'paralelo':
            Fij = ff.fator_paralelo_dict(A, B, kwargs['X'], kwargs['Y'])
      elif tipo.lower() == 'perpendicular':
            Fij = ff.fator_perpendicular_dict(A, B, kwargs['X'], kwargs['Y'], kwargs['Z'])
      return Fij

def ff_interno(df):
      """
      :param df: DataFrame with nodes
      :return: Array with all the Form Factors for a square cavity
      """
      lump = df
      N = len(df) # número de nós
      M = np.zeros((N,N))
      for i in range(N):
            for j in range(N):
                  #A = np.dot(lump.iloc[i,7], lump.iloc[j,7])
                  #print(A)
                  if i == j :
                        M[i,j] = 0

                  elif lump.iloc[i,7] == lump.iloc[j,7]:
                        if lump.iloc[i,7] == 'i':
                              M[i,j] = ff.fator_paralelo_dict(A=lump.iloc[i], B=lump.iloc[j], X=lump.iloc[i,5], Y=lump.iloc[i,6])
                        elif lump.iloc[i,7] == 'j':
                              M[i, j] = ff.fator_paralelo_dict(A=lump.iloc[i], B=lump.iloc[j], X=lump.iloc[i, 4],
                                                               Y=lump.iloc[i, 6])
                        else:
                              M[i, j] = ff.fator_paralelo_dict(A=lump.iloc[i], B=lump.iloc[j], X=lump.iloc[i, 4],
                                                               Y=lump.iloc[i, 5])
                  else:
                        if lump.iloc[i,7] == 'i' and lump.iloc[j,7] == 'k':
                              X = lump.iloc[i,5]
                              Y = lump.iloc[i,6]
                              Z = lump.iloc[j,4]
                              M[i,j] = ff.fator_perpendicular_dict(A=lump.iloc[i], B = lump.iloc[j], X=X, Y=Y, Z=Z)

                        elif lump.iloc[i,7] == 'i' and lump.iloc[j,7] == 'j':
                              X = lump.iloc[i,6]
                              Y = lump.iloc[i,5]
                              Z = lump.iloc[j,4]
                              M[i, j] = ff.fator_perpendicular_dict(A=lump.iloc[i], B=lump.iloc[j], X=X, Y=Y, Z=Z)

                        elif lump.iloc[i,7] == 'j' and lump.iloc[j,7] == 'i':
                              X = lump.iloc[i,6]
                              Y = lump.iloc[i,4]
                              Z = lump.iloc[j,5]
                              M[i, j] = ff.fator_perpendicular_dict(A=lump.iloc[i], B=lump.iloc[j], X=X, Y=Y, Z=Z)

                        elif lump.iloc[i,7] == 'j' and lump.iloc[j,7] == 'k':
                              X = lump.iloc[i,4]
                              Y = lump.iloc[i,6]
                              Z = lump.iloc[j,5]
                              M[i, j] = ff.fator_perpendicular_dict(A=lump.iloc[i], B=lump.iloc[j], X=X, Y=Y, Z=Z)

                        elif lump.iloc[i,7] == 'k' and lump.iloc[j,7] == 'i':
                              X = lump.iloc[i,5]
                              Y = lump.iloc[i,4]
                              Z = lump.iloc[j,6]
                              M[i, j] = ff.fator_perpendicular_dict(A=lump.iloc[i], B=lump.iloc[j], X=X, Y=Y, Z=Z)

                        elif lump.iloc[i,7] == 'k' and lump.iloc[j,7] == 'j':
                              X = lump.iloc[i,4]
                              Y = lump.iloc[i,5]
                              Z = lump.iloc[j,6]
                              M[i, j] = ff.fator_perpendicular_dict(A=lump.iloc[i], B=lump.iloc[j], X=X, Y=Y, Z=Z)

      return M


""" Fator de Gebhart """

def fator_gebhart(df, F):
      """
      :param df: Dataframe with nodes
      :param F: Array with all Form Factors
      :return B: Array with Gebhart factors
      """
      F = ff_interno(lump)
      N = len(F)
      M = np.zeros((N,N))
      for i in range(N):
            for k in range(N):
                  if i == k:
                        d = 1
                  else:
                        d = 0
                  M[i,k] = (1 - lump.iloc[k,8])*F[i,k] - d
      b = []
      for j in range(N):
            b.append(-lump.iloc[j,8]*F[j])
      B = []
      for i in range(N):
            B.append(np.linalg.solve(M,b[i]))

      return B



if __name__ == '__main__':
      print(ff_interno(lump))
      #print(help(ff_interno))
      print(fator_gebhart(lump, ff_interno(lump)))


""" Conexão entre nós """

""" Matriz de condução """

K = {'1': [2,4,5,6],
     '2': [1,3,5,6],
     '3': [2,4,5,6],
     '4': [1,2,5,6],
     '5': [1,2,3,4],
     '6': [1,2,3,4]
}

print(K)
N = 6
T = {}

for x in range(1, N):
    T["T{0}".format(x)] = 50
print(T)

M_K = np.zeros((N,N))
A = []
for i in range(1, len(K)):
      A.append(K['{0}'.format(i)])

print(lump.loc(0))

















