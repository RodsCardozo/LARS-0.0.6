import pandas as pd

A1 = {'nome': 1,
      'x': 5.0,
      'y': 0.0,
      'z': 0.0,
      'Lx': 0.2,
      'Ly': 10,
      'Lz': 30,
      'n': [1,0,0],
      'emissividade': 0.99,
      'densidade': 1000.0
      }
print(A1['nome'])
lump1 = pd.DataFrame.from_dict(A1)

import numpy as np
a2 = {'nome': int(input('Nome: ')),
      'x': float(input('x: ')),
      'y': float(input('y: ')),
      'z': float(input('z: ')),
      'Lx': float(input('Lx: ')),
      'Ly': float(input('Ly: ')),
      'Lz': float(input('Lz: ')),
      'n': np.array(input('n: ')),
      'emissividade': float(input('Emissividade: ')),
      'densidade': float(input('Densidade: '))
      }

print(a2)