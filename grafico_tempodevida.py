from Propagador_Orbital import propagador_orbital
from Calor_Incidente import calor_incidente
from Periodo_Orbital import periodo_orbital
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Plots
import os, sys

def resource_path(relative_path):

    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

df = pd.read_csv('dados_ECI.csv', sep=',')
df.plot('ite', 'r')
plt.show()