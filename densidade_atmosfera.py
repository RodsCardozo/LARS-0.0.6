import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from datetime import datetime
from datetime import timedelta
from nrlmsise00 import msise_model
def rho(data, altitude, latitude, longitude):
    densidade = msise_model(data, altitude, latitude, longitude, 150, 150, 4, lst=16)
    rho = densidade[0][5] * 1000
    return rho

h = np.linspace(0,1000,50)
print(h)
#rho = rho(datetime(2023,2,24,10,0,0,0), h, 0,0)
densidade = []
for i in range(len(h)):
    densidade.append(rho(datetime(2023,2,24,10,0,0,0), h[i], 0,0))

plt.figure(figsize=(6, 5))
plt.plot(h, densidade, 'g')
plt.yscale('log')
plt.grid()
plt.xlabel('Altitude [km]')
plt.ylabel('Densidade [kg/m³]')
plt.show()


if __name__ == '__main__':
    densidade = rho(datetime(2023,2,24,10,0,0,0), 1.2, -26,-28)
    print(densidade)