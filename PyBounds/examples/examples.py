import numpy as np
import pandas as pd
import requests
import io

import PyBounds.pyBounds as pyBounds

# import data: p_t, q_t
url = "https://raw.githubusercontent.com/JMSLab/PyBounds/b9ab4a7718e36a862cd4fc362e085c39149b05bb/PyBounds/roberts_schlenker_2013.csv?token=GHSAT0AAAAAABSFNAYKQWY7QYRYMLLYARG4YRLLEFA"
download = requests.get(url).content

data = pd.read_csv(io.StringIO(download.decode('utf-8')))
df = pd.DataFrame(data)
p_t = df['logPD0']
q_t = df['logD0']

if __name__ == '__main__':
    kinf = pyBounds.Bounds(p_t, q_t, k=np.inf, maxB=0.1)
    kinf.plot()
    k2 = pyBounds.Bounds(p_t, q_t, k=2, maxB=0.04)
    k2.plot()