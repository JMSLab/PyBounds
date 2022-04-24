import numpy as np
import pandas as pd
import requests
import io

import PyBounds.pyBounds as pyBounds

# import data: p_t, q_t
url = "https://raw.githubusercontent.com/JMSLab/PyBounds/8f0ee4fc5ab97b9afd116d3d48df232c0b455f51/PyBounds/examples/roberts_schlenker_2013.csv"
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
