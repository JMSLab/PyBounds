import numpy as np
import pandas as pd
import requests
import io
import pyBounds

# import data: p_t, q_t
url = "https://raw.githubusercontent.com/JMSLab/PyBounds/master/PyBounds/roberts_schlenker_2013.csv?token=GHSAT0AAAAAABSA4OUBEMYZWDWOTESJ3ISGYRI737Q"
download = requests.get(url).content

data = pd.read_csv(io.StringIO(download.decode('utf-8')))
df = pd.DataFrame(data)
p_t = df.iloc[:, 4]
q_t = df.iloc[:, 3]

if __name__ == '__main__':
    kinf = pyBounds.Bounds(p_t, q_t, k=np.inf, maxB=0.1)
    kinf.plot()
    k2 = pyBounds.Bounds(p_t, q_t, k=2, maxB=0.04)
    k2.plot()
