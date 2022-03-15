import numpy as np
import pandas as pd
import requests
import io

import PyBounds.pyBounds as pyBounds

# import data: p_t, q_t
url = "https://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@raw.githubusercontent.com/JMSLab/PyBounds/8e7bd52752a6fde6ba2ec4a2c3aeb6a4e1aebf02/PyBounds/tests/roberts_schlenker_2013.csv?token=GHSAT0AAAAAABSA4OUB3OKSG5NFCC7X74CEYRZB36Q"
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
