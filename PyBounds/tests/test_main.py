from unittest import main, TestCase
import numpy as np
import pandas as pd
import requests
import io

from .. import pyBounds

url = "https://raw.githubusercontent.com/JMSLab/PyBounds/master/PyBounds/roberts_schlenker_2013.csv?token=GHSAT0AAAAAABSA4OUBEMYZWDWOTESJ3ISGYRI737Q"
download = requests.get(url).content

data = pd.read_csv(io.StringIO(download.decode('utf-8')))
df = pd.DataFrame(data)
p_t = df.iloc[:, 4]
q_t = df.iloc[:, 3]

class Test(TestCase):
    def test_add(self):
        kinf = pyBounds.Bounds(p_t, q_t, k=np.inf, maxB=0.1)
        self.assertAlmostEqual(kinf.Theta_hat_k_B(0.07)[0], -0.122187908332218, delta = 0.000001)
        self.assertAlmostEqual(kinf.Theta_hat_k_B(0.07)[1], 0.134092308480860, delta = 0.000001)
        k2 = pyBounds.Bounds(p_t, q_t, k=2, maxB=0.04)
        self.assertAlmostEqual(k2.Theta_hat_k_B(0.03)[0], -0.123518421564038, delta = 0.000001)
        self.assertAlmostEqual(k2.Theta_hat_k_B(0.03)[1], 0.106798443278028, delta = 0.000001)
        return


if __name__ == '__main__':
    main()

