from unittest import main, TestCase
import numpy as np
import pandas as pd
import requests
import io

import PyBounds.pyBounds as pyBounds

url = "https://raw.githubusercontent.com/JMSLab/PyBounds/master/PyBounds/examples/roberts_schlenker_2013.csv?token=GHSAT0AAAAAABSA4OUAMVHX6ZDUVHMSQRGGYTEQC5Q"
download = requests.get(url).content

data = pd.read_csv(io.StringIO(download.decode('utf-8')))
df = pd.DataFrame(data)
p_t = df['logPD0']
q_t = df['logD0']


class Test(TestCase):
    def test_add(self):
        kinf = pyBounds.Bounds(p_t, q_t, k=np.inf, maxB=0.1)
        self.assertAlmostEqual(kinf.Theta_hat_k_B(0.07)[0], -0.122187908332218, delta=0.000001)
        self.assertAlmostEqual(kinf.Theta_hat_k_B(0.07)[1], 0.134092308480860, delta=0.000001)
        k2 = pyBounds.Bounds(p_t, q_t, k=2, maxB=0.04)
        self.assertAlmostEqual(k2.Theta_hat_k_B(0.03)[0], -0.123518421564038, delta=0.000001)
        self.assertAlmostEqual(k2.Theta_hat_k_B(0.03)[1], 0.106798443278028, delta=0.000001)
        with self.assertRaises(Exception):
            pyBounds.Bounds(p_t, q_t, k=-1, maxB=0.04)
        with self.assertRaises(Exception):
            kneg = pyBounds.Bounds(p_t, q_t, k=9, maxB=0)
            kneg.plot()
        with self.assertRaises(Exception):
            kneg = pyBounds.Bounds(p_t, q_t, k=9, maxB=0)
            kneg.intervals()

        # Checking for k close to 2
        k_2 = pyBounds.Bounds(p_t, q_t, k=2, maxB=0.04)
        k_under = pyBounds.Bounds(p_t, q_t, k=1.99, maxB=0.04)
        interval_2 = k_2.intervals()
        interval_under = k_under.intervals()
        for i in range(len(interval_2)):
            self.assertAlmostEqual(interval_2[i][0], interval_under[i][0], delta=0.001)
            self.assertAlmostEqual(interval_2[i][1], interval_under[i][1], delta=0.001)

        # Checking for k close to infinity
        k_inf = pyBounds.Bounds(p_t, q_t, k=np.inf, maxB=0.1)
        k_under = pyBounds.Bounds(p_t, q_t, k=2 ** 7, maxB=0.1)
        interval_inf = k_inf.intervals()
        interval_under = k_under.intervals()
        for i in range(len(interval_inf)):
            self.assertAlmostEqual(interval_inf[i][0], interval_under[i][0], delta=0.01)
            self.assertAlmostEqual(interval_inf[i][1], interval_under[i][1], delta=0.01)

        # Checking for T = 2
        p_t2 = [1, 2]
        q_t2 = [1, 2]
        k7 = pyBounds.Bounds(p_t2, q_t2, k=7, maxB=0.1)
        self.assertAlmostEqual(k7.underline_theta(0.08), 0.92, delta=0.000001)
        self.assertAlmostEqual(k7.overline_theta(0.08), 1.08, delta=0.000001)
        k9 = pyBounds.Bounds(p_t2, q_t2, k=9, maxB=0.1)
        self.assertAlmostEqual(k9.underline_theta(0.08), 0.92, delta=0.000001)
        self.assertAlmostEqual(k9.overline_theta(0.08), 1.08, delta=0.000001)

        # testing k-mean
        self.assertAlmostEqual(k7.kmean(2, p_t2), np.sqrt(5/2), delta=0.000001)
        self.assertAlmostEqual(k7.kmean(5, p_t2), (33 / 2) ** (1/5), delta=0.000001)

        # testing delta_epsilon
        self.assertAlmostEqual(k7.delta_epsilon(theta=4), -3, delta=0.000001)
        self.assertAlmostEqual(k7.delta_epsilon(theta=10), -9, delta=0.000001)

        # testing M_hat_k
        self.assertAlmostEqual(k7.M_hat_k(theta=4), 3, delta=0.000001)
        self.assertAlmostEqual(k9.M_hat_k(theta=10), 9, delta=0.000001)

        # testing B_tilde
        self.assertAlmostEqual(k7.B_tilde(), 0, delta=0.000001)
        self.assertAlmostEqual(k9.B_tilde(), 0, delta=0.000001)

        # testing underline_B
        self.assertAlmostEqual(k7.underline_B(), 0, delta=0.000001)
        self.assertAlmostEqual(k9.underline_B(), 0, delta=0.000001)

        p_t3 = [6, 7]
        q_t3 = [10, 12]
        k7 = pyBounds.Bounds(p_t3, q_t3, k=7, maxB=0.1)
        self.assertAlmostEqual(k7.underline_theta(0.08), 1.92, delta=0.000001)
        self.assertAlmostEqual(k7.overline_theta(0.08), 2.08, delta=0.000001)
        k9 = pyBounds.Bounds(p_t3, q_t3, k=9, maxB=0.1)
        self.assertAlmostEqual(k9.underline_theta(0.08), 1.92, delta=0.000001)
        self.assertAlmostEqual(k9.overline_theta(0.08), 2.08, delta=0.000001)

        # testing kmean
        self.assertAlmostEqual(k7.kmean(2, p_t3), np.sqrt(85 / 2), delta=0.000001)
        self.assertAlmostEqual(k7.kmean(2, q_t3), np.sqrt(244 / 2), delta=0.000001)
        self.assertAlmostEqual(k7.kmean(5, p_t3), (24583 / 2) ** (1/5), delta=0.000001)
        self.assertAlmostEqual(k7.kmean(5, q_t3), (348832 / 2) ** (1/5), delta=0.000001)

        # testing delta_epsilon
        self.assertAlmostEqual(k7.delta_epsilon(7), -5, delta=0.000001)
        self.assertAlmostEqual(k7.delta_epsilon(13), -11, delta=0.000001)

        # testing M_hat_k
        self.assertAlmostEqual(k7.M_hat_k(theta=7), 5, delta=0.000001)
        self.assertAlmostEqual(k9.M_hat_k(theta=13), 11, delta=0.000001)

        # testing B_tilde
        self.assertAlmostEqual(k7.B_tilde(), 0, delta=0.000001)
        self.assertAlmostEqual(k9.B_tilde(), 0, delta=0.000001)

        # testing underline_B
        self.assertAlmostEqual(k7.underline_B(), 0, delta=0.000001)
        self.assertAlmostEqual(k9.underline_B(), 0, delta=0.000001)

        return


if __name__ == '__main__':
    main()
