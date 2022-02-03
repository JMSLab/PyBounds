"""Main file of the package"""
#from .log import init_logger

#logger = init_logger('main.log')
# NEEDS TESTING!!!!!


import pandas as pd
import numpy as np
import requests
import io

# import data: p_t, q_t
url = "https://raw.githubusercontent.com/JMSLab/PyBounds/870170743f155b3a9da680ec64e0d5079e266601/PyBounds/roberts_schlenker_2013.csv?token=GHSAT0AAAAAABQYAQ7N7L2PIVTKPDHHG7WWYQEWC6A"
download = requests.get(url).content

data = pd.read_csv(io.StringIO(download.decode('utf-8')))
df = pd.DataFrame(data)
p_t = df.iloc[:, 0]
q_t = df.iloc[:, 1]


class Bounds:
    def __init__(self, p_t, q_t, k = np.inf, B = None, ngridpoints = 1000):
        self.p = p_t
        self.q = q_t
        if k < 1:
            raise Exception("Requires k>=1")
        else:
            self.k = k
        self.B = B
        self.ngridpoints = ngridpoints
    # implementation of k-mean

    @staticmethod
    def kmean(k, vector):
        """implementation of the k-mean function
            :param k: the "k" in k-mean
            :param vector: list of elements of which to take the k-mean of
            :return: outputs the k-mean of all the elements in the list
            """
        value = 0
        length = len(vector)
        if k == np.inf:
            for i in vector:
                value = max(value, i)
            return value
        else:
            for i in vector:
                value = value + i ** k
            return ((1. / length) * value) ** (1 / k)

    # calculate Delta_epsilon
    def delta_epsilon(self, theta):
        """
            implementation of the Delta_epsilon function as described in Section 2
            :param theta: given slope
            :return: Delta q_t - theta Delta p_t
            """
        return np.subtract(np.diff(self.q), theta * np.diff(self.p))

    def B_tilde(self):
        """
        calculates B_tilde used in underline{B} for the case k = infty; introduced in Proposition 1
        :return: B_tilde operation in the case of k=infty as defined in Proposition 1
        """
        return

    def underline_B(self):
        """
        case k = infty uses B_tilde as defined above
        case k in (1, infty): implement the definitions of underline{B} as stated in Proposition 1 and Proposition 2

        :return: returns the underline_{B_k} operation as defined in Propositions 1 and 2
        """
        if self.k == np.inf:
            return
        else:
            return

    def Theta_hat_k_B(self, value):
        """
        :param value: value provided
        :return: output [underline{theta}_k, overline{theta}_k] for a given value
            """
        if self.k == np.inf:
            # calculate \underline{\theta_k}
            # calculate \overline{\theta_k}
            # ^^ as defined in Prop. 1
            return
        else:
            # calculate \underline{\theta_k}, \overline{\theta_k} as defined in Prop. 2
            return

    def intervals(self):
        if self.B is None:
            self.B = self.underline_B()
        """
        :param B: upper bound B
        :param ngridpoints: finiteness of grid that user provides
        :return: for each i in the range [underline_B, B], output [underline{theta}_k, overline{theta}_k]
        """
        array_of_intervals = []
        i = self.underline_B()
        while i < self.B:
            array_of_intervals.append(self.Theta_hat_k_B(i))
            i += (self.B - self.underline_B()) / self.ngridpoints
        return array_of_intervals

    def plot(self):
        """
        :return: outputs a plot of B vs. theta similar to Figure 3
        """
        # returns a plot of B vs. \theta similar to Figure 3
        return

def add(a, b):
    """Add two number.
    This return the Addition of these 2 numbers.
    You will never want to use this function because it's just an example.
    :param a: First Number
    :type a: float
    :param b: Second Number
    :type b: float
    :return: The addition of the two number
    :rtype: float
    """
    #logger.debug('add(%s, %s)' % (a, b))
    return a + b


if __name__ == '__main__':
    bounds_data = Bounds(p_t, q_t)
    print(bounds_data.delta_epsilon(5))
