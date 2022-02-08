"""Main file of the package"""
# from .log import init_logger

# logger = init_logger('main.log')
# NEEDS TESTING!!!!!

import numpy as np


class Bounds:
    def __init__(self, p_t, q_t, k=np.inf, maxB=None, ngridpoints=1000):
        self.p = p_t
        self.q = q_t
        if k < 1:
            raise Exception("Requires k>=1")
        else:
            self.k = k
        if maxB is None:
            self.maxB = 2 * self.underline_B()
        else:
            self.maxB = maxB
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
        case k = infinity uses B_tilde as defined above
        case k in (1, infinity): implement the definitions of underline{B} as stated in Proposition 1 and Proposition 2

        :return: returns the underline_{B_k} operation as defined in Propositions 1 and 2
        """
        if self.k == np.inf:
            return 0
            # temporary placeholder to ensure program runs
        else:
            return 0
            # temporary placeholder to ensure program runs

    def Theta_hat_k_B(self, B):
        """
        :param B: upper bound B provided
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
        """
        :return: for each i in the range [underline_B, maxB], output [underline{theta}_k, overline{theta}_k]
        """
        array_of_intervals = []
        i = self.underline_B()
        while i < self.maxB:
            array_of_intervals.append(self.Theta_hat_k_B(i))
            i += (self.maxB - self.underline_B()) / self.ngridpoints
        return array_of_intervals

    def plot(self):
        """
        :return: outputs a plot of B vs. theta similar to Figure 3
        """
        # returns a plot of B vs. \theta similar to Figure 3
        return
