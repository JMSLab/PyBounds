"""Main file of the package"""
# from .log import init_logger

# logger = init_logger('main.log')
# NEEDS TESTING!!!!!

import numpy as np
import scipy as scipy
from scipy.optimize import minimize
from matplotlib import pyplot as plt


class Bounds:
    def __init__(self, p_t, q_t, k=np.inf, maxB=None, ngridpoints=1000):
        self.p = p_t
        self.q = q_t
        self.ngridpoints = ngridpoints
        if k < 1:
            raise Exception("Requires k>=1")
        else:
            self.k = k
        if maxB is None:
            self.maxB = 2 * self.underline_B()
        else:
            self.maxB = maxB

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
            return np.max(vector)
        else:
            for i in vector:
                value = value + i ** k
            return ((1 / length) * value) ** (1 / k)

    # calculate Delta_epsilon
    def delta_epsilon(self, theta):
        """
            implementation of the Delta_epsilon function as described in Section 2
            :param theta: given slope
            :return: Delta q_t - theta Delta p_t
            """
        return np.subtract(np.diff(self.q), theta * np.diff(self.p))

    def M_hat_k(self, theta):
        """
        implementation of tilde{M_k}
        used in underline_B
        :return: the k-mean of delta-epsilon
        """
        return self.kmean(self.k, np.abs(self.delta_epsilon(theta)))

    def B_tilde(self):
        """
        calculates B_tilde used in underline{B} for the case k = infty; introduced in Proposition 1
        :return: B_tilde operation in the case of k=infty as defined in Proposition 1
        """
        def f(B):
            return self.overline_theta(B) - self.underline_theta(B)
        return scipy.optimize.root_scalar(f, x0=0, x1=0.001).root

    def underline_theta(self, B):
        if self.k == np.inf:
            theta_underline_values = []
            for i in range(1, len(self.p)):
                if self.p[i] - self.p[i - 1] != 0:
                    delta_q_t = self.q[i] - self.q[i - 1]
                    delta_p_t = self.p[i] - self.p[i - 1]
                    theta_underline_values.append(delta_q_t / delta_p_t - B / (abs(delta_p_t)))
            return max(theta_underline_values)

        else:
            def f(theta):
                return self.M_hat_k(theta) - B

            sol = scipy.optimize.root_scalar(f, x0=-2, x1=-1)
            return sol.root

    def overline_theta(self, B):
        if self.k == np.inf:
            theta_overline_values = []
            for i in range(1, len(self.p)):
                if self.p[i] - self.p[i - 1] != 0:
                    delta_q_t = self.q[i] - self.q[i - 1]
                    delta_p_t = self.p[i] - self.p[i - 1]
                    theta_overline_values.append(delta_q_t / delta_p_t + B / (abs(delta_p_t)))
            return min(theta_overline_values)
        else:
            def f(theta):
                return self.M_hat_k(theta) - B

            sol = scipy.optimize.root_scalar(f, x0=2, x1=1)
            return sol.root

    def underline_B(self):
        """
        case k = infinity uses B_tilde as defined above
        case k in (1, infinity): implement the definitions of underline{B} as stated in Proposition 1 and Proposition 2
        :return: returns the underline_{B_k} operation as defined in Propositions 1 and 2
        """
        if self.k == np.inf:
            delta_q = []
            for i in range(1, len(self.p)):
                if self.p[i] - self.p[i - 1] == 0:
                    delta_q.append(self.q[i] - self.q[i - 1])
            if not delta_q:
                return self.B_tilde()
            else:
                return max(max(delta_q), self.B_tilde())
        else:
            # print(minimize(self.M_tilde_k, 0))
            return minimize(self.M_hat_k, np.array([0])).fun

    def Theta_hat_k_B(self, B):
        """
        :param B: upper bound B provided
        :return: output [underline{theta}_k, overline{theta}_k] for a given value
            """
        # print(self.underline_theta(B))
        return [self.underline_theta(B), self.overline_theta(B)]

    def intervals(self):
        """
        :return: for each i in the range [underline_B, maxB], output [underline{theta}_k, overline{theta}_k]
        """
        array_of_intervals = []
        i = self.underline_B()
        while i <= self.maxB:
            array_of_intervals.append(self.Theta_hat_k_B(i))
            i += (self.maxB - self.underline_B()) / self.ngridpoints
        # print(len(array_of_intervals))
        return array_of_intervals

    def plot(self):
        """
        :return: outputs a plot of B vs. theta similar to Figure 3
        """
        intervals = self.intervals()
        lower_bounds = []
        upper_bounds = []
        x = []
        i = self.underline_B()
        while i <= self.maxB:
            x.append(i)
            lower_bounds.append(self.underline_theta(i))
            upper_bounds.append(self.overline_theta(i))
            i = i + (self.maxB-self.underline_B())/self.ngridpoints
        plt.plot(x, lower_bounds, linewidth=.3)
        plt.xlim([0, self.maxB])
        plt.plot(x, upper_bounds, linewidth=.3)
        plt.ylim(top=0)
        #plt.axhline(y=0, color='r', linestyle='-')
        plt.savefig('Bounds Interval Plot k=' + str(self.k) + ".png")
        plt.show()
