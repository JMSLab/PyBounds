"""Main file of the package"""
# from .log import init_logger

# logger = init_logger('main.log')

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

    @staticmethod
    def kmean(k, vector):
        value = 0
        length = len(vector)
        if k == np.inf:
            return np.max(vector)
        else:
            for i in vector:
                value = value + i ** k
            return ((1 / length) * value) ** (1 / k)

    def delta_epsilon(self, theta):
        return np.subtract(np.diff(self.q), theta * np.diff(self.p))

    def M_hat_k(self, theta):
        return self.kmean(self.k, np.abs(self.delta_epsilon(theta)))

    def B_tilde(self):
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
            return minimize(self.M_hat_k, np.array([0])).fun

    def Theta_hat_k_B(self, B):
        return [self.underline_theta(B), self.overline_theta(B)]

    def intervals(self):
        array_of_intervals = []
        i = self.underline_B()
        while i <= self.maxB:
            array_of_intervals.append(self.Theta_hat_k_B(i))
            i += (self.maxB - self.underline_B()) / self.ngridpoints
        return array_of_intervals

    def plot(self):
        lower_bounds = []
        upper_bounds = []
        x = []
        i = self.underline_B()
        while i <= self.maxB:
            x.append(i)
            lower_bounds.append(self.underline_theta(i))
            upper_bounds.append(self.overline_theta(i))
            i = i + (self.maxB - self.underline_B()) / self.ngridpoints
        fig, (ax, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={
            "width_ratios": [self.underline_B(), self.maxB - self.underline_B()]})
        ax.plot(x, lower_bounds, linewidth=.3)
        ax2.plot(x, lower_bounds, linewidth=.3)
        ax.plot(x, upper_bounds, linewidth=.3)
        ax2.plot(x, upper_bounds, linewidth=.3)

        ax.set_xlim(0, self.underline_B(), auto=True)
        ax2.set_xlim(self.underline_B(), self.maxB, auto=True)

        ax.spines['bottom'].set_linestyle((0, (8, 5)))
        ax.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax.yaxis.tick_left()
        ax2.yaxis.set_ticks_position('none')
        plt.subplots_adjust(wspace=0)

        fig.text(0.5, 0.04, 'Bound B', ha='center', va='center')
        ax.set_ylabel('Slope $\\theta$')
        plt.savefig('Bounds Interval Plot k=' + str(self.k) + ".png")
        plt.show()
