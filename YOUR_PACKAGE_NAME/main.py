"""Main file of the package"""
#from .log import init_logger

#logger = init_logger('main.log')
#NEEDS TESTING!!!!!


import pandas as pd
import numpy as np
# import data: p_t, q_t
data = pd.read_csv('roberts_schlenker_2013.csv') # change dataset if needed to test
df = pd.DataFrame(data)
p_t = df.iloc[:, 0]
q_t = df.iloc[:, 1]

#determines/reads k
k = input("Enter k value: ")
if k == "":
    #default is -1, which represents the infinity norm
    k = -1
else:
    k = int(k)

#implementation of k-mean
def kmean(lst1, k):
    """implementation of the k-mean function
    :param lst1: list of elements of which to take the k-mean of
    :param k: the "k" in k-mean
    :return: outputs the k-mean of all the elements in the list
    """
    value = 0
    length = len(lst1)
    if k == -1:
        for i in lst1:
            value = max(value, i)
        return value
    else:
        for i in lst1:
            value = value+i^k
        return ((1./length)* value)**(1/k)

#print(kmean(p_t,k))

#print(np.abs(p_t))

#calculate Delta_epsilon
def delta_epsilon(theta, p_t, q_t):
    """
    implementation of the Delta_epsilon function as described on pg. 5
    :param theta: given slope
    :param p_t: first set of data that user inputs/the p_t variable on pg. 5
    :param q_t: second set of data that user inputs/the q_t variable on pg. 5
    :return: Delta q_t - theta Delta p_t
    """
    return np.subtract(np.diff(q_t), theta * np.diff(p_t))
#print(delta_epsilon(2.,p_t,q_t))

def B_tilde(p_t, q_t):
    """
    calculates B_tilde used in underline{B} for the case k = infty; introduced in Proposition 1
    :param p_t: first set of data that user provides
    :param q_t: second set of data that user provides
    :return: B_tilde operation in the case of k=infty as defined on pg. 6
    """
    return

def underline_B(k, p_t, q_t):

    """
    case k = infty uses B_tilde as defined above
    case k in (1, infty): implement the definitions of underline{B} as stated in Proposition 1 and Proposition 2

    :param k: Value k>=1 that user provides
    :param p_t: first set of data that user provides
    :param q_t: second set of data that user provides
    :return: returns the underline_{B_k} operation as defined on pg. 6-7
    """
    if k == -1:
        return
    else:
        return

#determines B
B = int(input("Enter B value: "))

grid_finiteness = int(input("Enter finiteness: "))

def intervals(k, p_t, q_t, B, grid_finiteness):
    """

    :param k: value k>=1 that user provides
    :param p_t: first set of data that user provides
    :param q_t: second set of data that user provides
    :param B: upper bound B
    :param grid_finiteness: finiteness of grid that user provides
    :return: for each i in the range [underline_B, B], output [underline{theta}_k, overline{theta}_k] as defined on pg. 6-7
    """
    if k == -1:
        for i in range(underline_B(k, p_t, q_t), B, (B-underline_B(k, p_t, q_t)) / grid_finiteness):
            #calculate \underline{\theta_k}
            #calculate \overline{\theta_k}
            #^^ as defined in Prop. 1
            return
    else:
        # calculate \underline{\theta_k}, \overline{\theta_k} as defined in Prop. 2
        return


def plot(k, B, p_t, q_t):
    """
    :param k: value k>=1 that user provides
    :param B: upper bound B that user provides
    :param p_t: first set of data that user provides
    :param q_t: second set of data that user provides
    :return: outputs a plot of B vs. theta similar to Figure 3
    """
    # returns a plot of B vs. \theta similar to Figure 3
    return

if __name__ == '__main__':
    print(underline_B(k, p_t, q_t))
    print(intervals(k, p_t, q_t, B, grid_finiteness))
    plot(k, B, p_t, q_t)