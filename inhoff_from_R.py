import math
import numpy as np

def theta(u, lambda_, lambdalen, h, x, delta2):
    m = lambdalen[0]
    sum = 0
    for i in range(m):
        sum += h[i] * math.atan(lambda_[i] * u[0]) + delta2[i] * lambda_[i] * u[0] / (1 + lambda_[i] ** 2 * u[0] ** 2)
    sum = 0.5 * sum - 0.5 * x[0] * u[0]
    return sum

def rho(u, lambda_, lambdalen, h, delta2):
    m = lambdalen[0]
    prod = 1
    for i in range(m):
        prod *= (1 + lambda_[i] ** 2 * u[0] ** 2) ** (0.25 * h[i]) * math.exp(0.5 * delta2[i] * lambda_[i] ** 2 * u[0] ** 2 / (1 + lambda_[i] ** 2 * u[0] ** 2))
    return prod

def imhoffunc(u, lambda_, lambdalen, h, x, delta2):
    res = math.sin(theta(u, lambda_, lambdalen, h, x, delta2)) / (u[0] * rho(u, lambda_, lambdalen, h, delta2))
    return res

def f(x, n, ex):
    xx = np.array([ex[0]])
    lambdalen = np.array([int(ex[1])])
    lambda_ = np.array(ex[2:2 + lambdalen[0]])
    h = np.array(ex[2 + lambdalen[0]:2 + 2 * lambdalen[0]])
    delta2 = np.array(ex[2 + 2 * lambdalen[0]:])
    u = np.array([x[i] for i in range(n)])
    for i in range(n):
        u[i] = x[i]
        x[i] = imhoffunc(u[i:i+1], lambda_, lambdalen, h, xx, delta2)

if __name__ == '__main__':
    print(imhoffunc(2,[]))