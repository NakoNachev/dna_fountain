import math

def ideal_soliton_distribution(K):
    """ this is the ideal solution distribution
    Currently not used"""
    return [1.0 / K if d == 1 else 1.0 / (d * (d - 1)) for d in range(1, K + 1)]

def robust_soliton_distribution(K: int, c: float, delta: float):
    """ this is the robust soliton distribution according to the paper
    Calculations are based on the provided formulas in the extra materials """
    s = c * math.sqrt(K * math.log(K / delta) ** 2)
    rho = [1 / K if d == 1 else 1 / (d * (d - 1)) for d in range(1, K + 1)]
    tau = [s / (K * d) if d <= math.floor(K / s) - 1 else 0 for d in range(1, K + 1)]
    index = math.floor(K / s) - 1
    if index <= K - 1:   # otherwise fails with index out of range
        tau[math.floor(K / s) - 1] = s * math.log(s / delta) / K

    #Combine to form mu(d)
    mu = [tau_d + rho_d for tau_d, rho_d in zip(tau, rho)]
    Z = sum(mu)
    normalized_mu = [mu_d / Z for mu_d in mu]

    return normalized_mu
