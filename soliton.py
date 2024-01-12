import math

def robust_soliton_distribution(K: int, c: float, delta: float):
    # Calculate parameter 's'
    s = c * math.sqrt(K * math.log(K / delta) ** 2)

    # Define tau(d)
    tau = [1 / K if d == 1 else 1 / (d * (d - 1)) for d in range(1, K + 1)]

    # Define rho(d)
    rho = [s / (K * d) if d <= math.floor(K / s) - 1 else 0 for d in range(1, K + 1)]
    rho[math.floor(K / s) - 1] = s * math.log(s / delta) / K

    # Combine to form mu(d)
    mu = [tau_d + rho_d for tau_d, rho_d in zip(tau, rho)]
    Z = sum(mu)
    normalized_mu = [mu_d / Z for mu_d in mu]

    return normalized_mu