import matplotlib.pyplot as plt

import math

def ideal_soliton_distribution(K):
    """ Generate the ideal soliton distribution for a given K. """
    return [1.0 / K if d == 1 else 1.0 / (d * (d - 1)) for d in range(1, K + 1)]

def robust_soliton_distribution(K, c, delta):
    """ Generate the robust soliton distribution for given K, c, and delta. """
    R = c * math.sqrt(K * math.log(K / delta))
    tau = [0] * K
    for i in range(1, K):
        if i < R:
            tau[i - 1] = 1 / (i * K)
        else:
            tau[i - 1] = math.log(R / delta) / K
            break

    # Combine ideal and tau to form robust distribution
    mu = [d + t for d, t in zip(ideal_soliton_distribution(K), tau)]
    Z = sum(mu)
    return [m / Z for m in mu]


def test_robust_soliton_distribution():
    K = 553
    c = 0.9
    delta = 0.001
    distribution = robust_soliton_distribution(K, c, delta)

    plt.bar(range(1, K + 1), distribution)
    plt.xlabel('Segment Size')
    plt.ylabel('Probability')
    plt.title('Robust Soliton Distribution')
    plt.show()

test_robust_soliton_distribution()
