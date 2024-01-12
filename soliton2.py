import numpy as np
EPSILON = 0.0001
ROBUST_FAILURE_PROBABILITY = 0.01
import math

# code reference: https://github.com/Spriteware/lt-codes-python/blob/master/distributions.py

def ideal_distribution(N):
    """ Create the ideal soliton distribution. 
    In practice, this distribution gives not the best results
    Cf. https://en.wikipedia.org/wiki/Soliton_distribution
    """

    probabilities = [1 / N]
    probabilities += [1 / (k * (k - 1)) for k in range(2, N+1)]
    probabilities_sum = sum(probabilities)

    assert probabilities_sum >= 1 - EPSILON and probabilities_sum <= 1 + EPSILON, "The ideal distribution should be standardized"
    return probabilities


def robust_distribution(N):
    """ Create the robust soliton distribution. 
    This fixes the problems of the ideal distribution
    Cf. https://en.wikipedia.org/wiki/Soliton_distribution
    """

    # The choice of M is not a part of the distribution ; it may be improved
    # We take the median and add +1 to avoid possible division by zero 
    M = N // 2 + 1 
    R = N / M

    extra_proba = [0] + [1 / (i * M) for i in range(1, M)]
    extra_proba += [math.log(R / ROBUST_FAILURE_PROBABILITY) / M]  # Spike at M
    extra_proba += [0 for k in range(M+1, N+1)]

    probabilities = np.add(extra_proba, ideal_distribution(N))
    probabilities /= np.sum(probabilities)
    probabilities_sum = np.sum(probabilities)
    print(f'probabilities len is {probabilities}')

    assert probabilities_sum >= 1 - EPSILON and probabilities_sum <= 1 + EPSILON, "The robust distribution should be standardized"
    return probabilities




def sample_elements(array_to_sample, probabilities):
    return np.random.choice(array_to_sample, p=probabilities)

# Example usage:
array_to_sample = np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, 5], [1, 2, 3, 4, 5], [1, 2, 3, 4, 5]])
N = len(array_to_sample)
print(f'N is {N}')

# Using ideal distribution
ideal_probabilities = ideal_distribution(N)
print(f'ideal {ideal_probabilities}')
sample_ideal = sample_elements(array_to_sample, ideal_probabilities)
print("Ideal Sample:", sample_ideal)

# Using robust distribution
robust_probabilities = robust_distribution(N)
print(f'robust_probabilities {robust_probabilities}')
sample_robust = sample_elements(array_to_sample, robust_probabilities)
print("Robust Sample:", sample_robust)
