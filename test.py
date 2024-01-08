class GaloisLFSR:
    def __init__(self, seed, polynomial):
        self.state = seed
        self.polynomial = polynomial

    def shift(self):
        feedback_bit = self.state & 1
        self.state >>= 1
        if feedback_bit:
            self.state ^= self.polynomial

    def generate_sequence(self, count):
        sequence = []
        for _ in range(count):
            sequence.append(self.state)
            self.shift()
        return sequence


def generate_seeds(seed, polynomial, count):
    lfsr = GaloisLFSR(seed, polynomial)
    return lfsr.generate_sequence(count)


# Example usage
initial_seed = 42
lfsr_polynomial = 0b10000000000000000000001011010111  # The given primitive polynomial
seed_count = 10

seeds = generate_seeds(initial_seed, lfsr_polynomial, seed_count)
print("Generated Seeds:", seeds)
