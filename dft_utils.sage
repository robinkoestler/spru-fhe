def root(index, exponent = None, precision: int = 53):
    # Calculate the primitive index-th root of unity. Using SageMath's zeta function for efficiency.
    if exponent:
        exponent %= index
        return ComplexField(precision).zeta(index) ** exponent
    return ComplexField(precision).zeta(index)

def br(binary_number, length): # bit-reversal of a binary number of length length
    return int('{:0{width}b}'.format(binary_number, width=length)[::-1], 2)

def br_lookup(log_N): # generates a look-up table for the bit-reversal permutation
    return [np.array([br(i, l) for i in range(2**l)], dtype=int) for l in range(log_N + 1)]
