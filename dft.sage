import numpy as np

load("dft_utils.sage")

class DFT: # a module to perform the DFT in the clear    
    def __init__(self, N) -> None:
        assert N & (N - 1) == 0, "N must be a power of 2"
        self.N = N
        self.log_N = ZZ(log(N, 2))
        self.br_lookup = br_lookup(self.log_N)
        self.precomp()
        
    def encode(self, v, input_order: str = "n"):
        """
        Encode a complex vector using the inverse Fast Discrete Fourier Transform (DFT) in O(log N) steps.
        
        This method implements Algorithm 6 from the paper, performing a forward DFT
        on a complex vector of length N/2 and returning a real vector of length N.
        The algorithm uses a butterfly structure with precomputed roots of unity indexed by 5**j
        for efficient computation.
        
        Parameters:
        -----------
        v : numpy.ndarray
            Input complex vector of length N/2 to be encoded
        input_order : str, optional
            Input ordering: "n" for normal order (default), "b" for bit-reversed order
            The fast variant of the iDFT bit-reverses the input automatically.
            Thus, by default, the input is assumed to be in normal order,
            and will get bit-reversed before applying the iDFT, to retain the normal order output.
            
        Returns:
        --------
        numpy.ndarray
            Real vector of length N containing the encoded DFT coefficients.
            The first N/2 elements are the real parts, the last N/2 are the imaginary parts.
            
        Algorithm Description:
        ---------------------
        1. Input validation and bit-reversal (if input_order="n")
        2. Iterative butterfly operations over log₂(N)-1 levels:
           - Reshape vector into (N/delta, delta/2) matrix
           - Split into even (a) and odd (b) components
           - Apply butterfly: a' = (a+b)/2, b' = (a-b)*ω/2
           - Update delta and root index for next level
        3. Reshape to 1D and separate real/imaginary parts
        
        Notes:
        ------
        - Requires N to be a power of 2
        - Precomputation must be called before encoding for efficiency
        """
                
        assert type(v) == np.ndarray, "Input vector must be a numpy array"
        assert len(v) == self.N // 2, "Input vector has the wrong length {}".format(len(v))
        assert hasattr(self, 'sequence'), "Precomputation not done, call precomp() first"

        if input_order == "n": # We need to bit-reverse manually
            v = v[self.br_lookup[-2]]
        else:
            assert input_order == "b", "Input order must be either 'n' or 'b'"
            
        delta, index = 2, 1
        for _ in range(self.log_N - 1):
            v = np.reshape(v, (self.N // delta, delta // 2))
            a, b = v[::2], v[1::2]
            # The line below must be done in one go, otherwise we overwrite incorrectly
            v[::2], v[1::2] = (a + b), (a - b) * self.sequence_inv[index]
            v /= 2
            delta *= 2
            index += 1
        v = np.reshape(v, self.N // 2)
        return np.append(np.real(v), np.imag(v))
    
    def decode(self, v, output_order: str = "n"): # Algorithm 5 in the paper
        """ 
        The inverse of the above encoding method, aka the DFT. 
        """        
        assert type(v) == np.ndarray, "Input vector must be a numpy array"
        assert len(v) == self.N, "Input vector has the wrong length {}".format(len(v))
        assert hasattr(self, 'sequence'), "Precomputation not done, call precomp() first"

        v = v[:self.N//2] + v[self.N//2:] * 1j # Undoing the real() and imag() parts 
        
        delta, index = 2 ** (self.log_N), self.log_N
        for _ in range(self.log_N - 2, -1, -1):
            v = np.reshape(v, (2 * self.N // delta, delta // 4))
            delta //= 2
            index -= 1
            # We do not need to multiply by 2 here, since we get it for free during e.g. a + ub below!
            a = v[::2]
            ub = self.sequence[index] * v[1::2]
            v[::2], v[1::2] = a + ub, a - ub
        
        if output_order == "n":
            v = v[self.br_lookup[-2]]
        else:
            assert output_order == "b", "Output order must be either 'n' or 'b'"
            
        return np.reshape(v, self.N//2)
    
    def precomp(self):
        """
        Precomputes the sequence(_inv) of roots of unity for the (inverse) DFT.
        Corresponds to all values of u resp. u^{-1} in the paper.
        """
        delta, list_roots = 2, []
        for k in range(self.log_N - 1, 0, -1):
            for i in range(self.N // 2 // delta):
                u = root(2*self.N, (5 ** self.br_lookup[k-1][i]) * self.N // (2**(k+1)))
                list_roots.append(u)
            delta *= 2
        sequence_inv, sequence = [[]], [[]]
        index, delta = 0, 2
        for _ in range(self.log_N - 1):
            step = self.N // 2 // delta
            A = [[list_roots[i] ** (-1)] for i in range(index, index + step)]
            B = [[list_roots[i]] for i in range(index, index + step)]
            sequence_inv.append(np.array(A))
            sequence.append(np.array(B))
            index += step
            delta *= 2
        self.sequence_inv = sequence_inv
        self.sequence = sequence
