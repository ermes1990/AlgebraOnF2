from LinearAlgebraF2 import *
from PolynomialF2 import *


class LFSR:
    # Initialize with polynomial and initial state.
    # Polynomial is reppresented as binary vector.
    # Example: X^5 + X^4 + X^3 + 1 -> [1, 0 ,0, 1, 1, 1]
    def __init__(self, polynomial, initial_state):
        self.dimension = len(initial_state)
        # Convert initial state to a 1 x n matrix.
        self.state = MatrixF2(1, self.dimension, [initial_state])

        # Create companion matrix.
        C = MatrixF2.identity(self.dimension - 1)
        C.insert_row(VectorF2([0] * (self.dimension - 1)), self.dimension)
        C.insert_column(polynomial[-2::-1], 0)
        self.companion_matrix = C

    # Clock LFSR by one.
    def Clock(self):
        res = self.state[0][-1]
        self.state = self.state * self.companion_matrix
        return res

    # Generate some stream.
    def Stream(self, stream_length):
        stream = []
        for i in range(stream_length):
            stream = stream + [self.Clock()]
        return stream


# Class simulating an LFSR with filter:
# Combining an LFSR and a boolean function.
class F_LFSR:
    def __init__(self, lfsr, boolean_function, filter_input):
        self.lfsr = lfsr
        self.lfsr_filter = boolean_function
        self.filter_input_index = filter_input

    def Clock(self):
        # Select filter input from the lfsr cells.
        filter_input = [self.lfsr.state[0][i] for i in self.filter_input_index]
        # Evaluate the filter with the correct input.
        res = self.lfsr_filter.Evaluate(filter_input)
        # Clock the lfsr.
        self.lfsr.Clock()
        # Return the evaluated value.
        return res

    # Generate some stream.
    def Stream(self, stream_length):
        stream = []
        for i in range(stream_length):
            stream = stream + [self.Clock()]
        return stream


def LFSR_test():
    formula = [{0, 1}, {2}, {3}, {4}]
    boolean_filter = PolynomialF2(formula, 5)

    # It represent an LFSR with 41 cells and polynomial X^41 + X^3 + 1
    poly = [0] * 42
    poly[41] = 1
    poly[3] = 1
    poly[0] = 1

    initial_state = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 1, 0, 1, 0, 0, 0, 0, 0, 0]

    lfsr = LFSR(poly, initial_state)
    filter_lfsr = F_LFSR(lfsr, boolean_filter, [0, 1, 5, 6, 40])



    print(filter_lfsr.Stream(25) == [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0])
