import copy


# Generate all vectors of weight less or equal than target
def low_weight(length, weight, weight_set):
    if weight > length:
        print("Impossible to have a weight bigger than length")
        return
    if weight == 0:
        return weight_set
    res = {element[:i] + "1" + element[i + 1:]
           for i in range(length) for element in weight_set}.union(weight_set)
    return low_weight(length, weight - 1, res)


# Class simulating Vectors over F2 and some operations.
class VectorF2(list):

    @staticmethod
    def low_weight_vectors(length, weight):
        res = low_weight(length, weight, {"0" * length})
        return [[int(i) for i in word] for word in res]

    # Sum of two vectors pointwise.
    def __add__(self, other):
        res = [i ^ j for (i, j) in zip(self, other)]
        return VectorF2(res)

    # Scalar product of two vectors.
    def __mul__(self, other):
        res = [i and j for (i, j) in zip(self, other)]
        return sum(res) % 2

    # Negation of vector: Vector XoR [1, ..., 1].
    def negation(self):
        return VectorF2([i ^ 1 for i in self])

    # Walsh Transform in recursive way.
    def walsh_transform(self, truth_table):
        res = []
        if len(truth_table) == 1:
            return truth_table[0]
        for i in range(len(truth_table) // 2):
            res.append([a + b for (a, b) in zip(truth_table[2 * i], truth_table[(2 * i) + 1])] +
                       [a - b for (a, b) in zip(truth_table[2 * i], truth_table[(2 * i) + 1])])
        return self.WalshTransform(res)

    # Walsh Spectrum
    def walsh_spectrum(self):
        # Convert truth table 0 -> 1, 1 -> -1
        return self.WalshTransform([[(-1) ** i] for i in self])

    # Shift vector on the right inserting 0.
    def left_shift(self, shift):
        for i in range(shift):
            self.append(0)
            self.pop(0)
        return self

    def right_shift(self, shift):
        for i in range(shift):
            self = [0] + self
            self.pop(len(self) - 1)
        return self


# Class representing matrices in F_2.
class MatrixF2(list):

    @staticmethod
    def identity(dimension):
        values = []
        for i in range(dimension):
            temp = [0] * dimension
            temp[i] = 1
            values.append(VectorF2(temp))
        return MatrixF2(dimension, dimension, values)

    # Matrix constructor.
    def __init__(self, rows, cols, values=[[]]):
        super(MatrixF2, self).__init__((VectorF2(v) for v in values))
        self.row_num = rows
        self.col_num = cols

    # Add two matrices component by component.
    def __add__(self, other):
        res = [i + j for (i, j) in zip(self, other)]
        return MatrixF2(self.row_num, self.col_num, res)

    # Get the column at the index index.
    def get_column(self, index):
        res = [i[index] for i in self]
        return VectorF2(res)

    # Get the row at the index index.
    def get_row(self, index):
        return self[index]

    # Transpose the matrix.
    def transpose(self):
        return MatrixF2(self.col_num, self.row_num, [self.get_column(i) for i in range(self.col_num)])

    # Insert column
    def insert_column(self, column, index):
        self.col_num += 1
        [row.insert(index, column[i]) for i, row in enumerate(self)]
        return self

    # Insert row at index index.
    def insert_row(self, row, index):
        self.row_num += 1
        self.insert(index, VectorF2(row))
        return self

    # Matrices multiplication.
    def __mul__(self, other):
        values = []
        row = []
        trans_other = other.transpose()
        for i in self:
            values.append(VectorF2([i * j for j in trans_other]))
        return MatrixF2(self.row_num, other.col_num, values)

    # Exponentiation of a matrix ( ** operator)
    def __pow__(self, power, modulo=None):
        if self.row_num != self.col_num:
            print("Power defined only for square matrices.")
            return
        if power == 0:
            return MatrixF2.identity(self.row_num)
        return self * self.__pow__(power - 1)

    # Reduce matrix over F2 to echelon form.
    def row_echelon(self):
        t = 0
        r = 0
        while (t < self.col_num - 1) and (r < self.row_num - 1):
            # Check if r^th row does not have a 1 in position t
            if self[r][t] == 0:
                for i in range(r, self.row_num):
                    if self[i][t] == 1:
                        # Swap if i^th row does have a 1 in position t.
                        self[r], self[i] = self[i], self[r]
                        break
            # If row r still has a 0 in position t, increment t
            if self[r][t] == 0:
                t = t + 1
            else:
                # If row r has a 1 in position t use it to reduce the following rows.
                for i in range(r + 1, self.row_num):
                    # If a row after the selected has a 1 in position t reduce using the r^th row.
                    if self[i][t] == 1:
                        self[i] = self[i] + self[r]
                r = r + 1
                t = t + 1
        return self

    # Solve A x = b.
    # Just find a solution (all free variable are put to 0).
    # TODO: find col_num - rank l.i. solutions.
    def solve_right(self, b, verbose=False):
        A = copy.deepcopy(self)
        # Initialize res with zero. (all free variable are put to zero)
        res = [0] * self.col_num
        A.insert_column(VectorF2(b), A.col_num)
        A.row_echelon()
        # Keep track of echelon steps.
        step = []
        for i in A:
            for j in range(A.col_num):
                if i[j] != 0:
                    step.append(j)
                    break
        rank = len(step)
        # Check if solvable
        if (step[-1] == A.col_num - 1):
            if verbose:
                print("System has no solution.")
            return VectorF2([])
        else:
            if verbose:
                print("Solution space has degree: ", self.col_num - rank)
        # Start from the bottom going up.
        # Index rank - 1 is the first row non null.
        for i in range(rank - 1, -1, -1):
            # looking for the first non 0 element in the row
            # Non free variable
            res[step[i]] = (
                    A[i][A.col_num - 1] ^
                    VectorF2(res[step[i] + 1:]) * VectorF2(A[i][step[i] + 1:])
            )
        return VectorF2(res)

    @staticmethod
    def shift_matrix(vector, col_num):
        first_row = VectorF2(vector + [0] * (col_num - len(vector)))
        rows = [copy.deepcopy(first_row).right_shift(i) for i in range(col_num - len(vector) + 1)]
        return MatrixF2(len(rows), col_num, rows)

    def __str__(self):
        res = ""
        for i in self:
            res = res + str(i) + "\n "
        return res


def LinearAlgebraF2_test():
    ## Just some print test.
    A = MatrixF2(2, 3, [[0, 1, 0], [1, 1, 0]])
    print("A:\n", A)

    B = MatrixF2(3, 3, [[1, 1, 0], [1, 0, 0], [0, 1, 0]])
    print("B:\n", B)

    print("A*B:\n", A * B)

