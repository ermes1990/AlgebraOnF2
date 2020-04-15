from math import comb


# Get the j-th vector of length n and weight w.
def get_vector_by_index(n, w, j):
    if n == 0:
        return ""
    if w == 0:
        return "0" * n
    limit = comb(n - 1, w - 1)
    if j <= limit - 1:
        return "1" + get_vector_by_index(n - 1, w - 1, j)
    else:
        return "0" + get_vector_by_index(n - 1, w, j - limit)


# Get the j-th vector of length n independent from weight.
# First weight 0, then weight 1 ...
def get_vector_by_abs_index(n, j):
    # Find weight.
    for weight in range(n + 1):
        binom = comb(n, weight)
        if j < binom:
            return get_vector_by_index(n, weight, j)
        j = j - binom


def get_abs_index(vector):
    n = len(vector)
    w = vector.count("1")
    res = 0
    for i in range(w):
        res = res + comb(n, i)
    for i in vector:
        if i == "0":
            res = res + comb(n-1, w-1)
            n = n - 1
        elif i == "1":
            n = n - 1
            w = w - 1
        if n == 0 or w == 0:
            return res


class MonomialF2(set):
    def __init__(self, monomial, order="DEGLEX"):
        self.order = order
        super(MonomialF2, self).__init__(monomial)

    def __str__(monomial):
        # Empty set is the constant monomial 1.
        if monomial == set():
            return "1"
        res = ""
        for i in sorted(list(monomial)):
            res = res + "X" + str(i) + "*"
        return res[:-1:]

    def __mul__(self, other):
        return self.union(other)

    def __lt__(self, other):
        # If equal other is not greater.
        if self == other:
            return False
        # USe DEGLEX order
        if self.order == "DEGLEX":
            # First look at deg.
            if len(self) < len(other):
                return True
            elif len(self) > len(other):
                return False
            # Sort and Look at lexicographic order.
            else:
                first = sorted(list(self))
                second = sorted(list(other))
                for i in range(len(first) - 1, -1, -1):
                    if first[i] < second[i]:
                        return True
                    elif first[i] > second[i]:
                        return False

    # Evaluate monomial
    # Variable is a tuple of 0 and 1
    def Evaluate(self, variable):
        # Empty set is the constant monomial 1
        res = 1
        # Monomial is a set of indexes
        for x in self:
            res = res * variable[x]
        return res

    @staticmethod
    # Look at monomial as binary representation:
    # X_0 * X_2 * X_3 in 5 variables is represented as: "100110"
    # These vectors can be sorted in a deg_lex way.
    def index_to_monomial(variable_num, index):
        return MonomialF2({i for i, val in enumerate(get_vector_by_abs_index(variable_num, index)) if val == "1"})

    def monomial_to_index(self, variable_num):
        if self == set():
            return 0
        if max(self) >= variable_num:
            print("the monomial cannot belong to this polinomial.")
        monomial_string = list("0" * variable_num)
        for i in self:
            monomial_string[i] = "1"
        return get_abs_index("".join(monomial_string))

# Class representing polynomial in F2 with more variables.
class PolynomialF2(list):
    # Coefficient:
    # The empty set() is the constant 1
    # Example: [{1,2}, {1}, {0}, set() ] -> x2x3 + x2 + +x1 + 1
    def __init__(self, coefficient, variable_num):
        self.variable_num = variable_num
        # Reduce module 2.
        unique_coefficient = set([frozenset(i) for i in coefficient])
        reduced_coefficient = [i for i in unique_coefficient if coefficient.count(i) % 2 == 1]
        self.degree = max([len(i) for i in reduced_coefficient])
        # Initialize polynomial.
        super(PolynomialF2, self).__init__([MonomialF2(i) for i in reduced_coefficient])

    # Convert the polynomial in a vector.
    def to_vector(self):
        # Empty vector
        res = [0] * sum([comb(self.variable_num, i) for i in range(self.degree + 1)])
        # Place a 1 for each monomial in the right place.
        active_index = [i.monomial_to_index(self.variable_num) for i in self]
        for monomial in self:
            res[monomial.monomial_to_index(self.variable_num)] = 1
        return res

    # TODO Not implemented
    @staticmethod
    def from_vector(vector, variable_num):
        monomial_list = []
        for i, val in enumerate(vector):
            if val == "1":
                monomial_list.append(MonomialF2.index_to_monomial(variable_num, i))
        return PolynomialF2(monomial_list, variable_num)

    # Add two polynomials in F2.
    def __add__(self, other):
        # Cast to list so '+' is used as list concatenation.
        right = set([frozenset(i) for i in self])
        left = set([frozenset(i) for i in other])
        res = right.union(left).difference(right.intersection(left))
        return PolynomialF2(list(res), max(self.variable_num, other.variable_num))

    # Polynomial multiplication in F2.
    def __mul__(self, other):
        temp = []
        # Multiply each monomial in self with each one in other.
        for i in self:
            for j in other:
                temp.append(i * j)
        # We work in characteristic 2, delete every 2 same monomials.
        return PolynomialF2(temp, max(self.variable_num, other.variable_num))

    # Convert the polynomial to a string to print
    def __str__(self):
        expression = ""
        for i in sorted(self, reverse=True):
            expression = expression + str(i) + " + "
        return expression[:-3:]

    # Evaluate polynomial
    def Evaluate(self, variable):
        res = 0
        for monomial in self:
            res = res ^ monomial.Evaluate(variable)
        return res

    # Get truth table
    def TruthTable(self):
        res = []
        for i in range(2 ** self.variable_num):
            res.append(self.Evaluate([int(j) for j in bin(i)[2:].zfill(self.variable_num)]))
        return res

    # Walsh Transform in recursive way
    def WalshTransform(self, truth_table):
        res = []
        if len(truth_table) == 1:
            return truth_table[0]
        for i in range(len(truth_table) // 2):
            res.append([a + b for (a, b) in zip(truth_table[2 * i], truth_table[(2 * i) + 1])] +
                       [a - b for (a, b) in zip(truth_table[2 * i], truth_table[(2 * i) + 1])])
        return self.WalshTransform(res)

    # Walsh Spectrum
    def WalshSpectrum(self):
        # Convert truth table 0 -> 1, 1 -> -1
        return self.WalshTransform([[(-1) ** i] for i in self.TruthTable()])


# Just some test for this file.
def PolynomialF2_print__eval_test():
    p = PolynomialF2([{1}, {1, 2, 3}, set()], 4)
    print("p(X):\n", p)
    print("p(1,0,1,1)", p.Evaluate((1, 0, 1, 1)))
    print("p(1,1,1,1)", p.Evaluate((0, 1, 1, 1)))
    print("p(1,1,1,1)", p.Evaluate((0, 1, 0, 1)))

    print("truthTable:", p.TruthTable())
    print("Walsh spectrum:", p.WalshSpectrum())


def PolynomialF2_sum_mul_test():
    a = PolynomialF2([{0}, {1}, {2}], 4)
    b = PolynomialF2([{0}, {1}, {3}], 4)
    c = PolynomialF2([{0}, {1}], 4)

    print("a: ", a)
    print("b: ", b)
    print("c: ", c)
    print("c + b: ", c + b)
    print("c * b:", c * b )
    print("a + c * b: ", a + (c * b))


def PolynomialF2_to_vector_test():
    a = PolynomialF2([{0}, {1}, {2}], 4)
    b = PolynomialF2([{0}, {1}, {3}], 4)
    c = PolynomialF2([{0}, {1}], 4)

    print("a: ", a.to_vector())
    print("c * b: ", (c * b).to_vector())
    print("a + c * b: ", (a + (c * b)).to_vector())
    print("a: ", a)
    print("a: ", PolynomialF2.from_vector("".join([str(i) for i in a.to_vector()]), 4))


#PolynomialF2_sum_mul_test()