# Algebraic attack
import copy
from PolynomialF2 import *
from LinearAlgebraF2 import *
from LFSR import *

# Mandatory assignment solve 3:
lfsr_poly = [0] * 32

# Update states.
lfsr_poly[31] = 1
lfsr_poly[26] = 1
lfsr_poly[5] = 1
lfsr_poly[3] = 1
lfsr_poly[0] = 1

# Initialize parameters from text.
lfsr = LFSR(lfsr_poly, [0] * 31)

variable_num = 31
key_stream = [1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1,
              0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1,
              1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1,
              1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1,
              0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
              0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1,
              1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0,
              1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1,
              0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0,
              1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0,
              1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0,
              0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0,
              0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
              0, 0, 0, 1, 1, 1]


# From the powers of companion matrix get the internal state at time i as linear combinations of the initial state.
# At time 0 we have s3, s2, s1, s0
# At time i: si  = M^i last column.
print("Generate linear relation with initial state.")
M = MatrixF2.identity(variable_num)
s = []
for i in range(len(key_stream) + lfsr.dimension - 1):
    s.append(M.get_column(M.col_num - 1))
    M *= lfsr.companion_matrix

# Convert states in linear polynomial.
# So that p[i] represent state s[i] as a linear polynomial in the initial state variables s[0]... s[30]
p = [PolynomialF2([{variable_num - 1 - j} for j, val in enumerate(i) if val == 1], variable_num) for i in s]

# Apply polynomial x0 + x5 + x1*x3 + x2*x4 using as variables the polynomials p.
print("Calculate the polynomials resulting applying the Filter to the linear polynomials representing the LFSR states.")
print("It will take some time...")
# p[i] represent state si in the variables s0,...,s30
# After i step:
# x0 -> i + 30
# x1 -> i + 25
# x2 -> i + 20
# x3 -> i + 11
# x4 -> i + 7
# x5 -> i + 0
row = [(p[i + 30] + p[i] + (p[i + 25] * p[i + 11]) + (p[i + 20] * p[i + 7])).to_vector() for i in range(len(key_stream))]

# Create linear system.
print("Create matrix associated to the linearized system.")
# Number of position needed to represent 4 variables monomial up to degree 2 included.
col_num = sum([comb(variable_num, i) for i in range(3)])
A = MatrixF2(len(key_stream), col_num, row)

# Solve the linear system
print("Solve the linear system.")
solution = A.solve_right(key_stream)
# I print just the interting part of the solution.
# No need to print first position representing the constant 1
# and the positions representing the monomials of degree bigger than 1.
print("Found the solution: \n", solution[1:32])

# Generate the stream using the solution to verify it is the correct solution.
print("Verify solution.")
sol = solution[1:32]
# Index counting from 0 to 30.
index_filter = [0, 5, 10, 19, 23, 30]

# Insert solution in from state s[30] to state s[0].
lfsr = LFSR(lfsr_poly, sol[::-1])
bool_filter = PolynomialF2([{0}, {5}, {1, 3}, {2, 4}], 6)
f_lfsr = F_LFSR(lfsr, bool_filter, index_filter)

print("Solution keystream:\n", f_lfsr.Stream(496) == key_stream)
