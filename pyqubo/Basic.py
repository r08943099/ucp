from pyqubo import Binary
x1, x2, x3, x4 = Binary("x1"), Binary("x2"), Binary("x3"), Binary("x4")
H = (4*x1 + 2*x2 + 7*x3 + x4)**2
model = H.compile()
qubo, offset = model.to_qubo()
print(qubo)
print(offset)