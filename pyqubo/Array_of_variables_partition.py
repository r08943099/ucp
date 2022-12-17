from pyqubo import Array
numbers = [4, 2, 7, 1] #The parameter for the variable
s = Array.create('s', shape = 4, vartype='SPIN')
# H = (4*s1 + 2*s2 + 7*s3 + s4)**2
H = sum(n*s for n,s in zip(numbers, s))**2
model = H.compile()
qubo, offset = model.to_qubo()
print(qubo)