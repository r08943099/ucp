from pyqubo import Binary, Array, Constraint, Placeholder
import dimod
from dwave.system import LeapHybridSampler
import dwavebinarycsp
import numpy as np

cellnum=400
positionnum=400
E_cp = set()
E_cc = set()
tuple_cp = tuple()
tuple_cc = tuple()
M = 10
x = Array.create('x', shape=(cellnum,positionnum), vartype = 'BINARY')
# print(x)
q = np.array([])
for i in x:
    for j in i:
        q = np.append(q, j)
#===========Construct H =============
M = 2*M
# # M = Placeholder('M')
#===Objective===
H = q[0]+q[1]
# HA=sum(i*(q[j]*q[j]) for (i, j) in E_cp)
# HB=sum(i*(q[j]*q[k]) for (i, j, k) in E_cc)
# #===Constraint===
# HC = sum(M*Constraint(((sum(x[i][j] for j in range(positionnum))-1))**2,label="cell%s in different positions" %i) for i in range(cellnum))
# # HD = sum(sum(sum(x[j][i]*x[k][i] for k in range(j+1,cellnum)) for j in range(cellnum)) for i in range(positionnum))
# HD = sum(M*Constraint(sum(sum(x[j][i]*x[k][i] for k in range(j+1,cellnum)) for j in range(cellnum)),label="number of cell for position%s" %i) for i in range(positionnum))



# HC=Constraint((sum(q[i] for i in range(positionnum))))
model = H.compile()
bqm = model.to_bqm()
# print(HC)


#solve
sampler = LeapHybridSampler(solver={'category': 'hybrid'})
sampleset = sampler.sample(bqm)
decoded_samples = model.decode_sampleset(sampleset)
best_sample = min(decoded_samples, key=lambda x: x.energy)
# print(best_sample.sample)
print(best_sample.energy)
print(best_sample.constraints(only_broken=True))
print(sampleset.info)