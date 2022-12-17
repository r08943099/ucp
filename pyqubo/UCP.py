from pyqubo import Binary, Array, Constraint, Placeholder
import neal
import dwavebinarycsp
import numpy as np
import time
st = time.time()

# f = open("input/2D_three.txt")
# f = open("input/cm151a.txt")
# f = open("input/cm138a.txt")
# f = open("input/cm150a.txt")
f = open("input/cm162a.txt")
line = f.readline()
cellnum=0
positionnum=0
padnum=0
E_cp = set()
E_cc = set()
tuple_cp = tuple()
tuple_cc = tuple()
while line:
    s = line.split(' ')
    if(s[0]=='CellNum:'): cellnum=int(s[1])
    elif(s[0]=='PositionNum:'): positionnum=int(s[1])
    elif(s[0]=='PadNum:'): padnum=int(int(s[1]))
    elif(s[0]=='cell_pad:'): 
        tuple_cp=(float(s[1]),int(s[2]))
        E_cp.add(tuple_cp)
    elif(s[0]=='cell_cell:'): 
        tuple_cc=(float(s[1]),int(s[2]),int(s[3]))
        E_cc.add(tuple_cc)
    elif(s[0]=='max_weight:'):
        M = float(s[1])    
    line = f.readline()
x = Array.create('x', shape=(cellnum,positionnum), vartype = 'BINARY')
# print(x)
# print(cellnum)
# print(positionnum)
q = np.array([])
for i in x:
    for j in i:
        q = np.append(q, j)
#===========Construct H =============
M = 2*M
# # M = Placeholder('M')
#===Objective===
HA=sum(i*(q[j]*q[j]) for (i, j) in E_cp)
HB=sum(i*(q[j]*q[k]) for (i, j, k) in E_cc)
#===Constraint===
HC = sum(M*Constraint(((sum(x[i][j] for j in range(positionnum))-1))**2,label="cell%s in different positions" %i) for i in range(cellnum))
# HD = sum(sum(sum(x[j][i]*x[k][i] for k in range(j+1,cellnum)) for j in range(cellnum)) for i in range(positionnum))
HD = sum(M*Constraint(sum(sum(x[j][i]*x[k][i] for k in range(j+1,cellnum)) for j in range(cellnum)),label="number of cell for position%s" %i) for i in range(positionnum))



# HC=Constraint((sum(q[i] for i in range(positionnum))))
H = HA+HB+HC+HD
model = H.compile()
bqm = model.to_bqm()
# print(HC)


#solve
sa = neal.SimulatedAnnealingSampler()
sampleset = sa.sample(bqm, num_reads=10, num_sweeps=10*cellnum**(4/3))
decoded_samples = model.decode_sampleset(sampleset)
best_sample = min(decoded_samples, key=lambda x: x.energy)
# print(best_sample.sample)
print(best_sample.energy)
print(best_sample.constraints(only_broken=True))

et = time.time()
# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')