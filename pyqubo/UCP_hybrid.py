from pyqubo import Binary, Array, Constraint, Placeholder
import dimod
from dwave.system import LeapHybridSampler
import dwavebinarycsp
import numpy as np  
import numba as nb

# f = open("input/2D_three.txt")
# f = open("input/UCP/benchmark_01.txt")
# f = open("input/UCP/benchmark_02.txt")
# f = open("input/UCP/benchmark_03.txt")
# f = open("input/UCP/benchmark_04.txt")
# f = open("input/UCP/benchmark_05.txt")
# f = open("input/UCP/benchmark_06.txt")
# f = open("input/UCP/benchmark_07.txt")
# f = open("input/UCP/benchmark_08.txt")
# f = open("input/UCP/benchmark_09.txt")
f = open("input/UCP/benchmark_10.txt")

line = f.readline()
cellnum=0
positionnum=0
padnum=0
M = 0
E_cp = set()
E_cc = set()
tuple_cp = tuple()
tuple_cc = tuple()
print("Start parse")
while line:
    s = line.split(' ')
    if(s[0]=='CellNum:'): cellnum=int(s[1])
    elif(s[0]=='PositionNum:'): positionnum=int(s[1])
    elif(s[0]=='PadNum:'): padnum=int(int(s[1]))
    elif(s[0]=='cell_pad:'): 
        tuple_cp=(float(s[1]),int(s[2]))
        E_cp.add(tuple_cp)
    # elif(s[0]=='cell_cell:'): 
    #     tuple_cc=(float(s[1]),int(s[2]),int(s[3]))
    #     E_cc.add(tuple_cc)
    elif(s[0]=='max_weight:'):
        M = float(s[1])  
    else:
        tuple_cc=(float(s[0]),int(s[1]),int(s[2]))
        E_cc.add(tuple_cc)       
    line = f.readline()
x = Array.create('x', shape=(cellnum,positionnum), vartype = 'BINARY')
# print(x)
q = np.array([])
for i in x:
    for j in i:
        q = np.append(q, j)
print("Start Construct H")        
#===========Construct H =============
M = 2*M
# # M = Placeholder('M')
#===Objective===
print("Start Construct HA,HB")  
HA=sum(i*(q[j]*q[j]) for (i, j) in E_cp)
HB=sum(i*(q[j]*q[k]) for (i, j, k) in E_cc)
#===Constraint===
print("Start Construct HC")  
HC = sum(M*Constraint(((sum(x[i][j] for j in range(positionnum))-1))**2,label="cell%s in different positions" %i) for i in range(cellnum))
print("Start Construct HD")
# HD = sum(M*Constraint(sum(sum(x[j][i]*x[k][i] for k in range(j+1,cellnum)) for j in range(cellnum)),label="number of cell for position%s" %i))
HD_tmp = list()
for i in range(positionnum):
    print(i)
    HD_tmp.append(M*Constraint(sum(sum(x[j][i]*x[k][i] for k in range(j+1,cellnum)) for j in range(cellnum)),label="number of cell for position%s" %i))
HD = sum(HD_tmp)

H = HA+HB+HC+HD
model = H.compile()
bqm = model.to_bqm()

#solve
print("Start solve")
sampler = LeapHybridSampler(solver={'category': 'hybrid'})
sampleset = sampler.sample(bqm)
decoded_samples = model.decode_sampleset(sampleset)
best_sample = min(decoded_samples, key=lambda x: x.energy)
# print(best_sample.sample)
print(best_sample.energy)
print(best_sample.constraints(only_broken=True))
print(sampleset.info)