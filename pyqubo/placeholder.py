from pyqubo import Binary, Constraint, Placeholder
import neal
a,b,c = Binary('a'),Binary('b'),Binary('c')
# M = Placeholder('M')
M = 5
H = 3*a + 2*b + 2*c + 2*a*c + M* Constraint((a+b+c-2)**2, label = 'a+b+c=2')
model = H.compile()
# bqm = model.to_bqm()
bqm = model.to_bqm(feed_dict={'M':6})

#solve
sa = neal.SimulatedAnnealingSampler()
sampleset = sa.sample(bqm, num_reads=10)
decoded_samples = model.decode_sampleset(sampleset)
best_sample = min(decoded_samples, key=lambda x: x.energy)
print(best_sample.sample)
print(best_sample.energy)
print(best_sample.constraints())
print(best_sample.constraints(only_broken=True))

