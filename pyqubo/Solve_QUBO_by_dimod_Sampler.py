from pyqubo import Binary
x1, x2 = Binary('x1'), Binary('x2')
H = (x1 + x2 - 1)**2
print(H)
# model = H.compile()
# bqm = model.to_bqm()
# print(bqm)

# import neal
# sa = neal.SimulatedAnnealingSampler()
# sampleset = sa.sample(bqm, num_reads=10)
# decoded_samples = model.decode_sampleset(sampleset)
# best_sample = min(decoded_samples, key=lambda x: x.energy)
# print(best_sample.sample)
# print(best_sample.energy)