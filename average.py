import numpy as np
import h5py
import glob
import os

os.makedirs("data/mean", exist_ok=True)

n = 25
λ = 1.5
seed = 47607
for β in [0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0]:
    print("Processing β={}".format(β))
    fns = []
    for i in range(1, 61):
        input_filename = "data/autocorr_n={}_λ={}_β={}_seed={}.h5".format(n, λ, β, seed + i)
        with h5py.File(input_filename, "r") as input:
            fns.append(np.asarray(input["data"]))
        # np.loadtxt(input_filename))
    fns = np.vstack(fns)
    fns = np.mean(fns, axis=0)
    output_filename = "data/mean/autocorr_n={}_λ={}_β={}.csv".format(n, λ, β)
    np.savetxt(output_filename, fns)

for β in [0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0]:
    print("Processing β={}".format(β))
    fns = []
    for filename in glob.glob("data/structure_n={}_λ={}_β={}_seed=*.h5".format(n, λ, β)):
        with h5py.File(filename, "r") as input:
            fns.append(np.expand_dims(input["data"], 0))
    fns = np.vstack(fns)
    fns = np.mean(fns, axis=0)
    output_filename = "data/mean/structure_n={}_λ={}_β={}.h5".format(n, λ, β)
    with h5py.File(output_filename, "w") as output:
        output["data"] = fns

# for t_w in [32, 128, 512, 2048, 8192, 32768]:
#     print("Processing t_w={}".format(t_w))
#     fns = []
#     for i in range(1, 601):
#         input_filename = "data/autocorr_n={}_β={}_t={}_seed={}.h5".format(n, β, t_w, seed + i)
#         with h5py.File(input_filename, "r") as input:
#             fns.append(np.asarray(input["autocorr"]))
#         # np.loadtxt(input_filename))
#     fns = np.vstack(fns)
#     fns = np.mean(fns, axis=0)
#     output_filename = "data/mean_autocorr_n={}_β={}_t={}.csv".format(n, β, t_w)
#     np.savetxt(output_filename, fns)
