import numpy as np
import h5py
import glob
import os

os.makedirs("data/mean", exist_ok=True)

def _mean_impl(prefix, n, λ, β, expand):
    fns = []
    for filename in glob.glob("data/{}_n={}_λ={}_β={}_seed=*.h5".format(prefix, n, λ, β)):
        with h5py.File(filename, "r") as input:
            f = np.asarray(input["data"])
        if expand:
            f = np.expand_dims(f, 0)
        fns.append(f)
    fns = np.vstack(fns)
    fns = np.mean(fns, axis=0)
    return fns

def mean_autocorr(n, λ, β):
    return _mean_impl("autocorr", n, λ, β, expand=False)

def mean_structure_factor(n, λ, β):
    return _mean_impl("structure", n, λ, β, expand=True)

n = 25
for λ in [2.5, 4.0]:
    for β in [0.05, 0.2, 0.6, 1.0, 1.4, 2.2, 3.0]:
        print("[*] Processing λ={}, β={}".format(λ, β))
        np.savetxt(
            "data/mean/autocorr_n={}_λ={}_β={}.csv".format(n, λ, β),
            mean_autocorr(n, λ, β)
        )
        with h5py.File("data/mean/structure_n={}_λ={}_β={}.h5".format(n, λ, β), "w") as output:
            output["data"] = mean_structure_factor(n, λ, β)

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
