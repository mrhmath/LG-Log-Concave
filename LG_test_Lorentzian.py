import pandas as pd
import sympy as sp
import numpy as np
import multiprocessing as mp
import os
import shutil
import time

from d2cube_generator import d2cube

# =========================
# HPC SAFE SETTINGS
# =========================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

mp.set_start_method("spawn", force=True)

# =========================
# SYMBOLS (must be global)
# =========================
t0, t1 = sp.symbols("t0 t1")
z, p, q = sp.symbols("z p q")

# =========================
# SLURM CPU DETECTION
# =========================
def get_slurm_cpus():
    return int(os.environ.get("SLURM_CPUS_PER_TASK", 1))

# =========================
# WORKER FUNCTIONS
# =========================
#%%% Check whether polynomial can be made positive

def evaluate_negative(f):
    return sp.expand(f.subs([(t0,-t0),(t1,-t1)]))

def is_positive(f):
    f = sp.Poly(f)
    if not (f-f.abs()).is_zero:
        return 0
    return 1

def make_positive(f):
    f0 = f
    if is_positive(f0) == 0:
        f0 = evaluate_negative(f0)
        if is_positive(f0) == 0:
            f0 = -f0
            if is_positive(f0) == 0:
                f0 = evaluate_negative(f0)
                if is_positive(f0) == 0:
                    return 0
    return f0

#%%% Normalization

def homogeneous(f):
    return sp.Poly(sp.Poly(f).homogenize(z),(t0,t1,z)) #(t0,t1,z) assigns the lexography

def multifactorial(v): #v is a vector
    mf = 1
    for d in v:
        mf *= sp.factorial(d)
    return mf

def normalize(f):
    homf = homogeneous(f)
    monomials = homf.monoms()
    newf = 0
    for m in monomials:
        newf += homf.coeff_monomial(m)*multifactorial(m)**-1*t0**m[0]*t1**m[1]*z**m[2]
    return sp.Poly(newf,(t0,t1,z))

#%%% Count positive eigenvalues of Hessian

def countpos(a):
    ct = 0
    for x in a:
        if x>0:
            ct += 1
    return ct

def H_count_positive_evals(f):
    Hess = sp.Matrix([[(f.as_expr()).diff(x).diff(y) for x in [t0,t1,z]] for y in [t0,t1,z]])
    npHess = sp.matrix2numpy(Hess,dtype=float)
    evals = np.linalg.eigvals(npHess)
    return np.sum(evals > 0)

#%%% Lorentzian test

def shift_vars(f): #f is a Laurant polynomial in t0 and t1, we shift f to polynomial in t0 and t1
    modif_poly = sp.Poly(f.subs([(t0,p),(t1,q)]).as_expr())
    shift_poly = ((modif_poly * p**sp.degree(modif_poly,p**-1) * q**sp.degree(modif_poly,q**-1)).as_expr()).subs([(p,t0),(q,t1)])
    return shift_poly

def is_Lorentz(p):  #Test if polynomial f is Lorenztian, constraint on at most one positive eigenvalue
    f = sp.sympify(p)
    f = shift_vars(f)
    f = make_positive(f)
    if f == 0:
        return "NP" #not positive
    f = normalize(f)
    degree_f = f.total_degree()
    if degree_f <=1:
        return "ND" #not enough degree
    if degree_f == 2:
        if H_count_positive_evals(f) > 1:
            return "N"
        return "Y"
    cube = d2cube(degree_f-2)
    for v in cube:
        fdeg2 = f.diff(*v)
        if H_count_positive_evals(fdeg2) > 1:
            return "N"
    return "Y"

def process_chunk(args):
    chunk_df, chunk_id, output_dir = args

    f_name = f"chunk_{chunk_id}.csv"
    file_path = os.path.join(output_dir, f_name)
    if os.path.exists(file_path):
        return file_path #skip over already computed chunks, in case the computation was interrupted

    chunk_df.loc[:, "Lorentzian"] = [
    is_Lorentz(p) if alt == 'Y' else np.nan
    for p, alt in zip(
        chunk_df["LG(t0,t1) Polynomial"],
        chunk_df["Alternating Knot"]
        )
    ]

    chunk_df.to_csv(file_path, index=False)
    return file_path

# =========================
# MAIN
# =========================
def run_computation(crossing):
    if crossing == 15:
        input_file = "LG-data_3-15c-alt-lc-unim.csv"
        chunks_dir = "chunks_output_3-15c-alt-lc-unim-Lorz"
        final_file = "LG-data_3-15c-alt-lc-unim-Lorz.csv"
        result_file = "LG_computation_results.txt"

        crossing_text = "3-15c"

    elif crossing == 16:
        input_file = "LG-data_16c-alt-lc-unim.csv"
        chunks_dir = "chunks_output_16c-alt-lc-unim-Lorz"
        final_file = "LG-data_16c-alt-lc-unim-Lorz.csv"
        result_file = "LG_computation_results.txt"

        crossing_text = "16c"

    else:
        print("Invalid input")
        return

    duration_text = f"Elapsed time to compute denomalized Lorentzian for {crossing_text} knots: "
    Lorz_test_text = f"Denomalized Lorentzian holds for {crossing_text} knots: "
    not_pos_text = f"Number of alternating {crossing_text} knots which have non-alternating LG: "
    low_deg_text = f"Number of alternating {crossing_text} knots which have degree 1 LG: "

    print("Beginning LG_test_Lorentzian:", crossing_text)

    os.makedirs(chunks_dir, exist_ok=True)

    n_cores = get_slurm_cpus()
    print("Using cores:", n_cores)

    chunk_size = 500

    ctx = mp.get_context("spawn")

    start_time = time.perf_counter()

    jobs = []
    with ctx.Pool(n_cores) as pool:
        for i, chunk in enumerate(pd.read_csv(input_file, chunksize=chunk_size)):
            jobs.append(pool.apply_async(process_chunk, ((chunk, i, chunks_dir),)))

        file_paths = [job.get() for job in jobs]

    # =========================
    # STREAM CONCATENATION
    # =========================
    if os.path.exists(final_file):
        os.remove(final_file)

    for i, fp in enumerate(sorted(file_paths)):
        chunk = pd.read_csv(fp)
        chunk.to_csv(final_file, mode='a', index=False, header=(i == 0))


    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(duration_text, elapsed_time)

    print(f"LG{crossing_text} + Lorenztian .csv printed to: {final_file}")

    if os.path.exists(chunks_dir):
        shutil.rmtree(chunks_dir)

    # =========================
    # TEST PROPERTIES
    # =========================
    final = pd.read_csv(final_file)

    Lorz_num_passed = final.loc[final["Lorentzian"] == "Y"].shape[0]
    num_not_pos = final.loc[final["Lorentzian"] == "NP"].shape[0]
    num_low_deg = final.loc[final["Lorentzian"] == "ND"].shape[0]
    num_alternating = final.loc[final["Alternating Knot"] == "Y"].shape[0]

    pass_Lorz = (Lorz_num_passed == num_alternating)

    with open(result_file, "a") as f:
        f.write(duration_text + str(elapsed_time)+"\n")
        f.write(Lorz_test_text + str(pass_Lorz)+"(" + str(Lorz_num_passed) + "/" + str(num_alternating) + " alternating knots)\n")
        f.write(not_pos_text + str(num_not_pos)+"\n")
        f.write(low_deg_text + str(num_low_deg)+"\n")

    print(duration_text + str(elapsed_time))
    print(Lorz_test_text + str(pass_Lorz)+"(" + str(Lorz_num_passed) + "/" + str(num_alternating) + " alternating knots)")
    print(not_pos_text + str(num_not_pos))
    print(low_deg_text + str(num_low_deg))

    # if os.path.exists("d2cubes"):
    #     shutil.rmtree("d2cubes")

    print(f"Cleaned temporary chunk files and tested denormalized Lorentzian properties. Done with {crossing_text}.")

def main():
    for c in [15,16]:
        run_computation(c)

if __name__ == "__main__":
    main()
