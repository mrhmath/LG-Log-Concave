import pandas as pd
import sympy as sp
import numpy as np
import multiprocessing as mp
import os
import shutil
import time

from numba import njit

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

# =========================
# SLURM CPU DETECTION
# =========================
def get_slurm_cpus():
    return int(os.environ.get("SLURM_CPUS_PER_TASK", 1))

# =========================
# WORKER FUNCTIONS
# =========================

def get_coefficient_matrix_laurent_fast(E):
    E = sp.expand(E)

    # Create a sparse coefficient dictionary: keys are (i, j), values are the coefficients
    coeff_dict = {}
    min_i = min_j = float('inf')
    max_i = max_j = float('-inf')

    for term in E.as_ordered_terms():
        term_dict = term.as_powers_dict()
        i = int(term_dict.get(t0, 0))
        j = int(term_dict.get(t1, 0))
        coeff = term / (t0**i * t1**j)

        coeff_dict[(i, j)] = float(coeff)

        # Update the range of exponents
        min_i = min(min_i, i)
        max_i = max(max_i, i)
        min_j = min(min_j, j)
        max_j = max(max_j, j)

    # Initialize the coefficient matrix
    A = np.zeros((max_i - min_i + 1, max_j - min_j + 1))

    # Fill in the coefficients
    for (i, j), coeff in coeff_dict.items():
        A[i - min_i, j - min_j] = coeff

    return A, (min_i, min_j)

@njit
def test_concave(T,i,j,l,m):
    if np.abs(T[i-l,j-m]) * np.abs(T[i+l,j+m]) <= T[i,j]**2:
        return 0
    else:
        return 1

@njit
def check_concave(AA, M, N):
    for l in range(-M, M+1):
        for m in range(-N, N+1):
            for i in range(M, 2*M):
                for j in range(N, 2*N):
                    if test_concave(AA, i, j, l, m) == 1:
                        return False
    return True

def is_concave(polynomial):
    E = sp.sympify(polynomial)
    A, offset = get_coefficient_matrix_laurent_fast(E)
    A = A.astype(np.int64)
    M,N = A.shape
    AA = np.zeros((3*M, 3*N),dtype=np.int64)
    AA[M:2*M, N:2*N] = A

    if check_concave(AA, M, N):
        return 'Y'
    else:
        return 'N'


def level_matrices(A):
    A_entries = np.unique(A).tolist()

    # Compute matrices for each entry
    matrices = []
    for a in A_entries:
        mask = A >= a
        B_a = np.where(mask, A, 0)
        matrices.append(B_a)

    return A_entries, matrices

def no_internal_zeros(A):
    # Step 1: support points
    points = [(i, j)
              for i, row in enumerate(A)
              for j, val in enumerate(row)
              if val != 0]

    if len(points) <= 1:
        return True

    # ---------- Convex Hull (Monotone Chain) ----------
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

    points_sorted = sorted(points)
    lower = []
    for p in points_sorted:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(points_sorted):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    hull = lower[:-1] + upper[:-1]  # vertices in CCW order

    # ---------- Point-in-Convex-Polygon Test ----------
    def point_in_convex_polygon(pt, poly):
        sign = None
        n = len(poly)
        for i in range(n):
            o = poly[i]
            a = poly[(i+1)%n]
            cp = cross(o, a, pt)
            if cp != 0:
                if sign is None:
                    sign = cp > 0
                elif (cp > 0) != sign:
                    return False
        return True  # inside or on boundary

    support_set = set(points)

    # ---------- Scan bounding box ----------
    min_i = min(p[0] for p in points)
    max_i = max(p[0] for p in points)
    min_j = min(p[1] for p in points)
    max_j = max(p[1] for p in points)

    for i in range(min_i, max_i+1):
        for j in range(min_j, max_j+1):
            if point_in_convex_polygon((i, j), hull):
                if (i, j) not in support_set:
                    return False  # internal zero found

    return True

def is_unimodal(polynomial):
    E = sp.sympify(polynomial)
    A, offset = get_coefficient_matrix_laurent_fast(E)
    A = abs(A.astype(int))
    A_entries, matrices = level_matrices(A)
    if all(no_internal_zeros(M) for M in matrices):
        return "Y"
    return "N"

def process_chunk(args):
    chunk_df, chunk_id, output_dir = args

    f_name = f"chunk_{chunk_id}.csv"
    file_path = os.path.join(output_dir, f_name)
    if os.path.exists(file_path):
        return file_path #skip over already computed chunks, in case the computation was interrupted

    chunk_df.loc[:, "Log-Concave"] = [
    is_concave(p) if alt == 'Y' else np.nan
    for p, alt in zip(
        chunk_df["LG(t0,t1) Polynomial"],
        chunk_df["Alternating Knot"]
        )
    ]

    chunk_df.loc[:, "Unimodal"] = [
        is_unimodal(p) if "Y" == alt else np.nan
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
        input_file = "LG-data_3-15c-alt.csv"
        chunks_dir = "chunks_output_3-15c-alt-lc-unim"
        final_file = "LG-data_3-15c-alt-lc-unim.csv"
        result_file = "LG_computation_results.txt"

        crossing_text = "3-15c"

    elif crossing == 16:
        input_file = "LG-data_16c-alt.csv"
        chunks_dir = "chunks_output_16c-lc-unim"
        final_file = "LG-data_16c-alt-lc-unim.csv"
        result_file = "LG_computation_results.txt"

        crossing_text = "16c"

    else:
        print("Invalid input")
        return


    duration_text = f"Elapsed time to compute log-concavity and unimodality for {crossing_text} knots: "
    lc_test_text = f"Log concavity holds for {crossing_text} knots: "
    unim_test_text = f"Unimodality holds for {crossing_text} knots: "

    print("Beginning LG_test_log-concavity_unimodality:", crossing_text)

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

    print(f"LG{crossing_text} + log-concave and unimodal csv printed to: {final_file}")

    if os.path.exists(chunks_dir):
        shutil.rmtree(chunks_dir)

    # =========================
    # TEST PROPERTIES
    # =========================
    final = pd.read_csv(final_file)

    lc_num_passed = final.loc[final["Log-Concave"] == "Y"].shape[0]
    unim_num_passed = final.loc[final["Unimodal"] == "Y"].shape[0]
    num_alternating = final.loc[final["Alternating Knot"] == "Y"].shape[0]

    pass_lc = lc_num_passed == num_alternating
    pass_unim = unim_num_passed == num_alternating

    with open(result_file, "a") as f:
        f.write(duration_text + str(elapsed_time)+"\n")
        f.write(lc_test_text + str(pass_lc)+"(" + str(lc_num_passed) + "/" + str(num_alternating) + " alternating knots)\n")
        f.write(unim_test_text + str(pass_unim)+"(" + str(unim_num_passed) + "/" + str(num_alternating) + " alternating knots)\n")

    print(lc_test_text + str(pass_lc)+"(" + str(lc_num_passed) + "/" + str(num_alternating) + " alternating knots)")
    print(unim_test_text + str(pass_unim)+"(" + str(unim_num_passed) + "/" + str(num_alternating) + " alternating knots)")

    print(f"Cleaned temporary chunk files and tested log-concave and unimodal properties. Done with {crossing_text}.")

def main():
    for c in [15,16]:
        run_computation(c)

if __name__ == "__main__":
    main()
