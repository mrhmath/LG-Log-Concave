import pandas as pd
import sympy as sp
import numpy as np
import multiprocessing as mp
import os
import shutil
# =========================
# HPC SAFE SETTINGS
# =========================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

mp.set_start_method("spawn", force=True)

# =========================
# SYMBOLS
# =========================
## V1(q,t) variables
q, t = sp.symbols('q t')
## LG(t0,t1) variables
t0, t1 = sp.symbols('t0 t1')

# =========================
# SLURM CPU DETECTION
# =========================
def get_slurm_cpus():
    return int(os.environ.get("SLURM_CPUS_PER_TASK", 1))

# =========================
# WORKER FUNCTIONS
# =========================
def V1toLG(polynomial): #performs the change of variables q -> (t0t1)^1/2, t -> (t0/t1)^1/2
    p = sp.sympify(polynomial)
    return sp.expand(p.subs([(q,t0**sp.Rational(1,2)*t1**sp.Rational(1,2)), (t,t0**(sp.Rational(1,2))*t1**(-sp.Rational(1,2)))]))

def process_chunk(args):
    chunk_df, chunk_id, output_dir = args

    f_name = f"chunk_{chunk_id}.csv"
    file_path = os.path.join(output_dir, f_name)
    if os.path.exists(file_path):
        return file_path #skip over already computed chunks, in case the computation was interrupted

    chunk_df["LG(t0,t1) Polynomial"] = [V1toLG(p) for p in chunk_df["V1 Polynomial"]]
    chunk_df.to_csv(file_path, index=False)
    return file_path

def run_computation(crossing):
    if crossing == 15:
        input_file = "V1-data_3-15c.csv"
        chunks_dir = "chunks_output_3-15c"
        mixed_file = "V1LG-data_3-15c.csv"
        final_file = "LG-data_3-15c.csv"

        crossing_text = "3-15c"

    elif crossing == 16:
        input_file = "V1-data_16c.csv"
        chunks_dir = "chunks_output_16c"
        mixed_file = "V1LG-data_16c.csv"
        final_file = "LG-data_16c.csv"

        crossing_text = "16c"

    else:
        print("Invalid input")
        return

    print("Beginning V1_to_LG:", crossing_text)

    os.makedirs(chunks_dir, exist_ok=True)

    n_cores = get_slurm_cpus()
    print("Using cores:", n_cores)

    chunk_size = 500

    ctx = mp.get_context("spawn")

    jobs = []
    with ctx.Pool(n_cores) as pool:
        for i, chunk in enumerate(pd.read_csv(input_file, chunksize=chunk_size)):
            jobs.append(pool.apply_async(process_chunk, ((chunk, i, chunks_dir),)))

        file_paths = [job.get() for job in jobs]

    # =========================
    # STREAM CONCATENATION
    # =========================
    if os.path.exists(mixed_file):
        os.remove(mixed_file)

    if os.path.exists(final_file):
        os.remove(final_file)

    for i, fp in enumerate(sorted(file_paths)):
        chunk = pd.read_csv(fp)
        chunk.to_csv(mixed_file, mode='a', index=False, header=(i == 0))

    final = pd.read_csv(mixed_file)[["SnapPy Name", "LG(t0,t1) Polynomial"]]
    final.to_csv(final_file, index=False)

    print(f"V1LG{crossing_text} .csv printed to: {mixed_file}")
    print(f"LG{crossing_text} only .csv printed to: {final_file}")

    if os.path.exists(chunks_dir):
        shutil.rmtree(chunks_dir)

    print(f"Cleaned temporary chunk files for V1 to LG. Done with {crossing_text}.")

def main():
    for c in [15,16]:
        run_computation(c)

if __name__ == "__main__":
    main()
