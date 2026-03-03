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

def alternating_name(name):
    if "n" not in name:
        return 'Y'
    return 'N'

def process_chunk(args):
    chunk_df, chunk_id, output_dir = args

    f_name = f"chunk_{chunk_id}.csv"
    file_path = os.path.join(output_dir, f_name)
    if os.path.exists(file_path):
        return file_path #skip over already computed chunks, in case the computation was interrupted

    chunk_df["Alternating Knot"] = [alternating_name(p) for p in chunk_df["SnapPy Name"]]
    chunk_df.to_csv(file_path, index=False)
    return file_path

def run_computation(crossing):
    if crossing == 15:
        input_file = "LG-data_3-15c.csv"
        chunks_dir = "chunks_output_3-15c-alt"
        prelim_file = "LG-data_3-15c-alt00.csv"
        final_file = "LG-data_3-15c-alt.csv"

        crossing_text = "3-15c"

    elif crossing == 16:
        input_file = "LG-data_16c.csv"
        chunks_dir = "chunks_output_16c-alt"
        final_file = "LG-data_16c-alt.csv"

        crossing_text = "16c"

    else:
        print("Invalid input")
        return

    print("Beginning LG_with_alternating:", crossing_text)

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
    if crossing == 16:
        if os.path.exists(final_file):
            os.remove(final_file)

        for i, fp in enumerate(sorted(file_paths)):
            chunk = pd.read_csv(fp)
            chunk.to_csv(final_file, mode='a', index=False, header=(i == 0))

    if crossing == 15:
        if os.path.exists(prelim_file):
            os.remove(prelim_file)

        if os.path.exists(final_file):
            os.remove(final_file)

        for i, fp in enumerate(sorted(file_paths)):
            chunk = pd.read_csv(fp)
            chunk.to_csv(prelim_file, mode='a', index=False, header=(i == 0))

        sm_alt_data = pd.read_csv("small_alternating.csv")
        final = pd.read_csv(prelim_file)
        final.loc[0:249,'Alternating Knot'] = sm_alt_data.loc[0:249,'Alternating']
        final.to_csv(final_file, index=False)

        if os.path.exists(prelim_file):
            os.remove(prelim_file)

    print(f"LG{crossing_text} + alternating knot .csv printed to: {final_file}")
    if os.path.exists(chunks_dir):
        shutil.rmtree(chunks_dir)

    print(f"Cleaned temporary chunk files for alternating. Done with {crossing_text}.")

def main():
    for c in [15,16]:
        run_computation(c)

if __name__ == "__main__":
    main()
