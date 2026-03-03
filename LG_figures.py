import pandas as pd
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import os

def get_coefficient_matrix_laurent_fast(E):
    # Define symbolic variables
    t0, t1 = sp.symbols('t0 t1')

    # Expand the expression
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

def visualize_coefficient_matrix(A, offset, title="Coefficient Matrix A"):
    min_i, min_j = offset

    extent = [min_i - 0.5, min_i + A.shape[0] - 0.5,
              min_j - 0.5, min_j + A.shape[1] - 0.5]

    fig = plt.figure(figsize=(8, 6))
    im = plt.imshow(A, extent=extent, cmap='coolwarm', interpolation='none', origin='lower')
    plt.colorbar(im, label="Coefficient value")

    plt.xticks(np.arange(min_j, min_j + A.shape[1]))
    plt.yticks(np.arange(min_i, min_i + A.shape[0]))
    plt.xlabel("Exponent of $t_0$")
    plt.ylabel("Exponent of $t_1$")
    plt.title(title)
    plt.grid(True, linestyle=':', alpha=0.5)
    # plt.show()
    return fig

def create_fold(fpath):
    if not os.path.exists(fpath):
        os.makedirs(fpath)
        print(f"Created folder: {fpath}")
    else:
        print(f"Folder already exists: {fpath}")

def main():
    fname = 'LG-data_3-15c-alt-lc-unim-Lorz.csv'
    LGdata15 = pd.read_csv(fname)
    LGdata15["Crossings"] = LGdata15["SnapPy Name"].str.split(r"[_an]").str[0]

    # generate LG diagrams of alternating knots with 3-10 crossings
    for c in range(3,11):
        figure_path = f"./figures/{c}c/"
        create_fold(figure_path)
        df = LGdata15.loc[(LGdata15["Crossings"]==str(c)) &
                          (LGdata15["Alternating Knot"]=="Y")
                          ]
        for p, name in zip(df["LG(t0,t1) Polynomial"], df["SnapPy Name"]):
            E = LGdata15.loc[LGdata15["SnapPy Name"] == name,'LG(t0,t1) Polynomial'].iloc[0]
            A, offset = get_coefficient_matrix_laurent_fast(E)
            fig = visualize_coefficient_matrix(np.abs(A), offset, title=f"LG(t0,t1) Polynomial Coefficients of {name}")
            fig.savefig(f"{figure_path}{name}.png",dpi=200)
            plt.close()

if __name__ == "__main__":
    main()
