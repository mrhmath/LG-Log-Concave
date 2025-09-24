# On Some Log-concavity Properties Of The Alexander--Conway And Links--Gould Invariants
**Authors:** Matthew Harper, Ben-Michael Kohli, Jiebo Song, Guillaume Tahar

This repository provides the official Python implementation of the paper \
[On Some Log-concavity Properties Of The Alexander--Conway And Links--Gould Invariants](https://arxiv.org/abs/2509.16868).

## Data Description

This repository includes several datasets derived from LG-explorer output and subsequent processing:

- *A_raw*: original data obtained directly from the [LG-explorer](https://people.smp.uq.edu.au/JonLinks/Links_Gould_Explorer.zip) 
- *A_new*: data produced by applying a change of variables to *A_raw* which unifies conventions across the dataset
- *A_figures*: matrix coefficients corresponding to each LG of a given knot.
- *A_whole*: a combined dataset that aggregates all of the above.
  
	Each dataset in A_whole is provided as Axx.csv files with the following columns:
	- key: Knot name or identifier.
	- original: the raw algebraic expression with $(x,y)$.
	- newExpr: the transformed expression after the change of variables $(t_0,t_1)$.
	- A: matrix of coefficients.
	- test: number of concavity conditions that do not hold.

- *knotinfo.csv*: a snapshot of the complete content of the KnotInfo and LinkInfo databases.
 
- *knotinfo_special_alternating.csv*: selected knots from knotinfo.csv (alternating=Y, positive=N)

## Notebooks

The repository contains Jupyter notebooks for data generation and processing:

- `A1_read_lg_newExpr.ipynb`
  - Generates the *A_new* dataset from *A_raw*.
- `A2_LG_Log_Concave.ipynb`
  - Checks the Log-Concavity property.
  - Generates both *A_whole* and *A_figures* from *A_new*.
- `Knotinfo.ipynb`
  - Selects the desired knot properties from the KnotInfo database.  
  - **Prerequisite:** Install the database first using `pip install database-knotinfo`.

The repository also contains a Python notebook for verifying the Lorenztian property of a polynomial in the variables $(t_0,t_1)$:
- `LG_Lorentzian.py`
  - Uses the *Lorentz* function to test all knots in *A_whole*.

## Citation
If you use this code, please cite the paper:
```bibtex
@misc{LG-Log-Concave,
      title={On some log-concavity properties of the Alexander-Conway and Links-Gould invariants}, 
      author={Harper, Matthew and Kohli, Ben-Michael and Song, Jiebo and Tahar, Guillaume},
      note={Preprint 2025, \href{https://arxiv.org/abs/math/2509.16868}{2509.16868}},
}
```
## Acknowledgments
The authors express their gratitude to David de Wit for his work implementing the LG-explorer, and to Jon Links for making it accessible on his webpage. 
MH was partially supported through the NSF-RTG grant #DMS-2135960. 
BMK was partially supported through the BJNSF grant IS24066.
GT was supported by the BJNSF grant IS23005.