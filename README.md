# On Some Log-concavity Properties Of The Alexander--Conway And Links--Gould Invariants
**Authors:** Matthew Harper, Ben-Michael Kohli, Jiebo Song, Guillaume Tahar

This repository provides the official Python implementation of the paper \
[On Some Log-concavity Properties Of The Alexander--Conway And Links--Gould Invariants](https://arxiv.org/abs/2509.16868).

The data files hosted and python scripts are mirrored at \
[Harvard Dataverse:10.7910/DVN/E3UFLP](https://arxiv.org/abs/math/2509.16868).

## Repository Description

We provide scripts for obtaining the Links--Gould polynomial (LG) in the variables $$(t_0,t_1)$$ for knots up to 16 crossings. We verify a version of the Fox-Stoimenow conjecture for this 2-variable polynomial on alternating knots. We also show that these polynomials are not denormalized Lorentzian.

The polynomials were computed from the dataset at [doi.org/10.7910/DVN/XE4TOF](https://doi.org/10.7910/DVN/XE4TOF). The .csv files `V-database_3-15c.db.tar.bz2` and `V1-data_16c.csv.tar.bz2` therein are necessary to run our python code and reproduce our dataset. These files must be provided.
   
## Scripts

With the exception of `d2cube_generator.py` and `LG_figures.py`, these scripts are designed to run on an HPC cluster managed by SLURM. However, the worker functions for these HPC scripts are modular and can also be executed locally.

**Script descriptions:**

- `d2cube_generator.py` \
Generates the set of all unordered combinations of $$t_0, t_1, z$$ (lexicographically) up to length 39. These combinations are used by `LG_test_Lorentzian.py`.


- `LG_figures.py` \
Generates heatmaps for alternating knots with 3-11 crossings using `LG-data_3-15c-alt-lc-unim-Lorz.csv`. Figures are sorted into folders by crossing number. Can also run from `LG-data_3-15c-alt.csv` by changing fname.


- `LG_test_log-concavity_unimodality.py` \
Tests Log-Concavity and Unimodality of LG for alternating knots and concatenates the results onto a new file. \
Input: `LG-data_3-15c-alt.csv` or `LG-data_16c-alt.csv` 
(from `LG_with_alternating.py`). \
Output: `LG-data_3-15c-alt-lc-unim.csv` or `LG-data_16c-alt-lc-unim.csv`, and a summary in `LG_computation_results.txt`.


- `LG_test_Lorentzian.py` \
Tests whether LG is denormalized Lorentzian for alternating knots and concatenates the results onto a new file. \
Input: `LG-data_3-15c-alt-lc-unim.csv` or `LG-data_16c-alt-lc-unim.csv` 
(from `LG_test_log-concavity_unimodality.py`). \
Output: `LG-data_3-15c-alt-lc-unim-Lorz.csv` or `LG-data_16c-alt-lc-unim-Lorz.csv` and a summary in `LG_computation_results.txt`.


- `LG_with_alternating.py` \
Determines whether a knot is alternating and concatenates the results onto a new file. \
Input: `LG-data_3-15c.csv` or `LG-data_16.csv` 
(from `V1_to_LG.py`) \
Output: `LG-data_3-15c-alt.csv` and `LG-data_16c-alt.csv`.


- `V1_to_LG.py` \
Converts the $$V_1(q,t)$$ polynomial to LG and concatenates the results onto a new file. \
Input: `V1-data_3-15c.csv` or `V1-data_16.csv` 
(from [doi.org/10.7910/DVN/XE4TOF](https://doi.org/10.7910/DVN/XE4TOF) and is not provided in the database)  \
Output: `V1LG-data_3-15c.csv` and `V1LG-data_16c.csv` (concatenated, not used in later scripts), \
	 	`LG-data_3-15c.csv` and `LG-data_16c.csv` (LG only).

## Script Execution Order

To reproduce the dataset and analyses, the scripts should be run in the 
following order:

1. `V1_to_LG.py`
   - Converts the $$V1(q,t)$$ polynomial to LG
   - Input: $$V_1$$ dataset files from [doi.org/10.7910/DVN/XE4TOF](https://doi.org/10.7910/DVN/XE4TOF)
   - Output: `LG-data_3-15c.csv` / `LG-data_16c.csv`

2. `LG_with_alternating.py`
   - Identifies alternating knots
   - Input: `LG-data_3-15c.csv` / `LG-data_16c.csv`
   - Output: `LG-data_3-15c-alt.csv` / `LG-data_16c-alt.csv`

3. `LG_test_log-concavity_unimodality.py`
   - Tests Log-Concavity and Unimodality
   - Input: `LG-data_3-15c-alt.csv` / `LG-data_16c-alt.csv`
   - Output: `LG-data_3-15c-alt-lc-unim.csv` / `LG-data_16c-alt-lc-unim.csv`
   - Summary written to `LG_computation_results.txt`

4. `d2cube_generator.py`
   - Generates unordered combinations of $$t_0, t_1, z$$
   - Used by `LG_test_Lorentzian.py`

5. `LG_test_Lorentzian.py`
   - Tests denormalized Lorentzian property
   - Input: `LG-data_3-15c-alt-lc-unim.csv` / `LG-data_16c-alt-lc-unim.csv`
   - Output: `LG-data_3-15c-alt-lc-unim-Lorz.csv` / `LG-data_16c-alt-lc-unim-Lorz.csv`
   - Summary written to `LG_computation_results.txt`

6. `LG_figures.py` (optional)
   - Generates heatmaps of alternating knots
   - Input: `LG-data_3-15c-alt-lc-unim-Lorz.csv`
   - Output: Figures sorted by crossing number


## Citation
If you use this dataset or scripts in your research, please cite both 
our and the original dataset as well as our publications/preprints.

```bibtex
@misc{Vn-Patterns,
  title = {Patterns of the {$V_2$}-polynomial of knots},
  author = {Garoufalidis, Stavros and Li, Shana Yunsheng},
  year = {2024},
  howpublished = {\href{https://arxiv.org/abs/math/2409.03557}{arXiv:2409.03557}},
}

@misc{Vn-Dataset,
  title = {Values of V_n-polynomials of knots},
  author = {Garoufalidis, Stavros and Li, Shana Yunsheng},
  year = {2024},
  howpublished = {\href{https://doi.org/10.7910/DVN/XE4TOF}{Harvard Dataverse:10.7910/DVN/XE4TOF}},
}


@misc{LG-Dataset,
      title = {Values of the Links--Gould polynomial}, 
      author = {Harper, Matthew and Kohli, Ben-Michael and Song, Jiebo and Tahar, Guillaume},
      note = {\href{https://doi.org/10.7910/DVN/E3UFLP}{Harvard Dataverse:10.7910/DVN/E3UFLP}},
}

@misc{LG-Log-Concave,
      title = {On some log-concavity properties of the Alexander--Conway and Links--Gould invariants}, 
      author = {Harper, Matthew and Kohli, Ben-Michael and Song, Jiebo and Tahar, Guillaume},
      note = {\href{https://arxiv.org/abs/math/2509.16868}{arXiv:2509.16868}},
}
```
## Acknowledgments
The authors express their gratitude to Stavros Garoufalidis and Shana Li for making their dataset available. The authors also thank David de Wit for his work implementing the LG-explorer which was used in a previous iteration of our investigations, and to Jon Links for making it accessible on his webpage.
MH was partially supported through the NSF-RTG grant #DMS-2135960. 
BMK was partially supported through the BJNSF grant IS24066.
GT was supported by the BJNSF grant IS23005.
