# Verma et al. 2024, Journal of Ionic Liquids, 4, 1, 100089, https://doi.org/10.1016/j.jil.2024.100089

# Python Script for Ionic Conductivity Calculation Using Dipole of the Simulation Box.

# This Python script calculates the ionic conductivity of pure ionic liquids using Einstein formalism, considering all correlation terms. The *.trr and *.tpr file from the GROMACS simulation are required for computing Ionic conductivity.

# To begin, run `transport_setup.sh` and adjust `n_species` according to the system composition. For pure ILs, set `n_species=2`; for mixtures of ILs and a solvent, use `n_species=3`.

# Once `transport_setup.sh` is complete, run `transport_analysis.py` and update `cation_list`, `anion_list`, `solvent_list`, and `charge_multiplier`.

# Output Units: delta or time in ps and product in Å^2

