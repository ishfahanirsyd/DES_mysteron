# MYSTERON Simulation Setup

This repository contains code for simulating supernovae and host galaxies by using Galaxy-driven model (Wiseman, 2022) on Mysteron.

---

## Directory Structure

- **Code (SIMS)**:  
  `$DESCODE`

- **Output**:  
  `$DESSIMS`
  (the output is in another Github repo: [GitHub](https://github.com/ishfahanirsyd/DES_output))

---

## Environment

Activate the conda environment:

`conda activate /priv/debass/software/dessn`

---

## shortcut

```bash
Base: $DESHOME 
$DESCODE       = $DESHOME/software/DES
$DESSIMS       = $DESHOME/SIMS
$AURA          = $DESSIMS/AURA
$config        = $DESSIMS/config
$efficiencies  = $DESSIMS/efficiencies
$mass_assembly = $DESSIMS/mass_assembly
$hostlib       = $DESSIMS/hostlib
$filters       = $DESSIMS/filters
```

--- 
# 1. SIMULATE STAR FORMATION HISTORY

How to run: `python  mass_assembly_quenched_mpi.py output_dir_name`

Example output dir: `/SFH_mpi/sfh_25_50`


## Set up:

```md
`--config config/config_rates.yaml`

Path to the configuration file (also defines the output base directory).

`--output test`

Name of the output directory.

`--dt 0.5`

Time step (in Myr) used to iterate and update the galaxy star formation rate.

`--early_step 25`

Seeding interval for galaxies in the early Universe (every 25 Myr).

`--late_step 50`

Seeding interval for galaxies in the late Universe (every 50 Myr).

`--n 100`

Number of galaxies seeded.
```


## nonMPI version
Dir: `$DESCODE/SIMS/AURA/simulations/scripts/mass_assembly_quenched.py`
## MPI version (mainly used)
Dir: `$DESCODE/SIMS/AURA/simulations/scripts/mass_assembly_quenched_mpi.py`
## Modified code: modified some steps from earlier version to run faster
Dir: `$DESCODE/SIMS/AURA/simulations/scripts/mass_assembly_adjusted.py`

---

**Final output used in the thesis:** `$mass_assembly/SFH_mpi/sfh_25_50`

----
# 2. SIMULATE HOST GALAXIES

Output from this part goes to: `$hostlib`


Example of output directory: `/output_hostlib_mpi/SFH_mpi/fixed/with_neb_av0-1`

When creating the output directory, a subdirectory named `SN_ages` must also be created inside it. This subdirectory is used to store the age distribution of supernovae.

## MPI version
Takes only one galaxy from step 1 at each time step.

Output:
- One HDF5 file per key (stored separately)
- One combined HDF5 file containing all keys

Dir: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched_mpi.py`

#### Running the code
```md
With nebular emission: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched_mpi.py –neb -o /output-directory`

Without nebular emission: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched_mpi.py -o /output-directory`
```

## nonMPI version
Takes all simulated galaxies from step 1 in each time steps

Dir: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched.py`

#### Running the code:
```md
With nebular emission: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched.py –neb -o output-directory`

Without nebular emission: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched.py -o output-directory`
```

## Modified version (mainly used)
Output:
- One HDF5 file per key (stored separately)

If this version is used, the output must be passed to **Step 3** before proceeding to **Step 4**.

Dir: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched_mpi_mod.py`

#### Running the code:
```md
With nebular emission: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched_mpi_mod.py –neb -o output-directory`

Without nebular emission: `$DESCODE/SIMS/AURA/simulations/scripts/make_hostlib_quenched_mpi_mod.py -o output-directory`
```

---

**Final output used in the thesis: `$hostlib/output_hostlib_mpi/SFH_mpi/fixed/with_neb_av0-1`**

---

# 3. COMBINE HOSTLIB
Used to combine separate files for each keys from step 2 

Dir: `$DESCODE/SIMS/AURA/simulations/scripts/combine_hostlib.py`

### Modify the input and output directory
- `dirname`: input directory
- `outfile`: output file

### Running the code:
`combine_hostlib.py`

---

** Final output used in the thesis: `$hostlib/output_hostlib_mpi/SFH_mpi/fixed/with_neb_av0-1/all_model_params_quench_BC03_z0.0005_1.10000_av0.00_1.00_rv_rand_full_age_dists_neb_U-2.00_res_1_beta_1.14_combined.h5`**

---

# 4. SUPERNOVAE SIMULATION
Dir: `$DESCODE/SIMS/AURA/simulations/scripts/run_sim_multi_thread_mod.py`

Running the code: `run_sim_multi_thread_mod.py grid_config_file`


## Grid Configuration

Grid configuration files define the **Rv parameter range** and the **age or mass step ranges** used in the supernova simulations.

**Directory:**
`$DESCODE/SIMS/AURA/simulations/config/grid_configs`

### Available Grid Configurations

- **Null model**
  - `test.yaml`

- **Rv–mass model**
  - `BS21_nostep.yaml`
  - `BS21_massstep.yaml`

- **Rv–galaxy age model**
  - `W21_nostep.yaml`
  - `W21_agestep.yaml`

- **Rv–SN age model**
  - `W21_linearage.yaml`

## Config
Used for define the property of simulated supernovae.

**Directory:**
`$DESCODE/SIMS/AURA/simulations/config`

### Available Configuration Files

- **Null model**
  - `baseline_model.yaml`

- **Rv–mass model**
  - `BS21_nostep.yaml`
  - `BS21_massstep.yaml`
  - `BS21_intrinmass.yaml`
  - `BS21_intrinage.yaml`

- **Rv–galaxy age model**
  - `W21_nostep.yaml`
  - `W21_intrinmass.yaml`
  - `W21_intrinage.yaml`
  - `W21_agestep.yaml`

- **Rv–SN age model**
  - `W21_Rv_linear_age.yaml`

---


#### Observational Data
**Directory:** `$DESCODE/SIMS/DES-OzDES_data`

**Data analysis and visualization notebooks:**
- `research_main_mass.ipynb` — mass step
- `research_main_color.ipynb` — color step
- `research_main_oii.ipynb` — [O II] step
- `compilation.ipynb` — data compilation

**Data files:**
- **DES 5-year sample**
  - `BBC1D_fixgamma0.FITRES`
- **OzDES with [OII]**
  - `ozdes_oII.csv`
- **Foundation sample**
  - `Foundation_Master_File.csv`
  - `FN_Host_Properties-3.csv`

---

#### Analysis
*Visualizes observed data compared to simulated data*

  
**Simulated cosmic star formation history (CSFH):**
- `SFH.ipynb`

**Simulated host galaxies:**
- `hostlib.ipynb`

**Observed vs simulated supernovae:**
- `SN.ipynb`  

**Miscellaneous plots from Master’s work:**
- `draft.ipynb`

**Combined observed datasets:**
- `df_merged.csv` — hosts without undetectable [O II]
- `df_merged_new.csv` — hosts including undetectable [O II]
