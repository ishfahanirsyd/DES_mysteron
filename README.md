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
Dir: `$AURA/simulations/scripts/mass_assembly_quenched.py`
## MPI version (mainly used)
Dir: `$AURA/simulations/scripts/mass_assembly_quenched_mpi.py`
## Modified code: modified some steps from earlier version to run faster
Dir: `$AURA/simulations/scripts/mass_assembly_adjusted.py`

---

**Final output used in the thesis: `$mass_assembly/SFH_mpi/sfh_25_50`**

----
# 2. SIMULATE HOST GALAXIES

Output from this part goes to: `$hostlib`


Example of output directory: `/output_hostlib_mpi/SFH_mpi/fixed/with_neb_av0-1`

When creating the output directory, a subdirectory named `SN_ages` must also be created inside it. This subdirectory is used to store the age distribution of supernovae.

## MPI version
(Takes only one galaxy from step 1 at each time step.)

Output:
- One HDF5 file per key (stored separately)
- One combined HDF5 file containing all keys

Dir: `$AURA$/simulations/scripts/make_hostlib_quenched_mpi.py`

#### Running the code
```md
- With nebular emission: `$AURA/simulations/scripts/make_hostlib_quenched_mpi.py –neb -o /output-directory`
- Without nebular emission: `$AURA/simulations/scripts/make_hostlib_quenched_mpi.py -o /output-directory`
```

## nonMPI version
(takes all simulated galaxies from step 1 in each time steps)

Dir: `$AURA/simulations/scripts/make_hostlib_quenched.py`

#### Running the code:
```md
- With nebular emission: `$AURA/simulations/scripts/make_hostlib_quenched.py –neb -o output-directory`
- Without nebular emission: `$AURA/simulations/scripts/make_hostlib_quenched.py -o output-directory`
```

## Modified version (mainly used)
Output:
- One HDF5 file per key (stored separately)

If this version is used, the output must be passed to **Step 3** before proceeding to **Step 4**.

Dir: `$AURA/simulations/scripts/make_hostlib_quenched_mpi_mod.py`

#### Running the code:
```md
- With nebular emission: `$AURA/simulations/scripts/make_hostlib_quenched_mpi_mod.py –neb -o output-directory`
- Without nebular emission: `$AURA/simulations/scripts/make_hostlib_quenched_mpi_mod.py -o output-directory`
```
---
**Final output used in the thesis: `$hostlib/output_hostlib_mpi/SFH_mpi/fixed/with_neb_av0-1`**
---
