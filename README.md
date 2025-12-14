# MYSTERON Simulation Setup

This repository contains code for simulating supernovae and host galaxies by using Galaxy-driven model (Wiseman, 2022) on Mysteron.

---

## Directory Structure

- **Code (SIMS)**:  
  `$DESCODE`

- **Output**:  
  `$DESSIMS`
  (the output is in another Github repo)

---

## Environment

Activate the conda environment:

`conda activate /priv/debass/software/dessn`

---

## shortcut

```bash
$DESHOME (base)
$DESSIMS       = $DESHOME/SIMS
$AURA          = $DESSIMS/AURA
$config        = $DESSIMS/config
$efficiencies  = $DESSIMS/efficiencies
$mass_assembly = $DESSIMS/mass_assembly
$hostlib       = $DESSIMS/hostlib
$filters       = $DESSIMS/filters
```

--- 
# Simulate Star Formation History

How to run: `python  mass_assembly_quenched_mpi.py output_dir_name`
Example output dir: `/SFH_mpi/sfh_25_50`

```bash
Set up:
--config config/config_rates.yaml
Path to the configuration file (also defines the output base directory).
--output test
Name of the output directory.
--dt 0.5
Time step (in Myr) used to iterate and update the galaxy star formation rate.
--early_step 25
Seeding interval for galaxies in the early Universe (every 25 Myr).
--late_step 50
Seeding interval for galaxies in the late Universe (every 50 Myr).
--n 100
Number of galaxies seeded.
```
--
### nonMPI version
Dir: $AURA/simulations/scripts/mass_assembly_quenched.py
### MPI version (mainly used)
Dir: $AURA/simulations/scripts/mass_assembly_quenched_mpi.py
### Modified code: modified some steps from earlier version to run faster
Dir: $AURA/simulations/scripts/mass_assembly_adjusted.py

--
Final output used in the thesis: $mass_assembly/SFH_mpi/sfh_25_50

---
