import glob
import pandas as pd
from tqdm import tqdm
import os

dirname= '/priv/debass/DES/SIMS/hostlib/output_hostlib_mpi/SFH_mpi/fixed/with_neb_av0-1.5'
outfile = '/priv/debass/DES/SIMS/hostlib/output_hostlib_mpi/SFH_mpi/fixed/with_neb_av0-1.5/all_model_params_quench_BC03_z0.0005_1.10000_av0.00_1.50_rv_rand_full_age_dists_neb_U-2.00_res_1_beta_1.14_combined.h5'

# Collect all dataframes into a list
df_list = []
for fn in tqdm(glob.glob(os.path.join(dirname,'*av0.00_1.50*.h5'))):
    df = pd.read_hdf(fn)
    df_list.append(df)


# Combine them efficiently
full_df = pd.concat(df_list, ignore_index=True)

# Try to cast all columns to float
for col in full_df.columns:
    try:
        full_df[col] = full_df[col].astype(float)
    except:
        print("Could not convert:", col)

# Save to HDF5
full_df.to_hdf(outfile, key='main')

# Rename .dat files
for fn in tqdm(glob.glob(os.path.join(dirname, 'SN_ages/*'))):
    new_fn = fn.replace('0.dat','0_combined.dat')
    os.rename(fn,new_fn)
