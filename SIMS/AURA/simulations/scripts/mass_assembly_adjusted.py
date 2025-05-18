import numpy as np
import pandas as pd
from astropy.cosmology import z_at_value
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(70,0.3)
import os
from scipy.special import erf
from astropy import units as u
import numpy as np
from tqdm import tqdm
import argparse
import yaml
import glob
import sys
sys.path.insert(0, os.path.abspath("../../.."))
from AURA.utils.utils import MyPool

def psi_Mz_alt(M, z):
    """Equation 8"""
    return ((M/1e10) ** 0.7) * np.exp(1.9 * z) / (np.exp(1.7 * (z - 2)) + np.exp(0.2 * (z - 2)))

def logMQ_z_alt_init(z):
    """Equation 10"""
    return 10.077 + 0.636 * z
    # return (13.077 + 0.636 * z) * (z <= 2) + (13.077 + 0.636 * 2) * (z > 2)

def pQ_Mz_alt_init(M, z):
    """Equation 9"""
    return 0.5*(1 - erf((np.log10(M) - logMQ_z_alt_init(z)) / 1.1))

def pmin(z):
    return 1 - ((z - 10) / 10) ** 2

def logMQ_z_alt(z):

    return (10.377 + 0.636*z)

def pQ_Mz_alt(M,z):

    return 0.5*(1 - erf((np.log10(M) - logMQ_z_alt(z)) / 1.5))

#def pQ_Mz_alt(M, z, Mq):
#     """Equation 9 with logMq fixed"""
    # return pmin(z) + (1 - pmin(z)) * pQ_Mz_alt_init(M, z)

    # return (1 - erf((np.log10(M) - (np.log10(Mq) - 0.85)) / 0.5)) #0.5*

def fml_t(t):
    """Equation 11"""
    return 0.046 * np.log(1 + t / 0.276)

def fml_dt(t, dt):
    """Equation 11"""
    return 0.046 * np.log((1 + (t + dt) / 0.276) / (1 + t / 0.276))

def draw_pQ(M, z):
    return np.random.uniform(0, 1, M.shape) > pQ_Mz_alt_init(M, z)

def draw_pQ_alt(M, z, f):
    return np.random.uniform(0, 1, M.shape) > pQ_Mz_alt(M, z) + f

#def draw_pQ_alt(M, z, f):
    #pq = np.minimum(1, pQ_Mz_alt(M, z) + f)
    #return np.random.uniform(0, 1, M.shape) < pq


def pQ_Mz_ft(M, z, isq, mqs, pq):    
    hold = draw_pQ(M[isq == False], z)
    mqs[isq == False] = M[isq == False] * hold

    isq[isq == False] = hold
    pq[isq == True] = pQ_Mz_alt(M[isq == True], z, mqs[isq == True])
    
    return isq, pq, mqs

def pmin_z(z):
    return 0


def pQ_Mz_ft2(M, z, isq, pq):
    frac = np.cumsum(isq) / np.arange(1, len(isq) + 1)
    isq[isq == False] = draw_pQ_alt(M[isq == False], z, frac[isq == False])
    pq[isq == True] = pmin_z(z)

    return isq, pq




# def sfr_Mz_alt(M, z, isq, mqs, pq):
#     isqs, pqs, mqs = pQ_Mz_ft(M, z, isq, mqs, pq)
#     return pqs * psi_Mz_alt(M, z), isqs, mqs, pqs

def sfr_Mz_alt(M, z, isq, pq):
    isq, pq = pQ_Mz_ft2(M, z, isq, pq)
    return pq * psi_Mz_alt(M, z) + 0.05 * np.random.choice([0,1],p=[0.95,0.05]) * psi_Mz_alt(M,z), isq, pq


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-dt','--dt',help='Time step (Myr)', default=0.5,type=float)
    parser.add_argument('-es','--early_step',help='Early Universe T_F step (Myr)', default=250,type=float)
    parser.add_argument('-ls','--late_step',help='Late Universe T_F step (Myr)', default=500, type=float)
    parser.add_argument('-ts','--tstart',help='Epoch to start the tracks (working backwards) (yr)', default = 0,type=float)
    parser.add_argument('-n','--n',help='Number of  galaxy',default=100,type=int)
    parser.add_argument('-t','--test',help='Is this a test?',action='store_true')
    args = parser.parse_args()
    return args

def script_worker(worker_args):
    args, tf, save_dir, zs, whole_fml_dt = [worker_args[i] for i in range(len(worker_args))]
    dt = args.dt
    N  = args.n
    #ages of the galaxies starting from 0(now) to tf (when galaxies are formed)
    ages = np.arange(0, tf + dt, dt)
    #mass of the galaxies (every galaxy starts with 1e6)
    m = 1e6 * np.ones(N)
    #quenching penalty of galaxies (initially set to 1 for all galaxies)
    pq = np.ones(N) 
    m_arr = []
    #galaxy is not quenched at first
    is_quenched = np.full(N, False)
    #basically is t, but specific for each galaxy, where ts is the time when the galaxy is formed minus the age of the galaxy. 
    # so when tf=ages, ts=0            
    ts = tf - ages
    # zs = whole_z[:len(ts)][::-1]
    # iterate over the ages (but using z)
    for i, z_t in tqdm(enumerate(zs), total = len(zs)):
        # pass
        # calculate mass created at each time step, and evaluate is the galaxy is quenched, also calculate the quenching penalty
        m_created, is_quenched, pq = sfr_Mz_alt(m, z_t, is_quenched, pq)
        m_created *= dt * 1e6
        #stacking the mass created at each time step
        if i == 0:
            m_formed = np.array([m_created.copy()])
        else:
            m_formed = np.vstack((m_formed, m_created.copy()))

        taus = ages[i::-1]
 
        # calculate the mass lost at each time step
        if i < 2:
            # Computing equation 11 and equation 12
            ml = fml_t(taus) * m_formed.T
            new_ml = np.sum(ml, axis = 1)
        else:
            new_ml = np.einsum('i, ij -> j', whole_fml_dt[i:0:-1], m_formed[:-1])
        
        m += m_created - new_ml
        m_arr.append(m.copy())
    
    m_formed = m_formed.T
    m_arr = np.array(m_arr).T
    final_age_weights = (1 - fml_t(taus)) * m_formed
    
    # df.to_hdf(os.path.join(save_dir,'SFHs_alt_ver2_%.1f_quenched.h5'%dt),key='m%3.0f'%tf)

    print("Saving to: ",os.path.join(save_dir,'SFHs_alt_%3.0f_%.1f_quenched_all_bursts.h5'%(tf,dt)))
    df = pd.DataFrame()
    #print(track)

    df = pd.DataFrame()
    # trace each galaxy, for each galaxy we trace the time, redshift, age, mass formed, final age weights and mass 
    # galaxy differ at when it was formed 
    for n in range(N):
        track = np.array([ts, zs, ages, m_formed[n], final_age_weights[n], m_arr[n]]).T

        df = pd.concat([df, pd.DataFrame(track,columns=['t','z','age','m_formed','final_age_weights','m_tot'],index=[n]*len(ages))])
    
    df.to_hdf(os.path.join(save_dir,'SFHs_alt_%3.0f_%.1f_quenched_all_bursts.h5'%(tf,dt)),key='main')

def main(args):
    #path to save the output
    save_dir = os.path.join(os.path.abspath("../../.."), 'mass_assembly/output_mass_assembly_adjusted')
    if args.test:
        save_dir = save_dir + '\\test_smaller_bursts'
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    #time step
    dt = args.dt
    #number of galaxies
    N = args.n

    if args.tstart ==0:
        #tstart = 0 means epoch to start the tracks (working backwards) is 0
        # making array of tf (time when galaxies are formed)
        tfs = np.concatenate([np.arange(1000,10000,args.late_step),np.arange(10000,13000,args.early_step)])
    elif args.tstart <10000:
        tfs = np.concatenate([np.arange(args.tstart,10000,args.late_step),np.arange(10000,13000,args.early_step)])
    else:
        tfs = np.arange(args.tstart,13000,args.early_step)
    #array of t, starting from 0 (now) to tf (when galaxies are formed)
    whole_t = np.arange(0, max(tfs) + dt, dt)
    #array of z
    whole_z = np.array([z_at_value(cosmo.lookback_time, t * u.Myr, zmin = 0) if t != 0 else 0 for t in whole_t])
    #calculate fractional mass lost since each epoch
    whole_fml_dt = fml_dt(whole_t, dt)
    #call worker function to pass the simulation into different processes, each process will run the simulation for each tf
    worker_args = [[args,tf,save_dir,whole_z[:int(tf//dt) + 1][::-1],whole_fml_dt] for tf in tfs]
    pool_size = 24
    pool = MyPool(processes=pool_size)
    for _ in tqdm(pool.imap_unordered(script_worker,worker_args),total=len(worker_args)):
        pass
    pool.close()
    pool.join()
    df =pd.DataFrame()
    for fn in glob.glob(os.path.join(save_dir,'SFHs_alt_*_%.1f_quenched_all_bursts.h5'%dt)):
        df=pd.read_hdf(fn)
        tf = os.path.split(fn)[-1].split('_')[2]
        df.to_hdf(os.path.join(save_dir,'SFHs_alt_%.1f_quenched_all_bursts.h5'%dt),key='%s'%tf)

if __name__=="__main__":
    main(parser())
