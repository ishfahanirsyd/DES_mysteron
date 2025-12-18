#!/usr/bin/env python

import numpy as np
import pandas as pd
from tqdm import tqdm
import glob
import sys, os
# CLi
# Commented out, as we will run this from other directories
#sys.path.insert(0, os.path.abspath("../../..")) 
from dust_extinction.parameter_averages import F19
from AURA.simulations.spectral_utils import load_spectrum, convert_escma_fluxes_to_griz_mags,interpolate_SFH,interpolate_SFH_pegase
from AURA.simulations.synspec import SynSpec, phi_t_pl
from AURA.utils.utils import MyPool
import argparse
from astropy.cosmology import FlatLambdaCDM
import warnings
from astropy.utils.exceptions import AstropyWarning
from tables import NaturalNameWarning
import pysynphot as S
S.locations.default_cdbs = os.path.join(os.environ["DESSIMS"], "grp", "redcat", "trds")
print("loc",S.locations.default_cdbs)
warnings.filterwarnings('ignore', category=NaturalNameWarning)

np.seterr(all='ignore')
warnings.simplefilter('ignore', category=AstropyWarning)

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
cosmo = FlatLambdaCDM(70,0.3)
bc03_flux_conv_factor = 3.12e7

# DTD parameters from W21
beta_x1hi = -1.68
norm_x1hi = 0.51E-13
beta_x1lo = -0.79
norm_x1lo = 1.19E-13
beta = -1.13
#beta=-1.5
dtd_norm = 2.08E-13

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',help='Filename for input SFHs', default='%s/SFH_mpi/sfh_25_50/SFHs_alt_0.5_quenched_all_bursts.h5'  % (os.environ["mass_assembly"]))
    parser.add_argument('-o','--output',help='Output directory',default='test')
    parser.add_argument('-z','--z',help='Redshift',default=0.5,type=str)
    parser.add_argument('-zl','--zlo',help='Redshift lower end',default=0.00,type=float)
    parser.add_argument('-zh','--zhi',help='Redshift upper end',default=1.1,type=float)
    parser.add_argument('-zs','--zstep',help='Redshift step',default=0.05,type=float)
    parser.add_argument('-al','--av_lo',help='Lowest Av',default=0,type=float)
    parser.add_argument('-ah','--av_hi',help='Highest Av',default=1.5,type=float)
    parser.add_argument('-na','--n_av',help='Av step',default=20,type=int)
    parser.add_argument('-at','--av_step_type',help='Av step type (lin or log)',default='lin')
    parser.add_argument('-u','--logU',help='Ionisation parameter',default=-2,type=float)
    parser.add_argument('-tr','--time_res',help='SFH time resolution',default=1,type=int)
    parser.add_argument('-t','--templates',help='Template library to use [BC03, PEGASE]',default='BC03',type=str)
    parser.add_argument('-tf','--templates_fn',help='Filename of templates',type=str,default='None')
    parser.add_argument('-ne','--neb',action='store_true')
    parser.add_argument('-b','--beta',help='Absolute value of the slope of the DTD',default=1.14,type=float)
    args, unknown = parser.parse_known_args() 
    return args

def sed_worker(worker_args):
    #get the parameters from the worker_args
    sfh_df,args,av_arr,z,tf,s,bc03_logt_float_array,counter= [worker_args[i] for i in range(8)]
    outputDir= os.path.join(os.environ["hostlib"],args.output)
    print('Starting %i'%counter)
    try:
        results = []
        # iterate through the SFHs for galaxies of different final masses
        # important! here only use one SFH 
        ## CLi We might want to explore options to use more than one
        i = np.random.randint(0,len(sfh_df.index.unique()))
        sfh_iter_df = sfh_df.loc[i]
        #print(sfh_iter_df)
        # get total mass of galaxy in each SFH (last row)
        mtot=sfh_iter_df['m_tot'].iloc[-1]
        #print('mtot', mtot)
        # get the age of the galaxy in each SFH (last row)
        age = sfh_iter_df['age'].iloc[-1]
        #print('Mass: ',np.log10(mtot),'age: ',age)
        # sSFR is the sum of the mass formed in the last 500 Myr divided by the total mass of the galaxy 
        # still not understand why this is divided by 250
        ssfr = np.sum(sfh_iter_df['m_formed'].iloc[-500:])/((250*1E+6)*mtot)
        sfr = ssfr*mtot
        sfh_iter_df['stellar_age'] = sfh_iter_df.age.values[::-1]
        ages = sfh_iter_df['stellar_age']/1000
        # calculating the DTD
        dtd_x1hi = phi_t_pl(ages,0.04,beta_x1hi,norm_x1hi)
        pred_rate_x1hi =np.sum(sfh_iter_df['m_formed']*dtd_x1hi)
        dtd_x1lo = phi_t_pl(ages,0.04,beta_x1lo,norm_x1lo)
        pred_rate_x1lo =np.sum(sfh_iter_df['m_formed']*dtd_x1lo)
        dtd_total =phi_t_pl(ages,0.04,-1*args.beta,dtd_norm)
        SN_age_dist = sfh_iter_df['m_formed']*dtd_total
        # get rate of SNe per unit mass formed
        pred_rate_total = np.sum(SN_age_dist) #R_G
        # weighting (in paper; volumetric rate of SNe) that then used to select reprentative galaxy to form SN Ia
        mwsa = np.average(sfh_iter_df['stellar_age'],weights=sfh_iter_df['m_formed']/mtot)
        # get the value of Rv based on the mass of the galaxy (Salim. 2018)
        if np.log10(mtot)<=9.5:
            mu_Rv = 2.61
        elif 9.5 <np.log10(mtot)<=10.5:
            mu_Rv = 2.99
            #avs_SBL =np.clip(np.random.normal(av_means_mhi(np.log10(mtot)),av_sigma(np.log10(mtot)),size=20),a_min=0,a_max=None)
        else:
            mu_Rv = 3.47
            #avs_SBL = np.clip(np.random.normal(av_means_mlo,av_sigma(np.log10(mtot)),size=20),a_min=0,a_max=None)
        #choose spectral template, which then used to get SED based on SFH, metallicity, IMF 
        if args.templates == 'BC03':
            sfh_coeffs_PW21 = interpolate_SFH(sfh_iter_df,mtot,bc03_logt_float_array)
            #print('sfhprint', len(sfh_coeffs_PW21))
            template=None
        elif args.templates =='PEGASE':
            if args.templates_fn =='None':
                templates = pd.read_hdf('/media/data3/wiseman/des/AURA/PEGASE/templates.h5',key='main')
            else:
                templates = pd.read_hdf(args.templates_fn,key='main')
            sfh_coeffs_PW21 = interpolate_SFH_pegase(sfh_iter_df,templates['time'],mtot,templates['m_star'])
        
        arr = np.zeros((len(ages),2))
        arr[:,0] = ages
        arr[:,1] = SN_age_dist
        #print(arr)
        np.savetxt('%s/all_model_params_quench_%s_z_%.5f_rv_rand_full_age_dists_neb_U%.2f_res_%i_beta_%.2f_%.1f.dat' % (outputDir,args.templates, z, args.logU, args.time_res, args.beta, tf),arr)

        for Av in av_arr:
            Rv = np.min([np.max([2.0,np.random.normal(mu_Rv,0.5)]),6.0])
            delta='None'
            #if args.templates =='PEGASE':
            #    sfh_coeffs_PW21 = None
            #    template = pd.read_hdf('/media/data3/wiseman/des/AURA/PEGASE/templates_analytic_orig_%i.h5' % tf,
            #                           key='main')
            galid_string = 'all_model_params_quench_%s_z_%.5f_rv_rand_full_age_dists_neb_U%.2f_res_%i_beta_%.2f_%.1f_%.3f_%.3f'%(args.templates,z,args.logU,args.time_res,args.beta,tf,Av,Rv)

            ## CLi
            U_R,fluxes,colours,EW_OII= s.calculate_model_fluxes_pw(z,sfh_coeffs_PW21,dust={'Av':Av,'Rv':Rv,'delta':'none','law':'CCM89'},
                                                    neb=args.neb,logU=args.logU,mtot=mtot,age=age,specsavename=galid_string,savespec=True,sfr=sfr)
            
            #U_R,fluxes,colours= s.calculate_model_fluxes_pw(z,sfh_coeffs_PW21,dust={'Av':Av,'Rv':Rv,'delta':'none','law':'CCM89'},
            #                                        neb=args.neb,logU=args.logU,mtot=mtot,age=age,specsavename=galid_string,savespec=True)
            
            # CLi moved from a list containing a single float to a float
            #obs_flux  = list(fluxes.values())#+cosmo.distmod(z).value
            desg,desr,desi,desz  = (fluxes[i][0] for i in fluxes.keys())
            #U,B,V,R,I,sdssu,sdssg,sdssr,sdssi,sdssz = (colours[i] for i in colours.keys())
            U,B,V,R,I,sdssu,sdssg,sdssr,sdssi,sdssz = (colours[i][0] for i in colours.keys())
            

            results.append([z, mtot, ssfr, mwsa, Av, Rv, delta, U_R[0], pred_rate_x1hi, pred_rate_x1lo, pred_rate_total, tf, desg, desr, desi, desz, U, B, V, R, I, sdssu, sdssg, sdssr, sdssi, sdssz, galid_string,EW_OII])

        df = pd.DataFrame(results,columns=['z','mass','ssfr','mean_age','Av','Rv','delta','U_R','pred_rate_x1_hi',
                                                'pred_rate_x1_lo','pred_rate_total','t_f',
                                                'm_g','m_r','m_i','m_z','U','B','V','R','I','sdssu','sdssg','sdssr','sdssi','sdssz','galid_spec','EW_OII'])


    except Exception as e:
        print(f"Error in worker {counter}: {e}")
        return
    print('Saving %i'%counter)
    df.to_hdf('%s/all_model_params_quench_%s_z%.5f_%.5f_av%.2f_%.2f_rv_rand_full_age_dists_neb_U%.2f_res_%i_beta_%.2f_%.5f_%i_sdss_u_r.h5' % 
              (outputDir, args.templates, args.zlo, args.zhi, av_arr[0], av_arr[-1],args.logU, args.time_res, args.beta, z, tf),key='main')
    print('Returning %i'%counter)
    return


def run(args):
    # DES filter objects

    filt_dir = os.environ["filters"]
    filt_obj_list = [
        load_spectrum(filt_dir+'/decam_g.dat'),
        load_spectrum(filt_dir+'/decam_r.dat'),
        load_spectrum(filt_dir+'/decam_i.dat'),
        load_spectrum(filt_dir+'/decam_z.dat'),
    ]

    nfilt = len(filt_obj_list)

#    aura_dir = '%s/templates/' % os.environ["DESSIMS"]
    aura_dir = '%s/' % os.environ["AURA"]
    dessims_dir= '%s/' % os.environ["DESSIMS"]
    outputDir= os.path.join(os.environ["hostlib"],args.output)
    #------------------------------------------------------------------------
    # BC03 SSPs as mc_spec Spectrum objects
    f1 = open(dessims_dir+'templates/bc03/bc03_logt_list.dat')
    if args.templates =='BC03':
        bc03_logt_list = [x.strip() for x in f1.readlines()]
        f1.close()
        bc03_logt_array = np.array(bc03_logt_list)
        ntemp = len(bc03_logt_array)
        bc03_logt_float_array =np.array([float(x) for x in (bc03_logt_array)])
        #print('bc03',bc03_logt_float_array)
        # CLi
        # bc03_dir = '%s/templates/bc03/bc03_ssp_templates_generated/' % os.environ["DESSIMS"]
        bc03_dir = '%s/AURA/bc03/bc03_ssp_templates_generated/' % os.environ["DESSIMS"]
        template_obj_list = []
        nLy_list = []
        for i in range(ntemp):
            bc03_fn = '%sbc03_chabrier_z02_%s.spec' % (bc03_dir, bc03_logt_list[i])
            new_template_spec =  load_spectrum(bc03_fn)
            template_obj_list.append(new_template_spec)

        s = SynSpec(template_obj_list = template_obj_list, neb=args.neb)
        neb=args.neb
    elif args.templates=='PEGASE':
        s = SynSpec(library='PEGASE',template_dir = '/media/data3/wiseman/des/AURA/PEGASE/templates/',neb=args.neb)

        neb=args.neb
    store = pd.HDFStore(args.input,'r')
    #ordered_keys = np.sort([int(x.strip('/')) for x in store.keys()])
    ordered_keys = np.sort([int(x.strip('/').replace('m', '')) for x in store.keys() if x.strip('/').replace('m', '').isdigit()])
    print("keys", ordered_keys)

    #z_array = [float(z) for z in args.z.split(',')]
    z_array = np.arange(args.zlo,args.zhi,args.zstep)
    print("z", z_array)
    if args.av_step_type == 'lin':
        av_arr = np.linspace(args.av_lo,args.av_hi,args.n_av)
    elif args.av_step_type == 'log':
        av_arr = np.logspace(args.av_lo,args.av_hi,args.n_av)

    for z in tqdm(z_array):
        print('Making hostlib for z=%.2f'%z)
        # call worker function to pass the simulation into different processes (only use one SFH)
        worker_args = []
        #for counter,tf in enumerate(ordered_keys): 
        for counter,tf in enumerate(ordered_keys[::-1][np.arange(0,len(ordered_keys),args.time_res)]):   # Iterate through the SFHs for galaxies of different final masses
            #with pd.HDFStore(os.path.join(save_dir, 'SFHs_alt_ver2_%.1f_quenched_m%3.0f.h5' % (dt, tf)), 'r') as store:
            sfh_df = store['%3.0f' % tf]

            #sfh_df = store['/SFHs_alt_ver2_0.5_quenched_m'+str(tf)+'.h5']
            #filter SFH file so that it only contains galaxies with redshift>z
            sfh_df = sfh_df[sfh_df['z']>z]

            # in parallel, only use one SFH
            if len(sfh_df)>0:
                worker_args.append([sfh_df,args,av_arr,z,tf,s,bc03_logt_float_array,counter])
        pool_size = 8
        #pool_size=1 # For testing
        pool = MyPool(processes=pool_size)
        print('Sending %i jobs'%len(worker_args))
        results =pool.map_async(sed_worker,worker_args)

        pool.close()
        pool.join()

    # Initialize empty DataFrame
    all_df = pd.DataFrame()

    # Loop over files
    for fn in glob.glob(os.path.join(
        outputDir,
        'all_model_params_quench_BC03_z*_*_av*_*_rv_rand_full_age_dists_neb_U-2.00_res_*_beta_1.14_*_*_sdss_u_r.h5')):
        print("Reading file:", fn)
        df = pd.read_hdf(fn, key='/main')  # read each file
        all_df = pd.concat([all_df, df], ignore_index=True)  # accumulate

    # Save all combined data to a single key
    all_df.to_hdf(
        '%s/all_model_params_quench_%s_z_rv_rand_full_age_dists_neb_U%.2f_res_%i_beta_%.2f.h5' % (
            outputDir, args.templates, args.logU, args.time_res, args.beta),
        key='main',
        mode='w'  # overwrite existing file
    )


    print("Done!")
if __name__=="__main__":
    args = parser()
    run(args)
