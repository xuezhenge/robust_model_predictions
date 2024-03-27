import os
import numpy as np
import pandas as pd
from pandas import read_csv
import time
import csv
import tqdm
import random
from scipy.integrate import odeint, ode,solve_ivp
from joblib import Parallel, delayed
import multiprocessing
import pdb
import matplotlib.pyplot as plt
import argparse
# import h5py
from sympy import symbols, Eq, solve
import math
import tqdm


parser = argparse.ArgumentParser()
parser.add_argument("--speciesid", type=int, 
    default=3, help="1, 2, 3, 4, 5")
parser.add_argument('--num_cores', type=int,
    default=3)
parser.add_argument('--K', type=int,
    default=50000000,help = "Carrying capacity")
parser.add_argument('--time_period', type=str,
     default="1980-2000", help="f5")
parser.add_argument('--pred_frac1', type=int,
     default=1, help="f5")
parser.add_argument('--pred_frac2', type=int,
     default=6, help="f5")
parser.add_argument('--start_idx', type=int,
    default=0)
parser.add_argument('--end_idx', type=int,
    default=15478)
parser.add_argument('--predation', type=str,
    default='P')
parser.add_argument('--dt', type=float,
    default=0.01)


args = parser.parse_args()
start_idx = args.start_idx
end_idx = args.end_idx
predation = args.predation
dt = args.dt

speciesid = args.speciesid
time_period = args.time_period
pred_frac1 = args.pred_frac1
pred_frac2 = args.pred_frac2
a_frac = pred_frac1*10**(-pred_frac2)

#processor
num_cores = args.num_cores 

# PARAMETER VALUES FOR LADYBIRD PREDATION
# attack rate
a_inst1 = 0.373
a_inst2 = 0.317
a_inst3 = 0.233
a_inst4 = 0.333
a_F = 0.353
a_M = 0.546
# handling time
Th_inst1 = 0.07708
Th_inst2 = 0.03596
Th_inst3 = 0.014125
Th_inst4 = 0.00725
Th_F = 0.010125
Th_M = 0.008291667
# transformation rate
Qp = 50
# Other parameters
K = args.K #carrying cacpacity of ladybird
theta = 0.5 # ratio of female ladybird to male ladybird

paras = pd.read_csv(f'../../outputs/species_parameters.csv')
# Function to get species parameter
def get_paras(paras,speciesid, stage):
    para = paras[(paras['speciesid'] == speciesid)&(paras['stage'] == stage)]
    m = para.m.item()
    CTmin = para.CTmin.item()
    TFmin = para.TFmin.item()
    CTmax = para.CTmax.item()
    TFmax = para.TFmax.item()
    Topt = para.Topt.item()
    q1 = para.q1.item()
    q2 = para.q2.item()
    b0 = para.b0.item()
    b1 = para.b1.item()
    b2 = para.b2.item()
    b3 = para.b3.item()
    return m,CTmin,TFmin,CTmax,TFmax,Topt,q1,q2,b0,b1,b2,b3 

#fecundity rate
mA_fec,CTminA_fec,TFminA_fec,CTmaxA_fec,TFmaxA_fec,ToptA_fec,q1L_fec,q2A_fec,b0A_fec,b1A_fec,b2A_fec,b3A_fec = get_paras(paras,speciesid = speciesid, stage='adult')

# Parameters for APHID instars developmental rate and mortality rate
mA_inst1,CTminA_inst1,TFminA_inst1,CTmaxA_inst1,TFmaxA_inst1,ToptA_inst1,q1A_inst1,q2A_inst1,b0A_inst1,b1A_inst1,b2A_inst1,b3A_inst1 = get_paras(paras,speciesid = speciesid, stage='instar1')
mA_inst2,CTminA_inst2,TFminA_inst2,CTmaxA_inst2,TFmaxA_inst2,ToptA_inst2,q1A_inst2,q2A_inst2,b0A_inst2,b1A_inst2,b2A_inst2,b3A_inst2 = get_paras(paras,speciesid = speciesid, stage='instar2')
mA_inst3,CTminA_inst3,TFminA_inst3,CTmaxA_inst3,TFmaxA_inst3,ToptA_inst3,q1A_inst3,q2A_inst3,b0A_inst3,b1A_inst3,b2A_inst3,b3A_inst3 = get_paras(paras,speciesid = speciesid, stage='instar3')
mA_inst4,CTminA_inst4,TFminA_inst4,CTmaxA_inst4,TFmaxA_inst4,ToptA_inst4,q1A_inst4,q2A_inst4,b0A_inst4,b1A_inst4,b2A_inst4,b3A_inst4 = get_paras(paras,speciesid = speciesid, stage='instar4')

# Parameters for LADYBIRD developmental rate and mortality rate
mL_egg,CTminL_egg,TFminL_egg,CTmaxL_egg,TFmaxL_egg,ToptL_egg,q1L_egg,q2L_egg,b0L_egg,b1L_egg,b2L_egg,b3L_egg = get_paras(paras,speciesid = 6, stage='egg')
mL_pupa,CTminL_pupa,TFminL_pupa,CTmaxL_pupa,TFmaxL_pupa,ToptL_pupa,q1L_pupa,q2L_pupa,b0L_pupa,b1L_pupa,b2L_pupa,b3L_pupa = get_paras(paras,speciesid = 6, stage='pupa')
mL_inst1,CTminL_inst1,TFminL_inst1,CTmaxL_inst1,TFmaxL_inst1,ToptL_inst1,q1L_inst1,q2L_inst1,b0L_inst1,b1L_inst1,b2L_inst1,b3L_inst1 = get_paras(paras,speciesid = 6, stage='instar1')
mL_inst2,CTminL_inst2,TFminL_inst2,CTmaxL_inst2,TFmaxL_inst2,ToptL_inst2,q1L_inst2,q2L_inst2,b0L_inst2,b1L_inst2,b2L_inst2,b3L_inst2 = get_paras(paras,speciesid = 6, stage='instar2')
mL_inst3,CTminL_inst3,TFminL_inst3,CTmaxL_inst3,TFmaxL_inst3,ToptL_inst3,q1L_inst3,q2L_inst3,b0L_inst3,b1L_inst3,b2L_inst3,b3L_inst3 = get_paras(paras,speciesid = 6, stage='instar3')
mL_inst4,CTminL_inst4,TFminL_inst4,CTmaxL_inst4,TFmaxL_inst4,ToptL_inst4,q1L_inst4,q2L_inst4,b0L_inst4,b1L_inst4,b2L_inst4,b3L_inst4 = get_paras(paras,speciesid = 6, stage='instar4')

# Parameters for LADYBIRD temperature dependent predation rate
mL_pred1,CTminL_pred1,TFminL_pred1,CTmaxL_pred1,TFmaxL_pred1,ToptL_pred1,q1L_pred1,q2L_pred1 = get_paras(paras,speciesid = 6, stage='pred1')[0:8]
mL_pred2,CTminL_pred2,TFminL_pred2,CTmaxL_pred2,TFmaxL_pred2,ToptL_pred2,q1L_pred2,q2L_pred2 = get_paras(paras,speciesid = 6, stage='pred2')[0:8]
mL_pred3,CTminL_pred3,TFminL_pred3,CTmaxL_pred3,TFmaxL_pred3,ToptL_pred3,q1L_pred3,q2L_pred3 = get_paras(paras,speciesid = 6, stage='pred3')[0:8]
mL_pred4,CTminL_pred4,TFminL_pred4,CTmaxL_pred4,TFmaxL_pred4,ToptL_pred4,q1L_pred4,q2L_pred4 = get_paras(paras,speciesid = 6, stage='pred4')[0:8]
mL_predF,CTminL_predF,TFminL_predF,CTmaxL_predF,TFmaxL_predF,ToptL_predF,q1L_predF,q2L_predF = get_paras(paras,speciesid = 6, stage='predF')[0:8]
mL_predM,CTminL_predM,TFminL_predM,CTmaxL_predM,TFmaxL_predM,ToptL_predM,q1L_predM,q2L_predM = get_paras(paras,speciesid = 6, stage='predM')[0:8]

CTmin_max = max(CTminA_fec,CTminA_inst1,CTminA_inst2,CTminA_inst3,CTminA_inst4,CTminL_egg,CTminL_inst1,CTminL_inst2,CTminL_inst3,CTminL_inst4,CTminL_pupa)
CTmin_min = min(CTminA_fec,CTminA_inst1,CTminA_inst2,CTminA_inst3,CTminA_inst4,CTminL_egg,CTminL_inst1,CTminL_inst2,CTminL_inst3,CTminL_inst4,CTminL_pupa)
CTmax_max = max(CTmaxA_fec,CTmaxA_inst1,CTmaxA_inst2,CTmaxA_inst3,CTmaxA_inst4,CTmaxL_egg,CTmaxL_inst1,CTmaxL_inst2,CTmaxL_inst3,CTmaxL_inst4,CTmaxL_pupa)

# Temperature performance curve
def thorneley_france(m,CTmin,TFmin,CTmax,TFmax,Topt,q1,q2,T):
    # Generic temperature-dependent function
    if (T>=(CTmin)) & (T<=(CTmax)):
        return m*(((T-(TFmin))**q1)*((TFmax-T)**q2))/(((Topt-TFmin)**q1)*((TFmax-Topt)**q2))
    else:
        return 0.0

# Mortality rate
def mor_rate(CTmin, CTmax,b0,b1,b2,b3,T):
    if (T>=(CTmin)) & (T<=(CTmax)):
        return b0 + b1*T + b2*T**2+ b3*T**3
    else:
        return 0.05

def fdpr(a,Th,Aden):
    #food_dependent_predation_rate
    a = a*a_frac
    return a/(1+a*Th*Aden)

def carring_capacity(Aden, K):
    return 1 - Aden/K

def fecL(Aden,a,Th,CTminL_predF,CTmaxL_predF,tdprF,T):
    # Fecundity rate of female ladybirds
    if (T>=(CTminL_predF)) & (T<=(CTmaxL_predF)):
        return tdprF*Aden*fdpr(a, Th, Aden)/Qp
    else:
        return 0

#all years daily temperature function
def get_temps(loc):
    temp_df = pd.read_csv(f'../../inputs/temperature_data/{time_period}/{loc}.csv')
    #Note: check temp_data
    temp_data = np.array(temp_df.Temp).reshape((20,365))
    return temp_data

def interpolated_temp(t,year,temp_data):
    temp = temp_data[year-1,:]
    temp_len = len(temp)
    t_tmp = int(np.floor(t))
    if t_tmp <= (temp_len-1):
        Temp_t = temp[t_tmp-1] + (temp[t_tmp] - temp[t_tmp-1])*(t - t_tmp)
    else:
        Temp_t = temp[temp_len-1]
    return Temp_t

def draw_multi_lines(years,x,pltt,folder,loc,xlabel,ylabel):
    t = np.arange(0,365*years,dt)
    plt.plot(t,x)
    plt.title(pltt + f': {loc}')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    fn = f'{loc}.png'
    out_file = os.path.join(folder, fn)
    if not os.path.isfile(out_file):
        plt.savefig(out_file,bbox_inches='tight')
    plt.close()

def draw_multi_scatters(x,y,folder,loc,xlabel,ylabel):
    plt.plot(x,y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('AL : loc = {}'.format(loc))
    fn = f'{loc}.png'.format(loc)
    out_file = os.path.join(folder, fn)
    if not os.path.isfile(out_file):
        plt.savefig(out_file,bbox_inches='tight')
    plt.close()

def writer_csv(rows, filename):
    with open(filename, "a+", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)

def batch(idx, fns,export_fns):
    tic = time.time()
    dataNP_folder, dataP_folder,plotAnp_folder,plotAp_folder, plotL_folder, plotAL_folder = export_fns
    #output file name
    fn = fns[idx]
    loc = fn.split(".csv", 1)[0]
    # lon,lat = loc.split("_", 1)
    # lon = float(lon)
    # lat = float(lat)
    temp_data = get_temps(loc)
    temp_fn = fn
    if predation == 'NP':
        out_file = os.path.join(dataNP_folder, temp_fn)
    else:
        out_file = os.path.join(dataP_folder, temp_fn)
    if os.path.exists(out_file):
        print('{} exits!!!'.format(temp_fn))
        return False

    #Model equations:
    def Solve_euler_model(year,var0,A_add,L_add,t_start,t_end,dt,predation):
        ts = np.arange(t_start,t_end,dt)
        n_t=len(ts)
        Aap = np.zeros([n_t]);A1 = np.zeros([n_t]); A2 = np.zeros([n_t]); A3 = np.zeros([n_t]); A4 = np.zeros([n_t])
        Legg = np.zeros([n_t]); L1 = np.zeros([n_t]); L2 = np.zeros([n_t]); L3 = np.zeros([n_t]); L4 = np.zeros([n_t]); Lpupa = np.zeros([n_t]); Lf = np.zeros([n_t]); Lm = np.zeros([n_t])
        Aden = np.zeros([n_t]); Lden = np.zeros([n_t]);Aborn = np.zeros([n_t]);Lborn = np.zeros([n_t])
        Aap[0],A1[0],A2[0],A3[0],A4[0],Legg[0],L1[0],L2[0],L3[0],L4[0],Lpupa[0],Lf[0],Lm[0] = var0
        Aden[0] = Aap[0] + A1[0] + A2[0] + A3[0] + A4[0]
        Lden[0] = Legg[0] + L1[0] + L2[0] + L3[0] + L4[0] + Lpupa[0] + Lf[0] + Lm[0]
        rate = np.zeros([n_t])
        num_change_A = 0; num_change_L = 0
        positive_dayA = 0; positive_dayL = 0
        for i in range(1, n_t):
            for j in np.arange(2):
                t = ts[i-1] #previous time step
                t_cur = ts[i] #current time step
                # Temperature at t:
                Temp_t = interpolated_temp(t,year,temp_data)
                if num_change_A == 0 and num_change_L == 0:
                    Aden_t, Aap_t, A1_t, A2_t, A3_t, A4_t = [A_add,A_add,0,0,0,0]
                    Lden_t, Legg_t, L1_t, L2_t, L3_t, L4_t, Lpupa_t, Lf_t, Lm_t = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                elif num_change_A == 1 and num_change_L == 0:
                    Aden_t = Aden[i-1];Aap_t = Aap[i-1]; A1_t = A1[i-1]; A2_t = A2[i-1]; A3_t = A3[i-1]; A4_t = A4[i-1]
                    if j == 0:
                        Lden_t, Legg_t, L1_t, L2_t, L3_t, L4_t, Lpupa_t, Lf_t, Lm_t = [L_add, 0, 0, 0, 0, 0, 0, L_add, 0]
                    else:
                        Lden_t, Legg_t, L1_t, L2_t, L3_t, L4_t, Lpupa_t, Lf_t, Lm_t = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                elif num_change_A == 1 and num_change_L == 1:
                    Aap_t = Aap[i-1]; A1_t = A1[i-1]; A2_t = A2[i-1]; A3_t = A3[i-1]; A4_t = A4[i-1]; Aden_t = Aden[i-1]
                    Legg_t = Legg[i-1]; L1_t = L1[i-1]; L2_t = L2[i-1]; L3_t = L3[i-1]; L4_t = L4[i-1]; Lpupa_t = Lpupa[i-1]; Lf_t = Lf[i-1]; Lm_t = Lm[i-1]; Lden_t = Lden[i-1]

                ## Temperature-dependent parameters of aphid
                #fecudity rate
                f_ap_t = thorneley_france(mA_fec,CTminA_fec,TFminA_fec,CTmaxA_fec,TFmaxA_fec,ToptA_fec,q1L_fec,q2A_fec, Temp_t)
                # development
                varphi_inst1_t = thorneley_france(mA_inst1,CTminA_inst1,TFminA_inst1,CTmaxA_inst1,TFmaxA_inst1,ToptA_inst1,q1A_inst1,q2A_inst1, Temp_t)
                varphi_inst2_t = thorneley_france(mA_inst2,CTminA_inst2,TFminA_inst2,CTmaxA_inst2,TFmaxA_inst2,ToptA_inst2,q1A_inst2,q2A_inst2, Temp_t)
                varphi_inst3_t = thorneley_france(mA_inst3,CTminA_inst3,TFminA_inst3,CTmaxA_inst3,TFmaxA_inst3,ToptA_inst3,q1A_inst3,q2A_inst3, Temp_t)
                varphi_inst4_t = thorneley_france(mA_inst4,CTminA_inst4,TFminA_inst4,CTmaxA_inst4,TFmaxA_inst4,ToptA_inst4,q1A_inst4,q2A_inst4, Temp_t)
                # mortality
                mu_inst1_t = mor_rate(CTminA_inst1, CTmaxA_inst1,b0A_inst1,b1A_inst1,b2A_inst1,b3A_inst1,Temp_t)
                mu_inst2_t = mor_rate(CTminA_inst2, CTmaxA_inst2,b0A_inst2,b1A_inst2,b2A_inst2,b3A_inst2,Temp_t)
                mu_inst3_t = mor_rate(CTminA_inst3, CTmaxA_inst3,b0A_inst3,b1A_inst3,b2A_inst3,b3A_inst3,Temp_t)
                mu_inst4_t = mor_rate(CTminA_inst4, CTmaxA_inst4,b0A_inst4,b1A_inst4,b2A_inst4,b3A_inst4,Temp_t)
                mu_ap_t = mor_rate(CTminA_fec, CTmaxA_fec,b0A_fec,b1A_fec,b2A_fec,b3A_fec,Temp_t)
                #carring capacity
                k_effect_t = carring_capacity(Aden_t,K)
                ## Temperature-dependent parameters of ladybird
                #Predartion rate as a function of temperature at maximal aphid density
                #temperature dependent predation rate
                tdprF = thorneley_france(mL_predF,CTminL_predF,TFminL_predF,CTmaxL_predF,TFmaxL_predF,ToptL_predF,q1L_predF,q2L_predF,Temp_t)/mL_predF
                tdprM = thorneley_france(mL_predM,CTminL_predM,TFminL_predM,CTmaxL_predM,TFmaxL_predM,ToptL_predM,q1L_predM,q2L_predM,Temp_t)/mL_predM
                tdpr_inst1 = thorneley_france(mL_pred1,CTminL_pred1,TFminL_pred1,CTmaxL_pred1,TFmaxL_pred1,ToptL_pred1,q1L_pred1,q2L_pred1,Temp_t)/mL_pred1
                tdpr_inst2 = thorneley_france(mL_pred2,CTminL_pred2,TFminL_pred2,CTmaxL_pred2,TFmaxL_pred2,ToptL_pred2,q1L_pred2,q2L_pred2,Temp_t)/mL_pred2
                tdpr_inst3 = thorneley_france(mL_pred3,CTminL_pred3,TFminL_pred3,CTmaxL_pred3,TFmaxL_pred3,ToptL_pred3,q1L_pred3,q2L_pred3,Temp_t)/mL_pred3
                tdpr_inst4 = thorneley_france(mL_pred4,CTminL_pred4,TFminL_pred4,CTmaxL_pred4,TFmaxL_pred4,ToptL_pred4,q1L_pred4,q2L_pred4,Temp_t)/mL_pred4
                # Fecundity rate of female ladybirds
                f_L_t =  fecL(Aden_t,a_F,Th_F,CTminL_predF,CTmaxL_predF,tdprF,Temp_t)
                #Stage-specific Development rates
                # #temperature dependent develoment rate fraction
                # tddr_t = thorneley_france(m_devL, TminL, TmaxL, Topt_devL, q1_devL, q2_devL, Temp_t)
                #Temperature-dependent development rates for egg and pupa
                delta_egg_t = thorneley_france(mL_egg,CTminL_egg,TFminL_egg,CTmaxL_egg,TFmaxL_egg,ToptL_egg,q1L_egg,q2L_egg,Temp_t)
                delta_pupa_t = thorneley_france(mL_pupa,CTminL_pupa,TFminL_pupa,CTmaxL_pupa,TFmaxL_pupa,ToptL_pupa,q1L_pupa,q2L_pupa,Temp_t)
                #Temperature-dependent developments rates at prey saturation
                delta_inst1_prey_saturation_t = thorneley_france(mL_inst1,CTminL_inst1,TFminL_inst1,CTmaxL_inst1,TFmaxL_inst1,ToptL_inst1,q1L_inst1,q2L_inst1,Temp_t)
                delta_inst2_prey_saturation_t = thorneley_france(mL_inst2,CTminL_inst2,TFminL_inst2,CTmaxL_inst2,TFmaxL_inst2,ToptL_inst2,q1L_inst2,q2L_inst2,Temp_t)
                delta_inst3_prey_saturation_t = thorneley_france(mL_inst3,CTminL_inst3,TFminL_inst3,CTmaxL_inst3,TFmaxL_inst3,ToptL_inst3,q1L_inst3,q2L_inst3,Temp_t)
                delta_inst4_prey_saturation_t = thorneley_france(mL_inst4,CTminL_inst4,TFminL_inst4,CTmaxL_inst4,TFmaxL_inst4,ToptL_inst4,q1L_inst4,q2L_inst4,Temp_t)
                # Mortality rate of various stages
                gamma_egg_t = mor_rate(CTminL_egg,CTmaxL_egg,b0L_egg,b1L_egg,b2L_egg,b3L_egg,Temp_t)
                gamma_inst1_t = mor_rate(CTminL_pupa,CTmaxL_pupa,b0L_pupa,b1L_pupa,b2L_pupa,b3L_pupa,Temp_t)
                gamma_inst2_t = mor_rate(CTminL_inst1,CTmaxL_inst1,b0L_inst1,b1L_inst1,b2L_inst1,b3L_inst1,Temp_t)
                gamma_inst3_t = mor_rate(CTminL_inst2,CTmaxL_inst2,b0L_inst2,b1L_inst2,b2L_inst2,b3L_inst2,Temp_t)
                gamma_inst4_t = mor_rate(CTminL_inst3,CTmaxL_inst3,b0L_inst3,b1L_inst3,b2L_inst3,b3L_inst3,Temp_t)
                gamma_pupa_t = mor_rate(CTminL_inst4,CTmaxL_inst4,b0L_inst4,b1L_inst4,b2L_inst4,b3L_inst4,Temp_t)
                gamma_f_t = gamma_pupa_t
                gamma_m_t = gamma_pupa_t
                # Predation rate for per aphid
                # Pt = P_perA_t*Aden_t
                P_perA_t = tdpr_inst1*L1_t*fdpr(a_inst1,Th_inst1,Aden_t) + tdpr_inst2*L2_t*fdpr(a_inst2,Th_inst2,Aden_t)+ tdpr_inst3*L3_t*fdpr(a_inst3,Th_inst3,Aden_t) + tdpr_inst4*L4_t*fdpr(a_inst4,Th_inst4,Aden_t) + tdprF*Lf_t*fdpr(a_F,Th_F,Aden_t) + tdprM*Lm_t*fdpr(a_M,Th_F,Aden_t)
                dA1_dt = f_ap_t*k_effect_t*Aap_t - mu_inst1_t*A1_t - varphi_inst1_t*A1_t - A1_t*P_perA_t
                dA2_dt = varphi_inst1_t*A1_t - mu_inst2_t*A2_t - varphi_inst2_t*A2_t - A2_t*P_perA_t
                dA3_dt = varphi_inst2_t*A2_t - mu_inst3_t*A3_t - varphi_inst2_t*A3_t - A3_t*P_perA_t
                dA4_dt = varphi_inst3_t*A3_t - mu_inst4_t*A4_t - varphi_inst4_t*A4_t - A4_t*P_perA_t
                dAap_dt = varphi_inst4_t*A4_t - mu_ap_t*Aap_t - Aap_t*P_perA_t

                #parameters for dL_dt
                delta_L1_t = delta_inst1_prey_saturation_t*tdpr_inst1*fdpr(a_inst1,Th_inst1,Aden_t)*Aden_t*Th_inst1
                delta_L2_t = delta_inst2_prey_saturation_t*tdpr_inst2*fdpr(a_inst2,Th_inst2,Aden_t)*Aden_t*Th_inst2
                delta_L3_t = delta_inst3_prey_saturation_t*tdpr_inst3*fdpr(a_inst3,Th_inst3,Aden_t)*Aden_t*Th_inst3
                delta_L4_t = delta_inst4_prey_saturation_t*tdpr_inst4*fdpr(a_inst4,Th_inst4,Aden_t)*Aden_t*Th_inst4

                rate[i] = delta_L3_t
                dLegg_dt = f_L_t*Lf_t - (gamma_egg_t + delta_egg_t)*Legg_t
                dL1_dt = delta_egg_t*Legg_t - delta_L1_t*L1_t - gamma_inst1_t*L1_t
                dL2_dt = delta_L1_t*L1_t - delta_L2_t*L2_t - gamma_inst2_t*L2_t
                dL3_dt = delta_L2_t*L2_t - delta_L3_t*L3_t - gamma_inst3_t*L3_t
                dL4_dt = delta_L3_t*L3_t - delta_L4_t*L4_t - gamma_inst4_t*L4_t
                dLpupa_dt = delta_L4_t*L4_t - delta_pupa_t*Lpupa_t - gamma_pupa_t*Lpupa_t
                dLf_dt = theta*delta_pupa_t*Lpupa_t - gamma_f_t*Lf_t
                dLm_dt = (1-theta)*delta_pupa_t*Lpupa_t - gamma_m_t*Lm_t

                Aap[i] = dt*dAap_dt + Aap_t
                A1[i] =dt*dA1_dt + A1_t
                A2[i] = dt*dA2_dt + A2_t
                A3[i] = dt*dA3_dt + A3_t
                A4[i] = dt*dA4_dt + A4_t
                Aden[i] = Aap[i] + A1[i] + A2[i] + A3[i] + A4[i]
                
                Legg[i] = dt*dLegg_dt + Legg_t
                L1[i] = dt*dL1_dt + L1_t
                L2[i] = dt*dL2_dt + L2_t
                L3[i] = dt*dL3_dt + L3_t
                L4[i] = dt*dL4_dt + L4_t
                Lpupa[i] = dt*dLpupa_dt + Lpupa_t
                Lf[i] = dt*dLf_dt + Lf_t
                Lm[i] = dt*dLm_dt + Lm_t
                Lden[i] = Legg[i] + L1[i] + L2[i] + L3[i] + L4[i] +  Lpupa[i] + Lf[i] + Lm[i]

                if Aap[i] < 0: Aap[i] = 0
                if A1[i] < 0: A1[i] = 0
                if A2[i] < 0: A2[i] = 0
                if A3[i] < 0: A3[i] = 0
                if A4[i] < 0: A4[i] = 0

                if Legg[i] < 0: Legg[i] = 0
                if L1[i] < 0: L1[i] = 0
                if L2[i] < 0: L2[i] = 0
                if L3[i] < 0: L3[i] = 0
                if L4[i] < 0: L4[i] = 0
                if Lpupa[i] < 0: Lpupa[i] = 0
                if Lf[i] < 0: Lf[i] = 0
                if Lm[i] < 0: Lm[i] = 0
                # set start condition
                if Aden[i] <= Aden_t and num_change_A == 0 and j == 0:
                    Aden[i] = 0; Aap[i] = 0; A1[i] = 0; A2[i] = 0; A3[i] = 0; A4[i] = 0
                    break
                if Aden[i] > Aden_t and num_change_A == 0 and j == 0:
                    positive_dayA = positive_dayA + 1
                    if positive_dayA == 5:
                        num_change_A = 1
                    break
                if num_change_A == 1 and num_change_L == 0 and Lden[i] > Lden_t and j == 0:
                    positive_dayL = positive_dayL + 1
                    if positive_dayL == 5:
                        num_change_L = 1
                    break
                else: continue
                # if num_change_A == 1 and num_change_L == 1:
                #     break
            # print(i,Aden[i],Lden[i])
            #End the simulation when one of them enter the overwintering period
            # if Temp_t < CTmin_max and num_change_A == 1 and num_change_L == 1 and i > 18250:
            #     break
            # if num_change_A == 1 and num_change_L == 0:
            #     A_end = Aden[i];L_end = 0
            # if num_change_A == 0 and num_change_L == 0:
            #     A_end = 0;L_end = 0
            # if A_end < 0.001: A_end = 0
            # if L_end < 0.001: L_end = 0
        outputs = [Aden, Lden, rate]
        return outputs

    #temp = temps(a,w)
    var0 = [0,0,0,0,0,0,0,0,0,0,0,0,0]
    A_add = 10000000
    if predation == 'P':
        L_add = 50000
    else:
        L_add = 0
    years = 20

    Adens_p = [];Ldens_p = []
    def yearly_output(A_add,L_add,year):
        Temps = temp_data[year-1,:]
        Temp_max = np.max(Temps)
        Temp_min = np.min(Temps)
        if Temp_max <= CTmin_min or Temp_min >= CTmax_max:
            Aden_p = np.zeros([int(365/dt)])
            Lden_p = np.zeros([int(365/dt)])
        else:
            outputs_p = Solve_euler_model(year,var0,A_add,L_add,t_start = 0, t_end = 365,dt=dt,predation = 'P')
            Aden_p = outputs_p[0]; Lden_p = outputs_p[1]; rate = outputs_p[2]
        Adens_p.extend(Aden_p)
        Ldens_p.extend(Lden_p)
        return Adens_p,Ldens_p

    for year in np.arange(years):
        year = year + 1
        yearly_output(A_add,L_add,year)

    Adens_p= np.array(Adens_p)
    Ldens_p= np.array(Ldens_p)
        
    AL = np.vstack([Adens_p,Ldens_p])
    AL = pd.DataFrame(np.transpose(AL))
    AL.columns = ['Aden_dt', 'Lden_dt']
    if not os.path.isfile(out_file):
        AL.to_csv(out_file,header=True)
    # if predation == 'NP':
    #     draw_multi_lines(years,Adens_p,'Anp',plotAnp_folder,loc,xlabel = 'Day', ylabel = 'Aphid Population Abundance')
    # else:
    #     draw_multi_lines(years,Adens_p,'Ap',plotAp_folder,loc,xlabel = 'Day', ylabel = 'Aphid Population Abundance')
    #     draw_multi_lines(years,Ldens_p,'L',plotL_folder,loc,xlabel = 'Day', ylabel = 'Ladybird Population Abundance')        
    #     draw_multi_scatters(Adens_p,Ldens_p,plotAL_folder,loc, xlabel = 'Aphid Population Abundance', ylabel = 'Ladybird Population Abundance')
    
    toc = time.time()
    print(temp_fn + " " + "Elapsed time: {}".format(toc - tic))
    return True

def folders(time_period):
    fold_dir = f"../../outputs/species{speciesid}/pred_frac{pred_frac1}_{pred_frac2}/exports_{time_period}_dt{dt}"
    export_fns = []
    for folder_name in ["data_np", "data_p","plotAnp", "plotAp", "plotL", "plotAL"]:
        folder = os.path.join(fold_dir, "{}".format(folder_name))
        if os.path.exists(folder) == False:
            os.makedirs(folder, exist_ok=True)
        export_fns += [folder]
    return export_fns

if __name__ == '__main__':
    export_fns = folders(time_period)
    num_idxs = 15478
    id_list = np.arange(num_idxs)[start_idx:end_idx]
    fns = os.listdir(f'../../inputs/temperature_data/{time_period}')
    for idx in np.arange(num_idxs):
        #batch(idx,fns,export_fns)
        processed_list = Parallel(n_jobs=num_cores)(delayed(batch)(idx,fns,export_fns) for idx in range(num_idxs))