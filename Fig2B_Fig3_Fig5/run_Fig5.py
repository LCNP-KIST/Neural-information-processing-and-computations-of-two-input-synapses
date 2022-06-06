#!/usr/bin/env python

# --Imports--
import sys
import os
import math
from time import time
import matplotlib.pyplot as plt
import numpy as np
import neuron
from neuron import h

import config_neuronmodel
from neuronmodel import *
from config_model_stim import *
from sim import *
from htools import *
import datetime
import multiprocessing as mp
import ctypes
import itertools


def pexec(args):
    ### seed for trial.
    trial_n = 0

    ### python seed for x1, x2 generation.
    random.seed(trial_n)

    ### x2 location arguments. use 'b' and 'c'.
    x2loc = 'b'

    ### arguments
    alphaidx = 0
    
    ### alpha parameter
    alphali = np.log2(np.array([i for i in range(1,11)]))
    alphali = alphali/np.max(alphali)
    alpha = alphali[alphaidx]
    
    ### Output directory
    path = 'data1/'
    identifier = '%s_%s'%(args[0],args[1])
    savepath = '%s%s'%(path,identifier)
    if not os.path.exists(path):
        os.mkdir(path)
    if not os.path.exists(savepath):
    	os.mkdir(savepath)
    
    ###  Synaptic mechanism  ###
    #BCM parameters
    P0_BCMTHRESHOLD = 0.003
    D0_BCMTHRESHOLD = 0.001
    SCOUNTTAU_BCMTHRESHOLD = 1500   # (ms) averaging time constant for postsynaptic spike count
    ALPHA_BCMTHRESHOLD = 62.5 # scaling constant, ALPHA=0 for non-BCM model
    PTAU = 20
    DTAU = 70
    TAU1 = 0.2
    TAU2 = 2.5
    WMAX = 1

    #background synapse frequency
    bg_freq = 1.0
    bg_interval = 1.0/f3*1000.0
    
    if args[0] == 0:
        P0_BCMTHRESHOLD = np.linspace(0, 0.006, num_r_glo.value)[args[1]]
    if args[0] == 1:
        D0_BCMTHRESHOLD = np.linspace(0, 0.002, num_r_glo.value)[args[1]]
    if args[0] == 2:
        SCOUNTTAU_BCMTHRESHOLD = np.linspace(0, 3000, num_r_glo.value)[args[1]]
    if args[0] == 3:
        ALPHA_BCMTHRESHOLD = np.linspace(12.5, 137.5, num_r_glo.value)[args[1]]
    if args[0] == 4:
        PTAU = np.linspace(0, 40, num_r_glo.value)[args[1]]
    if args[0] == 5:
        DTAU = np.linspace(0, 140, num_r_glo.value)[args[1]]
    if args[0] == 6:
        TAU1 = np.linspace(0.04, 0.44, num_r_glo.value)[args[1]]
    if args[0] == 7:
        TAU2 = np.linspace(0, 5, num_r_glo.value)[args[1]]
    if args[0] == 8:
        WMAX = np.linspace(0.4, 2.4, num_r_glo.value)[args[1]]
    
    
    SCOUNT0_BCMTHRESHOLD = 1/ALPHA_BCMTHRESHOLD
    BCMSLOW_BCMTHRESHOLD = 0   # Cliff's slowing of BCM rate, 0=no change, 0.5=50% slowing, 0.2=20% slowing etc
    
    SCEN = 11
    NO_REPS = 1
    RESET = True
    DT=0.1 # ms, set the integration step
    POST_AMP = 0.3# nA, amplitude of current injection to trigger the AP/bAP
    WARM_UP=1000 # ms
    AP_DELAY = 0 #ms, AP needs 4ms after stimulation to become initiated
    
    my_rawdata = {}
    
    #inputs
    freq = params['Input']['freq'].value
    wee = params['Input']['wee'].value
    wee_strong = params['Input']['wee_strong'].value
    
    # Synapses
    scen = SCEN
    delta_t = params['STDP']['delta_t'].value
    
    sim_params = params['sim']
    
    source = False
    
    cell = Neuron()
    if args[0] == 9:
        cell.ena = np.linspace(0, 120, num_r_glo.value)[args[1]]
    if args[0] == 10:
        cell.ek = np.linspace(0, -160, num_r_glo.value)[args[1]]
    if args[0] == 11:
        cell.eca = np.linspace(24, 264, num_r_glo.value)[args[1]]
    if args[0] == 12:
        cell.gna = np.linspace(0.0018, 0.0198, num_r_glo.value)[args[1]]
    if args[0] == 13:
        cell.gk = np.linspace(0.002, 0.022, num_r_glo.value)[args[1]]
    if args[0] == 14:
        cell.gk_kap = np.linspace(0.0058, 0.0638, num_r_glo.value)[args[1]]
    if args[0] == 15:
        cell.slope = np.linspace(1, 11, num_r_glo.value)[args[1]]
    if args[0] == 16:
        cell.gsca = np.linspace(0.3, 3.3, num_r_glo.value)[args[1]]
    if args[0] == 17:
        cell.git2 = np.linspace(0.001, 0.011, num_r_glo.value)[args[1]]
    if args[0] == 18:
        cell.gbar_kca = np.linspace(0.5, 5.5, num_r_glo.value)[args[1]]
    if args[0] == 19:
        cell.gna_ais = np.linspace(0.06, 0.66, num_r_glo.value)[args[1]]
    if args[0] == 20:
        cell.gna_ais_shifted = np.linspace(0.06, 0.66, num_r_glo.value)[args[1]]
        
    sim = Simulation(cell,sim_params)
    sim.dt = DT
    sim.v_init = -70
    
    ### Total simulation time, in ms
    total_time = 1000*30
    
    ### lists of synapse-related objects
    synapses=[]
    exstims=[]
    gli=[]
    ili=[]
    wli=[]
    netcons=[]
    psli=[]
    
    nsteps = int(total_time/DT)+1
    
    
    ### basic hidden state parameters (transition rates)
    r_on_slow = 16.7/1000.0
    r_on = r_on_slow
    r_off = r_on*2.0
    q_on1 = 0.05
    q_off1 = 0.001
    
    ### generate coupled states x1 and x2
    # print("creating x states...")
    x1_state = np.zeros(nsteps)
    x2_state = np.zeros(nsteps)
    generate_coupled_states(x1_state,x2_state,r_on,r_off,alpha,DT)
    # print("done.")
    
    ### pre-synaptic spikes
    computed_spikes=[]
    tvecli=[]
    vvli=[]
    
    if (scen == 11):
    
        bcm = h.BCMthreshold(cell.soma(0.5))
        bcm.p0 = P0_BCMTHRESHOLD
        bcm.d0 = D0_BCMTHRESHOLD
        bcm.scount0 = SCOUNT0_BCMTHRESHOLD
        bcm.scounttau = SCOUNTTAU_BCMTHRESHOLD
        bcm.alpha = ALPHA_BCMTHRESHOLD
    
        for si in range(1):
            syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(cell.branches[0](0.1))
            ## synapse parameters taken from dg model -- might need reconfig
            syn.tau1 = TAU1
            syn.tau2 = TAU2
            syn.e = 0 
            syn.ptau = PTAU
            syn.dtau = DTAU
            syn.start = 0
            syn.tau2 = WMAX
    
            vv1 = h.Vector(np.array([-1.0 for j in range(nsteps)]))
            tvec1 = h.Vector(np.arange(0,total_time+DT,DT))
            event_times = compute_spike_train(x1_state,q_on1,q_off1,DT)
            computed_spikes.append(event_times)
            for et in event_times:
                ti = int(et / DT)
                try:
                    vv1[ti] = 1.0
                    vv1[ti+1] = 1.0
                    vv1[ti+2] = 1.0
                except IndexError:
                    pass
            spik1 = h.strain()
            vv1.play(spik1._ref_x,tvec1,False)
            tvecli.append(tvec1)
            vvli.append(vv1)
    
            syn_g=h.Vector().record(syn._ref_g)
            syn_i=h.Vector().record(syn._ref_i)
            syn_w=h.Vector().record(syn._ref_wtrack)
    
            h.setpointer(bcm._ref_d, "d", syn)
            h.setpointer(bcm._ref_p, "p", syn)
            weight = wee
            nc = h.NetCon(spik1,syn,0,0,weight)
            post_spikes = h.Vector()
            nc.record(post_spikes)
            synapses.append(syn)
            exstims.append(spik1)
            gli.append(syn_g)
            ili.append(syn_i)
            wli.append(syn_w)
            netcons.append(nc)
            psli.append(post_spikes)
    
        for si in range(1):
            if x2loc == 'b':
                syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(cell.branches[0](0.11))
            elif x2loc == 'c':
                syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(cell.branches[0](0.09))
            else:
                raise Exception("Error in  x2loc")
            ## synapse parameters taken from dg model -- might need reconfig
            syn.tau1 = TAU1
            syn.tau2 = TAU2
            syn.e = 0 
            syn.ptau = PTAU
            syn.dtau = DTAU
            syn.start = 0
            syn.tau2 = WMAX
    
            vv1 = h.Vector(np.array([-1.0 for j in range(nsteps)]))
            tvec1 = h.Vector(np.arange(0,total_time+DT,DT))
            event_times = compute_spike_train(x2_state,q_on1,q_off1,DT)
            computed_spikes.append(event_times)
            for et in event_times:
                ti = int(et / DT)
                try:
                    vv1[ti] = 1.0
                    vv1[ti+1] = 1.0
                    vv1[ti+2] = 1.0
                except IndexError:
                    pass
            spik1 = h.strain()
            vv1.play(spik1._ref_x,tvec1,False)
            tvecli.append(tvec1)
            vvli.append(vv1)
    
            syn_g=h.Vector().record(syn._ref_g)
            syn_i=h.Vector().record(syn._ref_i)
            syn_w=h.Vector().record(syn._ref_wtrack)
    
            h.setpointer(bcm._ref_d, "d", syn)
            h.setpointer(bcm._ref_p, "p", syn)
            weight = wee
            nc = h.NetCon(spik1,syn,0,0,weight)
            post_spikes = h.Vector()
            nc.record(post_spikes)
            synapses.append(syn)
            exstims.append(spik1)
            gli.append(syn_g)
            ili.append(syn_i)
            wli.append(syn_w)
            netcons.append(nc)
            psli.append(post_spikes)
    
        locs=[0.1*i for i in range(1,11)]
        for si in range(10):
            syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(cell.branches[0](locs[si]))
            ## synapse parameters taken from dg model -- might need reconfig
            syn.tau1 = 0.2
            syn.tau2 = 2.5
            syn.e = 0 
            syn.ptau = 20
            syn.dtau = 70
            syn.start = 0
            syn.wMax = 1

            exstim = h.NetStim()
            exstim.interval = bg_interval
            exstim.number = 1e9
            exstim.start = 0
            exstim.noise = 1

            syn_g=h.Vector().record(syn._ref_g)
            syn_i=h.Vector().record(syn._ref_i)
            syn_w=h.Vector().record(syn._ref_wtrack)

            h.setpointer(bcm._ref_d, "d", syn)
            h.setpointer(bcm._ref_p, "p", syn)
            weight = wee
            nc = h.NetCon(exstim,syn,0,0,weight)
            post_spikes = h.Vector()
            nc.record(post_spikes)
            synapses.append(syn)
            exstims.append(exstim)
            gli.append(syn_g)
            ili.append(syn_i)
            wli.append(syn_w)
            netcons.append(nc)
            psli.append(post_spikes)

        for si in range(10):
            syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(cell.branches[1](locs[si]))
            ## synapse parameters taken from dg model -- might need reconfig
            syn.tau1 = 0.2
            syn.tau2 = 2.5
            syn.e = 0 
            syn.ptau = 20
            syn.dtau = 70
            syn.start = 0
            syn.wMax = 1

            exstim = h.NetStim()
            exstim.interval = bg_interval
            exstim.number = 1e9
            exstim.start = 0
            exstim.noise = 1

            syn_g=h.Vector().record(syn._ref_g)
            syn_i=h.Vector().record(syn._ref_i)
            syn_w=h.Vector().record(syn._ref_wtrack)

            h.setpointer(bcm._ref_d, "d", syn)
            h.setpointer(bcm._ref_p, "p", syn)
            weight = wee
            nc = h.NetCon(exstim,syn,0,0,weight)
            post_spikes = h.Vector()
            nc.record(post_spikes)
            synapses.append(syn)
            exstims.append(exstim)
            gli.append(syn_g)
            ili.append(syn_i)
            wli.append(syn_w)
            netcons.append(nc)
            psli.append(post_spikes)
        
        
        for si in range(10):
            syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(cell.basal_branches[0](locs[si]))
            ## synapse parameters taken from dg model -- might need reconfig
            syn.tau1 = 0.2
            syn.tau2 = 2.5
            syn.e = 0 
            syn.ptau = 20
            syn.dtau = 70
            syn.start = 0
            syn.wMax = 1

            exstim = h.NetStim()
            exstim.interval = bg_interval
            exstim.number = 1e9
            exstim.start = 0
            exstim.noise = 1

            syn_g=h.Vector().record(syn._ref_g)
            syn_i=h.Vector().record(syn._ref_i)
            syn_w=h.Vector().record(syn._ref_wtrack)

            h.setpointer(bcm._ref_d, "d", syn)
            h.setpointer(bcm._ref_p, "p", syn)
            weight = wee
            nc = h.NetCon(exstim,syn,0,0,weight)
            post_spikes = h.Vector()
            nc.record(post_spikes)
            synapses.append(syn)
            exstims.append(exstim)
            gli.append(syn_g)
            ili.append(syn_i)
            wli.append(syn_w)
            netcons.append(nc)
            psli.append(post_spikes)
        

        for si in range(10):
            syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(cell.basal_branches[1](locs[si]))
            ## synapse parameters taken from dg model -- might need reconfig
            syn.tau1 = 0.2
            syn.tau2 = 2.5
            syn.e = 0 
            syn.ptau = 20
            syn.dtau = 70
            syn.start = 0
            syn.wMax = 1

            exstim = h.NetStim()
            exstim.interval = bg_interval
            exstim.number = 1e9
            exstim.start = 0
            exstim.noise = 1

            syn_g=h.Vector().record(syn._ref_g)
            syn_i=h.Vector().record(syn._ref_i)
            syn_w=h.Vector().record(syn._ref_wtrack)

            h.setpointer(bcm._ref_d, "d", syn)
            h.setpointer(bcm._ref_p, "p", syn)
            weight = wee
            nc = h.NetCon(exstim,syn,0,0,weight)
            post_spikes = h.Vector()
            nc.record(post_spikes)
            synapses.append(syn)
            exstims.append(exstim)
            gli.append(syn_g)
            ili.append(syn_i)
            wli.append(syn_w)
            netcons.append(nc)
            psli.append(post_spikes)
    
        nsyn = len(synapses)

        exstim.seed(trial_n) ## This sets the seed for all background inputs (exstims 3~42)
        # print(nsyn, 'synapses made')
    
    
    sim.sim_time = total_time
    
    # recording
    trec = h.Vector()
    trec.record(h._ref_t)
    vrec = h.Vector()
    vrec.record(cell.soma(0.5)._ref_v)
    aisrec = h.Vector()
    aisrec.record(cell.ais(0.5)._ref_v)
    vdrec = h.Vector()
    vdrec.record(cell.branches[0](0.5)._ref_v)
    vbrec = h.Vector()
    vbrec.record(cell.basal_main(0.5)._ref_v)
    vorec = h.Vector()
    vorec.record(cell.oblique_branch(0.9)._ref_v)
    #BCM record
    bcm_d = h.Vector().record(bcm._ref_d)
    bcm_p = h.Vector().record(bcm._ref_p)
    bcm_alpha = h.Vector().record(bcm._ref_alpha)
    bcm_alpha_scount = h.Vector().record(bcm._ref_alpha_scount)
    if (scen == 2) or (scen == 4) or (scen == 11):
        currentrec = h.Vector()
        currentrec.record(syn._ref_i)
    
    spike_detector = h.NetCon(cell.soma(0.5)._ref_v, None, sec=cell.soma)
    post_spikes = h.Vector()
    spike_detector.record(post_spikes)
    
    # record state vars
    vit2m, vit2h, vscam, vscah, vkcan, vna3dendm, vna3dendh = [h.Vector() for x in range(7)]
    vit2m.record(cell.branches[0](0.5).it2._ref_m)
    vit2h.record(cell.branches[0](0.5).it2._ref_h)
    vscam.record(cell.branches[0](0.5).sca._ref_m)
    vscah.record(cell.branches[0](0.5).sca._ref_h)
    vkcan.record(cell.branches[0](0.5).kca._ref_n)
    vna3dendm.record(cell.branches[0](0.5).na3dend._ref_m)
    vna3dendh.record(cell.branches[0](0.5).na3dend._ref_h)
    
    # run simulation
    # print("Running simulation...")
    sim.go()
    # print("Done. Writing data...")
    
    t = np.array(trec)
    v = np.array(vrec)
    ais = np.array(aisrec)
    vd = np.array(vdrec)
    vb = np.array(vbrec)
    vo = np.array(vorec)
    
    it2m = np.array(vit2m)
    it2h = np.array(vit2h)
    scam = np.array(vscam)
    scah = np.array(vscah)
    kcan = np.array(vkcan)
    na3dendm = np.array(vna3dendm)
    na3dendh = np.array(vna3dendh)
    
    
    # delete
    del(cell)
    del(sim)
    
    del(trec);del(vrec)
    
    my_rawdata['t'] = t
    my_rawdata['v'] = v
    my_rawdata['vd'] = vd
    my_rawdata['vb'] = vb
    my_rawdata['vo'] = vo
    my_rawdata['it2m'] = it2m
    my_rawdata['it2h'] = it2h
    my_rawdata['scam'] = scam
    my_rawdata['scah'] = scah
    my_rawdata['kcan'] = kcan
    my_rawdata['na3dendm'] = na3dendm
    my_rawdata['na3dendh'] = na3dendh
    
    
    rawdata = {'raw_data': my_rawdata}
    
    Ndat = len(t)
    
    ### Data write
    
    # neuron variables
    g=open('%s/vars.txt'%savepath,'w')
    
    data_arrays = [v,ais,vd,vb,it2m,it2h,scam,scah,na3dendm,na3dendh, bcm_alpha_scount, bcm_p, bcm_d]
    
    for i in range(Ndat):
        tt = t[i]
        lin = "%s"%tt
        for datar in data_arrays:
            lin += "  %s"%datar[i]
        lin += '\n'
        g.write(lin)
    g.close()
    
    # synapse weights
    # g=open('%s/weights.txt'%savepath,'w')
    
    # for i in range(Ndat):
    #     tt = t[i]
    #     lin = "%s"%tt
    #     for w in wli:
    #         lin += "  %s"%w[i]
    #     lin += '\n'
    #     g.write(lin)
    
    # g.close()
    
    #presynaptic spike times
    for i in range(nsyn):
        posts=psli[i]
        g = open('%s/pre_%s.txt'%(savepath,i+1),'w')
        for ps in posts:
            g.write('%s\n'%ps)
    
        g.close()
    
    # x1 x2
    g=open("%s/x_state.txt"%savepath,'w')
    
    for i in range(Ndat):
        try: g.write("%s  %s  %s\n"%(t[i],x1_state[i],x2_state[i]))
        except IndexError:
            break
    
    g.close()

def parallelinit(num_P_c, num_r_c):
    global num_P_glo, num_r_glo
    num_P_glo = num_P_c
    num_r_glo = num_r_c

if __name__ == '__main__': 
    
    t1 = time()
    
    mp.freeze_support()

    num_P = 21 # Number of parameters
    num_r = 11 # Number of data points
    
    num_P_glo = mp.Value(ctypes.c_int, num_P)
    num_r_glo = mp.Value(ctypes.c_int, num_r)

    paramlist = list(itertools.product(range(num_P), range(num_r)))
    
    pool = mp.Pool(20, initializer=parallelinit, initargs=(num_P_glo, num_r_glo)) # Number of threads
    
    results = pool.map(pexec, paramlist)
    pool.close()
    
    t2 = time()
    
    print(t2-t1)
