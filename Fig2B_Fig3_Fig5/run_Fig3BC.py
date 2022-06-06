#!/usr/bin/env python

# --Imports--
import sys
import os
import math
from time import time
#import matplotlib.pyplot as plt
import numpy as np
import neuron
from neuron import h

from neuronmodel import *
from config_model_stim import *
from sim import *
from htools import *
import datetime

sst = datetime.datetime.now()
print('Figure:', sys.argv)
print("Starting time %s/%s/%s  %s:%s:%s"%(sst.day,sst.month,sst.year,sst.hour,sst.minute,sst.second))


### arguments
try:
    nsyn1 = int(sys.argv[1])
    nsyn2 = int(sys.argv[2])
    alpha = float(sys.argv[3])
    trial = int(sys.argv[4])
except IndexError:
    print('Enter 4 arguments: n_1, n_2, alpha, seed number.')
    sys.exit()

### set random seed
random.seed(trial)

### Output directory
path = 'data_Fig3/'
identifier = '%s_%s_%s_%s'%(nsyn1,nsyn2,alpha,trial)
savepath = '%s%s'%(path,identifier)
if not os.path.exists(path):
        os.mkdir(path)
if not os.path.exists(savepath):
	os.mkdir(savepath)
if os.path.isfile('%s/mi.txt'%savepath):
    print('mi.txt found, terminating.')
    sys.exit()


###  Synaptic mechanism  ###
#BCM parameters
P0_BCMTHRESHOLD=0.003
D0_BCMTHRESHOLD=0.001
SCOUNT0_BCMTHRESHOLD = 1/62.5
SCOUNTTAU_BCMTHRESHOLD = 1500   # (ms) averaging time constant for postsynaptic spike count
print("Time constant for the weighted spike count = %e min\n"%(SCOUNTTAU_BCMTHRESHOLD/1000/60))
ALPHA_BCMTHRESHOLD = 62.5*0.5               # scaling constant, ALPHA=0 for non-BCM model MODIFIED NSYN3
BCMSLOW_BCMTHRESHOLD = 0                # Cliff's slowing of BCM rate, 0=no change, 0.5=50% slowing, 0.2=20% slowing etc
print("ALPHA_BCMTHRESHOLD = %e\n"%ALPHA_BCMTHRESHOLD)
PTAU = 20
DTAU = 70
gamma = 1.0

DT=0.2 # ms, set the integration step

my_rawdata = {}

#inputs
freq = params['Input']['freq'].value
wee = params['Input']['wee'].value
wee_strong = params['Input']['wee_strong'].value

#background frequency in Hz
bgf = 1.0


# Other
delta_t = params['STDP']['delta_t'].value
sim_params = params['sim']
source = False

cell = Neuron()
sim = Simulation(cell,sim_params)
sim.dt = DT
sim.v_init = -70

### Total simulation time, in ms
total_time = 1000*40

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
q1 = 100.0
q2 = 100.0
r_on_slow = 16.7/1000.0
r_on = r_on_slow
r_off = r_on*2.0
q_on1 = q1/1000.0
q_off1 = 0.001
q_on2 = q2/1000.0
q_off2 = 0.001

### generate coupled states x1 and x2
print("creating x states...")
x1_state = np.zeros(nsteps)
x2_state = np.zeros(nsteps)
generate_coupled_states(x1_state,x2_state,r_on,r_off,alpha,DT)
print("done.")

### pre-synaptic spikes
computed_spikes=[]
tvecli=[]
vvli=[]

branch1segs = []
for seg in cell.branches[0]:
    branch1segs.append(seg)

branch2segs = []
for seg in cell.branches[1]:
    branch2segs.append(seg)

allapicalsegs = branch1segs + branch2segs
for seg in cell.apical:
    allapicalsegs.append(seg)
for seg in cell.apical_prox:
    allapicalsegs.append(seg)

allbasalsegs = []
for seg in cell.basal_main:
    allbasalsegs.append(seg)
for seg in cell.basal_branches[0]:
    allbasalsegs.append(seg)
for seg in cell.basal_branches[1]:
    allbasalsegs.append(seg)

targetsegs = []
for seg in cell.apical:
    targetsegs.append(seg)

dsec_list = []
dsec_list.append(cell.branches[0](1.0))
dsec_list.append(cell.branches[0](1.0/3.0))
dsec_list.append(cell.apical(0.025*30))
dsec_list.append(cell.apical(0.025*10))
dsec_list.append(cell.apical_prox(0.0))
dsec_list.append(cell.basal_main(0))
dsec_list.append(cell.basal_branches[0](1.0/3.0))

if True: #Initialize synapses and hidden state inputs

    bcm = h.BCMthreshold3(cell.soma(0.5))
    bcm.bcm_scale = gamma
    bcm.p0=P0_BCMTHRESHOLD
    bcm.d0=D0_BCMTHRESHOLD
    bcm.scount0 = SCOUNT0_BCMTHRESHOLD
    bcm.scounttau = SCOUNTTAU_BCMTHRESHOLD
    bcm.alpha = ALPHA_BCMTHRESHOLD

    #Inputs from Hidden State x1
    for si in range(nsyn1):
        syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(targetsegs[60+2*si])

        ## synapse parameters taken from dg model -- might need reconfig
        syn.tau1 = 0.2
        syn.tau2 = 2.5
        syn.e = 0 
        syn.ptau = 20
        syn.dtau = 70
        syn.start = 0
        syn.wMax = 1

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
        #syn_i=h.Vector().record(syn._ref_i)
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
        #ili.append(syn_i)
        wli.append(syn_w)
        netcons.append(nc)
        psli.append(post_spikes)

    # Inputs from hidden state x2
    for si in range(nsyn2):
        syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(targetsegs[61+2*si])

        ## synapse parameters taken from dg model -- might need reconfig
        syn.tau1 = 0.2
        syn.tau2 = 2.5
        syn.e = 0 
        syn.ptau = 20
        syn.dtau = 70
        syn.start = 0
        syn.wMax = 1

        vv1 = h.Vector(np.array([-1.0 for j in range(nsteps)]))
        tvec1 = h.Vector(np.arange(0,total_time+DT,DT))
        event_times = compute_spike_train(x2_state,q_on2,q_off2,DT)
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
        #syn_i=h.Vector().record(syn._ref_i)
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
        #ili.append(syn_i)
        wli.append(syn_w)
        netcons.append(nc)
        psli.append(post_spikes)

    # Background synapses in apical dendrite
    for si in range(71):
        randi = random.randrange(len(allapicalsegs))
        randomseg = allapicalsegs.pop(randi)
        syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(randomseg)
        ## synapse parameters taken from dg model -- might need reconfig
        syn.tau1 = 0.2
        syn.tau2 = 2.5
        syn.e = 0 
        syn.ptau = 20
        syn.dtau = 70
        syn.start = 0
        syn.wMax = 1

        exstim = h.NetStim()
        exstim.interval = 1000.0/bgf
        exstim.number = 1e9
        exstim.start = 0
        exstim.noise = 1

        #syn_g=h.Vector().record(syn._ref_g)
        #syn_i=h.Vector().record(syn._ref_i)
        syn_w=h.Vector().record(syn._ref_wtrack)

        h.setpointer(bcm._ref_d, "d", syn)
        h.setpointer(bcm._ref_p, "p", syn)
        weight = wee
        nc = h.NetCon(exstim,syn,0,0,weight)
        post_spikes = h.Vector()
        nc.record(post_spikes)
        synapses.append(syn)
        exstims.append(exstim)
        #gli.append(syn_g)
        #ili.append(syn_i)
        wli.append(syn_w)
        netcons.append(nc)
        psli.append(post_spikes)

    # Background synapses in basal dendrite
    for si in range(29):
        randi = random.randrange(len(allbasalsegs))
        randomseg = allbasalsegs.pop(randi)
        syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(randomseg)
        ## synapse parameters taken from dg model -- might need reconfig
        syn.tau1 = 0.2
        syn.tau2 = 2.5
        syn.e = 0 
        syn.ptau = 20
        syn.dtau = 70
        syn.start = 0
        syn.wMax = 1

        exstim = h.NetStim()
        exstim.interval = 1000.0/bgf
        exstim.number = 1e9
        exstim.start = 0
        exstim.noise = 1

        #syn_g=h.Vector().record(syn._ref_g)
        #syn_i=h.Vector().record(syn._ref_i)
        syn_w=h.Vector().record(syn._ref_wtrack)

        h.setpointer(bcm._ref_d, "d", syn)
        h.setpointer(bcm._ref_p, "p", syn)
        weight = wee
        nc = h.NetCon(exstim,syn,0,0,weight)
        post_spikes = h.Vector()
        nc.record(post_spikes)
        synapses.append(syn)
        exstims.append(exstim)
        #gli.append(syn_g)
        #ili.append(syn_i)
        wli.append(syn_w)
        netcons.append(nc)
        psli.append(post_spikes)

    nsyn = len(synapses)
    print(nsyn, 'synapses made')


sim.sim_time = total_time

# recording
trec = h.Vector()
trec.record(h._ref_t)
vrec = h.Vector()
vrec.record(cell.soma(0.5)._ref_v)
#aisrec = h.Vector()
#aisrec.record(cell.ais(0.5)._ref_v)
#vdrec = h.Vector()
#vdrec.record(cell.branches[0](0.5)._ref_v)
#vbrec = h.Vector()
#vbrec.record(cell.basal_main(0.5)._ref_v)
#vorec = h.Vector()
#vorec.record(cell.oblique_branch(0.9)._ref_v)
#BCM record
asgrec = h.Vector()
asgrec.record(bcm._ref_asg)
asrec = h.Vector()
asrec.record(bcm._ref_alpha_scount)
prec = h.Vector()
prec.record(bcm._ref_p)
drec = h.Vector()
drec.record(bcm._ref_d)

spike_detector = h.NetCon(cell.soma(0.5)._ref_v, None, sec=cell.soma)
post_spikes = h.Vector()
spike_detector.record(post_spikes)

# run simulation
print("Running simulation...")
sim.go()
print("Done. Writing data...")

t = np.array(trec)
v = np.array(vrec)
pp = np.array(prec)
dd = np.array(drec)

asg = np.array(asgrec)
alpha_scount = np.array(asrec)

# delete
del(cell)
del(sim)

del(trec);del(vrec)

my_rawdata['t'] = t
my_rawdata['v'] = v

rawdata = {'raw_data': my_rawdata}

Ndat = len(t)

### Data write

# neuron variables

only_mi = False 
if not only_mi:
    g=open('%s/vars.txt'%savepath,'w')
    
    #data_arrays = [v,dv1,dv2,dv3,ais,vd,vb,it2m,it2h,scam,scah,na3dendm,na3dendh]
    data_arrays = [v]
    
    for i in range(Ndat):
        if i%5 !=0: continue
        tt = t[i]
        lin = "%s"%tt
        for datar in data_arrays:
            lin += "  %s"%datar[i]
        lin += '\n'
        g.write(lin)
    g.close()
    
    # synapse weights
    g=open('%s/weights.txt'%savepath,'w')
    
    for i in range(Ndat):
        if i%5 !=0: continue
        tt = t[i]
        lin = "%s"%tt
        for w in wli:
            lin += "  %s"%w[i]
        lin += '\n'
        g.write(lin)
    
    g.close()
    
    # synapse conductances (just one)
    g=open('%s/syng.txt'%savepath,'w')
    
    for i in range(Ndat):
        if i%5 !=0: continue
        tt = t[i]
        lin = "%s"%tt
        for sg in gli:
            lin += "  %s"%sg[i]
        lin += '\n'
        g.write(lin)
    
    g.close()
    
    #presynaptic spike times
    for i in range(2):
        posts=psli[i]
        g = open('%s/pre_%s.txt'%(savepath,i+1),'w')
        for ps in posts:
            g.write('%s\n'%ps)
    
        g.close()
    
    #dendrite membrane potential
    #for i in range(len(dend_list)):
    #    g=open("%s/%s.tst"%(savepath,i+1), 'w')
    #    dvec = dend_list[i]
    #    for ijk in range(Ndat):
    #        g.write("%s  %s\n"%(t[ijk], dvec[ijk]))
    #    g.close()
    
    # x1 x2
    g=open("%s/x_state.txt"%savepath,'w')
    
    print(Ndat, len(x1_state))
    for i in range(Ndat):
        try: g.write("%s  %s  %s\n"%(t[i],x1_state[i],x2_state[i]))
        except IndexError:
            break
    
    g.close()
    
    g=open("%s/asg.txt"%savepath,'w')
    for i in range(Ndat):
        try: g.write("%s  %s  %s  %s  %s\n"%(t[i],alpha_scount[i],asg[i],pp[i],dd[i]))
        except IndexError:
            break
    g.close()

import micalc2
#total_time = 1000*60*5

g=open("%s/mi.txt"%savepath,'w')
for ii in range(3):
    starti = int(10000*ii/0.2)
    endi = int((10000+10000*ii)/0.2)
    MI1 = micalc2.get_mi_slice(v[starti:endi], x1_state[starti:endi], 0.2)
    MI2 = micalc2.get_mi_slice(v[starti:endi], x2_state[starti:endi], 0.2)
    g.write("%s, %s\n"%(MI1,MI2))
g.close()


sst = datetime.datetime.now()
print("Stop at %s/%s/%s  %s:%s:%s"%(sst.day,sst.month,sst.year,sst.hour,sst.minute,sst.second))
