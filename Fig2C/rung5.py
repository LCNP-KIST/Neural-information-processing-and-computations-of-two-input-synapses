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

#from neuronmodel import *
from config_model_stim import *
from sim import *
from htools import *
from datetime import datetime

###  Synaptic mechanism  ###
#BCM parameters
P0_BCMTHRESHOLD=0.002
D0_BCMTHRESHOLD=0.001
SCOUNT0_BCMTHRESHOLD = 1/62.5
SCOUNTTAU_BCMTHRESHOLD = 1500   # (ms) averaging time constant for postsynaptic spike count
print("Time constant for the weighted spike count = %e min\n"%(SCOUNTTAU_BCMTHRESHOLD/1000/60))
ALPHA_BCMTHRESHOLD = 62.5               # scaling constant, ALPHA=0 for non-BCM model
BCMSLOW_BCMTHRESHOLD = 0                # Cliff's slowing of BCM rate, 0=no change, 0.5=50% slowing, 0.2=20% slowing etc
print("ALPHA_BCMTHRESHOLD = %e\n"%ALPHA_BCMTHRESHOLD)
PTAU = 20
DTAU = 70


path = 'data/'

SCEN = 11
NO_REPS = 1
RESET = True
DT=0.2 # ms, set the integration step
POST_AMP = 0.3# nA, amplitude of current injection to trigger the AP/bAP
WARM_UP=1000 # ms
AP_DELAY = 0 #ms, AP needs 4ms after stimulation to become initiated

print(sys.argv)

try:
    nsyn1 = int(sys.argv[1])
    qall =  float(sys.argv[2])
    dcaapw = float(sys.argv[3])    
    stimw = float(sys.argv[4])
    trial = int(sys.argv[5])
except IndexError:
    print('needs arguments')
    sys.exit()

bgw = 1.0
alpha = 1.0

random.seed(trial)

identifier = '%s_%s_%s_%s_%s_%s'%(nsyn1,qall,dcaapw,stimw,bgw,trial)
savepath = '%s%s'%(path,identifier)
if os.path.exists("%s/vars.txt"%savepath):
        print("VARS found. terminating.")
        sys.exit()
if not os.path.exists(path):
        os.mkdir(path)
if not os.path.exists(savepath):
	os.mkdir(savepath)

def _get_current_trace(freq,delta_t,t_stop,pre=False,test=True) :
    trace = np.zeros(int(t_stop/DT))
    for i in range(NO_REPS) :
        if(pre) :
            start_t = (0 + i* (1000.0/freq) + WARM_UP)
        else :
            start_t = (0 + delta_t - AP_DELAY + i* (1000.0/freq) + WARM_UP)
        end_t = (start_t+2)
        if(test) :
            print('start_t=%g, end_t=%g (t_stop=%g, len(trace)=%f)' % (start_t,end_t,t_stop,len(trace)))
        trace[int(start_t/DT):int(end_t/DT)] = POST_AMP
    return trace

def create_syn(syn_type,pos, thresh):
    syn = h.ExpSynSTDP(pos)
    syn.thresh = thresh
    return syn


my_rawdata = {}

#inputs
freq = params['Input']['freq'].value
wee_modulation = 3
wee = params['Input']['wee'].value*wee_modulation
wee_strong = params['Input']['wee_strong'].value*wee_modulation

# Synapses
scen = SCEN
delta_t = params['STDP']['delta_t'].value

sim_params = params['sim']

source = False

#cell = Neuron()

#Define Gidon cell instead - by JH Woo
h.load_file('xor.hoc')
h.biophys()
class Cell():
    pass
cell = Cell()
cell.soma = h.soma
cell.E = -74.0
cell.apical_prox = h.apic[0]
cell.apical = h.apic[28]
cell.oblique_branch = h.apic[1]
cell.branches = [h.apic[35], h.apic[56]]
cell.second_branches = [h.apic[80], h.apic[71], h.apic[55], h.apic[51]]
cell.basal_main = h.dend[86]
cell.basal_branches = [h.dend[87], h.dend[94]]
cell.ais = h.axon

h.apic[76].insert('dCaAP')
dcaap_loc = 0.01
#dcaap_loc = 0.214
dcaap_w = dcaapw
seg_dcaap = h.apic[76](dcaap_loc)
seg_dcaap.dCaAP.vth = -36
seg_dcaap.dCaAP.D = 0.3
seg_dcaap.dCaAP.vrest = -75
seg_dcaap.dCaAP.refract_period = 50
seg_dcaap.dCaAP.sigma_diff = 21
seg_dcaap.dCaAP.tauA = 3
seg_dcaap.dCaAP.tauB = 0.4
seg_dcaap.dCaAP.w = dcaap_w


sim = Simulation(cell,sim_params)
sim.dt = DT
sim.v_init = -70
# total_time = WARM_UP+NO_REPS*(1000.0/freq)+100
total_time = 1000*30

synapses=[]
exstims=[]
gli=[]
ili=[]
wli=[]
netcons=[]
psli=[]

nsteps = int(total_time/DT)+1

x1_state = np.zeros(nsteps)
x2_state = np.zeros(nsteps)
r_on_slow = 16.7/1000.0
r_on = r_on_slow
r_off = r_on*2.0

print("creating x states...")
generate_coupled_states(x1_state,x2_state,r_on,r_off,alpha,DT)
print("done.")

q_on1 = 0.05
q_off1 = 0.001

q_on1 = qall / 1000.0
q_on2 = qall / 1000.0

computed_spikes=[]
tvecli=[]
vvli=[]

targetsegs = []
for i in [72,76]:
    #print('adding segs from apic',i)
    for seg in h.apic[i]:
        targetsegs.append(seg)

print(len(targetsegs),'targetsegs')


all_dend = []
for sec in h.dend:
    for seg in sec:
        all_dend.append(seg)
for sec in h.apic:
    for seg in sec:
        all_dend.append(seg)

all_dend.remove(seg_dcaap)


recdends = []
recdendsecs = []

#from soma to dcaap location
counter = 0
for di in [0, 6, 28, 34, 56, 60]:
    for seg in h.apic[di]:
        recdends.append(seg)
        recdendsecs.append(h.apic[di])

        #print('recording apic:',seg)

print(len(recdends),'recdends')

if (scen == 11):

    for si in range(nsyn1):
        syn = h.Exp2Syn(targetsegs[si])
        all_dend.remove(targetsegs[si])
        ## synapse parameters taken from dg model -- might need reconfig
        syn.tau1 = 0.2
        syn.tau2 = 2.5
        syn.e = 0 

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

        weight = wee*stimw
        nc = h.NetCon(spik1,syn,0,0,weight)
        post_spikes = h.Vector()
        nc.record(post_spikes)
        synapses.append(syn)
        exstims.append(spik1)
        gli.append(syn_g)
        ili.append(syn_i)
        netcons.append(nc)
        psli.append(post_spikes)
        print("X1 placed in ",syn.get_segment())

    for si in range(100):
        randi = random.randrange(len(all_dend))
        randomseg = all_dend.pop(randi)

        syn = h.Exp2Syn(randomseg)
        ## synapse parameters taken from dg model -- might need reconfig
        syn.tau1 = 0.2
        syn.tau2 = 2.5
        syn.e = 0 

        exstim = h.NetStim()
        exstim.interval = 1000
        exstim.number = 1e9
        exstim.start = 0
        exstim.noise = 1


        weight = wee*bgw
        nc = h.NetCon(exstim,syn,0,0,weight)
        post_spikes = h.Vector()
        nc.record(post_spikes)
        synapses.append(syn)
        exstims.append(exstim)
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

drecvecs = []
drecs = []

dend_posts = []
dend_detectors = []
for di, seg in enumerate(recdends):
    drec = h.Vector()
    drec.record(seg._ref_v)
    drecvecs.append(drec)
    drecs.append(drec)

    spike_detector = h.NetCon(seg._ref_v, None,sec = recdendsecs[di])
    post_spikes = h.Vector()
    spike_detector.record(post_spikes)

    dend_posts.append(post_spikes)
    dend_detectors.append(spike_detector)

#BCM record

spike_detector = h.NetCon(cell.soma(0.5)._ref_v, None, sec=cell.soma)
post_spikes = h.Vector()
spike_detector.record(post_spikes)

# record state vars
#vit2m, vit2h, vscam, vscah, vkcan, vna3dendm, vna3dendh = [h.Vector() for x in range(7)]
#vit2m.record(cell.branches[0](0.5).it2._ref_m)
#vit2h.record(cell.branches[0](0.5).it2._ref_h)
#vscam.record(cell.branches[0](0.5).sca._ref_m)
#vscah.record(cell.branches[0](0.5).sca._ref_h)
#vkcan.record(cell.branches[0](0.5).kca._ref_n)
#vna3dendm.record(cell.branches[0](0.5).na3dend._ref_m)
#vna3dendh.record(cell.branches[0](0.5).na3dend._ref_h)

# run simulation
print("Running simulation")
sim.go()
print("Done.")

t = np.array(trec)
v = np.array(vrec)
#ais = np.array(aisrec)
#vd = np.array(vdrec)
#vb = np.array(vbrec)
#vo = np.array(vorec)
#vdcaap = np.array(vdcaaprec)

dvs = []

for ddi, dd in enumerate(drecs):
    print("conv drecs...", ddi,len(drecs))
    dvs.append(np.array(dd))
print("...done.")

#it2m = np.array(vit2m)
#it2h = np.array(vit2h)
#scam = np.array(vscam)
#scah = np.array(vscah)
#kcan = np.array(vkcan)
#na3dendm = np.array(vna3dendm)
#na3dendh = np.array(vna3dendh)


#plt.plot(t,v)
#plt.plot(t,vd)
#plt.show()


# delete
del(cell)
del(sim)

del(trec);del(vrec)

#my_rawdata['t'] = t
#my_rawdata['v'] = v
#my_rawdata['vd'] = vd
#my_rawdata['vb'] = vb
#my_rawdata['vo'] = vo
#my_rawdata['it2m'] = it2m
#my_rawdata['it2h'] = it2h
#my_rawdata['scam'] = scam
#my_rawdata['scah'] = scah
#my_rawdata['kcan'] = kcan
#my_rawdata['na3dendm'] = na3dendm
#my_rawdata['na3dendh'] = na3dendh


rawdata = {'raw_data': my_rawdata}

Ndat = len(t)

print('print vars.')
g=open('%s/vars.txt'%savepath,'w')

#data_arrays = [v,ais,vd,vb,it2m,it2h,scam,scah,na3dendm,na3dendh, bcm_alpha_scount, bcm_p, bcm_d]
#data_arrays = [v,ais,vd,vb,bcm_alpha_scount, bcm_p, bcm_d, vdcaap]
data_arrays = [v] #, vdcaap]

for i in range(Ndat):
    tt = t[i]
    lin = "%s"%tt
    for datar in data_arrays[:3]:
        lin += "  %s"%datar[i]
    lin += '\n'
    g.write(lin)
g.close()

if trial == 0:
    print('print dv.')
    g =open("%s/dv.txt"%savepath,'w')
    for i in range(Ndat):
        tt = t[i]
        lin = "%s"%tt
        for datar in dvs:
            lin += "  %s"%datar[i]
        lin += '\n'
        g.write(lin)
    g.close()

print('print pre.')
for i in range(nsyn):
    posts=psli[i]
    g = open('%s/pre_%s.txt'%(savepath,i+1),'w')
    for ps in posts:
        g.write('%s\n'%ps)

    g.close()

print('print x1.')
g=open("%s/x_state.txt"%savepath,'w')

print(Ndat, len(x1_state))
for i in range(Ndat):
    try: g.write("%s  %s  %s\n"%(t[i],x1_state[i],x2_state[i]))
    except IndexError:
        break

g.close()

g = open('%s/post_soma.txt'%savepath,'w')
for ps in post_spikes:
    g.write('%s\n'%ps)

g.close()

for i, post_spikes in enumerate(dend_posts):
    g = open('%s/post_d%s.txt'%(savepath,i),'w')
    for ps in post_spikes:
        g.write('%s\n'%ps)
    
    g.close()



import micalc2
#total_time = 1000*60*5

g=open("%s/mi.txt"%savepath,'w')
starti = int(0)
endi = int(30000/0.2)
MI1 = micalc2.get_mi_slice(v[starti:endi], x1_state[starti:endi], 0.2)
MI2 = micalc2.get_mi_slice(v[starti:endi], x2_state[starti:endi], 0.2)
g.write("%s, %s\n"%(MI1,MI2))
for ii in range(len(dvs)):
    print('calculating MI', ii)
    starti = int(0)
    endi = int(30000/0.2)
    MI1 = micalc2.get_mi_slice(dvs[ii][starti:endi], x1_state[starti:endi], 0.2)
    MI2 = micalc2.get_mi_slice(dvs[ii][starti:endi], x2_state[starti:endi], 0.2)
    g.write("%s, %s\n"%(MI1,MI2))
g.close()

