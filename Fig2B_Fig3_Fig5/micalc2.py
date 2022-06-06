from glob import glob
from math import exp, log, log2
from htools import *
import numpy as np

def dLdt(L,I,r_on,r_off,theta):
    return r_on *(1.0+exp(-L)) -r_off*(1.0+exp(L)) + I - theta


def get_post_spikes(t,v, thresh):
    N=len(v)
    post_spikes = []
    last_v = v[0]
    for i in range(N):
        current_v = v[i]
        if current_v > thresh and last_v <= thresh:
            post_spikes.append(t[i])
        last_v = current_v
    return post_spikes

def get_mi_spikes(post_spikes,x_state,t,dt,r_on = 16.7/1000.0,r_off = 2.0*16.7/1000.0):
    ## Calculate  MI between given state x and post spikes in v
    if len(t)!= len(x_state):
        raise Exception()

    N=len(t)

    post_spikes = np.array(post_spikes)
    if t[0]!=0:
        post_spikes = post_spikes - t[0]
        t = t - t[0]

    ## get ests
    time_on = 0.0
    time_off = 0.0
    for i in range(N):
        if x_state[i]==0:
            time_off += dt
        elif x_state[i] ==1:
            time_on += dt
        else: raise Exception()
#     print(time_on, time_off, time_on + time_off)

    est_qon = 0.0
    est_qoff = 0.0

    for spike_time in post_spikes:
        ti = int(spike_time / dt)
#        print('post spike:',spike_time,ti)
        state_at_t = x_state[ti]
        if state_at_t ==1: 
            est_qon += 1.0
        elif state_at_t == 0:
            est_qoff += 1.0
        else: raise Exception()
    if est_qon ==0:
        est_qon +=1
    if est_qoff ==0:
        est_qoff +=1

    est_qon = est_qon/time_on
    est_qoff = est_qoff/time_off

#    print(est_qon, est_qoff)
    est_w = log(est_qon/est_qoff)
    theta = est_qon - est_qoff


    #WARNING: L is integrated with a simple euler method
    tmp_spikes = []
    for tspik in post_spikes:
        tmp_spikes.append(tspik)
    ttt=0.0
    ti = 0
    L = 0
    Lt = []
    Lt.append(L)

    while ti < N:
        dL = dLdt(L,0,r_on,r_off,theta) * dt
        L = L + dL
        try:
            if ttt>tmp_spikes[0]:
                #print('spike detect: %s, %s'%(ttt,tmp_spikes[0]))
                tmp_spikes.pop(0)
                L += est_w
        #         print(vvv,dLdt(L,vvv,r_on,r_off,0),dL)
        except IndexError:
            pass

        ttt += dt
        ti += 1
        Lt.append(L)

    Lt = np.array(Lt)

    p_1 = 1.0/(1.0 + np.exp(-Lt))
    p_0 = 1.0-p_1

    ave_x = 0.0
    ave_hxy = 0.0
    
    for i in range(N):
        x = float(x_state[i])
    
        L = Lt[i]
        p1 = p_1[i]
        p0 = p_0[i]
    
        ave_x += float(x)
        ave_hxy += x*log2(p1) + (1.0-x)*log2(1.0-p1)
    
    ave_x = ave_x / float(N)
    ave_hxy = -ave_hxy /float(N)
    hxx = - ave_x*log2(ave_x) - (1.0 - ave_x)*log2(1.0-ave_x)
    
    MI = hxx - ave_hxy
    
    #     print("Hxx1 = ",hxx1)
    #     print("Hx1y = ",ave_hx1y)
    #     print("Hxx2 = ",hxx2)
    #     print("Hx2y = ",ave_hx2y)
    #     print("MI1 = %s"%MI1)
    #     print("MI2 = %s"%MI2)
    #     print("F1 = %s"%(MI1/preMI1))
    #     print("F2 = %s"%(MI2/preMI2))
    #     print("pre12MI = %s"%pre12MI)
    
    return MI

def get_mi_slice(v,x_state,dt,r_on = 16.7/1000.0,r_off = 2.0*16.7/1000.0):
    ## Calculate  MI between given state x and post spikes in v
    if len(v)!= len(x_state):
        raise Exception()

    N=len(v)
    t = np.array([dt*i for i in range(N)])

    #get post spikes
    post_spikes = []
    last_v = v[0]
    for i in range(N):
        current_v = v[i]
        if current_v > 0.0 and last_v <=0.0:
            post_spikes.append(t[i])
        last_v = current_v

    ## get ests
    time_on = 0.0
    time_off = 0.0
    for i in range(N):
        if x_state[i]==0:
            time_off += dt
        elif x_state[i] ==1:
            time_on += dt
        else: raise Exception()
#     print(time_on, time_off, time_on + time_off)

    est_qon = 0.0
    est_qoff = 0.0

    for spike_time in post_spikes:
        ti = int(spike_time / dt)
#        print('post spike:',spike_time,ti)
        state_at_t = x_state[ti]
        if state_at_t ==1: 
            est_qon += 1.0
        elif state_at_t == 0:
            est_qoff += 1.0
        else: raise Exception()

    est_qon = est_qon/time_on
    est_qoff = est_qoff/time_off

    if est_qon == 0: return 0
    est_w = log(est_qon/est_qoff)
    theta = est_qon - est_qoff


    #WARNING: L is integrated with a simple euler method
    tmp_spikes = []
    for tspik in post_spikes:
        tmp_spikes.append(tspik)
    ttt=0.0
    ti = 0
    L = 0
    Lt = []
    Lt.append(L)

    while ti < N:
        vvv = v[ti]
        dL = dLdt(L,0,r_on,r_off,theta) * dt
        L = L + dL
        try:
            if ttt>tmp_spikes[0]:
                #print('spike detect: %s, %s'%(ttt,tmp_spikes[0]))
                tmp_spikes.pop(0)
                L += est_w
        #         print(vvv,dLdt(L,vvv,r_on,r_off,0),dL)
        except IndexError:
            pass

        ttt += dt
        ti += 1
        Lt.append(L)

    Lt = np.array(Lt)

    p_1 = 1.0/(1.0 + np.exp(-Lt))
    p_0 = 1.0-p_1

    ave_x = 0.0
    ave_hxy = 0.0
    
    for i in range(N):
        x = float(x_state[i])
    
        L = Lt[i]
        p1 = p_1[i]
        p0 = p_0[i]
    
        ave_x += float(x)
        ave_hxy += x*log2(p1) + (1.0-x)*log2(1.0-p1)
    
    ave_x = ave_x / float(N)
    ave_hxy = -ave_hxy /float(N)
    hxx = - ave_x*log2(ave_x) - (1.0 - ave_x)*log2(1.0-ave_x)
    
    MI = hxx - ave_hxy
    
    #     print("Hxx1 = ",hxx1)
    #     print("Hx1y = ",ave_hx1y)
    #     print("Hxx2 = ",hxx2)
    #     print("Hx2y = ",ave_hx2y)
    #     print("MI1 = %s"%MI1)
    #     print("MI2 = %s"%MI2)
    #     print("F1 = %s"%(MI1/preMI1))
    #     print("F2 = %s"%(MI2/preMI2))
    #     print("pre12MI = %s"%pre12MI)
    
    return max(0, MI)



#g=open('data2/mi.txt','w')
#data_all = []
#
#for nsyn in [2,12,22]:
#    data = []
#    for alphaidx in range(37):
#        pres = []
#        mihi = []
#        milo = []
#        for trial in range(10):
#            if len(glob('data2/%s_%s_%s/x_state.txt'%(nsyn,alphaidx,trial)))>0:
#                pre12mi, premi1, premi2, mi1, mi2 = get_mi('data2/%s_%s_%s/'%(nsyn,alphaidx,trial), nsyn)
#                pres.append(pre12mi)
#                mihi.append(max(mi1,mi2))
#                milo.append(min(mi1,mi2))
#                print(nsyn,alphaidx,trial,pre12mi,premi1,premi2,mi1,mi2)
#                g.write("%s  %s  %s  %s  %s  %s  %s  %s\n"%(nsyn,alphaidx,trial,pre12mi,premi1,premi1,mi1,mi2))
#
#        data.append([pres,mihi,milo])
#    data_all.append(data)
#
#g.close()
#                        
#g=open('data2b/mi.txt','w')
#datab_all = []
#
#for nsyn in [2,12,22]:
#    data = []
#    for alphaidx in range(37):
#        pres = []
#        mihi = []
#        milo = []
#        for trial in range(10):
#            if len(glob('data2b/%s_%s_%s/x_state.txt'%(nsyn,alphaidx,trial)))>0:
#                pre12mi, premi1, premi2, mi1, mi2 = get_mi('data2b/%s_%s_%s/'%(nsyn,alphaidx,trial), nsyn)
#                pres.append(pre12mi)
#                mihi.append(max(mi1,mi2))
#                milo.append(min(mi1,mi2))
#                print(nsyn,alphaidx,trial,pre12mi,premi1,premi2,mi1,mi2)
#                g.write("%s  %s  %s  %s  %s  %s  %s  %s\n"%(nsyn,alphaidx,trial,pre12mi,premi1,premi1,mi1,mi2))
#        data.append([pres,mihi,milo])
#    datab_all.append(data)
#
#g.close()
