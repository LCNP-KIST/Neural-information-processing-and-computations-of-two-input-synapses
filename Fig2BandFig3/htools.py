import random
from math import log2
import numpy as np
from glob import glob
from math import exp, log, log2


def generate_periodic(x1_state,onlen,offlen,dt,offset=0):
    #start with off.
    N = len(x1_state)
    current_ti = 0
    current_x1 = 0
    oncount = 0.0
    offcount = 0.0

    if offset >0:
        if offset < offlen:
            offcount = offset
        elif onlen+offlen > offset >= offlen:
            oncount = offset - offlen
            current_x1 = 1
        else:
            raise Exception("offset should be under onlen + offlen")


    while current_ti < N:
        x1_state[current_ti] = current_x1

        if current_x1 == 0:

            if offcount >= offlen:
                #transition
                offcount = 0
                current_x1 = 1
            else:
                offcount += dt

        elif current_x1 == 1:
            if oncount >= onlen:
                #transition
                oncount = 0
                current_x1 = 0
            else:
                oncount += dt

        current_ti += 1


def generate_ind_states(x1_state,x2_state,r_on,r_off,dt):
    N = len(x1_state)
    current_ti = 0
    current_x1 = 0
    current_x2 = 0
    pon = r_on*dt
    poff = r_off*dt

    ## col1: transition from 0, 0
    p00_00 = ( (1.0-pon)**2 ) 
    p00_01 = p00_00 + (1.0-pon)*pon 
    p00_10 = p00_01 + (1.0-pon)*pon
    p00_11 = p00_10 + pon*pon
    ## col2: transition from 0, 1
    p01_00 = (1.0-pon)*poff
    p01_01 = p01_00 + (1.0 - pon)*(1.0 - poff)
    p01_10 = p01_01 + pon * poff
    p01_11 = p01_10 + pon * (1.0 - poff)
    ## col3: transition from 1, 0
    p10_00 = (1.0-pon)*poff
    p10_01 = p10_00 + pon*poff
    p10_10 = p10_01 + (1.0 - pon) * (1.0 - poff)
    p10_11 = p10_10 + pon * (1.0 - poff) 
    ## col1: transition from 1, 1
    p11_00 = poff**2 
    p11_01 = p11_00 + (1.0-poff)*poff 
    p11_10 = p11_01 + (1.0-poff)*poff
    p11_11 = p11_10 + (1.0 - poff)**2 

#    ## col1: transition from 0, 0
#    print(p00_00)
#    print(p00_01)
#    print(p00_10)
#    print(p00_11)
#    ## col2: transition from 0, 1
#    print(p01_00)
#    print(p01_01)
#    print(p01_10)
#    print(p01_11)
#    ## col3: transition from 1, 0
#    print(p10_00)
#    print(p10_01)
#    print(p10_10)
#    print(p10_11)
#    ## col1: transition from 1, 1
#    print(p11_00)
#    print(p11_01)
#    print(p11_10)
#    print(p11_11)


    while current_ti < N:
        n = random.random()
        if current_x1 == 0 and current_x2 == 0:
            if n < p00_00:
                x1, x2 = 0,0
            elif n < p00_01:
                x1, x2 = 0,1
            elif n < p00_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 0 and current_x2 == 1:
            if n < p01_00:
                x1, x2 = 0,0
            elif n < p01_01:
                x1, x2 = 0,1
            elif n < p01_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 1 and current_x2 == 0:
            if n < p10_00:
                x1, x2 = 0,0
            elif n < p10_01:
                x1, x2 = 0,1
            elif n < p10_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 1 and current_x2 == 1:
            if n < p11_00:
                x1, x2 = 0,0
            elif n < p11_01:
                x1, x2 = 0,1
            elif n < p11_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        else:
            raise Exception()

        x1_state[current_ti] = x1
        x2_state[current_ti] = x2
        current_x1 = x1
        current_x2 = x2
        current_ti += 1

def generate_coupled_states(x1_state,x2_state,r_on,r_off,alpha,dt):
    N = len(x1_state)
    current_ti = 0
    current_x1 = 0
    current_x2 = 0
    pon = r_on*dt
    poff = r_off*dt

    p00_00 = ( (1.0-pon)**2 ) * (1.0 + (1/(1.0-pon)-1.0 )*alpha)
    p00_01 = p00_00+ (1.0-pon)*pon *(1-alpha)
    p00_10 = p00_01+ (1.0-pon)*pon *(1-alpha )
    p00_11 = p00_10+ pon*pon * (1.0 + alpha*(1.0 - pon)/pon)
    ## col2: transition from 0, 1

    p01_00 = (1.0-pon)*poff
    p01_01 = p01_00 + (1.0 - pon)*(1.0 - poff)
    p01_10 = p01_01 + pon * poff
    p01_11 = p01_10 + pon * (1.0 - poff)
    ## col3: transition from 1, 0
    p10_00 = (1.0-pon)*poff 
    p10_01 = p10_00 + pon*poff 
    p10_10 = p10_01 + (1.0 - pon) * (1.0 - poff) 
    p10_11 = p10_10 + pon * (1.0 - poff) 
    ## col1: transition from 1, 1
    p11_00 = poff**2  * (1.0 + alpha*(1.0 - poff)/poff)
    p11_01 = p11_00 + (1.0-poff)*poff  *(1.0-alpha)
    p11_10 = p11_01 + (1.0-poff)*poff *(1.0-alpha)
    p11_11 = p11_10 + (1.0 - poff)**2 * (1.0 + (1.0/(1.0-poff)-1.0)*alpha)

    while current_ti < N:
        n = random.random()
        if current_x1 == 0 and current_x2 == 0:
            if n < p00_00:
                x1, x2 = 0,0
            elif n < p00_01:
                x1, x2 = 0,1
            elif n < p00_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 0 and current_x2 == 1:
            if n < p01_00:
                x1, x2 = 0,0
            elif n < p01_01:
                x1, x2 = 0,1
            elif n < p01_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 1 and current_x2 == 0:
            if n < p10_00:
                x1, x2 = 0,0
            elif n < p10_01:
                x1, x2 = 0,1
            elif n < p10_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 1 and current_x2 == 1:
            if n < p11_00:
                x1, x2 = 0,0
            elif n < p11_01:
                x1, x2 = 0,1
            elif n < p11_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        else:
            raise Exception()

        x1_state[current_ti] = x1
        x2_state[current_ti] = x2
        current_x1 = x1
        current_x2 = x2
        current_ti += 1

def generate_coupled_states_old(x1_state,x2_state,r_on,r_off,alpha,dt):
    N = len(x1_state)
    current_ti = 0
    current_x1 = 0
    current_x2 = 0
    pon = r_on*dt
    poff = r_off*dt

    ## col1: transition from 0, 0
    p00_00 = ( (1.0-pon)**2 ) / (1.0 - alpha*pon)
    p00_01 = p00_00 + (1.0-pon)*pon *( (1.0+pon)/(2.0*pon) - alpha/2.0 - (1.0-pon)/(2.0*pon*(1.0-alpha*pon)) )
    p00_10 = p00_01 + (1.0-pon)*pon *( (1.0+pon)/(2.0*pon) - alpha/2.0 - (1.0-pon)/(2.0*pon*(1.0-alpha*pon))  )
    p00_11 = p00_10 + pon*pon * (1.0 + alpha*(1.0 - pon)/pon)
    ## col2: transition from 0, 1
    p01_00 = (1.0-pon)*poff * (1.0 + alpha*(poff/(pon+poff))/( (1.0-pon)*poff) - alpha)
    p01_01 = p01_00 + (1.0 - pon)*(1.0 - poff) * (1.0-alpha)
    p01_10 = p01_01 + pon * poff * (1.0-alpha)
    p01_11 = p01_10 + pon * (1.0 - poff) * (1.0 + alpha*(pon/(pon+poff))/( (1.0-poff)*pon) - alpha)
    ## col3: transition from 1, 0
    p10_00 = (1.0-pon)*poff * (1.0 + alpha*(poff/(pon+poff))/( (1.0-pon)*poff) - alpha)
    p10_01 = p10_00 + pon*poff * ( 1.0-alpha)
    p10_10 = p10_01 + (1.0 - pon) * (1.0 - poff) * ( 1.0-alpha)
    p10_11 = p10_10 + pon * (1.0 - poff)  *  (1.0 + alpha*(pon/(pon+poff))/( (1.0-poff)*pon) - alpha)
    ## col1: transition from 1, 1
    p11_00 = poff**2  * (1.0 + alpha*(1.0 - poff)/poff)
    p11_01 = p11_00 + (1.0-poff)*poff  *( (1.0+poff)/(2.0*poff) - alpha/2.0 - (1.0-poff)/(2.0*poff*(1.0-alpha*poff)) )
    p11_10 = p11_01 + (1.0-poff)*poff *( (1.0+poff)/(2.0*poff) - alpha/2.0 - (1.0-poff)/(2.0*poff*(1.0-alpha*poff)) )
    p11_11 = p11_10 + (1.0 - poff)**2 / (1.0 - alpha*poff)

#    ## col1: transition from 0, 0
#    print(p00_00)
#    print(p00_01)
#    print(p00_10)
#    print(p00_11)
#    ## col2: transition from 0, 1
#    print(p01_00)
#    print(p01_01)
#    print(p01_10)
#    print(p01_11)
#    ## col3: transition from 1, 0
#    print(p10_00)
#    print(p10_01)
#    print(p10_10)
#    print(p10_11)
#    ## col1: transition from 1, 1
#    print(p11_00)
#    print(p11_01)
#    print(p11_10)
#    print(p11_11)


    while current_ti < N:
        n = random.random()
        if current_x1 == 0 and current_x2 == 0:
            if n < p00_00:
                x1, x2 = 0,0
            elif n < p00_01:
                x1, x2 = 0,1
            elif n < p00_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 0 and current_x2 == 1:
            if n < p01_00:
                x1, x2 = 0,0
            elif n < p01_01:
                x1, x2 = 0,1
            elif n < p01_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 1 and current_x2 == 0:
            if n < p10_00:
                x1, x2 = 0,0
            elif n < p10_01:
                x1, x2 = 0,1
            elif n < p10_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        elif current_x1 == 1 and current_x2 == 1:
            if n < p11_00:
                x1, x2 = 0,0
            elif n < p11_01:
                x1, x2 = 0,1
            elif n < p11_10:
                x1, x2 = 1,0
            else:
                x1, x2 = 1,1
        else:
            raise Exception()

        x1_state[current_ti] = x1
        x2_state[current_ti] = x2
        current_x1 = x1
        current_x2 = x2
        current_ti += 1

def random_qs(mu_q,sigma_q):
    q_on = np.random.normal(mu_q,sigma_q)
    q_off = np.random.normal(mu_q,sigma_q)
    
    if q_on > q_off > 0: 
        return q_on, q_off
    else:
        return random_qs(mu_q,sigma_q)
    
def compute_spike_train(states,q_on,q_off,dt):
    t=0
    event_times = []
    for s in states:
        if s == 1:
            p = q_on*dt
            if random.random() < p:
                event_times.append(t)
        elif s == 0:
            p = q_off*dt
            if random.random() < p:
                event_times.append(t)
        else:
            raise Exception()
        t+=dt
    return event_times


def generate_state(x1_state,r_on,r_off,dt):
    N = len(x1_state)

    current_ti = 0
    current_x1 = 0
    while current_ti < N:
        
        #Get a random probability value from the uniform distribution's PDF
        n = random.random()
    
        if current_x1 == 0:
            if n < r_on*dt:
                x1 = 1
            else:
                x1 = 0
        elif current_x1 == 1:
            if n < r_off*dt:
                x1 = 0
            else:
                x1 = 1

        x1_state[current_ti] = x1
        current_x1 = x1
        current_ti += 1

def generate_dependent_state(x2_state,x1_state,r_on,r_off,dt,nu):
    N = len(x1_state)
    if len(x2_state)!=N: raise Exception()

    current_ti = 0
    current_x2 = 0
    while current_ti < N:
        
        #Get a random probability value from the uniform distribution's PDF
        n = random.random()
    
        x1 = x1_state[current_ti]
        if x1 == 0:
            if current_x2 == 0:
                if n < r_on*dt/nu:
                    x2 = 1
                else:
                    x2 = 0
            elif current_x2 == 1:
                if n < r_off*dt*nu:
                    x2 = 0
                else:
                    x2 = 1
        else:
            if current_x2 == 0:
                if n < r_on*dt*nu:
                    x2 = 1
                else:
                    x2 = 0
            elif current_x2 == 1:
                if n < r_off*dt/nu:
                    x2 = 0
                else:
                    x2 = 1
        x2_state[current_ti] = x2
        current_x2 = x2
        current_ti += 1


def calc_MI_btw_states(x1_state,x2_state,trec=0,debug=False):

    nsteps=len(x1_state)
    ## calculate x1, x2 MI
    pp_1_0=0.0
    pp_1_1=0.0
    pp_2_0=0.0
    pp_2_1=0.0
    pp_00=0.0
    pp_10=0.0
    pp_01=0.0
    pp_11=0.0
    for i in range(trec,nsteps):
        x1 = x1_state[i]
        x2 = x2_state[i]
    
        if x1 == 0:
            pp_1_0 +=1.0
        else:
            pp_1_1 +=1.0
    
        if x2 == 0:
            pp_2_0 +=1.0
        else:
            pp_2_1 +=1.0
    
        if x1 ==0:
            if x2 ==0: pp_00 += 1.0
            else: pp_01 += 1.0
        else:
            if x2 ==0: pp_10 += 1.0
            else: pp_11 += 1.0
    nn = float(nsteps-trec)
    pp_1_0 = pp_1_0/nn
    pp_1_1 = pp_1_1/nn
    pp_2_0 = pp_2_0/nn
    pp_2_1 = pp_2_1/nn
    pp_00 = pp_00/nn
    pp_01 = pp_01/nn
    pp_10 = pp_10/nn
    pp_11 = pp_11/nn
    if debug: print("PRECHECK: ",pp_1_0,pp_1_1, pp_1_0+pp_1_1)
    if debug: print("PRECHECK: ",pp_2_0,pp_2_1, pp_2_0+pp_2_1)
    if debug: print("PRECHECK: ",pp_00,pp_01,pp_10,pp_11, pp_00+pp_01+pp_10+pp_11)
    pre12MI = 0.0
    if pp_00 > 0:
        pre12MI += pp_00*log2( pp_00/(pp_1_0*pp_2_0))
    if pp_10 > 0:
        pre12MI += pp_10*log2( pp_10/(pp_1_1*pp_2_0)) 
    if pp_01 > 0:
        pre12MI += pp_01*log2( pp_01/(pp_1_0*pp_2_1)) 
    if pp_11 > 0:
        pre12MI += pp_11*log2( pp_11/(pp_1_1*pp_2_1)) 
    if debug: print('pre12MI:', pre12MI)
    return pre12MI



# MI calculation functions
def dLdt(L,I,r_on,r_off,theta):
    return r_on *(1.0+exp(-L)) -r_off*(1.0+exp(L)) + I - theta


def get_mi(fpath,nsyn,dt=0.1):

#     fpath='%s/%s_%s_%s/'%(path,nsyn,alphaidx,trial)

    f=open("%svars.txt"%fpath)
    t=[]
    v=[]
    vd=[]
    for l in f:
        stuff=l.split()
        t.append(float(stuff[0]))
        v.append(float(stuff[1]))
        vd.append(float(stuff[2]))

    f.close()
    nsyn = 2

#     f=open("%sweights.txt"%fpath)

#     weights=[[] for i in range(nsyn)]
#     for l in f:
#         stuff=l.split()
#         for wi in range(1,nsyn+1):
#             weights[wi-1].append(float(stuff[wi]))
#     f.close()

    f=open("%sx_state.txt"%fpath)

    x1=[]
    x2=[]
    for l in f:
        stuff=l.split()
        x1.append(float(stuff[1]))
        x2.append(float(stuff[2]))
    while len(t)>len(x1):
        x1.append(float(stuff[1]))
        x2.append(float(stuff[2]))
    f.close()
    x1_state=np.array(x1)
    x2_state=np.array(x2)


    pre_spike_lis=[]
    for si in range(nsyn):
        pre_spikes=[]
        f=open('%spre_%s.txt'%(fpath,si+1))
        for l in f:
            pre_spikes.append(float(l))
        f.close()
        pre_spike_lis.append(pre_spikes)

    post_spikes = []
    last_v = v[0]
    for i in range(len(t)):
        current_v = v[i]
        if current_v > 0.0 and last_v <=0.0:
            post_spikes.append(t[i])
        last_v = current_v
        
    q_on1 = 0.05
    q_off1 = 0.001
    r_on_slow = 16.7/1000.0
    r_on = r_on_slow
    r_off = r_on*2.0
    theta1 = q_on1 - q_off1
    theta2 = q_on1 - q_off1
    nsteps = len(t)
    opt_w = [log(q_on1/q_off1), log(q_on1/q_off1)]


    ######
    #  pre-spike MI
    ######



    #WARNING: L is integrated with a simple euler method
    pre_spike_lis_save = []
    for pre_spike_li in pre_spike_lis:
        pre_spike_lis_save.append(pre_spike_li.copy())

    ttt=0.0
    ti = 0
    L1 = 0
    Lt1 = []
    Lt1.append(L1)
    L2 = 0
    Lt2 = []
    Lt2.append(L2)

    t_cut = 5
    trec = int(t_cut/dt)
    tcut = []
    L1cut = []
    L2cut = []
    x1cut = []
    x2cut = []

#     print('thetas: ',theta1,theta2)

    while ti < nsteps:
        dL1 = dLdt(L1,0,r_on,r_off,theta1) * dt
        L1 = L1 + dL1
        dL2 = dLdt(L2,0,r_on,r_off,theta2) * dt
        L2 = L2 + dL2
        try:
            for si in [0]:
                if ttt>pre_spike_lis[si][0]:
                    #print('spike detect 1 : %s, %s, %s'%(ttt,si,pre_spike_lis[si][0]))
                    pre_spike_lis[si].pop(0)
                    L1 += opt_w[si]
                    #print(dLdt(L1,0,r_on,r_off,theta1),L1,dL1)
        except IndexError:
            pass
        try:
            for si in [1]:
                if ttt>pre_spike_lis[si][0]:
                    #print('spike detect 2 : %s, %s, %s'%(ttt,si,pre_spike_lis[si][0]))
                    pre_spike_lis[si].pop(0)
                    L2 += opt_w[si]
                    #print(dLdt(L2,0,r_on,r_off,theta2),L2,dL2)
        except IndexError:
            pass
        if ti>=trec:
            tcut.append(ttt)
            L1cut.append(L1)
            L2cut.append(L2)
            x1cut.append(x1_state[ti])
            x2cut.append(x2_state[ti])
        ttt += dt
        ti += 1
        Lt1.append(L1)
        Lt2.append(L2)

    tcut = np.array(tcut)
    L1cut = np.array(L1cut)
    L2cut = np.array(L2cut)
    x1cut = np.array(x1cut)
    x2cut = np.array(x2cut)

    Lt1 = np.array(Lt1)
    Lt2 = np.array(Lt2)

    p_11 = 1.0/(1.0 + np.exp(-Lt1))
    p_01 = 1.0-p_11
    p_12 = 1.0/(1.0 + np.exp(-Lt2))
    p_02 = 1.0-p_12


    ave_x1 = 0.0
    ave_hx1y = 0.0
    ave_x2 = 0.0
    ave_hx2y = 0.0
    hxy_test = []
    for i in range(trec,nsteps):
    #     print(i)
        x1 = float(x1_state[i])
        x2 = float(x2_state[i])

        L1 = Lt1[i]
        p11 = p_11[i]
        p01 = p_01[i]
        L2 = Lt2[i]
        p12 = p_12[i]
        p02 = p_02[i]
    #     if i %1000==0: print(i,x,L,p1,p0)

        ave_x1 += x1
        ave_hx1y += x1*log2(p11) + (1.0-x1)*log2(1.0-p11)
        ave_x2 += x2
        ave_hx2y += x2*log2(p12) + (1.0-x2)*log2(1.0-p12)
        #print(i,x1,x2,ave_x1,ave_x2)

    ave_x1 = ave_x1 / float(nsteps-trec)
    ave_hx1y = -ave_hx1y /float(nsteps-trec)
#     print('ave_x1=',ave_x1)
    hxx1 = - ave_x1*log2(ave_x1) - (1.0 - ave_x1)*log2(1.0-ave_x1)

    MI1 = hxx1 - ave_hx1y

    ave_x2 = ave_x2 / float(nsteps-trec)
    ave_hx2y = -ave_hx2y /float(nsteps-trec)
#     print('ave_x2=',ave_x2)
    hxx2 = - ave_x2*log2(ave_x2) - (1.0 - ave_x2)*log2(1.0-ave_x2)

    MI2 = hxx2 - ave_hx2y

    ## calculate x1, x2 MI
    pre12MI=calc_MI_btw_states(x1_state,x2_state,trec=trec)
#     print('pre12MI:', pre12MI)


    #g.write("preHxx = %s\n"%hxx)
    #g.write("preHxy = %s\n"%ave_hxy)
    #g.write("preMI = %s\n"%MI)

#     print("1.preHxx = ",hxx1)
#     print("1.preHxy = ",ave_hx1y)
#     print("1.preMI = %s"%MI1)
#     print("2.preHxx = ",hxx2)
#     print("2.preHxy = ",ave_hx2y)
#     print("2.preMI = %s"%MI2)
    # ## compute MI

    preMI1=MI1
    preMI2=MI2

#     print("computing MI")
    trec = int(t_cut/dt)

    ## get ests 1
    time_on = 0.0
    time_off = 0.0
    for i in range(trec,nsteps):
        if x1_state[i]==0:
            time_off += dt
        elif x1_state[i] ==1:
            time_on += dt
        else: raise Exception()
#     print(time_on, time_off, time_on + time_off)


    est_qon = 0.0
    est_qoff = 0.0

    for spike_time in post_spikes:
        ti = int(spike_time / dt)
        if (ti < trec): continue
        state_at_t = x1_state[ti]
        if state_at_t ==1: 
            est_qon += 1.0
        elif state_at_t == 0:
            est_qoff += 1.0
        else: raise Exception()

    est_qon = est_qon/time_on
    est_qoff = est_qoff/time_off

    est_w = log(est_qon/est_qoff)
    theta = est_qon - est_qoff

    est_qon1 = est_qon
    est_qoff1 = est_qoff
    est_w1 = est_w
    theta1 = theta

    ## get ests 2
    time_on = 0.0
    time_off = 0.0
    for i in range(trec,nsteps):
        if x2_state[i]==0:
            time_off += dt
        elif x2_state[i] ==1:
            time_on += dt
        else: raise Exception()
#     print(time_on, time_off, time_on + time_off)


    est_qon = 0.0
    est_qoff = 0.0

    for spike_time in post_spikes:
        ti = int(spike_time / dt)
        if (ti < trec): continue
        state_at_t = x2_state[ti]
        if state_at_t ==1: 
            est_qon += 1.0
        elif state_at_t == 0:
            est_qoff += 1.0
        else: raise Exception()

    est_qon = est_qon/time_on
    est_qoff = est_qoff/time_off

    est_w = log(est_qon/est_qoff)
    theta = est_qon - est_qoff

    est_qon2 = est_qon
    est_qoff2 = est_qoff
    est_w2 = est_w
    theta2 = theta


    #WARNING: L is integrated with a simple euler method
    ## MI 1 
    theta = theta1
    tmp_spikes = []
    for tspik in post_spikes:
        tmp_spikes.append(tspik)
    ttt=0.0
    ti = 0
    L1 = 0
    L2 = 0
    Lt1 = []
    Lt1.append(L1)
    Lt2 = []
    Lt2.append(L2)


    tcut = []
    vcut = []
    L1cut = []
    L2cut = []
    x1cut = []
    x2cut = []

    while ti < nsteps:
        vvv = v[ti]
        dL1 = dLdt(L1,0,r_on,r_off,theta1) * dt
        L1 = L1 + dL1
        dL2 = dLdt(L2,0,r_on,r_off,theta2) * dt
        L2 = L2 + dL2
        try:
            if ttt>tmp_spikes[0]:
                #print('spike detect: %s, %s'%(ttt,tmp_spikes[0]))
                tmp_spikes.pop(0)
                L1 += est_w1
                L2 += est_w2
        #         print(vvv,dLdt(L,vvv,r_on,r_off,0),dL)
        except IndexError:
            pass
        if ti>=trec:
            tcut.append(ttt)
            vcut.append(v[ti])
            L1cut.append(L1)
            L2cut.append(L2)
            x1cut.append(x1_state[ti])
            x2cut.append(x2_state[ti])
        ttt += dt
        ti += 1
        Lt1.append(L1)
        Lt2.append(L2)

    tcut = np.array(tcut)
    vcut = np.array(vcut)
    L1cut = np.array(L1cut)
    L2cut = np.array(L2cut)
    x1cut = np.array(x1cut)
    x2cut = np.array(x2cut)

    Lt1 = np.array(Lt1)
    Lt2 = np.array(Lt2)

    p_11 = 1.0/(1.0 + np.exp(-Lt1))
    p_01 = 1.0-p_11
    p_12 = 1.0/(1.0 + np.exp(-Lt2))
    p_02 = 1.0-p_12



    ave_x1 = 0.0
    ave_hx1y = 0.0
    ave_x2 = 0.0
    ave_hx2y = 0.0
    hxy_test = []

    for i in range(trec,nsteps):
    #     print(i)
        x1 = float(x1_state[i])
        x2 = float(x2_state[i])

        L1 = Lt1[i]
        p11 = p_11[i]
        p01 = p_01[i]
        L2 = Lt2[i]
        p12 = p_12[i]
        p02 = p_02[i]
    #    if i %1000==0: print(i,x1,x2,L1,L2,p11,p01,p12,p02)

        ave_x1 += float(x1)
        ave_hx1y += x1*log2(p11) + (1.0-x1)*log2(1.0-p11)
        ave_x2 += float(x2)
        ave_hx2y += x2*log2(p12) + (1.0-x2)*log2(1.0-p12)

    ave_x1 = ave_x1 / float(nsteps-trec)
    ave_hx1y = -ave_hx1y /float(nsteps-trec)
    hxx1 = - ave_x1*log2(ave_x1) - (1.0 - ave_x1)*log2(1.0-ave_x1)

    MI1 = hxx1 - ave_hx1y

    ave_x2 = ave_x2 / float(nsteps-trec)
    ave_hx2y = -ave_hx2y /float(nsteps-trec)
    hxx2 = - ave_x2*log2(ave_x2) - (1.0 - ave_x2)*log2(1.0-ave_x2)

    MI2 = hxx2 - ave_hx2y

#     print("Hxx1 = ",hxx1)
#     print("Hx1y = ",ave_hx1y)
#     print("Hxx2 = ",hxx2)
#     print("Hx2y = ",ave_hx2y)
#     print("MI1 = %s"%MI1)
#     print("MI2 = %s"%MI2)
#     print("F1 = %s"%(MI1/preMI1))
#     print("F2 = %s"%(MI2/preMI2))
#     print("pre12MI = %s"%pre12MI)

    return pre12MI, preMI1, preMI2, MI1, MI2

def get_mi_sliced(fpath,tint):
    ## Calculate pre1MI, pre2MI, pre12MI (MI between hidden states), MI1 (MI between hidden state 1 and post spikes), and MI2.
    ## calculate all values over 30 second intervals and return lists.
    nsyn = 2

    ## get time and voltage vectors
    f=open("%svars.txt"%fpath)
    t=[]
    v=[]
    vd=[]
    for l in f:
        stuff=l.split()
        t.append(float(stuff[0]))
        v.append(float(stuff[1]))
        vd.append(float(stuff[2]))

    f.close()

#     f=open("%sweights.txt"%fpath)

#     weights=[[] for i in range(nsyn)]
#     for l in f:
#         stuff=l.split()
#         for wi in range(1,nsyn+1):
#             weights[wi-1].append(float(stuff[wi]))
#     f.close()


#   get hidden state arrays
    f=open("%sx_state.txt"%fpath)

    x1=[]
    x2=[]
    for l in f:
        stuff=l.split()
        x1.append(float(stuff[1]))
        x2.append(float(stuff[2]))
    while len(t)>len(x1):
        x1.append(float(stuff[1]))
        x2.append(float(stuff[2]))
    f.close()
    x1_state=np.array(x1)
    x2_state=np.array(x2)


## get presynaptic spikes
    pre_spike_lis=[]
    for si in range(nsyn):
        pre_spikes=[]
        f=open('%spre_%s.txt'%(fpath,si+1))
        for l in f:
            pre_spikes.append(float(l))
        f.close()
        pre_spike_lis.append(pre_spikes)

##get post spikes
    post_spikes = []
    last_v = v[0]
    for i in range(len(t)):
        current_v = v[i]
        if current_v > 0.0 and last_v <=0.0:
            post_spikes.append(t[i])
        last_v = current_v
        
#presynaptic parameters!!
    q_on1 = 0.05
    q_off1 = 0.001
    r_on_slow = 16.7/1000.0
    r_on = r_on_slow
    r_off = r_on*2.0
    dt=0.1
    theta1 = q_on1 - q_off1
    theta2 = q_on1 - q_off1
    nsteps = len(t)
    opt_w = [log(q_on1/q_off1), log(q_on1/q_off1)]


    ######
    #  pre-spike MI
    ######



    #WARNING: L is integrated with a simple euler method
    pre_spike_lis_save = []
    for pre_spike_li in pre_spike_lis:
        pre_spike_lis_save.append(pre_spike_li.copy())

    ttt=0.0
    ti = 0
    L1 = 0
    Lt1 = []
    Lt1.append(L1)
    L2 = 0
    Lt2 = []
    Lt2.append(L2)

    t_cut = 5
    trec = 0
    tcut = []
    L1cut = []
    L2cut = []
    x1cut = []
    x2cut = []

#     print('thetas: ',theta1,theta2)

    while ti < nsteps:
        dL1 = dLdt(L1,0,r_on,r_off,theta1) * dt
        L1 = L1 + dL1
        dL2 = dLdt(L2,0,r_on,r_off,theta2) * dt
        L2 = L2 + dL2
        try:
            for si in [0]:
                if ttt>pre_spike_lis[si][0]:
                    #print('spike detect 1 : %s, %s, %s'%(ttt,si,pre_spike_lis[si][0]))
                    pre_spike_lis[si].pop(0)
                    L1 += opt_w[si]
                    #print(dLdt(L1,0,r_on,r_off,theta1),L1,dL1)
        except IndexError:
            pass
        try:
            for si in [1]:
                if ttt>pre_spike_lis[si][0]:
                    #print('spike detect 2 : %s, %s, %s'%(ttt,si,pre_spike_lis[si][0]))
                    pre_spike_lis[si].pop(0)
                    L2 += opt_w[si]
                    #print(dLdt(L2,0,r_on,r_off,theta2),L2,dL2)
        except IndexError:
            pass
        if ti>=trec:
            tcut.append(ttt)
            L1cut.append(L1)
            L2cut.append(L2)
            x1cut.append(x1_state[ti])
            x2cut.append(x2_state[ti])
        ttt += dt
        ti += 1
        Lt1.append(L1)
        Lt2.append(L2)

    tcut = np.array(tcut)
    L1cut = np.array(L1cut)
    L2cut = np.array(L2cut)
    x1cut = np.array(x1cut)
    x2cut = np.array(x2cut)

    Lt1 = np.array(Lt1)
    Lt2 = np.array(Lt2)

    p_11 = 1.0/(1.0 + np.exp(-Lt1))
    p_01 = 1.0-p_11
    p_12 = 1.0/(1.0 + np.exp(-Lt2))
    p_02 = 1.0-p_12

    t_stop = t[-1]
    n_intervals = int(t_stop/tint)
    preMI1li=[]
    preMI2li=[]
    
    print('tstop = %s, n_intervals = %s'%(t_stop,n_intervals))
    # Calculate pre MI's for each interval
    for i_int in range(n_intervals):
        int_start_t = i_int*tint
        int_stop_t = (i_int+1.0)*tint
        #print('pre %s, %s to %s'%(i_int,int_start_t,int_stop_t))
        int_start_i = int(int_start_t/dt)
        int_stop_i  = int(int_stop_t/dt)
        #print('pre indices from',int_start_i, 'to',int_stop_i)

        ave_x1 = 0.0
        ave_hx1y = 0.0
        ave_x2 = 0.0
        ave_hx2y = 0.0
        
        for i in range(int_start_i,int_stop_i):
        #     print(i)
            x1 = float(x1_state[i])
            x2 = float(x2_state[i])
    
            L1 = Lt1[i]
            p11 = p_11[i]
            p01 = p_01[i]
            L2 = Lt2[i]
            p12 = p_12[i]
            p02 = p_02[i]
        #     if i %1000==0: print(i,x,L,p1,p0)
    
            ave_x1 += x1
            ave_hx1y += x1*log2(p11) + (1.0-x1)*log2(1.0-p11)
            ave_x2 += x2
            ave_hx2y += x2*log2(p12) + (1.0-x2)*log2(1.0-p12)
            #print(i,x1,x2,ave_x1,ave_x2)
    
        ave_x1 = ave_x1 / float(int_stop_i-int_start_i)
        ave_hx1y = -ave_hx1y /float(int_stop_i-int_start_i)
        hxx1 = - ave_x1*log2(ave_x1) - (1.0 - ave_x1)*log2(1.0-ave_x1)
    
        MI1 = hxx1 - ave_hx1y
    
        ave_x2 = ave_x2 / float(int_stop_i-int_start_i)
        ave_hx2y = -ave_hx2y /float(int_stop_i-int_start_i)
        hxx2 = - ave_x2*log2(ave_x2) - (1.0 - ave_x2)*log2(1.0-ave_x2)

        MI2 = hxx2 - ave_hx2y

        preMI1li.append(MI1)
        preMI2li.append(MI2)

    ## calculate x1, x2 MI
    pre12MI=calc_MI_btw_states(x1_state,x2_state,trec=trec)
#     print('pre12MI:', pre12MI)


    #g.write("preHxx = %s\n"%hxx)
    #g.write("preHxy = %s\n"%ave_hxy)
    #g.write("preMI = %s\n"%MI)

#     print("1.preHxx = ",hxx1)
#     print("1.preHxy = ",ave_hx1y)
#     print("1.preMI = %s"%MI1)
#     print("2.preHxx = ",hxx2)
#     print("2.preHxy = ",ave_hx2y)
#     print("2.preMI = %s"%MI2)
    # ## compute MI


#     print("computing MI")
    trec = int(t_cut/dt)

    ## get ests 1
    time_on = 0.0
    time_off = 0.0
    for i in range(trec,nsteps):
        if x1_state[i]==0:
            time_off += dt
        elif x1_state[i] ==1:
            time_on += dt
        else: raise Exception()
#     print(time_on, time_off, time_on + time_off)


    est_qon = 0.0
    est_qoff = 0.0

    for spike_time in post_spikes:
        ti = int(spike_time / dt)
        if (ti < trec): continue
        state_at_t = x1_state[ti]
        if state_at_t ==1: 
            est_qon += 1.0
        elif state_at_t == 0:
            est_qoff += 1.0
        else: raise Exception()

    est_qon = est_qon/time_on
    est_qoff = est_qoff/time_off

    est_w = log(est_qon/est_qoff)
    theta = est_qon - est_qoff

    est_qon1 = est_qon
    est_qoff1 = est_qoff
    est_w1 = est_w
    theta1 = theta

    ## get ests 2
    time_on = 0.0
    time_off = 0.0
    for i in range(trec,nsteps):
        if x2_state[i]==0:
            time_off += dt
        elif x2_state[i] ==1:
            time_on += dt
        else: raise Exception()
#     print(time_on, time_off, time_on + time_off)


    est_qon = 0.0
    est_qoff = 0.0

    for spike_time in post_spikes:
        ti = int(spike_time / dt)
        if (ti < trec): continue
        state_at_t = x2_state[ti]
        if state_at_t ==1: 
            est_qon += 1.0
        elif state_at_t == 0:
            est_qoff += 1.0
        else: raise Exception()

    est_qon = est_qon/time_on
    est_qoff = est_qoff/time_off

    est_w = log(est_qon/est_qoff)
    theta = est_qon - est_qoff

    est_qon2 = est_qon
    est_qoff2 = est_qoff
    est_w2 = est_w
    theta2 = theta


    #WARNING: L is integrated with a simple euler method
    ## MI 1 
    theta = theta1
    tmp_spikes = []
    for tspik in post_spikes:
        tmp_spikes.append(tspik)
    ttt=0.0
    ti = 0
    L1 = 0
    L2 = 0
    Lt1 = []
    Lt1.append(L1)
    Lt2 = []
    Lt2.append(L2)


    tcut = []
    vcut = []
    L1cut = []
    L2cut = []
    x1cut = []
    x2cut = []

    while ti < nsteps:
        vvv = v[ti]
        dL1 = dLdt(L1,0,r_on,r_off,theta1) * dt
        L1 = L1 + dL1
        dL2 = dLdt(L2,0,r_on,r_off,theta2) * dt
        L2 = L2 + dL2
        try:
            if ttt>tmp_spikes[0]:
                #print('spike detect: %s, %s'%(ttt,tmp_spikes[0]))
                tmp_spikes.pop(0)
                L1 += est_w1
                L2 += est_w2
        #         print(vvv,dLdt(L,vvv,r_on,r_off,0),dL)
        except IndexError:
            pass
        if ti>=trec:
            tcut.append(ttt)
            vcut.append(v[ti])
            L1cut.append(L1)
            L2cut.append(L2)
            x1cut.append(x1_state[ti])
            x2cut.append(x2_state[ti])
        ttt += dt
        ti += 1
        Lt1.append(L1)
        Lt2.append(L2)

    tcut = np.array(tcut)
    vcut = np.array(vcut)
    L1cut = np.array(L1cut)
    L2cut = np.array(L2cut)
    x1cut = np.array(x1cut)
    x2cut = np.array(x2cut)

    Lt1 = np.array(Lt1)
    Lt2 = np.array(Lt2)

    p_11 = 1.0/(1.0 + np.exp(-Lt1))
    p_01 = 1.0-p_11
    p_12 = 1.0/(1.0 + np.exp(-Lt2))
    p_02 = 1.0-p_12

    MI1li=[]
    MI2li=[]


    for i_int in range(n_intervals):
        int_start_t = i_int*tint
        int_stop_t = (i_int+1.0)*tint
        #print('post %s, %s to %s'%(i_int,int_start_t,int_stop_t))
        int_start_i = int(int_start_t/dt)
        int_stop_i  = int(int_stop_t/dt)
        #print('post from',int_start_i, 'to',int_stop_i)

        ave_x1 = 0.0
        ave_hx1y = 0.0
        ave_x2 = 0.0
        ave_hx2y = 0.0
    
        for i in range(int_start_i,int_stop_i):
        #     print(i)
            x1 = float(x1_state[i])
            x2 = float(x2_state[i])
    
            L1 = Lt1[i]
            p11 = p_11[i]
            p01 = p_01[i]
            L2 = Lt2[i]
            p12 = p_12[i]
            p02 = p_02[i]
        #    if i %1000==0: print(i,x1,x2,L1,L2,p11,p01,p12,p02)
    
            ave_x1 += float(x1)
            ave_hx1y += x1*log2(p11) + (1.0-x1)*log2(1.0-p11)
            ave_x2 += float(x2)
            ave_hx2y += x2*log2(p12) + (1.0-x2)*log2(1.0-p12)
    
        ave_x1 = ave_x1 / float(int_stop_i-int_start_i)
        ave_hx1y = -ave_hx1y /float(int_stop_i-int_start_i)
        hxx1 = - ave_x1*log2(ave_x1) - (1.0 - ave_x1)*log2(1.0-ave_x1)
    
        MI1 = hxx1 - ave_hx1y
    
        ave_x2 = ave_x2 / float(int_stop_i-int_start_i)
        ave_hx2y = -ave_hx2y /float(int_stop_i-int_start_i)
        hxx2 = - ave_x2*log2(ave_x2) - (1.0 - ave_x2)*log2(1.0-ave_x2)
    
        MI2 = hxx2 - ave_hx2y

        MI1li.append(MI1)
        MI2li.append(MI2)
    
    #     print("Hxx1 = ",hxx1)
    #     print("Hx1y = ",ave_hx1y)
    #     print("Hxx2 = ",hxx2)
    #     print("Hx2y = ",ave_hx2y)
    #     print("MI1 = %s"%MI1)
    #     print("MI2 = %s"%MI2)
    #     print("F1 = %s"%(MI1/preMI1))
    #     print("F2 = %s"%(MI2/preMI2))
    #     print("pre12MI = %s"%pre12MI)
    
    return pre12MI, preMI1li, preMI2li, MI1li, MI2li
