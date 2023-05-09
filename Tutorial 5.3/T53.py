import numpy as np
import matplotlib.pyplot as plt

C = 1E-9 # 1 nF
R = 10E6 # 10 mOhm  
E = -70E-3 # -70 mV
Vth = -54E-3 # -54 mV
V_reset = -80E-3 # -80 mV
E_12_rev = -70E-3 # -70 mV
E_21_rev = -70E-3 # -70 mV
G_12 = 1E-6 # 1 uS
G_21 = 1E-6 # 1 uS
tau_syn = 10E-3 # 10 ms


pR = 1
D1 = 1
D2 = 1

dt = 0.1E-3 
t = np.arange(0, 6, dt)

Iapp1 = np.zeros(len(t))
Iapp1[:] = 2E-9
Iapp1[0:1000] += 3E-9
Iapp2 = np.zeros(len(t))
Iapp2[:] = 2E-9
Iapp2[30000:31000] += 3E-9

sigma = 0 
n = np.random.randn(len(t))

Vm1 = np.zeros(len(t))
Vm2 = np.zeros(len(t))

s1 = np.zeros(len(t))
s2 = np.zeros(len(t))

########
#Part A
########

##Part ii

for i in range(1, len(t)):
    s1[i] = s1[i-1] + dt *(-s1[i-1]/tau_syn)
    s2[i] = s2[i-1] + dt *(-s2[i-1]/tau_syn)
    Vm1[i] = Vm1[i-1] + ((E-Vm1[i-1])/R + G_21*s2[i]*(E_21_rev-Vm1[i-1]) + Iapp1[i-1] + sigma*n[i-1]) * dt/C
    Vm2[i] = Vm2[i-1] + ((E-Vm2[i-1])/R + G_12*s1[i]*(E_12_rev-Vm2[i-1]) + Iapp2[i-1] + sigma*n[i-1]) * dt/C

    if Vm1[i] > Vth:
        Vm1[i] = V_reset
        s1[i] = s1[i] + pR*D1*(1-s1[i])
    if Vm2[i] > Vth:
        Vm2[i] = V_reset
        s2[i] = s2[i] + pR*D2*(1-s2[i])
    

plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t, Vm1, label='Vm1', color='orange')
plt.ylabel('Membrane Potential (V)')
plt.ylim(-.085, -.05)
plt.title('Part ii: Vm1 and Vm2')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t, Vm2, label='Vm2', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.ylim(-.085, -.05)
plt.figure(1).set_size_inches(13, 8)
plt.legend()
plt.savefig('Part ii: Vm1 and Vm2.png')

plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t, s1, label='s1', color='orange')
plt.ylabel('Synaptic Conductance')
plt.title('Part ii: s1 and s2')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t, s2, label='s2', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Synaptic Conductance')
plt.figure(2).set_size_inches(13, 8)
plt.legend() 
plt.savefig('Part ii: s1 and s2.png')
plt.show()

##Part iii

sigma = 50E-12/np.sqrt(dt)

Iapp1[:] = 2E-9
Iapp2[:] = 2E-9

for i in range(1, len(t)):
    s1[i] = s1[i-1] + dt *(-s1[i-1]/tau_syn)
    s2[i] = s2[i-1] + dt *(-s2[i-1]/tau_syn)
    # multiply sigma by a randn value
    sigma_1 = sigma * np.random.randn()
    sigma_2 = sigma * np.random.randn()
    Vm1[i] = Vm1[i-1] + ((E-Vm1[i-1])/R + G_21*s2[i]*(E_21_rev-Vm1[i-1]) + Iapp1[i-1] + sigma_1*n[i-1]) * dt/C
    Vm2[i] = Vm2[i-1] + ((E-Vm2[i-1])/R + G_12*s1[i]*(E_12_rev-Vm2[i-1]) + Iapp2[i-1] + sigma_2*n[i-1]) * dt/C

    if Vm1[i] > Vth:
        Vm1[i] = V_reset
        s1[i] = s1[i] + pR*D1*(1-s1[i])
    if Vm2[i] > Vth:
        Vm2[i] = V_reset
        s2[i] = s2[i] + pR*D2*(1-s2[i])


plt.figure(3)
plt.subplot(2,1,1)
plt.plot(t, Vm1, label='Vm1', color='orange')
plt.ylabel('Membrane Potential (V)')
plt.title('Part iii: Vm1 and Vm2')
plt.subplot(2,1,2)
plt.plot(t, Vm2, label='Vm2', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.figure(3).set_size_inches(13, 8)
plt.legend()
plt.savefig('Part iii: Vm1 and Vm2.png')

plt.figure(4)
plt.subplot(2,1,1)
plt.plot(t, s1, label='s1', color='orange')
plt.ylabel('Synaptic Conductance')
plt.title('Part iii: s1 and s2')
plt.subplot(2,1,2)
plt.plot(t, s2, label='s2', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Synaptic Conductance')
plt.figure(4).set_size_inches(13, 8)
plt.legend()
plt.savefig('Part iii: s1 and s2.png')
plt.show()


## Part iv


states = np.zeros(len(t))
states[0] = 1
spike_times = np.zeros(len(t))

for i in range(1, len(t)):
    s1[i] = s1[i-1] + dt *(-s1[i-1]/tau_syn)
    s2[i] = s2[i-1] + dt *(-s2[i-1]/tau_syn)
    sigma_1 = sigma * np.random.randn()
    sigma_2 = sigma * np.random.randn()
    Vm1[i] = Vm1[i-1] + ((E-Vm1[i-1])/R + G_21*s2[i]*(E_21_rev-Vm1[i-1]) + Iapp1[i-1] + sigma_1*n[i-1]) * dt/C
    Vm2[i] = Vm2[i-1] + ((E-Vm2[i-1])/R + G_12*s1[i]*(E_12_rev-Vm2[i-1]) + Iapp2[i-1] + sigma_2*n[i-1]) * dt/C

    if Vm1[i] > Vth:
        Vm1[i] = V_reset
        s1[i] = s1[i] + pR*D1*(1-s1[i])
        states[i] = 2
        spike_times[i] = t[i]
    if Vm2[i] > Vth:
        Vm2[i] = V_reset
        s2[i] = s2[i] + pR*D2*(1-s2[i])
        states[i] = 1
        spike_times[i] = t[i]

spike_times = np.where(spike_times != 0)[0]
diff_spike_times = np.diff(spike_times)


plt.figure(5)
plt.subplot(2,1,1)
plt.plot(t, Vm1, label='Vm1', color='orange')
plt.ylabel('Membrane Potential (V)')
plt.title('Part iv: Vm1 and Vm2')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t, Vm2, label='Vm2', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.legend()
plt.figure(5).set_size_inches(12, 8)
plt.savefig('part iv Vm1 and Vm2.png')

plt.figure(6)
plt.subplot(3,1,1)
plt.plot(t, states, label='states', color='purple')
plt.title('States: sigma = 50pA')
plt.ylabel('State')
plt.legend()
plt.subplot(3,1,2)
plt.hist(diff_spike_times[::2], bins=100, color='purple')
plt.title('Histogram of one state switching times (even entries)')
plt.ylabel('Frequency')
plt.subplot(3,1,3)
plt.hist(diff_spike_times[1::2], bins=100, color='purple')
plt.title('Histogram of other state switching times (odd entries)')
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.figure(6).set_size_inches(12, 10)
plt.savefig('part iv.png')
plt.show()


#######
#Part B
#######

## Part v

pR = 0.2
tau_D = 0.2
D1 = np.zeros(len(t))
D2 = np.zeros(len(t))

sigma = 0

Iapp1 = np.zeros(len(t))
Iapp1[:] = 2E-9
Iapp1[0:1000] += 3E-9
Iapp2 = np.zeros(len(t))
Iapp2[:] = 2E-9
Iapp2[30000:31000] += 3E-9

for i in range(1, len(t)):
    D1[i] = D1[i-1] + dt * ((1-D1[i-1])/tau_D)
    D2[i] = D2[i-1] + dt * ((1-D2[i-1])/tau_D)
    s1[i] = s1[i-1] + dt *(-s1[i-1]/tau_syn)
    s2[i] = s2[i-1] + dt *(-s2[i-1]/tau_syn)
    Vm1[i] = Vm1[i-1] + ((E-Vm1[i-1])/R + G_21*s2[i]*(E_21_rev-Vm1[i-1]) + Iapp1[i-1] + sigma*n[i-1]) * dt/C
    Vm2[i] = Vm2[i-1] + ((E-Vm2[i-1])/R + G_12*s1[i]*(E_12_rev-Vm2[i-1]) + Iapp2[i-1] + sigma*n[i-1]) * dt/C

    if Vm1[i] > Vth:
        Vm1[i] = V_reset
        s1[i] = s1[i] + pR*D1[i]*(1-s1[i])
        D1[i] = D1[i] * (1-pR)
    if Vm2[i] > Vth:
        Vm2[i] = V_reset
        s2[i] = s2[i] + pR*D2[i]*(1-s2[i])
        D2[i] = D2[i] * (1-pR)

plt.figure(7)
plt.subplot(2,1,1)
plt.plot(t, Vm1, label='Vm', color='orange')
plt.ylabel('Membrane Potential (V)')
plt.title('Part v: Vm1 and Vm2')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t, Vm2, label='Vm', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.legend() 
plt.figure(7).set_size_inches(12, 8)
plt.savefig('part_v_Vm.png')

plt.figure(8)
plt.subplot(2,1,1)
plt.plot(t, s1, label='s1', color='orange')
plt.ylabel('Synaptic Conductance')
plt.title('Part v: s1 and s2')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t, s2, label='s2', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Synaptic Conductance')
plt.legend()
plt.figure(8).set_size_inches(12, 8)
plt.savefig('part_v_s.png')

plt.figure(9)
plt.subplot(2,1,1)
plt.plot(t, D1, label='D1', color='orange')
plt.ylabel('Depression')
plt.title('Part v: D1 and D2')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t, D2, label='D2', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Depression')
plt.legend()
plt.figure(9).set_size_inches(12, 8)
plt.savefig('part_v_D.png')
plt.show()


## Part vi

sigma = 5E-12/np.sqrt(dt)

states = np.zeros(len(t))
states[0] = 1
spike_times = np.zeros(len(t))

for i in range(1, len(t)):
    D1[i] = D1[i-1] + dt * ((1-D1[i-1])/tau_D)
    D2[i] = D2[i-1] + dt * ((1-D2[i-1])/tau_D)
    s1[i] = s1[i-1] + dt *(-s1[i-1]/tau_syn)
    s2[i] = s2[i-1] + dt *(-s2[i-1]/tau_syn)
    sigma_1 = sigma * np.random.randn()
    sigma_2 = sigma * np.random.randn()
    Vm1[i] = Vm1[i-1] + ((E-Vm1[i-1])/R + G_21*s2[i]*(E_21_rev-Vm1[i-1]) + Iapp1[i-1] + sigma_1*n[i-1]) * dt/C
    Vm2[i] = Vm2[i-1] + ((E-Vm2[i-1])/R + G_12*s1[i]*(E_12_rev-Vm2[i-1]) + Iapp2[i-1] + sigma_2*n[i-1]) * dt/C

    if Vm1[i] > Vth:
        Vm1[i] = V_reset
        s1[i] = s1[i] + pR*D1[i]*(1-s1[i])
        D1[i] = D1[i] * (1-pR)
        states[i] = 2
        spike_times[i] = t[i]
    if Vm2[i] > Vth:
        Vm2[i] = V_reset
        s2[i] = s2[i] + pR*D2[i]*(1-s2[i])
        D2[i] = D2[i] * (1-pR)
        states[i] = 1
        spike_times[i] = t[i]

spike_times = np.where(spike_times != 0)[0]
diff_spike_times = np.diff(spike_times)

plt.figure(10)
plt.subplot(2,1,1)
plt.plot(t, Vm1, label='Vm', color='orange')
plt.ylabel('Membrane Potential (V)')
plt.title('Part vi: Vm1 and Vm2')
plt.legend()
plt.subplot(2,1,2)  
plt.plot(t, Vm2, label='Vm', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.legend()
plt.figure(10).set_size_inches(12, 8)
plt.savefig('part_vi_Vm.png')

plt.figure(11)
plt.subplot(3,1,1)
plt.plot(t, states, label='states', color = 'purple')
plt.title('States: sigma = 5 pA')
plt.ylabel('States')
plt.legend()
plt.subplot(3,1,2)
plt.hist(diff_spike_times[0::2], bins=100, color='purple')
plt.title('Histogram of one state switching times (even entries)')
plt.ylabel('Frequency')
plt.subplot(3,1,3)
plt.hist(diff_spike_times[1::2], bins=100, color='purple')
plt.title('Histogram of other state switching times (odd entries)')
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.figure(11).set_size_inches(12, 10)
plt.savefig('part_vi.png')
plt.show()