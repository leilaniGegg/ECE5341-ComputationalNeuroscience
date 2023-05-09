import numpy as np
import matplotlib.pyplot as plt

G_leak = 30E-9 # 30 nS      # leak conductance
G_na_max = 12E-6 # 12 uS    # max Na conductance
G_k_max = 3.6E-6 # 3.6 uS   # max K conductance
E_na = 45E-3 # 45 mV        # Na reversal potential
E_k = -82E-3 # -82 mV       # K reversal potential
E_leak = -60E-3 # -60 mV    # leak reversal potential
Cm = 100E-12 # 100 pF       # membrane capacitance

dt = 0.01E-7 # 0.01 ms     
t = np.arange(0, 0.35, dt)  

Vm = np.zeros(len(t))       # membrane potential
Vm[0] = E_leak
m = np.zeros(len(t))        # Na activation
m[0] = 0
h = np.zeros(len(t))        # Na inactivation
h[0] = 0
n = np.zeros(len(t))        # K activation
n[0] = 0
Iapp = 0

########
##Part a
########
for i in range(1, len(t)):
    # membrane potential
    Vm[i] = Vm[i-1] + dt*(G_leak*(E_leak - Vm[i-1]) + G_na_max*m[i-1]**3*h[i-1]*(E_na - Vm[i-1]) + G_k_max*n[i-1]**4*(E_k - Vm[i-1]) + Iapp) / Cm
    # Na activation
    alpha_m = (10**5*(-Vm[i]-0.045)) / (np.exp(100*(-Vm[i] - 0.045)) - 1)
    beta_m = 4*10**3*np.exp((-Vm[i]-0.070)/0.018)
    m[i] = m[i-1] + dt*(alpha_m*(1-m[i-1]) - beta_m*m[i-1])
    # Na inactivation
    alpha_h = 70*np.exp(50*(-Vm[i]-0.070))
    beta_h = 10**3 / (1+np.exp(100*(-Vm[i] - 0.040)))
    h[i] = h[i-1] + dt*(alpha_h*(1-h[i-1]) - beta_h*h[i-1])
    # K activation
    #deal with the case where Vm[i-1] = -0.060
    if Vm[i] == -0.060:
        Vm[i] = -0.06001
    alpha_n = (10**4*(-Vm[i]-0.060)) / (np.exp(100*(-Vm[i] - 0.060)) - 1)
    beta_n = 125*np.exp((-Vm[i]-0.070)/0.080)
    n[i] = n[i-1] + dt*(alpha_n*(1-n[i-1]) - beta_n*n[i-1])

plt.figure(1)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title('Membrane Potential vs Time Part a')
plt.savefig('Part a.png')

########
##Part b
########
Iapp = np.zeros(len(t))
Iapp[10000:20000] = 0.22E-9

for i in range(1, len(t)):
    # membrane potential
    Vm[i] = Vm[i-1] + dt*(G_leak*(E_leak - Vm[i-1]) + G_na_max*m[i-1]**3*h[i-1]*(E_na - Vm[i-1]) + G_k_max*n[i-1]**4*(E_k - Vm[i-1]) + Iapp[i-1]) / Cm
    # Na activation
    alpha_m = (10**5*(-Vm[i]-0.045)) / (np.exp(100*(-Vm[i] - 0.045)) - 1)
    beta_m = 4*10**3*np.exp((-Vm[i]-0.070)/0.018)
    m[i] = m[i-1] + dt*(alpha_m*(1-m[i-1]) - beta_m*m[i-1])
    # Na inactivation
    alpha_h = 70*np.exp(50*(-Vm[i]-0.070))
    beta_h = 10**3 / (1+np.exp(100*(-Vm[i] - 0.040)))
    h[i] = h[i-1] + dt*(alpha_h*(1-h[i-1]) - beta_h*h[i-1])
    # K activation
    #deal with the case where Vm[i-1] = -0.060
    if Vm[i] == -0.060:
        Vm[i] = -0.06001
    alpha_n = (10**4*(-Vm[i]-0.060)) / (np.exp(100*(-Vm[i] - 0.060)) - 1)
    beta_n = 125*np.exp((-Vm[i]-0.070)/0.080)
    n[i] = n[i-1] + dt*(alpha_n*(1-n[i-1]) - beta_n*n[i-1])  

plt.figure(2)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title('Membrane Potential vs Time Part b')
plt.savefig('Membrane Potential Part b.png')

plt.figure(3)
plt.plot(t, Iapp)
plt.xlabel('Time (s)')
plt.ylabel('Applied current (A)')
plt.title('Applied Current vs Time Part b')
plt.savefig('Applied Current Part b.png')

#######
#Part c
#######

##First Delay
delay = 10E-3 # 10 ms

pulse_duration = 0.005  # duration of each pulse in seconds (5ms)
pulse_amplitude = .22E-9  # amplitude of each pulse in nanoampere
n_pulses = 10  # number of pulses
Iapp = np.zeros(len(t))
for i in range(0, n_pulses):
    Iapp[int((delay + i*(pulse_duration + delay))*(1/dt)):int((delay + i*(pulse_duration + delay) + pulse_duration)*(1/dt))] = pulse_amplitude


for i in range(1, len(t)):
    # membrane potential
    Vm[i] = Vm[i-1] + dt*(G_leak*(E_leak - Vm[i-1]) + G_na_max*m[i-1]**3*h[i-1]*(E_na - Vm[i-1]) + G_k_max*n[i-1]**4*(E_k - Vm[i-1]) + Iapp[i-1]) / Cm
    # Na activation
    alpha_m = (10**5*(-Vm[i]-0.045)) / (np.exp(100*(-Vm[i] - 0.045)) - 1)
    beta_m = 4*10**3*np.exp((-Vm[i]-0.070)/0.018)
    m[i] = m[i-1] + dt*(alpha_m*(1-m[i-1]) - beta_m*m[i-1])
    # Na inactivation
    alpha_h = 70*np.exp(50*(-Vm[i]-0.070))
    beta_h = 10**3 / (1+np.exp(100*(-Vm[i] - 0.040)))
    h[i] = h[i-1] + dt*(alpha_h*(1-h[i-1]) - beta_h*h[i-1])
    # K activation
    #deal with the case where Vm[i-1] = -0.060
    if Vm[i] == -0.060:
        Vm[i] = -0.06001
    alpha_n = (10**4*(-Vm[i]-0.060)) / (np.exp(100*(-Vm[i] - 0.060)) - 1)
    beta_n = 125*np.exp((-Vm[i]-0.070)/0.080)
    n[i] = n[i-1] + dt*(alpha_n*(1-n[i-1]) - beta_n*n[i-1])

plt.figure(4)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title('Membrane potential with delay of 10 ms')
plt.savefig('Membrane Potential PartC delay10.png')

plt.figure(5)
plt.plot(t, Iapp)
plt.xlabel('Time (s)')
plt.ylabel('Applied current (A)')
plt.title('Applied current with delay of 10 ms')
plt.savefig('Applied Current PartC delay10.png')

###Second Delay
delay = 20E-3 # 20 ms

pulse_duration = 0.005  # duration of each pulse in seconds (5ms)
pulse_amplitude = .22E-9  # amplitude of each pulse in nanoampere
n_pulses = 10  # number of pulses
Iapp = np.zeros(len(t))
for i in range(0, n_pulses):
    Iapp[int((delay + i*(pulse_duration + delay))*(1/dt)):int((delay + i*(pulse_duration + delay) + pulse_duration)*(1/dt))] = pulse_amplitude


for i in range(1, len(t)):
    # membrane potential
    Vm[i] = Vm[i-1] + dt*(G_leak*(E_leak - Vm[i-1]) + G_na_max*m[i-1]**3*h[i-1]*(E_na - Vm[i-1]) + G_k_max*n[i-1]**4*(E_k - Vm[i-1]) + Iapp[i-1]) / Cm
    # Na activation
    alpha_m = (10**5*(-Vm[i]-0.045)) / (np.exp(100*(-Vm[i] - 0.045)) - 1)
    beta_m = 4*10**3*np.exp((-Vm[i]-0.070)/0.018)
    m[i] = m[i-1] + dt*(alpha_m*(1-m[i-1]) - beta_m*m[i-1])
    # Na inactivation
    alpha_h = 70*np.exp(50*(-Vm[i]-0.070))
    beta_h = 10**3 / (1+np.exp(100*(-Vm[i] - 0.040)))
    h[i] = h[i-1] + dt*(alpha_h*(1-h[i-1]) - beta_h*h[i-1])
    # K activation
    #deal with the case where Vm[i-1] = -0.060
    if Vm[i] == -0.060:
        Vm[i] = -0.06001
    alpha_n = (10**4*(-Vm[i]-0.060)) / (np.exp(100*(-Vm[i] - 0.060)) - 1)
    beta_n = 125*np.exp((-Vm[i]-0.070)/0.080)
    n[i] = n[i-1] + dt*(alpha_n*(1-n[i-1]) - beta_n*n[i-1])

plt.figure(6)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title('Membrane potential with delay of 20 ms')
plt.savefig('Membrane Potential PartC delay20.png')

plt.figure(7)
plt.plot(t, Iapp)
plt.xlabel('Time (s)')
plt.ylabel('Applied current (A)')
plt.title('Applied current with delay of 20 ms')
plt.savefig('Applied Current PartC delay20.png')


##########
###Part d
##########
Iapp = np.full(len(t), 0.6E-9)
Vm[0] = -0.065
m[0] = 0.05
h[0] = 0.5
n[0] = 0.35

delay = 20E-3 # 20 ms
pulse_duration = 0.005  # duration of each pulse in seconds (5ms)
pulse_amplitude = 0  # bring current to 0
n_pulses = 10  # number of pulses
for i in range(0, n_pulses):
    Iapp[int((delay + i*(pulse_duration + delay))*100000):int((delay + i*(pulse_duration + delay) + pulse_duration)*100000)] = pulse_amplitude

for i in range(1, len(t)):
    # membrane potential
    Vm[i] = Vm[i-1] + dt*(G_leak*(E_leak - Vm[i-1]) + G_na_max*m[i-1]**3*h[i-1]*(E_na - Vm[i-1]) + G_k_max*n[i-1]**4*(E_k - Vm[i-1]) + Iapp[i-1]) / Cm
    # Na activation
    alpha_m = (10**5*(-Vm[i]-0.045)) / (np.exp(100*(-Vm[i] - 0.045)) - 1)
    beta_m = 4*10**3*np.exp((-Vm[i]-0.070)/0.018)
    m[i] = m[i-1] + dt*(alpha_m*(1-m[i-1]) - beta_m*m[i-1])
    # Na inactivation
    alpha_h = 70*np.exp(50*(-Vm[i]-0.070))
    beta_h = 10**3 / (1+np.exp(100*(-Vm[i] - 0.040)))
    h[i] = h[i-1] + dt*(alpha_h*(1-h[i-1]) - beta_h*h[i-1])
    # K activation
    #deal with the case where Vm[i-1] = -0.060
    if Vm[i] == -0.060:
        Vm[i] = -0.06001
    alpha_n = (10**4*(-Vm[i]-0.060)) / (np.exp(100*(-Vm[i] - 0.060)) - 1)
    beta_n = 125*np.exp((-Vm[i]-0.070)/0.080)
    n[i] = n[i-1] + dt*(alpha_n*(1-n[i-1]) - beta_n*n[i-1])

plt.figure(8)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title('Membrane potential vs Time Part d')
plt.savefig('Membrane Potential PartD.png')

plt.figure(9)
plt.plot(t, Iapp)
plt.xlabel('Time (s)')
plt.ylabel('Applied current (A)')
plt.title('Applied current vs Time Part d')
plt.savefig('Applied Current PartD.png')

plt.figure(10)
plt.plot(t, m, label='m')
plt.plot(t, h, label='h')
plt.plot(t, n, label='n')
plt.xlabel('Time (s)')
plt.ylabel('Gating variable')
plt.title('Gating variables vs Time Part d')
plt.legend()
plt.savefig('Gating Variables PartD.png')


##########
###Part e
##########
Vm[0] = -0.065
m[0] = 0.05
h[0] = 0.5
n[0] = 0.35

Iapp = np.full(len(t), 0.65E-9)
Iapp[10000:10500] = 1E-9

for i in range(1, len(t)):
    # membrane potential
    Vm[i] = Vm[i-1] + dt*(G_leak*(E_leak - Vm[i-1]) + G_na_max*m[i-1]**3*h[i-1]*(E_na - Vm[i-1]) + G_k_max*n[i-1]**4*(E_k - Vm[i-1]) + Iapp[i-1]) / Cm
    # Na activation
    alpha_m = (10**5*(-Vm[i]-0.045)) / (np.exp(100*(-Vm[i] - 0.045)) - 1)
    beta_m = 4*10**3*np.exp((-Vm[i]-0.070)/0.018)
    m[i] = m[i-1] + dt*(alpha_m*(1-m[i-1]) - beta_m*m[i-1])
    # Na inactivation
    alpha_h = 70*np.exp(50*(-Vm[i]-0.070))
    beta_h = 10**3 / (1+np.exp(100*(-Vm[i] - 0.040)))
    h[i] = h[i-1] + dt*(alpha_h*(1-h[i-1]) - beta_h*h[i-1])
    # K activation
    #deal with the case where Vm[i-1] = -0.060
    if Vm[i] == -0.060:
        Vm[i] = -0.06001
    alpha_n = (10**4*(-Vm[i]-0.060)) / (np.exp(100*(-Vm[i] - 0.060)) - 1)
    beta_n = 125*np.exp((-Vm[i]-0.070)/0.080)
    n[i] = n[i-1] + dt*(alpha_n*(1-n[i-1]) - beta_n*n[i-1])

plt.figure(11)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title('Membrane potential vs Time Part e')
plt.savefig('Membrane Potential PartE.png')

plt.figure(12)
plt.plot(t, Iapp)
plt.xlabel('Time (s)')
plt.ylabel('Applied current (A)')
plt.title('Applied current vs Time Part e')
plt.savefig('Applied Current PartE.png')

plt.figure(13)
plt.plot(t, m, label='m')
plt.plot(t, h, label='h')
plt.plot(t, n, label='n')
plt.xlabel('Time (s)')
plt.ylabel('Gating variable')
plt.title('Gating variables vs Time Part e')
plt.legend()
plt.savefig('Gating Variables PartE.png')


##########
###Part f
##########
Vm[0] = -0.065
m[0] = 0
h[0] = 0
n[0] = 0

Iapp = np.full(len(t), 0.7E-9)
Iapp[10000:10500] = 1E-9

for i in range(1, len(t)):
    # membrane potential
    Vm[i] = Vm[i-1] + dt*(G_leak*(E_leak - Vm[i-1]) + G_na_max*m[i-1]**3*h[i-1]*(E_na - Vm[i-1]) + G_k_max*n[i-1]**4*(E_k - Vm[i-1]) + Iapp[i-1]) / Cm
    # Na activation
    alpha_m = (10**5*(-Vm[i]-0.045)) / (np.exp(100*(-Vm[i] - 0.045)) - 1)
    beta_m = 4*10**3*np.exp((-Vm[i]-0.070)/0.018)
    m[i] = m[i-1] + dt*(alpha_m*(1-m[i-1]) - beta_m*m[i-1])
    # Na inactivation
    alpha_h = 70*np.exp(50*(-Vm[i]-0.070))
    beta_h = 10**3 / (1+np.exp(100*(-Vm[i] - 0.040)))
    h[i] = h[i-1] + dt*(alpha_h*(1-h[i-1]) - beta_h*h[i-1])
    # K activation
    #deal with the case where Vm[i-1] = -0.060
    if Vm[i] == -0.060:
        Vm[i] = -0.06001
    alpha_n = (10**4*(-Vm[i]-0.060)) / (np.exp(100*(-Vm[i] - 0.060)) - 1)
    beta_n = 125*np.exp((-Vm[i]-0.070)/0.080)
    n[i] = n[i-1] + dt*(alpha_n*(1-n[i-1]) - beta_n*n[i-1])

plt.figure(14)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title('Membrane potential vs Time Part f')
plt.savefig('Membrane Potential PartF.png')

plt.figure(15)
plt.plot(t, Iapp)  
plt.xlabel('Time (s)')
plt.ylabel('Applied current (A)')
plt.title('Applied current vs Time Part f')
plt.savefig('Applied Current PartF.png')

plt.figure(16)
plt.plot(t, m, label='m')
plt.plot(t, h, label='h')
plt.plot(t, n, label='n')
plt.xlabel('Time (s)')
plt.ylabel('Gating variable')
plt.title('Gating variables vs Time Part f')
plt.legend()
plt.savefig('Gating Variables PartF.png')
plt.show()