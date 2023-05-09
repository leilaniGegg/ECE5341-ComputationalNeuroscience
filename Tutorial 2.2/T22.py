import numpy as np
import matplotlib.pyplot as plt
import math

EL = -70E-3 #-70 mV
Rm = 100E6
Cm = .1E-9
tau = Rm*Cm
dt = 0.1E-3
t = np.arange(0, 2, dt)
V_peak = 50E-3 #50 mV

I_app = np.arange(1E-10, 6E-10, 2.5E-11) #100-600 pA, step 25 pA

Vm = np.zeros(len(t))
Vm[0] = EL

global firing_rates_q1
global firing_rates_q2
global firing_rates_q3

global membrane_potential_q1
global membrane_potential_q2
global membrane_potential_q3

#############
## Question 1
#############
# 1.1) Plot mean firing rate as a function of input current
# 1.2) Plot mean membrane potential as a function of input current
# 1.3) Plot mean membrane potential as a function of firing rate
Vth = -50E-3 # -50 mV
V_reset = -65E-3 # -65 mV
tau_ref = 2.5E-3 #2.5 ms
firing_rates_q1 = np.zeros(len(I_app))
membrane_potential_q1 = np.zeros(len(I_app))
for i in range(0, len(I_app)):
    numspike = 0
    spike_time = -tau_ref-.1 #some value less than -tau_ref
    for n in range(1, len(t)):
        if(t[n] > spike_time + tau_ref):
            Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[i]*dt/Cm
            if Vm[n] >= Vth:
                Vm[n] = V_reset
                numspike += 1
                spike_time = t[n]
        else:
            Vm[n] = V_reset
                
    membrane_potential_q1[i] = np.average(Vm)
    firing_rates_q1[i] = numspike/2
   


I_app_2 = [2.2E-10, 6E-10]  #100 pA, 600 pA
for i in range(0, len(I_app_2)):
    spike_time = -tau_ref-.1 #some value less than -tau_ref
    for n in range(1, len(t)): 
        if(t[n] > spike_time + tau_ref):
            Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app_2[i]*dt/Cm
            if Vm[n] >= Vth:
                Vm[n] = V_reset
                Vm[n-1] = V_peak
                spike_time = t[n]
        else:
            Vm[n] = V_reset
    plt.figure(1)
    plt.plot(t[0:1000], Vm[0:1000])
    plt.legend(['220 pA', '600 pA'])
    plt.xlabel('Time (s)')
    plt.ylabel('Membrane Potential (V)')
    plt.title('Q1 Membrane Potential as a Function of Time')
    plt.savefig('Q1 Membrane Potential as a Function of Time.png')
    

               

#############
## Question 2
#############
# 2.1) Plot mean firing rate as a function of input current
# 2.2) Plot mean membrane potential as a function of input current
# 2.3) Plot mean membrane potential as a function of firing rate
Vth = np.zeros(len(Vm))
Vth[0] = -50E-3 # -50 mV
tau_Vth = 1E-3 #1 ms
Vth_max = 200E-3 #200 mV

firing_rates_q2 = np.zeros(len(I_app))
membrane_potential_q2 = np.zeros(len(I_app))
for i in range(0, len(I_app)):
    numspike = 0
    for n in range(1, len(t)): 
        Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[i]*dt/Cm
        Vth[n] = Vth[n-1] + ((Vth[0] - Vth[n-1])/tau_Vth) * dt
        if Vm[n] > Vth[n]:
            Vm[n] = V_reset
            Vth[n] = Vth_max
            numspike += 1
                
    firing_rates_q2[i] = numspike/2
    membrane_potential_q2[i] = np.average(Vm)
    
for i in range(0, len(I_app_2)):
    for n in range(1, len(t)): 
        Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app_2[i]*dt/Cm
        Vth[n] = Vth[n-1] + ((Vth[0] - Vth[n-1])/tau_Vth) * dt
        if Vm[n] > Vth[n]:
            Vm[n] = V_reset
            Vm[n-1] = V_peak
            Vth[n] = Vth_max

   
    plt.figure(2)
    plt.plot(t[0:1000], Vm[0:1000])
    plt.legend(['220 pA', '600 pA'])
    plt.xlabel('Time (s)')
    plt.ylabel('Membrane Potential (V)')
    plt.title('Q2 Membrane Potential as a Function of Time')
    plt.savefig('Q2 Membrane Potential as a Function of Time.png')


plt.figure(4)
plt.plot(I_app, firing_rates_q1)
plt.plot(I_app, firing_rates_q2)
plt.legend(['Q1', 'Q2'])
plt.xlabel('Input Current (A)')
plt.ylabel('Mean Firing Rate (Hz)')
plt.title('Mean Firing Rate as a Function of Input Current')
plt.savefig('Mean Firing Rate as a Function of Input Current.png')

plt.figure(5)
plt.plot(I_app, membrane_potential_q1)
plt.plot(I_app, membrane_potential_q2)
plt.legend(['Q1', 'Q2'])
plt.xlabel('Input Current (A)')
plt.ylabel('Membrane Potential (V)')
plt.title('Mean Membrane Potential as a Function of Input Current')
plt.savefig('Mean Membrane Potential as a Function of Input Current.png')

plt.figure(6)
plt.plot(firing_rates_q1, membrane_potential_q1)
plt.plot(firing_rates_q2, membrane_potential_q2)
plt.legend(['Q1', 'Q2'])
plt.xlabel('Mean Firing Rate (Hz)')
plt.ylabel('Membrane Potential (V)')
plt.title('Mean Membrane Potential as a Function of Firing Rate')
plt.savefig('Mean Membrane Potential as a Function of Firing Rate.png')
plt.show()