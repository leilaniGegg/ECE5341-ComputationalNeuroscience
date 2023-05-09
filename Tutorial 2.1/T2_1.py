import numpy as np
import matplotlib.pyplot as plt
import math


#Tutorial 2.1
#initialize constants and Vm
EL = -70E-3 #leak potential
print(EL)
Rm = 5E6 #resistance
Cm = 2E-9 #membrane capacitance
tau = Rm*Cm
Vth = -50E-3 # spike threshold for neuron
V_reset = -65E-3 # reset potential
GL = Cm/tau   #leak conductance

dt = 0.1E-3
t = np.arange(0, 2, dt) #time vector
Vm = np.zeros(len(t)) # vector of membrane potential
Vm[0] = EL
I = GL*(Vth-EL)
I_0 = I; 
I_app = I_0*np.ones(len(t)) #vector of applied current

#Part 1a
for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[n-1]*dt/Cm
    if Vm[n] > Vth:
        Vm[n] = V_reset

plt.figure(1)
plt.plot(t,Vm)
plt.xlabel('Time')
plt.ylabel('Voltage (mV)')
plt.title('LIF Model with minimum applied current I_0=' + str(I_0))
plt.savefig('Part_1a', 'jpg')

#Part 1b
#I_0 less than minimum applied current
I_0 = I - .4E-9
I_app = I_0*np.ones(len(t))
for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[n-1]*dt/Cm
    if Vm[n] > Vth:
        Vm[n] = V_reset

plt.figure(2)
plt.plot(t[1:300],Vm[1:300])
plt.xlabel('Time')
plt.ylabel('Voltage (mV)')
plt.title('LIF Model with I_0=' + str(I_0))
plt.savefig('Part_1b_less_than', 'jpg')

#The minimum applied current needed for the neuron to produce spikes
#is 4.0e-9

#I_0 greater than minimum applied current
I_0 = I + .4E-9
I_app = I_0*np.ones(len(t))
for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[n-1]*dt/Cm
    if Vm[n] > Vth:
        Vm[n] = V_reset
plt.figure(3)
plt.plot(t[1:600],Vm[1:600])
plt.xlabel('Time')
plt.ylabel('Voltage (mV)')
plt.title('LIF Model with I_0=' + str(I_0))
plt.savefig('Part_1b_greater_than', 'jpg')

#Part 1c
I_app_currents = [I+.2E-9, I+.3E-9, I+.4E-9, I+.5E-9, I+.6E-9, I+.7E-9, I+.8E-9, I+.9E-9, I+1E-9, I+1.1E-9]
firing_rates = np.zeros(len(I_app_currents))
theoretical_FR = np.zeros(len(I_app_currents))
for i in range(0, len(I_app_currents)):
    #count spikes and convert to a firing rate
    #firing rate = spike coumt / time interval
    num_spikes = 0
    I_0 = I_app_currents[i]
    I_app = I_0*np.ones(len(t))
    for n in range(1, len(t)):
        Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[n-1]*dt/Cm
        if Vm[n] > Vth:
            Vm[n] = V_reset
            num_spikes = num_spikes + 1

    firing_rates[i] = num_spikes/2
    #theoretical calculation
    if((I_app_currents[i]*Rm + EL - V_reset > 0) & (I_app_currents[i]*Rm + EL - Vth > 0)):
        theoretical_FR[i] = 1/(tau*np.log(I_app_currents[i]*Rm + EL - V_reset) - tau*np.log(I_app_currents[i]*Rm + EL - Vth))
 
#Part 1d
plt.figure(4)
plt.plot(I_app_currents, firing_rates,'+', linestyle='solid')
plt.plot(I_app_currents, theoretical_FR, '+', linestyle='solid')
plt.legend(['Trials', 'Theoretical'])
plt.xlabel('I_{app} Values')
plt.ylabel('Firing Rates (Hz)')
plt.title('Firing Rates: Trials vs. Theoretical')
plt.savefig('Firing_Rates_Trials_vs_Theoretical', 'jpg')

#Part 2a & 2b
plt.figure(5)
#Orignal firing_rates from part c
plt.plot(I_app_currents, firing_rates, '*', linestyle='solid')
plt.xlabel('I_{app} Values')
plt.ylabel('Firing Rates (Hz)')
plt.title('Firing Rates with Noise')
plt.savefig('Firing_Rates_with_Noise', 'jpg')

#sigma1
firing_rates1 = np.zeros(len(I_app_currents))
sigma1 = 0.01
for i in range(0, len(I_app_currents)):
    num_spikes = 0
    I_0 = I_app_currents[i]
    I_app = I_0*np.ones(len(t))
    for n in range(1, len(t)):
        noise = np.random.randn(1)*sigma1*math.sqrt(dt)
        Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[n-1]*dt/Cm + noise
        if Vm[n] > Vth:
            Vm[n] = V_reset
            num_spikes = num_spikes + 1
    firing_rates1[i] = num_spikes/2

plt.plot(I_app_currents, firing_rates1, '*', linestyle='solid')

#sigma2
firing_rates2 = np.zeros(len(I_app_currents))
sigma2 = 0.05
for i in range(0, len(I_app_currents)):
    num_spikes = 0
    I_0 = I_app_currents[i]
    I_app = I_0*np.ones(len(t))
    for n in range(1, len(t)):
        noise = np.random.randn(1)*sigma2*np.sqrt(dt)
        Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[n-1]*dt/Cm + noise
        if Vm[n] > Vth:
            Vm[n] = V_reset
            num_spikes = num_spikes + 1
    firing_rates2[i] = num_spikes/2

plt.plot(I_app_currents, firing_rates2, '*', linestyle='solid')

#sigma3
firing_rates3 = np.zeros(len(I_app_currents))
sigma3 = 0.1
for i in range(0, len(I_app_currents)):
    num_spikes = 0
    I_0 = I_app_currents[i]
    I_app = I_0*np.ones(len(t))
    for n in range(1, len(t)):
        noise = np.random.randn(1)*sigma3*np.sqrt(dt)
        Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[n-1]*dt/Cm + noise
        if Vm[n] > Vth:
            Vm[n] = V_reset
            num_spikes = num_spikes + 1
     
    firing_rates3[i] = num_spikes/2

plt.plot(I_app_currents, firing_rates3,'*', linestyle='solid')
plt.legend(['sigma = 0', 'sigma=0.04', 'sigma=0.07', 'sigma=0.12'])
plt.xlabel('I_{app} Values')
plt.ylabel('Firing Rates (Hz)')
plt.title('Firing Rates with Noise')
plt.savefig('Firing_Rates_with_Noise', 'jpg')
plt.show()

#As sigma is increased and more noise is added, the firing rate increases
#on average but also becomes bendy and chaotic after a certain point.

#Part c
#if dt is scaled by .10 at the beginning before any calculations are made
#then the results are not significantly different.