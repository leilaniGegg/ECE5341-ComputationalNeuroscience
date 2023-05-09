import math
import numpy as np
import matplotlib.pyplot as plt

EL = -75E-3
Vth = -50E-3
V_reset = -80E-3
Rm = 100E6
Cm = 100E-12   # 100 pF
tau = Rm*Cm
Ek = -80E-3
dGsra = 1E-9
tau_sra = 200E-3

dt = .1E-3
t = np.arange(0, 1.5, dt)
Vm = np.zeros(len(t))
Vm[0] = EL
Gsra = np.zeros(len(Vm))
Gsra[0] = 0


I_app = np.zeros(len(t))
I_app[5000:10000] = 5E-10


## Q1 part a

for n in range(1, len(t)):
    Gsra[n] = Gsra[n-1] - Gsra[n-1]*dt/tau_sra
    Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[n]*dt/Cm + Gsra[n]*(Ek - Vm[n-1])*dt/Cm    

    if Vm[n] > Vth:
        Vm[n] = V_reset
        Gsra[n] += dGsra
        

# plot Vm as a function of time
plt.figure(1)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.title('Membrane Potential vs Time')
plt.savefig('Q1a_Vm.png')

#plot Gsra as a function of time
plt.figure(2)
plt.plot(t, Gsra)
plt.xlabel('Time (s)')
plt.ylabel('Gsra (S)')
plt.title('Gsra vs Time')
plt.savefig('Q1a_Gsra.png')

#plot I_app as a function of time
plt.figure(3)
plt.plot(t, I_app)
plt.xlabel('Time (s)')
plt.ylabel('I_app (A)')
plt.title('I_app vs Time')
plt.savefig('Q1a_Iapp.png')

# Q1 part b
I_app = np.arange(1E-10, 5E-10, 2E-11) #100-500 pA, step 20 pA
t = np.arange(0, 5, dt)
Vm = np.zeros(len(t))
Vm[0] = EL
Gsra = np.zeros(len(Vm))
Gsra[0] = 0

firing_rate = np.zeros(len(I_app))
initial_isi = np.zeros(len(I_app))

for i in range(0, len(I_app)):
    num_spikes = 0
    for n in range(1, len(t)):
        Gsra[n] = Gsra[n-1] - Gsra[n-1]*dt/tau_sra
        Vm[n] = Vm[n-1] + (EL - Vm[n-1])*dt/tau + I_app[i]*dt/Cm + Gsra[n]*(Ek - Vm[n-1])*dt/Cm    

        if Vm[n] > Vth:
            #calculate first interspike interval for each I_app
            if num_spikes == 0:
                initial_isi[i] = (n*dt)
            if num_spikes == 1:
                initial_isi[i] = 1/(n*dt - initial_isi[i])
            Vm[n] = V_reset
            num_spikes += 1
            Gsra[n] += dGsra

    firing_rate[i] = num_spikes/5   # inverse of ISI is firing rate, so just plot firing rate


plt.figure(4)
plt.plot(I_app, firing_rate)  #1/ISI
plt.plot(I_app, initial_isi)
plt.legend(['1/(ISI)', 'Initial ISI'])
plt.xlabel('I_app (A)')
plt.ylabel('Spike Rate (Hz)')
plt.title('Spike Rate vs I_app')
plt.savefig('Q1b.png')


## Q2 part a
dth = 2E-3
GL = 1E-8  # 10 nS
a = 2E-9 # 2 nS
b = 2E-11 # 2 pA

t = np.arange(0, 1.5, dt)
Vm = np.zeros(len(t))
Vm[0] = EL
Isra = np.zeros(len(Vm))
Isra[0] = 0

I_app = np.zeros(len(t))
I_app[5000:10000] = 5E-10

for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth)/dth)) - Isra[n-1]*dt/Cm + I_app[n]*dt/Cm
    Isra[n] = a*dt/tau_sra*(Vm[n-1] - EL) - Isra[n-1]*dt/tau_sra + Isra[n-1]
    
    if Vm[n] > Vth:
        Vm[n] = V_reset
        Isra[n] += b

#plot Vm as a function of time
plt.figure(5)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.legend(['Vm', 'Isra'])
plt.title('Membrane Potential vs Time')
plt.savefig('Q2a_Vm.png')

# plot Isra as a function of time
plt.figure(6)
plt.plot(t, Isra)
plt.xlabel('Time (s)')
plt.ylabel('Isra (A)')
plt.title('Isra vs Time')
plt.savefig('Q2a_Isra.png')

## Q2 part b
I_app = np.arange(1E-10, 5E-10, 2E-11) #100-500 pA, step 20 pA
t = np.arange(0, 5, dt)
Vm = np.zeros(len(t))
Vm[0] = EL
Isra = np.zeros(len(Vm))
Isra[0] = 0

firing_rate_q2 = np.zeros(len(I_app))
initial_isi_q2 = np.zeros(len(I_app))

for i in range(0, len(I_app)):
    num_spikes = 0
    for n in range(1, len(t)):
        Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth)/dth)) - Isra[n-1]*dt/Cm + I_app[i]*dt/Cm
        Isra[n] = a*dt/tau_sra*(Vm[n-1] - EL) - Isra[n-1]*dt/tau_sra + Isra[n-1]  

        if Vm[n] > Vth:
            #calculate first interspike interval for each I_app
            if num_spikes == 0:
                initial_isi_q2[i] = (n*dt)
            if num_spikes == 1:
                initial_isi_q2[i] = 1/(n*dt - initial_isi_q2[i])
            Vm[n] = V_reset
            num_spikes += 1
            Isra[n] += b

    firing_rate_q2[i] = num_spikes/5   # inverse of ISI is firing rate, so just plot firing rate

plt.figure(7)
plt.plot(I_app, firing_rate_q2)  #1/ISI
plt.plot(I_app, initial_isi_q2)
plt.legend(['1/(ISI)', 'Initial ISI'])
plt.xlabel('I_app (A)')
plt.ylabel('Spike Rate (Hz)')
plt.title('Spike Rate vs I_app')
plt.savefig('Q2b.png')
plt.show()



## Q3 Challenge
I_app = np.arange(0, 10E-10, 1E-11) #0-100 pA, step 10 pA
Vth = np.zeros(len(Vm))
Vth[0] = -50E-3 # -50 mV
tau_Vth = 2E-3 #2 ms
Vth_max = 200E-3 #200 mV

Gsra = np.zeros(len(t))
Gsra[0] = 0
tau_sra = .5E-3 #500 us
firing_rate_q3 = np.zeros(len(I_app))
mean_membrane_potential = np.zeros(len(I_app))

for i in range(0, len(I_app)):
    num_spikes = 0
    for n in range(1, len(t)):
        Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth[n])/dth)) - Gsra[n-1]*(Ek-Vm[n])+ I_app[i]*dt/Cm
        Vth[n] = Vth[n-1] + ((Vth[0] - Vth[n-1])/tau_Vth) * dt
        Gss = a*np.tan(0.25*math.pi*min((Vm[n]-Ek)/(Vth[0]+Vth_max-Ek), 1.999))
        Gsra[n] = Gsra[n-1] + dt*(Gss - Gsra[n-1])/tau_sra 

        if Vm[n] > Vth[n]:
            #calculate first interspike interval for each I_app
            Vm[n] = V_reset
            num_spikes += 1
            Gsra[n] += b
            Vth[n] = Vth_max

    firing_rate_q3[i] = num_spikes/5   # inverse of ISI is firing rate, so just plot firing rate
    mean_membrane_potential[i] = np.average(Vm)

plt.figure(8)
plt.plot(I_app, firing_rate_q3)
plt.xlabel('I_app (A)')
plt.ylabel('Spike Rate (Hz)')
plt.title('Spike Rate vs I_app')

plt.figure(9)
plt.plot(t[0:2000], Vm[0:2000])
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.title('Membrane Potential vs Time')

plt.figure(10)
plt.plot(I_app, mean_membrane_potential)
plt.xlabel('I_app (A)')
plt.ylabel('Mean Membrane Potential (V)')
plt.title('Mean Membrane Potential vs I_app')

plt.figure(11)
plt.plot(firing_rate_q3, mean_membrane_potential)
plt.xlabel('Spike Rate (Hz)')
plt.ylabel('Mean Membrane Potential (V)')
plt.title('Mean Membrane Potential vs Spike Rate')
plt.show()