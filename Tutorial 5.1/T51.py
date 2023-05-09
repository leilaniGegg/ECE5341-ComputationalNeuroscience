import numpy as np
import matplotlib.pyplot as plt

#Part a
dt = 0.1E-3
t = np.arange(0, 4, dt)

#Part b
presynaptic_fr = np.zeros(len(t))

presynaptic_fr[0:10000] = 20
presynaptic_fr[10000:20000] = 100
presynaptic_fr[20000:30000] = 10
presynaptic_fr[30000:40000] = 50

#Part c
presynaptic_spikes = np.random.rand(len(t)) < presynaptic_fr*dt


#Part d: Synaptic Conductance
dG = 1E-9
G_syn = np.zeros(len(t))
for i in range(1, len(t)):
    if presynaptic_spikes[i]:
        G_syn[i] = G_syn[i-1] + dG
    else:
        G_syn[i] = G_syn[i-1] + dt*(-G_syn[i-1]/.1)


#Part e
p_0 = 0.5
D = np.zeros(len(t))
D[0] = 1
for i in range(1, len(t)):
    if presynaptic_spikes[i]:
        D[i] = D[i-1] - p_0*dt     
    else:
        D[i] = D[i-1] + dt*((1-D[i-1]/.25))

plt.figure(2)
plt.plot(t, D)
plt.xlabel('Time (s)')
plt.ylabel('Synaptic depression')
plt.title('Synaptic depression over time')
plt.gcf().set_size_inches(10, 5)
plt.savefig('Synaptic depression.png')

#Part f: Synaptic Conductance with Depression
dG_max = 5E-9
p_0 = 0.5
G_syn_2 = np.zeros(len(t))
for i in range(1, len(t)):
    if presynaptic_spikes[i]:
        dG = dG_max * p_0 * D[i-1]
        G_syn_2[i] = G_syn_2[i-1] + dG
    else:
        G_syn_2[i] = G_syn_2[i-1] + dt*(-G_syn_2[i-1]/.1)



#Part g
p_0 = 0.2
F = np.zeros(len(t))
F[0] = 1
f_fac = 0.25  
F_max = 1/p_0
for i in range(1, len(t)):
    if presynaptic_spikes[i]:
        F[i] = F[i-1] + f_fac*(F_max-F[i-1])
    else:
        F[i] = F[i-1] + dt*((1-F[i-1]/0.25))
    
plt.figure(4)
plt.plot(t, F)
plt.xlabel('Time (s)')
plt.ylabel('Synaptic facilitation')
plt.title('Synaptic facilitation over time')
plt.gcf().set_size_inches(10, 5)
plt.savefig('Synaptic facilitation.png')


#Part h
D_2 = np.zeros(len(t))
D_2[0] = 1
p_0 = 0.2
for i in range(1, len(t)):
    if presynaptic_spikes[i]:
        D_2[i] = D_2[i-1] + p_0*D[i-1]*F[i-1]           
    else:
        D_2[i] = D_2[i-1] + dt*((1-D_2[i-1]/.25))

plt.figure(5)
plt.plot(t, D_2)
plt.xlabel('Time (s)')
plt.ylabel('Synaptic depression')
plt.title('Synaptic depression over time Part h')
plt.gcf().set_size_inches(10, 5)
plt.savefig('Synaptic depression part h.png')


#Part i: Synaptic Conductance with Depression and Facilitation
G_syn_3 = np.zeros(len(t))
dG_max = 4E-9
p_0 = 0.2
for i in range(1, len(t)):
    if presynaptic_spikes[i]:
        dG = dG_max * p_0 * D_2[i-1] * F[i-1]
        G_syn_3[i] = G_syn_3[i-1] + dG
    else:
        G_syn_3[i] = G_syn_3[i-1] + dt*(-G_syn_3[i-1]/.1)


plt.figure(8)
plt.subplot(4, 1, 1)
plt.plot(t, G_syn)
plt.ylabel('Synaptic conductance (S)')
plt.title('Synaptic conductance (d)')
plt.subplot(4, 1, 2)
plt.plot(t, G_syn_2)
plt.ylabel('Synaptic conductance (S)')
plt.title('Synaptic conductance w/ depression (f)')
plt.subplot(4, 1, 3)
plt.plot(t, G_syn_3)
plt.ylabel('Synaptic conductance (S)')
plt.title('Synaptic conductance w/ depression and facilitation (i)')
plt.subplot(4, 1, 4)
plt.plot(t, presynaptic_spikes)
plt.xlabel('Time (s)')
plt.ylabel('Presynaptic spikes')
plt.title('Presynaptic spikes')
plt.gcf().set_size_inches(10, 15)
plt.savefig('Synaptic conductances(d,f,i).png')
plt.show()