import numpy as np
import matplotlib.pyplot as plt

from PIR import PIR


Iapp_base = np.arange(-200E-12, 200E-12, 20E-12)
step_sizes = np.arange(0, 100E-12, 20E-12)

num_spikes = np.zeros((len(Iapp_base), len(step_sizes)))
min_ISI = np.zeros((len(Iapp_base), len(step_sizes)))

global Vm, t

for Iapp in Iapp_base:
    for step in step_sizes:
        Vm, t, temp = PIR(Iapp, step)
        spikes = np.where(np.diff(np.sign(Vm)) == -2)[0]   #where Vm crosses 0 
        if len(spikes) > 1:
            print('Iapp = ' + str(Iapp) + 'pA, step = ' + str(step) + 'pA')
            print('Number of spikes = ' + str(len(spikes)))
            print('Minimum interspike interval = ' + str(np.min(np.diff(spikes))))
            num_spikes[Iapp_base == Iapp, step_sizes == step] = len(spikes)
            min_ISI[Iapp_base == Iapp, step_sizes == step] = np.min(np.diff(spikes)) * 0.01E-3



#use aspect='auto' since 'equal' will distort the image
plt.matshow(num_spikes, cmap='plasma', interpolation='nearest', aspect='auto', origin='lower', extent=[-200, 200, 0, 100], fignum=1)
plt.colorbar(label='Number of spikes')
plt.xlabel('Baseline current (pA)')
plt.ylabel('Current step (pA)')
plt.title('Number of spikes')
plt.savefig('Num_Spikes.png', dpi=300)

plt.matshow(min_ISI, cmap= 'plasma', interpolation='nearest', aspect='auto', origin='lower', extent=[-200, 200, 0, 100], fignum=2)
plt.colorbar(label='Minimum ISI (ms)')
plt.xlabel('Baseline current (pA)')
plt.ylabel('Current step (pA)')
plt.title('Minimum ISI')
plt.savefig('Min_ISI.png', dpi=300)
plt.show()

#Different types of behavior
Vm, t, Iapp = PIR(0, 100E-12)

plt.figure(3)
plt.plot(t, Iapp)
plt.xlabel('Time (ms)')
plt.ylabel('Current (pA)')
plt.title('Current vs time: A')

plt.figure(4)
plt.plot(t, Vm)
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.title('Membrane potential vs time: A')


Vm, t, Iapp = PIR(-200E-12, 400E-12)

plt.figure(5)
plt.plot(t, Iapp)
plt.xlabel('Time (ms)')
plt.ylabel('Current (pA)')
plt.title('Current vs time: B')

plt.figure(6)
plt.plot(t, Vm)
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.title('Membrane potential vs time B')


Vm, t, Iapp = PIR(100E-12, 50E-12)

plt.figure(7)
plt.plot(t, Iapp)
plt.xlabel('Time (ms)')
plt.ylabel('Current (pA)')
plt.title('Current vs time: C')

plt.figure(8)
plt.plot(t, Vm)
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.title('Membrane potential vs time: C')
plt.show()