import numpy as np
import matplotlib.pyplot as plt


#########
##Functions
#########

def expandbin(initial_V, bin_width_old, bin_width_new):
    new_V_size = int(len(initial_V)*bin_width_old/bin_width_new)
    new_V = np.zeros(new_V_size)

    for i in range(0, new_V_size):
        new_V[i] = np.mean(initial_V[i*int(bin_width_new/bin_width_old):(i+1)*int(bin_width_new/bin_width_old)])

    return new_V

def STA(I_app, spikes, dt, tminus=75E-3, tplus=25E-3):
    
    nminus = int(tminus/dt) #number of time bins before a spike
    nplus = int(tplus/dt)   #number of time bins after a spike
    tcorr = np.arange(-nminus*dt, nplus*dt, dt)
    STA = np.zeros(len(tcorr))

    #create a vector of the time bins that contain spikes
    spike_times = np.where(spikes == 1)[0]
    for i in range(0, len(spike_times)):
        begin = 0 if spike_times[i] - nminus < 0 else spike_times[i] - nminus
        #define the window
        end = len(spike_times) if spike_times[i] + nplus > len(I_app) else spike_times[i] + nplus
        #append to the STA vector the values of I_app[begin:end] but keep STA the same size
        for n in range (begin, end):
            STA[n - spike_times[i] + nminus] += I_app[n]/len(spikes)

    STA = STA/len(spike_times)
    return STA, tcorr


#########
##Parameters
#########
EL = -60E-3 #-60 mV
Vth = -50E-3 #-50 mV
Vreset = -80E-3 #-80 mV
Vmax = 40E-3 #40 mV
dth = 2E-3 #2 mV
GL = 8E-9 #8 nS
Cm = 100E-12 #100 pF
a = 10E-9 #10 nS
b = 5E-10 # 0.5 nA
tau_sra = 50E-3 #50 ms


I_app_temp = np.random.uniform(-5E-10, 5E-10, 40000)    #40000*5ms random values between -5 and 5 nA

dt = 2E-5 #0.02 ms
t = np.arange(0,200,dt) #200 seconds
I_app = np.repeat(I_app_temp, 250) #250 timestep (250*0.02) = 5ms

Vm = np.zeros(len(t))
Vm[0] = EL

Isra = np.zeros(len(Vm))
Isra[0] = 0

spikes = np.zeros(len(t))

for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth)/dth)) - Isra[n-1]*dt/Cm + I_app[n]*dt/Cm
    Isra[n] = a*dt/tau_sra*(Vm[n-1] - EL) - Isra[n-1]*dt/tau_sra + Isra[n-1]
    if Vm[n] > Vmax:
        Vm[n] = Vreset
        Isra[n] += b
        spikes[n] = 1

# Plot Vm as a function of time
plt.figure(1)
plt.plot(t, Vm)
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.title('Membrane Potential as a function of time')
plt.savefig('T31_Vm.png')



new_dt = 1E-3 #1 ms
t = np.arange(0,200,new_dt) #200 seconds

#downsampled stimulus
DS_I_app = expandbin(I_app, dt, new_dt)

#downsampled spikes
DS_spikes = expandbin(spikes, dt, new_dt)
#change all the reduced spike values back to 1.
DS_spikes[DS_spikes != 0] = 1

sta, tcorr = STA(DS_I_app, DS_spikes, new_dt)

#Plot I_app as a function of time
plt.figure(2)
plt.plot(t, DS_I_app)
plt.xlabel('Time (s)')
plt.ylabel('Applied Current (A)')
plt.title('Applied Current as a function of time')
plt.savefig('T31_I_app.png')

plt.figure(3)
plt.plot(tcorr, sta)
plt.xlabel('Time Window (s)')
plt.ylabel('Spike Triggered Average')
plt.title('Spike Triggered Average')
plt.savefig('T31_STA.png')
plt.show()