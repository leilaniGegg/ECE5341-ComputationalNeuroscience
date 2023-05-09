import numpy as np
import matplotlib.pyplot as plt

EL = -70E-3 
Vth = -50E-3
Vreset = -80E-3
Vmax = 10E-3
dth = 2E-3
GL = 10E-9
Cm = 100E-12
a = 2E-9
b = 0
tau_sra = 150E-3

dt = .01E-3
t = np.arange(0, 100, dt)
Vm = np.zeros(len(t))
Vm[0] = EL

Isra = np.zeros(len(t))
Isra[0] = 0
spikes = np.zeros(len(t))

sigma = 50E-12 #50 pA
I_app = np.zeros(len(t))
for i in range(len(t)):
    I_app[i] = np.random.normal(0, sigma/np.sqrt(dt))   #mean, std dev

num_spikes = 0


#########
##Part 1a
#########
for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth)/dth)) - Isra[n-1]*dt/Cm + I_app[n]*dt/Cm
    Isra[n] = a*dt/tau_sra*(Vm[n-1] - EL) - Isra[n-1]*dt/tau_sra + Isra[n-1]
    if Vm[n] > Vmax:
        Vm[n] = Vreset
        Isra[n] += b
        spikes[n] = t[n]
        num_spikes += 1
        
spikes = np.where(spikes != 0)
isi = np.zeros(len(spikes))
isi = np.diff(spikes)

isi = [i * dt for i in isi]
#plot histogram of ISI with 25 bins for the values of ISI
plt.figure(1)
plt.hist(isi, 25)
plt.xlabel('ISI (s)')
plt.ylabel('Count')
plt.title('Histogram of ISI for b = 0')
plt.savefig('T32_Hist_a.png')
plt.show()

CV = np.std(isi)/np.mean(isi)
print("The coefficient of variation is: ", CV)

#multiply each element of spikes by dt to get the spike times in seconds
spikes = [i * dt for i in spikes]

#1a.iii
num_spikes = np.zeros(1000)
for i in range(np.size(spikes)):
    x = int(np.fix(spikes[0][i]/.1))
    num_spikes[x] += 1

mean = np.mean(num_spikes)
variance = np.var(num_spikes)
fano_factor = variance/mean

print("The fano factor is: ", fano_factor)


window_size = np.arange(10E-3, 1, 10E-3)
fano_factor = np.zeros(len(window_size))
for i in range(len(window_size)):
    num_spikes = np.zeros(int(100/window_size[i]))
    for j in range(np.size(spikes)):
        x = int(spikes[0][j]/window_size[i])
        #handle the case where the spike time is exactly equal to the window size
        try:
            num_spikes[x] += 1
        except:
            pass
    
    mean = np.mean(num_spikes)
    variance = np.var(num_spikes)
    fano_factor[i] = variance/mean


plt.figure(2)
plt.scatter(window_size, fano_factor)
plt.xlabel('Window Size (s)')
plt.ylabel('Fano Factor Part 1a')
plt.title('Fano Factor vs. Window Size for Part 1a')
plt.savefig('T32_Fano_a.png')
plt.show()


#########
##Part 1b
#########

b = 1E-9
spikes = np.zeros(len(t))
num_spikes = 0
isi = np.zeros(len(t))

for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth)/dth)) - Isra[n-1]*dt/Cm + I_app[n]*dt/Cm
    Isra[n] = a*dt/tau_sra*(Vm[n-1] - EL) - Isra[n-1]*dt/tau_sra + Isra[n-1]
    if Vm[n] > Vmax:
        Vm[n] = Vreset
        Isra[n] += b
        spikes[n] = t[n]
        num_spikes += 1

spikes = np.where(spikes != 0)
isi = np.diff(spikes)
isi = [i * dt for i in isi]
#plot histogram of ISI with 25 bins for the values of ISI
plt.figure(3)
plt.hist(isi, 25)
plt.xlabel('ISI (s)')
plt.ylabel('Count')
plt.title('Histogram of ISI for b = 1 nA')
plt.savefig('T32_Hist_b.png')
plt.show()


CV = np.std(isi)/np.mean(isi)
print("The coefficient of variation is: ", CV)

#multiply each element of spikes by dt to get the spike times in seconds
spikes = [i * dt for i in spikes]

#1b.iii
num_spikes = np.zeros(1000)

for i in range(np.size(spikes)):
    x = int(np.fix(spikes[0][i]/.1))
    num_spikes[x] += 1

mean = np.mean(num_spikes)
variance = np.var(num_spikes)
fano_factor = variance/mean

print("The fano factor is: ", fano_factor)

window_size = np.arange(10E-3, 1, 10E-3)
fano_factor = np.zeros(len(window_size))
for i in range(len(window_size)):
    num_spikes = np.zeros(int(100/window_size[i]))
    for j in range(np.size(spikes)):
        x = int(spikes[0][j]/window_size[i])
        #handle the case where the spike time is exactly equal to the window size
        try:
            num_spikes[x] += 1
        except:
            pass
    mean = np.mean(num_spikes)
    variance = np.var(num_spikes)
    fano_factor[i] = variance/mean


plt.figure(4)
plt.scatter(window_size, fano_factor)
plt.xlabel('Window Size (s)')
plt.ylabel('Fano Factor')
plt.title('Fano Factor vs. Window Size for Part 1b')
plt.savefig('T32_Fano_b.png')
plt.show()

#########
##Part c
#########

#0 nA added to the input current
b = 0
spikes = np.zeros(len(t))
sigma = 20E-12 
I_app = np.zeros(len(t))
for i in range(len(t)):
    I_app[i] = np.random.normal(0, sigma/np.sqrt(dt))   #mean, std dev

num_spikes = 0
isi = np.zeros(len(t))

for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth)/dth)) - Isra[n-1]*dt/Cm + I_app[n]*dt/Cm
    Isra[n] = a*dt/tau_sra*(Vm[n-1] - EL) - Isra[n-1]*dt/tau_sra + Isra[n-1]
    if Vm[n] > Vmax:
        Vm[n] = Vreset
        Isra[n] += b
        spikes[n] = t[n]
        num_spikes += 1
        
spikes = np.where(spikes != 0)
isi = np.diff(spikes)
#Multiply by dt to get the ISI in seconds
isi = [i * dt for i in isi]
#plot histogram of ISI with 25 bins for the values of ISI
plt.figure(5)
plt.hist(isi, 25)
plt.xlabel('ISI (s)')
plt.ylabel('Count')
plt.title('Histogram of ISI: 0 nA added')
plt.savefig('T32_Hist_c0.png')
plt.show()


CV = np.std(isi)/np.mean(isi)
print("The coefficient of variation is: ", CV)

spikes = [i * dt for i in spikes]

#1c.iii
num_spikes = np.zeros(1000)
for i in range(np.size(spikes)):
    x = int(np.fix(spikes[0][i]/.1))
    num_spikes[x] += 1

mean = np.mean(num_spikes)
variance = np.var(num_spikes)
fano_factor = variance/mean

print("The fano factor is: ", fano_factor)


##0.1 nA added to the input current
spikes = np.zeros(len(t))
I_app = np.zeros(len(t))
for i in range(len(t)):
    I_app[i] = np.random.normal(0, sigma/np.sqrt(dt)) + .1E-9   #add .1 nA 

num_spikes = 0
isi = np.zeros(len(t))

for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth)/dth)) - Isra[n-1]*dt/Cm + I_app[n]*dt/Cm
    Isra[n] = a*dt/tau_sra*(Vm[n-1] - EL) - Isra[n-1]*dt/tau_sra + Isra[n-1]
    if Vm[n] > Vmax:
        Vm[n] = Vreset
        Isra[n] += b
        spikes[n] = t[n]
        num_spikes += 1
        
spikes = np.where(spikes != 0)
isi = np.diff(spikes)
isi = [i * dt for i in isi]
#plot histogram of ISI with 25 bins for the values of ISI
plt.figure(6)
plt.hist(isi, 25)
plt.xlabel('ISI (s)')
plt.ylabel('Count')
plt.title('Histogram of ISI: 0.1 nA added')
plt.savefig('T32_Hist_c.1.png')
plt.show()

CV = np.std(isi)/np.mean(isi)
print("The coefficient of variation is: ", CV)

spikes = [i * dt for i in spikes]

#1c.iii
num_spikes = np.zeros(1000)
for i in range(np.size(spikes)):
    x = int(np.fix(spikes[0][i]/.1))
    num_spikes[x] += 1

mean = np.mean(num_spikes)
variance = np.var(num_spikes)
fano_factor = variance/mean

print("The fano factor is: ", fano_factor)


##0.2 nA added to the input current
spikes = np.zeros(len(t))
I_app = np.zeros(len(t))
for i in range(len(t)):
    I_app[i] = np.random.normal(0, sigma/np.sqrt(dt)) + .2E-9   #mean, std dev

num_spikes = 0
isi = np.zeros(len(t))

for n in range(1, len(t)):
    Vm[n] = Vm[n-1] + GL*dt/Cm*(EL-Vm[n-1] + dth*np.exp((Vm[n-1]-Vth)/dth)) - Isra[n-1]*dt/Cm + I_app[n]*dt/Cm
    Isra[n] = a*dt/tau_sra*(Vm[n-1] - EL) - Isra[n-1]*dt/tau_sra + Isra[n-1]
    if Vm[n] > Vmax:
        Vm[n] = Vreset
        Isra[n] += b
        spikes[n] = t[n]
        num_spikes += 1
        
spikes = np.where(spikes != 0)
isi = np.diff(spikes)
isi = [i * dt for i in isi]
#plot histogram of ISI with 25 bins for the values of ISI
plt.figure(7)
plt.hist(isi, 25)
plt.xlabel('ISI (s)')
plt.ylabel('Count')
plt.title('Histogram of ISI: 0.2 nA added')
plt.savefig('T32_Hist_c.2.png')
plt.show()

CV = np.std(isi)/np.mean(isi)
print("The coefficient of variation is: ", CV)

spikes = [i * dt for i in spikes]

#1c.iii
num_spikes = np.zeros(1000)
for i in range(np.size(spikes)):
    x = int(np.fix(spikes[0][i]/.1))
    num_spikes[x] += 1

mean = np.mean(num_spikes)
variance = np.var(num_spikes)
fano_factor = variance/mean

print("The fano factor is: ", fano_factor)