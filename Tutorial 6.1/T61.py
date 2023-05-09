import numpy as np
import matplotlib.pyplot as plt

D = 1
alpha_0 = 0.5
sigma = 0.5
WEE = 8
p_r = 1
tau_s = 2E-3
tau_r = 10E-3
r_max = 100
r_0 = 0.1

#####
## 1
#####

##Part a
ds = .1E-3
s1 = np.arange(0, 1/WEE, ds)
f = np.zeros(len(s1))

for i in range(len(s1)):
    S = WEE * s1[i]
    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0

dr = .1E-3
r = np.arange(0, r_max, dr)
s2 = np.zeros(len(r))

for i in range(len(r)):
    s2[i] = (alpha_0 * D * p_r * r[i] * tau_s) / (1 + alpha_0 * D * p_r * r[i] * tau_s)

plt.figure(1)
plt.plot(s1, f, label='f(WEEs1)')
plt.xlabel('synaptic input')
plt.ylabel('firing rate curve')
plt.plot(s2, r, label="s2(r)")
plt.xlabel('synaptic input')
plt.ylabel('firing rate')
plt.title('1a. firing rate vs s')
plt.legend()
plt.savefig('1a. firing rate vs s.png')

##Part b
dt = .1E-3
t = np.arange(0, 20, dt)

s_in = 0.05
f = np.ones(len(t))
s = np.zeros(len(t))
r = np.zeros(len(t))

for i in range(1, len(t)):
    #if in the first 50 ms of the simulation 
    if i >= 100000 and i <= 100500:
        S = WEE*s[i-1] + s_in
    else:
        S = WEE*s[i-1]

    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0
    
    r[i] = r[i-1] + dt * (-r[i-1] + f[i]) / tau_r

    s[i] = s[i-1] + dt * (-s[i-1]/tau_s + alpha_0 * D * p_r * r[i-1] * (1-s[i-1]))
    

plt.figure(2)
plt.plot(t, r, label='firing rate')
plt.xlabel('time')
plt.ylabel('firing rate')
plt.title('1b. firing rate vs time')
plt.legend()
plt.savefig('1b. firing rate vs time.png')

#####
## 2
#####

tau_D = 250E-3
p_r = 0.2
WEE = 60

##Part a

ds = .1E-3
s1 = np.arange(0, 1/WEE, ds)
f = np.zeros(len(s1))

for i in range(len(s1)):
    S = WEE * s1[i]
    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0

dr = .1E-3
r = np.arange(0, r_max, dr)
s2 = np.zeros(len(r))
D = np.zeros(len(r))

for i in range(len(r)):
    D[i] = 1/(1 + p_r*r[i-1]*tau_D)
    s2[i] = (alpha_0 * D[i-1] * p_r * r[i-1] * tau_s) / (1 + alpha_0 * D[i-1] * p_r * r[i-1] * tau_s)

plt.figure(3)
plt.plot(s1, f, label='f(WEEs1)')
plt.xlabel('synaptic input')
plt.ylabel('firing rate curve')
plt.plot(s2, r, label="s2(r)")
plt.xlabel('synaptic input')
plt.ylabel('firing rate')
plt.title('2a. firing rate vs s')
plt.legend()
plt.savefig('2a. firing rate vs s.png')

##Part b

dt = .1E-3
t = np.arange(0, 20, dt)
s_in = 0.002

s = np.zeros(len(t))
f = np.ones(len(t))
r = np.zeros(len(t))
D = np.zeros(len(r))

for i in range(1, len(t)):
    #if in the first 50 ms of the simulation 
    if t[i] >= 10 and t[i] <= 12:
        S = s_in * WEE
    else:
        S = WEE*s[i-1]

    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0
    
    r[i] = r[i-1] + dt * (-r[i-1] + f[i]) / tau_r
    D[i] = D[i-1] + dt * ((1-D[i-1])/tau_D - p_r*D[i-1]*r[i-1])
    s[i] = s[i-1] + dt * (-s[i-1]/tau_s + alpha_0 * D[i-1] * p_r * r[i-1] * (1-s[i-1]))
    

plt.figure(4)
plt.plot(t, r, label='firing rate')
plt.xlabel('time')
plt.ylabel('firing rate r(t) ')
plt.title('2b. firing rate vs time')
plt.legend()
plt.savefig('2b. firing rate vs time.png')


#####
## 3
#####

p_r = 0.5
WEE = 35
r_0 = -0.1

##Part a

ds = .1E-3
s1 = np.arange(0, 1/WEE, ds)
f = np.zeros(len(s1))

for i in range(len(s1)):
    S = WEE * s1[i]
    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0

dr = .1E-3
r = np.arange(0, r_max, dr)
s2 = np.zeros(len(r))
D = np.zeros(len(r))

for i in range(len(r)):
    D[i] = 1/(1 + p_r*r[i-1]*tau_D)
    s2[i] = (alpha_0 * D[i-1] * p_r * r[i] * tau_s) / (1 + alpha_0 * D[i-1] * p_r * r[i] * tau_s)

plt.figure(5)
plt.plot(s1, f, label='f(WEEs1)')
plt.xlabel('synaptic input')
plt.ylabel('firing rate curve')
plt.plot(s2, r, label="s2(r)")
plt.xlabel('synaptic input')
plt.ylabel('firing rate')
plt.title('3a. firing rate vs s')
plt.legend()
plt.savefig('3a. firing rate vs s.png')

##Part b
s_in = 0.002

dt = .1E-3
t = np.arange(0, 20, dt)

f = np.ones(len(t))
s = np.zeros(len(t))
r = np.zeros(len(t))
D = np.zeros(len(r))

for i in range(1, len(t)):
    #if in the first 50 ms of the simulation 
    if i >= 100000 and i <= 100500:
        S = WEE*s_in
    else:
        S = WEE*s[i-1]

    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0
    
    r[i] = r[i-1] + dt * (-r[i-1] + f[i]) / tau_r
    D[i] = D[i-1] + dt * ((1-D[i-1])/tau_D - p_r*D[i-1]*r[i-1])
    s[i] = s[i-1] + dt * (-s[i-1]/tau_s + alpha_0 * D[i-1] * p_r * r[i-1] * (1-s[i-1]))
    

plt.figure(6)
plt.plot(t, r, label='firing rate')
plt.xlabel('time')
plt.ylabel('firing rate')
plt.title('3b. firing rate vs time')
plt.legend()
plt.savefig('3b. firing rate vs time.png')


##Part c

f = np.ones(len(t))
s = np.zeros(len(t))
r = np.zeros(len(t))
D = np.zeros(len(r))

r[0] = 9
D[0] = 1/(1 + p_r*r[0]*tau_D)
s[0] = s[0] + dt * (-s[0]/tau_s + alpha_0 * D[0] * p_r * r[0] * (1-s[0]))

for i in range(1, len(t)):
    #if in the first 50 ms of the simulation 
    if i >= 100000 and i <= 100500:
        S = WEE*s_in
    else:
        S = WEE*s[i-1]

    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0
    
    r[i] = r[i-1] + dt * (-r[i-1] + f[i]) / tau_r
    D[i] = D[i-1] + dt * ((1-D[i-1])/tau_D - p_r*D[i-1]*r[i-1])
    s[i] = s[i-1] + dt * (-s[i-1]/tau_s + alpha_0 * D[i-1] * p_r * r[i-1] * (1-s[i-1]))
    

plt.figure(7)
plt.plot(t, r, label='firing rate')
plt.xlabel('time')
plt.ylabel('firing rate')
plt.title('3c. firing rate vs time')
plt.legend()
plt.savefig('3c. firing rate vs time.png')


#####
## 4
#####

tau_D = 125E-3
alpha_0 = 0.25
p_r = 1

#Part a
ds = .1E-3
s1 = np.arange(0, 1/WEE, ds)
f = np.zeros(len(s1))

for i in range(len(s1)):
    S = WEE * s1[i]
    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0

dr = .1E-3
r = np.arange(0, r_max, dr)
s2 = np.zeros(len(r))
D = np.zeros(len(r))

for i in range(len(r)):
    D[i] = 1/(1 + p_r*r[i-1]*tau_D)
    s2[i] = (alpha_0 * D[i-1] * p_r * r[i] * tau_s) / (1 + alpha_0 * D[i-1] * p_r * r[i] * tau_s)

plt.figure(8)
plt.plot(s1, f, label='f(WEEs1)')
plt.xlabel('synaptic input')
plt.ylabel('firing rate curve')
plt.plot(s2, r, label="s2(r)")
plt.xlabel('synaptic input')
plt.ylabel('firing rate')
plt.title('4a. firing rate vs s')
plt.legend()
plt.savefig('4a. firing rate vs s.png')

#Part b
s_in = 0.002
dt = .1E-3
t = np.arange(0, 20, dt)

f = np.ones(len(t))
s = np.zeros(len(t))
r = np.zeros(len(t))
D = np.zeros(len(r))

for i in range(1, len(t)):
    if i >= 100000 and i <= 100500:
        S = WEE* s_in
    else:
        S = WEE*s[i-1]

    if S > 0:
        f[i] = r_0 + r_max * (S**1.2 / (S**1.2 + sigma**1.2))
    else:
        f[i] = r_0
    
    r[i] = r[i-1] + dt * (-r[i-1] + f[i]) / tau_r
    D[i] = D[i-1] + dt * ((1-D[i-1])/tau_D - p_r*D[i-1]*r[i-1])
    s[i] = s[i-1] + dt * (-s[i-1]/tau_s + alpha_0 * D[i-1] * p_r * r[i-1] * (1-s[i-1]))
    

plt.figure(9)
plt.plot(t, r, label='firing rate')
plt.xlabel('time')
plt.ylabel('firing rate')
plt.title('4b. firing rate vs time')
plt.legend()
plt.savefig('4b. firing rate vs time.png')
plt.show()