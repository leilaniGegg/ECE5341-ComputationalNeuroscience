import numpy as np
import matplotlib.pyplot as plt

#define PIR function
def PIR(Iapp_base, step):
    GL = 10E-9 #10 nS
    G_Na_max = 3.6E-6 #3.6 uS
    G_K_max = 1.6E-6 #1.6 uS
    G_T_max = 0.22E-6 #0.22 uS
    E_Na = 55E-3 #55 mV
    E_k = -90E-3 #-90 mV
    E_Ca = 120E-3 #120 mV
    EL = -70E-3 #-70 mV
    Cm = 100E-12 #100 pF

    dt = 0.01E-3 #0.01 ms
    t = np.arange(0, 0.75, dt)

    Vm = np.zeros(len(t))       # membrane potential
    Vm[0] = EL
    h = np.zeros(len(t))        # Na inactivation
    n = np.zeros(len(t))        # K activation
    h_T = np.zeros(len(t))      
    m = np.zeros(len(t))        # Na activation
    m_T = np.zeros(len(t))    

    Iapp = np.zeros(len(t))
    Iapp[0:25000] = Iapp_base
    Iapp[25000:50000] = Iapp_base + step
    Iapp[50000:75000] = Iapp_base

    for i in range(1, len(t)):
        Vm[i] = Vm[i-1] + dt*(GL*(EL-Vm[i-1]) +G_Na_max*m[i-1]**3*h[i-1]*(E_Na-Vm[i-1]) + G_K_max*n[i-1]**4*(E_k-Vm[i-1]) + G_T_max*m_T[i-1]**2*h_T[i-1]*(E_Ca-Vm[i-1]) + Iapp[i-1]) / Cm 
        alpha_h = 350*np.exp(-50*(Vm[i]+0.058))
        beta_h = 5000/(1+np.exp(-100*(Vm[i]+0.028)))
        h[i] = dt*(alpha_h*(1-h[i-1])-beta_h*h[i-1]) + h[i-1]

        alpha_n = (5*10**4*(Vm[i] + 0.034))/(1-np.exp(-100*(Vm[i] + 0.034)))
        beta_n = 625 * np.exp(-12.5*(Vm[i] + 0.044))
        n[i] = dt*(alpha_n*(1-n[i-1])-beta_n*n[i-1]) + n[i-1]

        h_T_inf = 1/(1+np.exp(500*(Vm[i]+0.076)))
        tau_h_T = 0
        if Vm[i] < -0.080:
            tau_h_T = 0.001 * np.exp(15*(Vm[i] + 0.467))
        else:
            tau_h_T = 0.028 + 0.001 * np.exp(-(Vm[i]+0.022)/0.0105)
        h_T[i] = dt*(h_T_inf - h_T[i-1])/tau_h_T + h_T[i-1]

        alpha_m = (10**5*(Vm[i] + 0.035))/(1-np.exp(-100*(Vm[i] + 0.035)))
        beta_m = 4000 * np.exp(-(Vm[i]+0.06)/0.018)
        m[i] = alpha_m / (alpha_m + beta_m)

        m_T[i] = 1/(1+np.exp(-(Vm[i]+0.052)/0.0074))

    return Vm, t, Iapp