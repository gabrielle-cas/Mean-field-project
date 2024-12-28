from brian2 import *

start_scope()

DT=0.1 # time step
defaultclock.dt = DT*ms
N1 = 400 # number of inhibitory neurons
N2 = 1600 # number of excitatory neurons 

TotTime=5000 #Simulation duration (ms)
duration = TotTime*ms


Cm = 200*pF
gL = 10*nS
#EL = -65*mV
#VT = -50.4*mV
#DeltaT = 2*mV
tauw = 500*ms
a =0.0*nS# 4*nS
#b = 0.08*nA
I = 0.*nA
Ee=0*mV
Ei=-80*mV
#Vcut = VT + 5 * DeltaT  # practical threshold condition
#N = 2000
#NRS=8000

# Synapse parameters
V_E = 0. * mV  # reversal potential for excitatory synapses
V_I = -70. * mV  # reversal potential for inhibitory synapses
tau_AMPA = 2.0 * ms  # AMPA synapse decay
tau_NMDA_rise = 2.0 * ms  # NMDA synapse rise
tau_NMDA_decay = 100.0 * ms  # NMDA synapse decay
tau_GABA = 10.0 * ms  # GABA synapse decay
alpha = 0.5 * kHz  # saturation of NMDA channels at high presynaptic firing rates
C = 1 * mmole  # extracellular magnesium concentration


eqs = """
dvm/dt=(gL*(EL-vm) + I - I_AMPA - I_NMDA - I_GABA -  I_AMPA_ext)/Cm : volt (unless refractory)

I_AMPA = s_AMPA * (vm - V_E) : amp
ds_AMPA / dt = - s_AMPA / tau_AMPA : siemens

I_GABA = s_GABA * (vm - V_I) : amp
ds_GABA / dt = - s_GABA / tau_GABA : siemens

I_AMPA_ext = s_AMPA_ext * (vm - V_E) : amp  
ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : siemens

I_NMDA = gEEN * s_NMDA_tot * (vm - V_E) / ( 1 + exp(-0.062 * vm/mvolt) ) : amp
s_NMDA_tot : 1
ds_NMDA / dt = - s_NMDA / tau_NMDA_decay + alpha * x * (1 - s_NMDA) : 1
dx / dt = - x / tau_NMDA_rise : 1


Vr:volt
VT:volt
EL:volt  
"""

# I_NMDA = gEEN * s_NMDA_tot * (vm - V_E) / ( 1 + exp(-0.062 * vm/mvolt) ) : amp
# s_NMDA_tot : 1

# ds_NMDA / dt = - s_NMDA / tau_NMDA_decay + alpha * x * (1 - s_NMDA) : 1
# dx / dt = - x / tau_NMDA_rise : 1


# Population 1 - Fast Spiking

G_inh = NeuronGroup(N1, model=eqs, threshold='vm > VT', reset="vm = Vr", refractory='5*ms', method='heun')
G_inh.vm = -60*mV#EL
G_inh.EL=-67*mV
G_inh.Vr = -65*mV #
G_inh.VT=-50.*mV



# Population 2 - Regular Spiking

G_exc = NeuronGroup(N2, model=eqs, threshold='vm > VT', reset="vm = Vr", refractory='5*ms', method='heun')
G_exc.vm = -60*mV#EL
G_exc.EL=-63*mV
G_exc.Vr = -65*mV 
G_exc.VT=-50.*mV


# external drive--------------------------------------------------------------------------

P_ed=PoissonGroup(N2, rates=10*Hz)

# Network-----------------------------------------------------------------------------

# connections-----------------------------------------------------------------------------
prbC=0.1 #0.05

# Synaptic conductances
gextE = 2.496 * nS  # external -> excitatory neurons (AMPA)
gextI = 1.944 * nS  # external -> inhibitory neurons (AMPA)
gEEA = 0.104 * nS   # excitatory -> excitatory neurons (AMPA)
gEIA = 0.081 * nS   # excitatory -> inhibitory neurons (AMPA)
gEEN = 0.327 * nS   # excitatory -> excitatory neurons (NMDA)
gEIN = 0.258 * nS  # excitatory -> inhibitory neurons (NMDA)
gIE = 4.375 * nS  # inhibitory -> excitatory neurons (GABA)
gII = 3.4055 * nS  # inhibitory -> inhibitory neurons (GABA)
 

S_AMPA_EE = Synapses(G_exc, G_exc, on_pre='s_AMPA_post+=1.55*gEEA',) 
S_AMPA_EE.connect('i!=j', p=prbC)

S_AMPA_EI = Synapses(G_exc, G_inh, on_pre='s_AMPA_post+=gEIA')
S_AMPA_EI.connect('i!=j',p=prbC)

S_NMDA_EE = Synapses(G_exc, G_exc, on_pre='x_post += 1.55')
S_NMDA_EE.connect('i!=j',p=prbC)

S_NMDA_EI = Synapses(G_exc, G_inh, on_pre='x_post += 1')
S_NMDA_EI.connect('i!=j', p=prbC)

S_GABA_II = Synapses(G_inh, G_inh, on_pre='s_GABA_post+=gII')
S_GABA_II.connect('i!=j',p=prbC)

S_GABA_IE = Synapses(G_inh, G_exc, on_pre='s_GABA_post+=gIE')
S_GABA_IE.connect('i!=j', p=prbC)

S_ed_in = Synapses(P_ed, G_inh, on_pre='s_AMPA_ext_post += gextI')
S_ed_in.connect(p=prbC)

S_ed_ex = Synapses(P_ed, G_exc, on_pre='s_AMPA_ext_post += gextE')
S_ed_ex.connect(p=prbC)



# Recording tools -------------------------------------------------------------------------------

M1G_inh = SpikeMonitor(G_inh)
FRG_inh = PopulationRateMonitor(G_inh)
M1G_exc = SpikeMonitor(G_exc)
FRG_exc = PopulationRateMonitor(G_exc)



# Run simulation -------------------------------------------------------------------------------

print('--##Start simulation##--')
run(duration, report='stdout', profile=True)
print(profiling_summary())
print('--##End simulation##--')


# Plots -------------------------------------------------------------------------------

# prepare raster plot
RasG_inh = array([M1G_inh.t/ms, [i+N2 for i in M1G_inh.i]])
RasG_exc = array([M1G_exc.t/ms, M1G_exc.i])



''' The following function is very important as it defines the time windows through which the rates will be calculated. 
In other words, out of an ensemble of discrete spiking events, we set bins to count them, and form a mean population 
firing rate that is time dependent. The size of the bins most likely influence the results we can obtain, and this is 
one of the questions you can ask'''


def bin_array(array, BIN, time_array):
    N0 = int(BIN/(time_array[1]-time_array[0]))
    N1 = int((time_array[-1]-time_array[0])/BIN)
    return array[:N0*N1].reshape((N1,N0)).mean(axis=1)

BIN=5 ## Size of the time windows in ms
time_array = arange(int(TotTime/DT))*DT



LfrG_exc=array(FRG_exc.rate/Hz)
TimBinned,popRateG_exc=bin_array(time_array, BIN, time_array),bin_array(LfrG_exc, BIN, time_array)

LfrG_inh=array(FRG_inh.rate/Hz)
TimBinned,popRateG_inh=bin_array(time_array, BIN, time_array),bin_array(LfrG_inh, BIN, time_array)



# create the figure

fig=figure(figsize=(8,8))
ax1=fig.add_subplot(211)
ax2=fig.add_subplot(212)


ax1.plot(RasG_inh[0], RasG_inh[1], ',r')
ax1.plot(RasG_exc[0], RasG_exc[1], ',g')

ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Neuron index')

ax2.plot(TimBinned,popRateG_inh, 'r')
ax2.plot(TimBinned,popRateG_exc, 'g')

ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Firing Rate (Hz)')


fig.tight_layout()

plt.show()


