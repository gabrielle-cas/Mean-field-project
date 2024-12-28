from brian2 import *

start_scope()

DT=0.1 # time step
defaultclock.dt = DT*ms
N1 = 2000 # number of inhibitory neurons
N2 = 8000 # number of excitatory neurons 

TotTime=5000 #Simulation duration (ms)
duration = TotTime*ms


C = 200*pF
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



eqs = """
dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)-GsynE*(vm-Ee)-GsynI*(vm-Ei)+I-w)/C : volt (unless refractory)
dw/dt=(a*(vm-EL)-w)/tauw : amp
dGsynI/dt = -GsynI/TsynI : siemens
dGsynE/dt = -GsynE/TsynE : siemens
TsynI:second
TsynE:second
Vr:volt
b:amp
DeltaT:volt
Vcut:volt
VT:volt
EL:volt
"""


# Population 1 - Fast Spiking

G_inh = NeuronGroup(N1, model=eqs, threshold='vm > Vcut', reset="vm = Vr; w += b", refractory='5*ms', method='heun')
G_inh.vm = -60*mV#EL
G_inh.EL=-67*mV
G_inh.w = a * (G_inh.vm - G_inh.EL)
G_inh.Vr = -65*mV #
G_inh.TsynI =5.0*ms
G_inh.TsynE =5.0*ms
G_inh.b=0*pA
G_inh.DeltaT=0.5*mV
G_inh.VT=-50.*mV
G_inh.Vcut=G_inh.VT + 5 * G_inh.DeltaT



# Population 2 - Regular Spiking

G_exc = NeuronGroup(N2, model=eqs, threshold='vm > Vcut', reset="vm = Vr; w += b", refractory='5*ms', method='heun')
G_exc.vm = -60*mV#EL
G_exc.EL=-63*mV
G_exc.w = a * (G_exc.vm - G_exc.EL)
G_exc.Vr = -65*mV 
G_exc.TsynI =5.0*ms
G_exc.TsynE =5.0*ms
G_exc.b=5*pA
G_exc.DeltaT=2*mV
G_exc.VT=-50.*mV
G_exc.Vcut=G_exc.VT + 5 * G_exc.DeltaT


# external drive--------------------------------------------------------------------------

P_ed=PoissonGroup(N2, rates=5*Hz)

# Network-----------------------------------------------------------------------------

# connections-----------------------------------------------------------------------------
#seed(0)
Qi=5.0*nS
Qe=1.5*nS

prbC=.05 #0.05
 

S_12 = Synapses(G_inh, G_exc, on_pre='GsynI_post+=Qi',) 
S_12.connect('i!=j', p=prbC)

S_11 = Synapses(G_inh, G_inh, on_pre='GsynI_post+=Qi')
S_11.connect('i!=j',p=prbC)

S_21 = Synapses(G_exc, G_inh, on_pre='GsynE_post+=Qe')
S_21.connect('i!=j',p=prbC)

S_22 = Synapses(G_exc, G_exc, on_pre='GsynE_post+=Qe')
S_22.connect('i!=j', p=prbC)

S_ed_in = Synapses(P_ed, G_inh, on_pre='GsynE_post+=Qe')
S_ed_in.connect(p=prbC)

S_ed_ex = Synapses(P_ed, G_exc, on_pre='GsynE_post+=Qe')
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


