from brian2 import *

start_scope()

DT=0.1 # time step
defaultclock.dt = DT*ms
N1 = 400 # number of inhibitory neurons
N2 = 1600 # number of excitatory neurons 

TotTime=5000 #Simulation duration (ms)
print('Time=5s')
duration = TotTime*ms

#(0.04*v**2+5*v+140-u-GsynE*(v-Ee)-GsynI*(v-Ei)-I)/Tn:1
eqs="""
dv/dt=(g_iz*(v-E_iz)**2-u-GsynE*(v-Ee)-GsynI*(v-Ei)-I)/Tv:1 (unless refractory)
du/dt=(a*(b*v-u))/Tu:1
dGsynI/dt = -GsynI/Tsyn : 1
dGsynE/dt = -GsynE/Tsyn : 1  
a:1
b:1
c:1
d:1
g_iz:1
E_iz:1
Tsyn:second
Tv:second
Tu:second
Ee:1
Ei:1
I:1
"""
#Pvar:1
#Qi:1
#Qe:1
#"""

# Population 1 - Fast Spiking

G_inh = NeuronGroup(N1, eqs, threshold='v > -30.', reset='v = c; u=u+d',refractory='5*ms', method='heun')

G_inh.v=-50
G_inh.u=0
G_inh.a=0.#0.1
G_inh.b=0.#0.2
G_inh.c=-55
G_inh.d=0.0
G_inh.Tsyn=5.*ms
G_inh.Tv=1*second
G_inh.Tu=1*second
G_inh.Ee=0.0
G_inh.Ei=-80.
G_inh.I=0.
G_inh.g_iz=0.04
G_inh.E_iz=-60


# Population 2 - Regular Spiking
G_exc = NeuronGroup(N2, eqs, threshold='v > -30.', reset='v = c; u=u+d',refractory='5*ms', method='heun')
G_exc.v=-50
G_exc.u=0#25
G_exc.a=0.1
G_exc.b=0.0#0.2
G_exc.c=-65
G_exc.d=0#15.0#8.0
G_exc.Tsyn=5.*ms
G_exc.Tv=1*second
G_exc.Tu=1*second
G_exc.Ee=0.0
G_exc.Ei=-80.
G_exc.I=0.0
G_exc.g_iz=0.01
G_exc.E_iz=-65

# external drive--------------------------------------------------------------------------

P_ed=PoissonGroup(N2, rates=10*Hz)

# Network-----------------------------------------------------------------------------

# connections-----------------------------------------------------------------------------
#seed(0)
Qi=5.0#e-6 #*nS
Qe=1.5#e-6 #*nS

prbC=.05 #0.05
 


S_12 = Synapses(G_inh, G_exc, on_pre='GsynI_post+=Qi') #'v_post -= 1.*mV')
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


PgroupE = NeuronGroup(1, 'P:1', method='heun')
		
PE=Synapses(G_exc, PgroupE, 'P_post = u_pre : 1 (summed)')
PE.connect(p=1)
P2mon = StateMonitor(PgroupE, 'P', record=0)


PgroupMuVE = NeuronGroup(1, 'Pv:1', method='heun')
		
PmuE=Synapses(G_exc, PgroupMuVE, 'Pv_post = v_pre : 1 (summed)')
PmuE.connect(p=1)
P2MuVemon = StateMonitor(PgroupMuVE, 'Pv', record=0)







# Recording tools -------------------------------------------------------------------------------
rec1=1
rec2=2

M1G_inh = SpikeMonitor(G_inh)
FRG_inh = PopulationRateMonitor(G_inh)
M1G_exc = SpikeMonitor(G_exc)
FRG_exc = PopulationRateMonitor(G_exc)

M2G1 = StateMonitor(G_inh, 'v', record=range(rec1))
M3G1 = StateMonitor(G_inh, 'u', record=range(rec1))
M4G1 = StateMonitor(G_inh, 'GsynE', record=range(rec1))
M5G1 = StateMonitor(G_inh, 'GsynI', record=range(rec1))

M2G2 = StateMonitor(G_inh, 'v', record=range(rec2))
M3G2 = StateMonitor(G_inh, 'u', record=range(rec2))
M4G2 = StateMonitor(G_inh, 'GsynE', record=range(rec2))
M5G2 = StateMonitor(G_inh, 'GsynI', record=range(rec2))

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



# Lt1G1=array(M2G1.t/ms)
# Lt2G1=array(M3G1.t/ms)
# Lt1G2=array(M2G2.t/ms)
# Lt2G2=array(M3G2.t/ms)
# LVG1=[]
# LwG1=[]
# LVG2=[]
# LwG2=[]

# LgseG1=[]
# LgsiG1=[]
# LgseG2=[]
# LgsiG2=[]

# for a in range(rec1):
#     LVG1.append(array(M2G1[a].v))
#     LwG1.append(array(M3G1[a].u))
#     LgseG1.append(array(M4G1[a].GsynE/nS))
#     LgsiG1.append(array(M5G1[a].GsynI/nS))

# for a in range(rec2):
#     LVG2.append(array(M2G2[a].v))
#     LwG2.append(array(M3G2[a].u))
#     LgseG2.append(array(M4G2[a].GsynE/nS))
#     LgsiG2.append(array(M5G2[a].GsynI/nS))

# #create the figure
# fig=plt.figure(figsize=(12,8))
# ax1=fig.add_subplot(221)
# ax2=fig.add_subplot(222)
# ax3=fig.add_subplot(223)
# ax4=fig.add_subplot(224)

# for a in range(len(LVG1)):
#     ax1.plot(Lt1G1, LVG1[a], color ='r')
#     ax1.plot(Lt1G1, LVG1[a], ls=(0, (2,6))) 
#     ax3.plot(Lt2G1, LwG1[a],color ='r')
#     ax3.plot(Lt2G1, LwG1[a],ls=(0, (2,6))) 


# for a in range(len(LgseG1)):
#     ax2.plot(Lt1G1, LgsiG1[a], color ='r')
#     ax2.plot(Lt1G1, LgsiG1[a], ls=(0, (2,6))) 
#     ax2.plot(Lt1G1, LgseG1[a], color ='g')
#     ax2.plot(Lt1G1, LgseG1[a], ls=(0, (2,6))) 

# for a in range(len(LVG2)):
#     ax1.plot(Lt1G2, LVG2[a],color ='g')
#     ax1.plot(Lt1G2, LVG2[a], ls=(0, (2,6))) 
#     ax3.plot(Lt2G2, LwG2[a],color ='g')
#     ax3.plot(Lt2G2, LwG2[a], ls=(0, (2,6)))

# for a in range(len(LgsiG2)):
#     ax4.plot(Lt1G1, LgsiG2[a], color ='r')
#     ax4.plot(Lt1G1, LgsiG2[a], ls=(0, (2,6))) 
#     ax4.plot(Lt1G1, LgseG2[a], color ='g')
#     ax4.plot(Lt1G1, LgseG2[a], ls=(0, (2,6))) 



# # prepare raster plot
# RasG_inh = array([M1G_inh.t/ms, [i+N2 for i in M1G_inh.i]])
# RasG_exc = array([M1G_exc.t/ms, M1G_exc.i])



# # prepare firing rate
# def bin_array(array, BIN, time_array):
#     N0 = int(BIN/(time_array[1]-time_array[0]))
#     N1 = int((time_array[-1]-time_array[0])/BIN)
#     return array[:N0*N1].reshape((N1,N0)).mean(axis=1)

# BIN=5
# time_array = arange(int(TotTime/DT))*DT



# LfrG_exc=array(FRG_exc.rate/Hz)
# TimBinned,popRateG_exc=bin_array(time_array, BIN, time_array),bin_array(LfrG_exc, BIN, time_array)

# LfrG_inh=array(FRG_inh.rate/Hz)
# TimBinned,popRateG_inh=bin_array(time_array, BIN, time_array),bin_array(LfrG_inh, BIN, time_array)



# # create the figure

# fig=figure(figsize=(8,12))
# ax1=fig.add_subplot(211)
# ax2=fig.add_subplot(212)


# ax1.plot(RasG_inh[0], RasG_inh[1], ',r')
# ax1.plot(RasG_exc[0], RasG_exc[1], ',g')

# ax1.set_xlabel('Time (ms)')
# ax1.set_ylabel('Neuron index')

# ax2.plot(TimBinned,popRateG_inh, 'r')
# ax2.plot(TimBinned,popRateG_exc, 'g')

# ax2.set_xlabel('Time (ms)')
# ax2.set_ylabel('Firing Rate (Hz)')

# TimBinned,Pu=bin_array(time_array, BIN, time_array),bin_array(P2mon[0].P, BIN, time_array)
# fig=plt.figure(figsize=(8,5))
# ax3=fig.add_subplot(111)
# ax2 = ax3.twinx()
# ax3.plot(TimBinned/1000,popRateG_inh, 'r')
# ax3.plot(TimBinned/1000,popRateG_exc, 'SteelBlue')
# ax2.plot(TimBinned/1000,(Pu/8000), 'orange')
# ax2.set_ylabel('mean u')
# #ax2.set_ylim(0.0, 0.045)
# ax3.set_xlabel('Time (s)')
# ax3.set_ylabel('population Firing Rate')


# TimBinned,Pv=bin_array(time_array, BIN, time_array),bin_array(P2MuVemon[0].Pv, BIN, time_array)
# fig=plt.figure(figsize=(12,5))
# plt.plot(TimBinned/1000,(Pv/8000), 'blue')
# #ax1=fig.add_subplot(211)
# #ax2=fig.add_subplot(212)
# #ax1.set_title('P neuron group')
# #ax1.plot(P1mon.P[0])
# #ax1.plot(P2mon.P[0])
# #ax2.set_title('P received by the neuron')
# #ax2.plot(G1mon.Pvar[0])
# #ax2.plot(G2mon.Pvar[0])

# fig.tight_layout()

# show()



