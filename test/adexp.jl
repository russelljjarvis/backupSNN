include("../src/SpikingNeuralNetworks.jl")
include("../src/unit.jl")
include("../src/plot.jl")

SNN = SpikingNeuralNetworks
adparam = SNN.ADEXParameter(;a = 6.050246708405076, b = 7.308480222357973,
    cm = 803.1019662706587,
    v_rest= -63.22881649139353,
    tau_m=19.73777028610565,
    tau_w=351.0551915202058,
    v_thresh=-39.232165554444265,
    delta_T=6.37124632135508,
    v_spike = -36.58205819488362,
    v_reset = -59.18792270568965,
    spike_delta = 16.33506432689027)

E = SNN.AD(;N = 1, param=adparam)
E.I = [795.57128906]

SNN.monitor(E, [:v,:I])
SNN.sim!(E, [], 0.1ms)
SNN.vecplot(E, :v) |> display

#SNN.sim!([Etest],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=2000ms)
#SNN.sim!([E],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=2000ms)






Ne = 1000;
Ni = 600#int(200.0*(10.0/8.0))
#E = SNN.IZ(;N = Ne, param = SNN.IZParameter(;a = 0.02, b = 0.2, c = -65, d = 8));#, C=1))
#I = SNN.IZ(;N = Ni, param = SNN.IZParameter(;a = 0.1, b = 0.2, c = -65, d = 2));#, C=1))

param = SNN.IZParameter_more(;a = 0.02, b = 0.2, c = -65.0, d = 8.0, C=100.0, vr=-65.0, k=4.0, vt=-50.0, vPeak=20.0, cellnumber=7)
E = SNN.IZ_more(;N = 1, param)
#E = SNN.IZ_more(;N = 1, param)
E.I = [190.0015]
SNN.monitor(E, [:v,:I])

SNN.sim!([E],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=2000ms)
SNN.vecplot(E, :v) |> display

param = SNN.IZParameter_more(;a = 0.02, b = 0.2, c = -65.0, d = 8.0, C=90.0, vr=-65.0, k=1.6, vt=-62.0, vPeak=20.0, cellnumber=6)
E = SNN.IZ_more(;N = Ne, param)

param = SNN.IZParameter_more(;a = 0.02, b = 0.2, c = -65.0, d = 8.0, C=190.0, vr=-75.0, k=1.6, vt=-62.0, vPeak=20.0, cellnumber=5)

I = SNN.IZ_more(;N = Ni, param)
I.I = [160.0015]

SNN.monitor(I, [:v,:I])
SNN.sim!([I],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=2000ms)
SNN.vecplot(I, :v) |> display



SNN.vecplot(I, :I) |> display

Etest =SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.02, b = 0.2, c = -65, d = 8))#, C=1))
Etest.I = [3.1]
SNN.monitor(Etest, [:v])
SNN.sim!([Etest],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=1000ms)
SNN.vecplot(Etest, :v) |> display

EE = SNN.SpikingSynapse(E, E, :v; σ = 0.25,  p = 0.25);
EI = SNN.SpikingSynapse(E, I, :v; σ = 0.25,  p = 0.25);
IE = SNN.SpikingSynapse(I, E, :v; σ = -1.0, p = 1.0);
II = SNN.SpikingSynapse(I, I, :v; σ = -1.0, p = 0.25);

#EE = SNN.SpikingSynapse(E, E, :v; σ = 1,  p = 1);
#EI = SNN.SpikingSynapse(E, I, :v; σ = 1,  p = 1);
#IE = SNN.SpikingSynapse(I, E, :v; σ = -1.0, p = 1);
#II = SNN.SpikingSynapse(I, I, :v; σ = -1.0, p = 1);
println(size(E.N),size(E.N))

P = [E, I];
C = [EE, EI, IE, II];
println("passed !!!")

SNN.monitor([E, I], [:fire])
for t = 1:2000
    E.I .= 5randn(Ne)
    I.I .= 2randn(Ni)
    SNN.sim!(P, C, 1ms)
end
SNN.raster(P) #|> display
#SNN.vecplot(E, :fire) |> display

#{'C': 114.25027724782701, 'k': 0.714430107979936, 'vr': -68.3530583578185, 'vt': -49.88514142246343, 'vPeak': 41.99002172075032, 'a': 0.01223913281254824, 'b': -1.9721425473513066, 'c': -49.74145056836865, 'd': 94.78925041952787}
Etest = SNN.IZ(;N = 1, param = SNN.IZParameter(;a =  0.01223913281254824, b = -1.9721425473513066, c = -49.74145056836865, d = 94.78925041952787))#,C=1.1425,vr=-68.3530583578185));
Etest.I = [168.275];
SNN.monitor(Etest, [:v]);
SNN.sim!([Etest],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=1000ms);

SNN.vecplot(Etest, :v)# |> display







Ne = 1000;
Ni = 600;
Nb = 50
E = SNN.IZ(;N = Ne, param = SNN.IZParameter(;a = 0.02, b = 0.2, c = -65, d = 8));#, C=1))
#E = SNN.IZ_more(;N = Ne, param = SNN.IZ_moreParameter(;a =  0.01223913281254824, b = -1.9721425473513066, c = -49.74145056836865, d = 94.78925041952787))#,C=1.1425,vr=-68.3530583578185));

I = SNN.IZ(;N = Ni, param = SNN.IZParameter(;a = 0.1, b = 0.2, c = -65, d = 2));
EE = SNN.SpikingSynapse(E, E, :v; σ = 0.25,  p = 0.25);
EI = SNN.SpikingSynapse(E, I, :v; σ = 0.25,  p = 0.25);
IE = SNN.SpikingSynapse(I, E, :v; σ = -1.0, p = 1.0);
II = SNN.SpikingSynapse(I, I, :v; σ = -1.0, p = 0.25);

E_new = SNN.IZ(;N = Nb, param = SNN.IZParameter(;a =  0.01223913281254824, b = -1.9721425473513066, c = -49.74145056836865, d = 94.78925041952787))#,C=1.1425,vr=-68.3530583578185));
#EE_new = SNN.SpikingSynapse(E_new,E_new, :v; σ = 0.25,  p = 0.25);
#P = [E_new,E, I];
#C = [EE_new,EE, EI, IE, II];

P = [E, I];
C = [EE, EI, IE, II]#,EE_new];

SNN.monitor([E, I], [:fire])
for t = 1:1000
    #E_new.I .= 2randn(Nb)
    E.I .= 0.0112482752(Ne)
    I.I .= 1.68275(Ni)
    SNN.sim!(P, C, 1ms)
end
SNN.raster(P) |> display

typeof(E)

Ne = 1000;      Ni = 200#int(200.0*(10.0/8.0))
Etest = SNN.IZ(;N = 1, param = SNN.IZParameter(;a =  0.07976319652352468, b = -1.545038853189859, c = -52.726764975633166, d = 144.76185252557173))#, C=0.9336,vr=-65.55690901865944))
Etest.I = [1100]
SNN.monitor(Etest, [:v])
SNN.sim!([Etest],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=1000ms)
#SNN.sim!([Etest],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=1000ms)

SNN.vecplot(Etest, :v) |> display

#{'C': 93.36810828872358, 'k': 0.7094612672234731, 'vr': -65.55690901865944, 'vt': -49.93900899948406, 'vPeak': 46.969034033727034, 'a': 0.07976319652352468, 'b': -1.545038853189859, 'c': -52.726764975633166, 'd': 144.76185252557173}
E = SNN.IZ(;N = Ne, param = SNN.IZParameter(;a =  0.07976319652352468, b = -1.545038853189859, c = -52.726764975633166, d = 144.76185252557173))#, C=0.9336,vr=-65.55690901865944))
I = SNN.IZ(;N = Ni, param = SNN.IZParameter(;a = 0.1, b = 0.2, c = -65, d = 2))
#I = SNN.IZ(;N = Ni, param = SNN.IZParameter(;a =  0.07976319652352468, b = -1.545038853189859, c = -52.726764975633166, d = 144.76185252557173))#, C=0.9336,vr=-65.55690901865944))

EE = SNN.SpikingSynapse(E, E, :v; σ = 0.25,  p = 0.5);
EI = SNN.SpikingSynapse(E, I, :v; σ = 0.25,  p = 0.5);
IE = SNN.SpikingSynapse(I, E, :v; σ = -1.0, p = 1.0);
II = SNN.SpikingSynapse(I, I, :v; σ = -1.0, p = 0.25);
P = [E, I]
C = [EE, EI, IE, II]

SNN.monitor([E, I], [:fire])
for t = 1:2000
    E.I .= 5randn(Ne)
    I.I .= 2randn(Ni)
    SNN.sim!(P, C, 1ms)
end
SNN.raster(P) |> display

println("\n olf mit \n")

Etest = SNN.IZ(;N = 1, param = SNN.IZParameter(;a =  0.015190425633023837, b = -1.989511454925138, c = -54.859511795883556, d = 129.7659026061511))#, C=1.99615880961252,vr=-59.78177430465574))
Etest.I = [400]
SNN.monitor(Etest, [:v])
#SNN.sim!([Etest],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=1000ms)
SNN.sim!([Etest],[0],dt = 0.25ms, simulation_duration = 2000ms, delay = 500ms,stimulus_duration=1000ms)

SNN.vecplot(Etest, :v) #|> display
#E = SNN.HH(;N = 1)
#E.I = [0.0752nA]
#SNN.monitor(E, [:v])
#SNN.sim!([E], []; dt = 0.015ms, delay=100ms, stimulus_duration=1000ms, duration = 1300ms)
#SNN.vecplot(E, :v) |> display

Ne = 1000;      Ni = 200#int(200.0*(10.0/8.0))
#{'C': 93.36810828872358, 'k': 0.7094612672234731, 'vr': -65.55690901865944, 'vt': -49.93900899948406, 'vPeak': 46.969034033727034, 'a': 0.07976319652352468, 'b': -1.545038853189859, 'c': -52.726764975633166, 'd': 144.76185252557173}
E = SNN.IZ(;N = Ne, param = SNN.IZParameter(;a =  0.015190425633023837, b = -1.989511454925138, c = -54.859511795883556, d = 129.7659026061511))#, C=1.99615880961252,vr=-59.78177430465574))
I = SNN.IZ(;N = Ni, param = SNN.IZParameter(;a = 0.1, b = 0.2, c = -65, d = 2))
#I = SNN.IZ(;N = Ni, param = SNN.IZParameter(;a =  0.015190425633023837, b = -1.989511454925138, c = -54.859511795883556, d = 129.7659026061511))#, C=1.99615880961252,vr=-59.78177430465574))

EE = SNN.SpikingSynapse(E, E, :v; σ = 0.25,  p = 0.5)
EI = SNN.SpikingSynapse(E, I, :v; σ = 0.25,  p = 0.5)
IE = SNN.SpikingSynapse(I, E, :v; σ = -1.0, p = 1.0)
II = SNN.SpikingSynapse(I, I, :v; σ = -1.0, p = 0.25)
P = [E, I]
C = [EE, EI, IE, II]

SNN.monitor([E, I], [:fire])
for t = 1:2000
    E.I .= 5randn(Ne)
    I.I .= 2randn(Ni)
    SNN.sim!(P, C, 1ms)
end
SNN.raster(P) |> display

#using Plots
#=
RS = SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.02, b = 0.2, c = -65, d = 8))
IB = SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.02, b = 0.2, c = -55, d = 4))
CH = SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.02, b = 0.2, c = -50, d = 2))
FS = SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.1, b = 0.2, c = -65, d = 2))
TC1 = SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.02, b = 0.25, c = -65, d = 0.05))
TC2 = SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.02, b = 0.25, c = -65, d = 0.05))
RZ = SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.1, b = 0.26, c = -65, d = 2))
LTS = SNN.IZ(;N = 1, param = SNN.IZParameter(;a = 0.1, b = 0.25, c = -65, d = 2))
P = [RS, IB, CH, FS, TC1, TC2, RZ, LTS]

SNN.monitor(P, [:v])
T = 2second
for t = 1:20000
    for p in [RS, IB, CH, FS, LTS]
        p.I = [10]
    end
    TC1.I = [(t < 0.2T) ? 0mV : 2mV]
    TC2.I = [(t < 0.2T) ? -30mV : 0mV]
    RZ.I =  [(0.5T < t < 0.6T) ? 10mV : 0mV]
    SNN.sim!(P, [], 0.1ms)
end
=#
