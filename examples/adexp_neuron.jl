using Plots
using SpikingNeuralNetworks
SNN.@load_units

adparam = SNN.ADEXParameter(;a = 6.050246708405076,
    b = 7.308480222357973,
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
SNN.sim!([E],[], dt=0.1ms, duration=2000ms)
SNN.vecplot(E, :v) |> display
