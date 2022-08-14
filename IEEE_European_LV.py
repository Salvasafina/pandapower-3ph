import pandas as pd
import pandapower.networks as nw
import matplotlib.pyplot as plt
import pandapower as pp
import pandapower.plotting as pplot
from pandapower.pf.runpp_3ph import runpp_3ph
from math import acos, sqrt, tan

# load IEEE European LV network from Pandapower
net = nw.ieee_european_lv_asymmetric()

# fix load to nominal values
net.asymmetric_load = net.asymmetric_load.iloc[0:0]
loads_df = pd.read_csv('European_LV_Test_Feeder_v2\European_LV_CSV\Loads.csv', header = 2)
for i in loads_df.index:
    p_load = loads_df.kW[i]/1000
    q_load = p_load*tan(acos(loads_df.PF[i]))/1000
    if loads_df.phases[i] == 'A':
        pp.create_asymmetric_load(net, loads_df.Bus[i], p_load,0,0, q_load,0,0 )
    elif loads_df.phases[i] == 'B':
        pp.create_asymmetric_load(net, loads_df.Bus[i], 0, p_load,0,0, q_load,0)
    elif loads_df.phases[i] == 'C':
        pp.create_asymmetric_load(net, loads_df.Bus[i], 0,0,p_load, 0,0,q_load)

# fictional loads added just to be plotted
pp.create_loads(net= net, buses = net.asymmetric_load.bus, p_mw = 0, q_mvar = 0)

# network diagram including loads
pplot.simple_plot(net, show_plot=False, bus_size= .4, plot_loads= True, load_size= .8)

# assing profiles to loads
loads_df['Profile_number'] = 'NaN'
num_loads = loads_df.shape[0]

load_profiles_df = pd.DataFrame() #empty df that will include all the load profiles
load_profiles_df['time'] = pd.read_csv(f'European_LV_Test_Feeder_v2\European_LV_CSV\Load Profiles\Load_profile_{1}.csv', usecols= [0])
for i in range(num_loads):
    load_profiles_df[f'profile_{i+1}'] = pd.read_csv(f'European_LV_Test_Feeder_v2\European_LV_CSV\Load Profiles\Load_profile_{i+1}.csv', usecols= [1])
    profile = loads_df.Yearly[i].split('_')[1]
    loads_df['Profile_number'][i] = profile # asign the corresponding load profile to the load i
print(load_profiles_df)

# --- plot load profiles --- #
load_profiles_df.plot(x='time', y=load_profiles_df.columns[1:])

# --- 3ph power flow considering profiles --- #
evol_res_bus = []
for i in range(1440):
    print(load_profiles_df['time'][i])
    for j in range(num_loads):
        net.asymmetric_load.scaling[j] = load_profiles_df[f'profile_{loads_df.Profile_number[j]}'][i]
    runpp_3ph(net)
    evol_res_bus.append(net.res_bus_3ph)

# number of elements by type present on the network
print(net)

# on peak case voltage results
res_snapshot = evol_res_bus[565][['vm_a_pu','vm_b_pu','vm_c_pu']]*416/sqrt(3)
with pd.option_context('display.max_rows', None, 'display.max_columns', None,'display.precision', 4):
    print(res_snapshot)

plt.show()
