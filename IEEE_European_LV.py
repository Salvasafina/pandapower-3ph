import pandas as pd
import pandapower.networks as nw
import matplotlib.pyplot as plt
import pandapower as pp
import pandapower.plotting as pplot
from pandapower.pf.runpp_3ph import runpp_3ph
from math import acos, tan
from functions_pp_3ph import *

# load IEEE European LV network from Pandapower
net = nw.ieee_european_lv_asymmetric()

# fix load to nominal values
net.asymmetric_load = net.asymmetric_load.iloc[0:0]
loads_df = pd.read_csv('European_LV_Test_Feeder_v2\European_LV_CSV\Loads.csv', header = 2)
for i in loads_df.index:
    p_load = loads_df.kW[i]/1000
    q_load = p_load*tan(acos(loads_df.PF[i]))
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


# --- 3ph power flow considering profiles --- #
evol_res_bus = []
evol_p_lines = []
impedances = get_zabc(net)
for i in range(560,570):
    time = load_profiles_df['time'][i]
    for j in range(num_loads):
        net.asymmetric_load.scaling[j] = load_profiles_df[f'profile_{loads_df.Profile_number[j]}'][i]
    print(f'Solving for time:{time}')
    runpp_3ph(net)
    p_rlcp = rlcp(net, impedances)
    evol_res_bus.append(net.res_bus_3ph)
    evol_p_lines.append(p_rlcp)


# number of elements by type present on the network
print(net)

# on peak case voltage and losses results
res_v_snapshot = evol_res_bus[5][['vm_a_pu','vm_b_pu','vm_c_pu']]
res_p_snaphot = evol_p_lines[5]
print(f'Total losses by phase at time 09:26:00: \n{res_p_snaphot.sum()}')
print(f'Total losses at time 09:26:00: \n{res_p_snaphot.sum().sum()}')
plt.show()

