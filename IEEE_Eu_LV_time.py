import pandas as pd
import pandapower as pp
import pandapower.networks as nw
from math import acos, tan
from pandapower.pf.runpp_3ph import runpp_3ph
from functions_pp_3ph import *
import pandapower.plotting as pplot
import matplotlib.pyplot as plt

# load IEEE European LV network from Pandapower
net = nw.ieee_european_lv_asymmetric()

# create dataframe with loads nominal values
loads_df = pd.read_csv('European_LV_Test_Feeder_v2\European_LV_CSV\Loads.csv', header = 2)

# fix load to nominal values
net.asymmetric_load = net.asymmetric_load.iloc[0:0]
for i in loads_df.index:
    p_load = loads_df.kW[i]/1000
    q_load = p_load*tan(acos(loads_df.PF[i]))
    if loads_df.phases[i] == 'A':
        pp.create_asymmetric_load(net, loads_df.Bus[i], p_load,0,0, q_load,0,0 )
    elif loads_df.phases[i] == 'B':
        pp.create_asymmetric_load(net, loads_df.Bus[i], 0, p_load,0,0, q_load,0)
    elif loads_df.phases[i] == 'C':
        pp.create_asymmetric_load(net, loads_df.Bus[i], 0,0,p_load, 0,0,q_load)


# add column 'Profile_number' to loads dataframe
loads_df['Profile_number'] = 'NaN'
num_loads = loads_df.shape[0]

# create dataframe with all load profiles
load_profiles_df = pd.DataFrame() #empty df that will include all the load profiles
load_profiles_df['time'] = pd.read_csv(f'European_LV_Test_Feeder_v2\European_LV_CSV\Load Profiles\Load_profile_{1}.csv', usecols= [0])
for i in range(num_loads):
    load_profiles_df[f'profile_{i+1}'] = pd.read_csv(f'European_LV_Test_Feeder_v2\European_LV_CSV\Load Profiles\Load_profile_{i+1}.csv', usecols= [1])
    # assing profiles to loads in the loads dataframe
    profile = loads_df.Yearly[i].split('_')[1]
    loads_df['Profile_number'][i] = profile # asign the corresponding load profile to the load i


# --- 3ph power flow considering profiles --- #
evol_res_line = []
evol_p_lines = []
evol_p_trafos = []
evol_losses_k = []

zabc_lines, zabc_trafos = get_zabc(net)
mat_incide = mat_incid(net)
mat_gamma =np.linalg.inv(mat_incide)

for i in range(0,1440,10):
    time = load_profiles_df['time'][i]
    # scaling all loads to the corresponding values at the current time
    for j in range(num_loads):
        net.asymmetric_load.scaling[j] = load_profiles_df[f'profile_{loads_df.Profile_number[j]}'][i]
    print(f'Solving for time:{time}')
    # running power flow
    runpp_3ph(net)
    i_long_b = complex_currents(net, zabc_lines)
    i_shunt = shunt_currents(net, mat_gamma, i_long_b)
    # calculating losses by RLCP and BCDLA
    res_rlcp_lines, res_rlcp_trafos = rlcp(net, zabc_lines, zabc_trafos)
    losses_k =  bcdla_3ph(net, zabc_lines, zabc_trafos, i_long_b, i_shunt, mat_gamma)

    evol_res_line.append(net.res_line_3ph)
    evol_p_lines.append(res_rlcp_lines)
    evol_p_trafos.append(res_rlcp_trafos)
    evol_losses_k.append(losses_k)

    net.res_line_3ph.to_excel(f'IEEE_lv_time/line_pp/line_pp_{i}.xlsx')
    res_rlcp_lines.to_excel(f'IEEE_lv_time/line_rlcp/line_rlcp_{i}.xlsx')
    res_rlcp_trafos.to_excel(f'IEEE_lv_time/trafo_rlcp/trafo_rlcp_{i}.xlsx')
    losses_k.to_excel(f'IEEE_lv_time/node_bcdla/nodes_bcdla_{i}.xlsx')


# fictional loads added just to be plotted
pp.create_loads(net= net, buses = net.asymmetric_load.bus, p_mw = 0, q_mvar = 0)

# network diagram including loads
pplot.simple_plot(net, show_plot=False, bus_size= .4, plot_loads= True, load_size= .8)
plt.show()
