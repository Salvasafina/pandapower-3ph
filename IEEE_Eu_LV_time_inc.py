import pandas as pd
import pandapower as pp
import pandapower.networks as nw
from math import acos, tan
from pandapower.pf.runpp_3ph import runpp_3ph
from functions_pp_3ph import *

#  ---- load IEEE European LV network from Pandapower ---- #
print('Loading network')
net = nw.ieee_european_lv_asymmetric()
# -------------------------------------------------------- #

#  ---- create dataframe with loads nominal values ---- #
loads_df = pd.read_csv('European_LV_Test_Feeder_v2/European_LV_CSV/Loads.csv', header = 2)
loads_df['Profile_number'] = 'NaN'   # add column 'Profile_number' to loads dataframe
num_loads = loads_df.shape[0]
# ----------------------------------------------------- #

#  ----  fix load to nominal values  ---- #
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
# --------------------------------------- #

# ---- create dataframe with all load profiles ---- #
print('Loading load profiles')
load_profiles_df = pd.DataFrame() #empty df that will include all the load profiles
load_profiles_df['time'] = pd.read_csv(f'European_LV_Test_Feeder_v2\European_LV_CSV\Load Profiles\Load_profile_{1}.csv', usecols= [0])
for i in range(num_loads):
    load_profiles_df[f'profile_{i+1}'] = pd.read_csv(f'European_LV_Test_Feeder_v2\European_LV_CSV\Load Profiles\Load_profile_{i+1}.csv', usecols= [1])
    # assing profiles to loads in the loads dataframe
    profile = loads_df.Yearly[i].split('_')[1]
    loads_df['Profile_number'][i] = profile # asign the corresponding load profile to the load i
# ------------------------------------------------- #


# ---- miscelaneus matrices needed for losses methods ---- #
print('Computing needed matrices')
zabc_lines, zabc_trafos = get_zabc(net)
yabc_lines = np.linalg.inv(zabc_lines)
rabc_lines = np.real(zabc_lines)
rabc_trafos = np.real(zabc_trafos)
mat_incide = mat_incid(net)
mat_gamma =np.linalg.inv(mat_incide)
# -------------------------------------------------------- #

# ----------- uncomment to apply load increase of 15% --------------#
net.asymmetric_load.p_c_mw = net.asymmetric_load.p_c_mw*1.15
net.asymmetric_load.q_c_mvar = net.asymmetric_load.q_c_mvar*1.15
# ------------------------------------------------------------------#


# ---- execution of 3ph PF for entire day and post procesing ---- # 
for i in range(1440): # (1440) to consider the entire day
    time = load_profiles_df['time'][i]
    print(f'Solving for time:{time}')
    
    for j in range(num_loads):  # scaling all loads to the corresponding values at the current time
        net.asymmetric_load.scaling[j] = load_profiles_df[f'profile_{loads_df.Profile_number[j]}'][i]
    
    runpp_3ph(net)  # running 3ph power flow

    # ---- computing currents to apply RLCP and BCDLA ---- #
    i_long_b = complex_currents(net, yabc_lines)
    i_shunt = shunt_currents(net, mat_gamma, i_long_b)
    # ---------------------------------------------------- #

    #  ---- computing losses by RLCP and BCDLA ---- #
    res_rlcp_lines, res_rlcp_trafos = rlcp(net, rabc_lines, rabc_trafos,i_long_b) # res_rlcp_line , res_rlcp_trafos are df for the current time
    losses_k =  bcdla_3ph(net, rabc_lines, rabc_trafos, i_long_b, i_shunt, mat_gamma) #losses_k is a df for the current time
    # --------------------------------------------- #

    # ---- saving current time results in csv files ---- #
    net.res_ext_grid_3ph.to_csv(f'IEEE_lv_time/load_increase/pp/ext_grid_pp/ext_grid_pp_{i}.csv')
    net.res_bus_3ph.to_csv(f'IEEE_lv_time/load_increase/pp/buses_pp/buses_pp_{i}.csv')
    net.res_line_3ph.to_csv(f'IEEE_lv_time/load_increase/pp/line_pp/line_pp_{i}.csv')
    net.res_trafo_3ph.to_csv(f'IEEE_lv_time/load_increase/pp/trafo_pp/trafo_pp_{i}.csv')
    net.res_asymmetric_load_3ph.to_csv(f'IEEE_lv_time/load_increase/pp/load_pp/load_pp_{i}.csv')

    res_rlcp_lines.to_csv(f'IEEE_lv_time/load_increase/rlcp/line_rlcp/line_rlcp_{i}.csv')
    res_rlcp_trafos.to_csv(f'IEEE_lv_time/load_increase/rlcp/trafo_rlcp/trafo_rlcp_{i}.csv')

    losses_k.to_csv(f'IEEE_lv_time/load_increase/bcdla/node_bcdla/nodes_bcdla_{i}.csv')
    # -------------------------------------------------- #
# --------------------------------------------------------------- #

