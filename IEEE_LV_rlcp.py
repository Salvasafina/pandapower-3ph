import pandapower.networks as nw
from pandapower.pf.runpp_3ph import runpp_3ph
from functions_pp_3ph import *

net = nw.ieee_european_lv_asymmetric()
#net.asymmetric_load.q_b_mvar = net.asymmetric_load.q_b_mvar*1000
#net.asymmetric_load.q_c_mvar = net.asymmetric_load.q_c_mvar*1000
runpp_3ph(net)


impedances = get_zabc(net)
p_rlcp = rlcp(net, impedances)
p_tot = p_rlcp.sum()
print(p_tot)
print(p_tot.sum())
