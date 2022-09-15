import pandapower.networks as nw
from functions_pp_3ph import *
from pandapower.pf.runpp_3ph import runpp_3ph

'''
Applicazione dei metodi RLCP e BCDLA per il calcolo e allocazione delle perdite sulla rete IEEE European Low
Voltage Test Feeder, al momento di picco di domanda alle 09:26:00.

I risultati calcolati da Pandapower, RLCP, BCDLA si salvano in file excel 

Perdite totali calcolate con RLCP e BCDLA coincidono, ma sono diverse da quelle ottenute da Pandapower dovuto
alla forma di calcolare le perdite nel trasformatore che alimenta il sistema. 
'''

# load IEEE European LV network from Pandapower
net = nw.ieee_european_lv_asymmetric()
net.asymmetric_load.q_b_mvar = net.asymmetric_load.q_b_mvar*1000
net.asymmetric_load.q_c_mvar = net.asymmetric_load.q_c_mvar*1000

# Preliminari per RLCP e BCDLAs
zabc_lines, zabc_trafos = get_zabc(net)
mat_incidenze = mat_incid(net)
mat_gamma =np.linalg.inv(mat_incidenze)

# Flusso di potenza trifase e calcolo di correnti shunt
runpp_3ph(net)
i_long_b = complex_currents(net, zabc_lines)
i_shunt = shunt_currents(net, mat_gamma, i_long_b)

# Perdite calcolate con RLCP
res_rlcp_lines, res_rlcp_trafos = rlcp(net, zabc_lines, zabc_trafos)
res_rlcp_lines.to_excel('IEEE_lv_peak/reslines566rlcp.xlsx')
res_rlcp_trafos.to_excel('IEEE_lv_peak/restrafos566rlcp.xlsx')
total_losses_rlcp = res_rlcp_lines.sum().sum() + res_rlcp_trafos.sum().sum()

print('Risultati perdite RLCP\nLines')
print(f'Total losses by phase at time 09:26:00: \n{res_rlcp_lines.sum()}')
print(f'Total line losses at time 09:26:00: \n{res_rlcp_lines.sum().sum()}\n')
print('Trafos')
print(f'Total losses by phase at time 09:26:00: \n{res_rlcp_trafos.sum()}')
print(f'Total trafo losses at time 09:26:00: \n{res_rlcp_trafos.sum().sum()}')
print(f'Total losses at time 09:26:00: \n{total_losses_rlcp}')


# Perdite calcolate con BCDLA
res_bcdla =  bcdla_3ph(net, zabc_lines, zabc_trafos, i_long_b, i_shunt, mat_gamma)
res_bcdla.to_excel('IEEE_lv_peak/ressist566bcdla.xlsx')

print('Risultati perdite BCDLA')
print(f'Losses allocated to nodes, by phase, at time 09:26:00: \n{res_bcdla}')
print(f'Total losses at time 09:26:00: \n{res_bcdla.sum().sum()}')
