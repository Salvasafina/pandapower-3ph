from functions_pp_3ph import *
from pandapower.pf.runpp_3ph import runpp_3ph
import pandapower.plotting as pplot
import matplotlib.pyplot as plt

'''
Prova flusso di potenza trifase con Pandapower e ripartizione delle perdite tramite metodo RLCP e BCDLA
a partire dei risultati di Pandapower. Si usa l'esempio dell'articolo 'The Computation of Neutral and Dirt
Currents and Power Losses', W.H Kersting.

Risultati diversi: qui il carico trifase al nodo 1 va considerato come una potenza assegnata, 
invece che una impedenza costante.

Pandapower fa il calcolo delle perdite per fase con il metodo classico, assegna valori negativi per fasi a e b,
l'applicazione di RLCP o BCDLA riporta alle stesse perdite totali con diversi risultati per fase, tutti positivi. 
'''

net = basic_network_3ph()
print(net)

# Preliminari per RLCP e BCDLA
zabc_lines, zabc_trafos = get_zabc(net)
mat_incidenze = mat_incid(net)
mat_gamma =np.linalg.inv(mat_incidenze)

# Risultati del flusso di potenza si salvano in file excel
runpp_3ph(net)
i_long_b = complex_currents(net, zabc_lines)
i_shunt = shunt_currents(net, mat_gamma, i_long_b)
net.res_line_3ph.to_excel('rete2nodi/reslinea.xlsx')
net.res_bus_3ph.to_excel('rete2nodi/resbus.xlsx')
net.res_asymmetric_load_3ph.to_excel('rete2nodi/resload.xlsx')

# Perdite calcolate con RLCP
res_rlcp_lines, res_rlcp_trafos = rlcp(net, zabc_lines, zabc_trafos)
print(f'Risultati perdite con RLCP: \n{res_rlcp_lines}')
print(f'Perdite totali: {res_rlcp_lines.sum().sum()}')

# Perdite calcolate con BCDLA
res_bcdla = bcdla_3ph(net, zabc_lines, zabc_trafos, i_long_b, i_shunt, mat_gamma)
print(f'Losses by BCDLA: \n{res_bcdla}')
print(f'Perdite totali: {res_bcdla.sum().sum()}')

#  ---- Schema della rete ----  #
ax = pplot.simple_plot(net, show_plot=False, plot_loads= True, load_size= 1.2)
sc = pplot.create_line_switch_collection(net, size = 0.01, distance_to_bus = 0.05)
pplot.draw_collections([sc], ax = ax)
plt.show()