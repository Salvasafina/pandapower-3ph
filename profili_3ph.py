from math import cos,sin,acos,pi
import numpy as np
import pandapower as pp
from funzioni_rete import *
import matplotlib.pyplot as plt
from pandapower.pf.runpp_3ph import runpp_3ph

#  ---- Preliminari e sistemazione dati ----  #
# Matrice per trasformata Fortescue
A = np.array([  [1,1,1],
                [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j],
                [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j]])
invA = 1/3*np.array([[1,1,1],
                    [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j],
                    [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j]])

# Matrice di impedenze di fase -dopo Kron-  [ohm/mi]
Zabc = np.array([[.3375 + 1.0478j, .1535+.3849j, .1559+.5017j],[.1535+.3849j, .3414 +1.0348j, .1580+.4236j],[.1559+.5017j, .1580+.4236j, .3465+1.0179j]])
Zabc = 3*Zabc #linea lung 3 mi

# Sistemo Zabc per considerare la linea trasposta
Zabc[0,0] = Zabc[1,1] = Zabc[2,2] = (Zabc[0,0] + Zabc[1,1] + Zabc[2,2])/3
Zabc[0,1] = Zabc[0,2] = Zabc[1,2] = (Zabc[0,1] + Zabc[0,2] + Zabc[1,2])/3
Zabc[1,0] = Zabc[2,0] = Zabc[2,1] = Zabc[0,1]
#print(f'Zabc\n {Zabc}\n')

#Matrice delle impedenze in sequenza 012
Z012 = invA@Zabc@A
#print(f'Z012\n {Z012}\n')

# Resistenza e reattanza in sequenza 012 (Z1 = Z2)
r0 = np.real(Z012[0,0])
x0 = np.imag(Z012[0,0])
c0 = 0
r1 = np.real(Z012[1,1])
x1 = np.imag(Z012[1,1])
c1 = 0

# Potenza del carico trifase (MW y Mvar)
pa = 1.5*0.88
pb = 1*0.95
pc = 2*.8
qa = 1.5*sin(acos(.88))
qb = 1*sin(acos(0.95))
qc = 2*sin(acos(.8))

#  ---- Costruzione rete e calcolo del flusso di carico con Pandapower  ----  #
#Costruzione rete
net = pp.create_empty_network(f_hz=60)      # Rete vuota
pp.create_bus(net, 12.47, 'bus0')           # Aggiungo y nodi
pp.create_bus(net, 12.47, 'bus1')
pp.create_ext_grid(net, 0, 1, 0, s_sc_max_mva=1000000000, rx_max = 0.1, x0x_max = 0.1, r0x0_max=0.1) # Rete esterna / nodo slack
pp.create_line_from_parameters(net, 0, 1, 1, r1, x1, c1, 1000000, name = 'linea01'
    , r0_ohm_per_km= r0, x0_ohm_per_km= x0, c0_nf_per_km= c0) # Aggiungo la linea (Si crea una linea di 1km con dati corrispondenti a una linea di 3 miglia)
pp.create_switch(net, bus = 0, element=0, et= 'l') # Aggiungo gli interruttori
pp.create_switch(net, bus = 1, element=0, et= 'l')
pp.create_asymmetric_load(net=net, bus = 1, p_a_mw=pa, q_a_mvar=qa, sn_mva=1.5) # Aggiungo il carico trifase
pp.create_asymmetric_load(net=net, bus = 1, p_b_mw=pb, q_b_mvar=qb, sn_mva=1)
pp.create_asymmetric_load(net=net, bus = 1, p_c_mw=pc, q_c_mvar=qc, sn_mva=2)

print(net)


# --- Carica profili di carico --- #
profili_df = carica_profili_df('ReteRurale_profili.mat')

# --- Soluzione del sistema ad ogni intervallo di tempo  --- #
tre_df = net.res_bus_3ph.copy()
for i in range(96):
    net.asymmetric_load.scaling = profili_df.ind[i]
    runpp_3ph(net)
    tre_df = tre_df.append(net.res_bus_3ph.loc[1],ignore_index=True)

# --- Grafico --- #
evol_tens_df = tre_df[['vm_a_pu','vm_b_pu','vm_c_pu']]
plt.plot(evol_tens_df)
plt.title('Evoluzione tensioni per fase')
plt.xlabel('Time step')
plt.legend(evol_tens_df)
plt.show()
