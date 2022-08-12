from math import cos,sin,acos,pi, sqrt
import numpy as np
import pandapower as pp
import pandapower.plotting as pplot
import matplotlib.pyplot as plt

#  ---- Preliminari e sistemazione dati ----  #
#Matrice per trasformata Fortescue
A = np.array([  [1,1,1],
                [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j],
                [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j]])
invA = 1/3*np.array([[1,1,1],
                    [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j],
                    [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j]])

# Matrice di impedenze di fase -dopo Kron-  [ohm/mile]
Zabc = np.array([[.3375 + 1.0478j, .1535+.3849j, .1559+.5017j],[.1535+.3849j, .3414 +1.0348j, .1580+.4236j],[.1559+.5017j, .1580+.4236j, .3465+1.0179j]])
Zabc = 3*Zabc #linea de 3 millas

# Sistemo Zabc per considerare la linea trasposta
Zabc[0,0] = Zabc[1,1] = Zabc[2,2] = (Zabc[0,0] + Zabc[1,1] + Zabc[2,2])/3
Zabc[0,1] = Zabc[0,2] = Zabc[1,2] = (Zabc[0,1] + Zabc[0,2] + Zabc[1,2])/3
Zabc[1,0] = Zabc[2,0] = Zabc[2,1] = Zabc[0,1]
print(f'Zabc\n {Zabc}\n')

#Matrice delle impedenze in sequenza 012
Z012 = invA@Zabc@A
print(f'Z012\n {Z012}\n')

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

#  ---- Costruzione rete con Pandapower  ----  #
#Costruzione rete
# Rete vuota
net = pp.create_empty_network(f_hz=60)
# Aggiungo y nodi
pp.create_bus(net, 12.47, 'bus0')
pp.create_bus(net, 12.47, 'bus1')
# Rete esterna / nodo slack
pp.create_ext_grid(net, 0, 1, 0, s_sc_max_mva=1000000000, rx_max = 0.1, x0x_max = 0.1, r0x0_max=0.1) 
 # Aggiungo la linea (Si crea una linea di 1km con dati corrispondenti a una linea di 3 miles)
pp.create_line_from_parameters(net, 0, 1, 1, r1, x1, c1, 1000000, name = 'linea01'
    , r0_ohm_per_km= r0, x0_ohm_per_km= x0, c0_nf_per_km= c0)
# Aggiungo gli interruttori
pp.create_switch(net, bus = 0, element=0, et= 'l')
pp.create_switch(net, bus = 1, element=0, et= 'l')
# Aggiungo il carico trifase
pp.create_asymmetric_load(net=net, bus = 1, p_a_mw=pa, q_a_mvar=qa, sn_mva=1.5)
pp.create_asymmetric_load(net=net, bus = 1, p_b_mw=pb, q_b_mvar=qb, sn_mva=1)
pp.create_asymmetric_load(net=net, bus = 1, p_c_mw=pc, q_c_mvar=qc, sn_mva=2)

# informazione della rete
print(net)
# ---  Esegue flusso di potenza trifase e salva i risultati --- # 
# Flusso di potenza (per trifase solo si ha Newton-Rhapson)
pp.pf.runpp_3ph.runpp_3ph(net, calculate_voltage_angles = True)

# Risultati del flusso di potenza si salvano in file excel
net.res_line_3ph.to_excel('reslinea.xlsx')
net.res_bus_3ph.to_excel('resbus.xlsx')
net.res_asymmetric_load_3ph.to_excel('resload.xlsx')

print(net.res_bus_3ph)

#  ---- Metodo RLCP per ripartizione delle perdite fra le fasi ----  #
Rabc = np.real(Zabc) #[ohm]
print(f'Rabc\n {Rabc}\n')

# Tensioni ai nodi in cartesiano [V]
V= 12470/sqrt(3)
V0 = V*np.array([[net.res_bus_3ph.vm_a_pu[0]*cos(net.res_bus_3ph.va_a_degree[0]*pi/180)  +  net.res_bus_3ph.vm_a_pu[0]*sin(net.res_bus_3ph.va_a_degree[0]*pi/180)*1j],
                 [net.res_bus_3ph.vm_b_pu[0]*cos(net.res_bus_3ph.va_b_degree[0]*pi/180)  +  net.res_bus_3ph.vm_b_pu[0]*sin(net.res_bus_3ph.va_b_degree[0]*pi/180)*1j],
                 [net.res_bus_3ph.vm_c_pu[0]*cos(net.res_bus_3ph.va_c_degree[0]*pi/180)  +  net.res_bus_3ph.vm_c_pu[0]*sin(net.res_bus_3ph.va_c_degree[0]*pi/180)*1j]])
V1 = V*np.array([[net.res_bus_3ph.vm_a_pu[1]*cos(net.res_bus_3ph.va_a_degree[1]*pi/180)  +  net.res_bus_3ph.vm_a_pu[1]*sin(net.res_bus_3ph.va_a_degree[1]*pi/180)*1j],
                 [net.res_bus_3ph.vm_b_pu[1]*cos(net.res_bus_3ph.va_b_degree[1]*pi/180)  +  net.res_bus_3ph.vm_b_pu[1]*sin(net.res_bus_3ph.va_b_degree[1]*pi/180)*1j],
                 [net.res_bus_3ph.vm_c_pu[1]*cos(net.res_bus_3ph.va_c_degree[1]*pi/180)  +  net.res_bus_3ph.vm_c_pu[1]*sin(net.res_bus_3ph.va_c_degree[1]*pi/180)*1j]])

# Ricalcolo le correnti per averli in cartesiano[A]
Iabc = np.linalg.inv(Zabc)@(V0-V1)

# Calcolo delle perdite per fase tramite RLCP [W]
Pabc = np.real(Iabc*(Rabc@np.conj(Iabc)))
print(f'Pabc\n {Pabc}\n')
print(f'Perdite totali: {sum(Pabc)}')

#  ---- Schema della rete ----  #
ax = pplot.simple_plot(net, show_plot=False)
sc = pplot.create_line_switch_collection(net, size = 0.01, distance_to_bus = 0.05)
pplot.draw_collections([sc], ax = ax)
plt.show()
