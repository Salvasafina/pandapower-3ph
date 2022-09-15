import pandas as pd
#pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
from math import cos,sin,pi, acos, sqrt
import numpy as np
import pandapower as pp

def get_zabc(net):
    zabc_lines = []
    zabc_trafos = []
    for trafo_index in net.trafo.index:
        z = (net.trafo.vn_lv_kv[trafo_index]**2/net.trafo.sn_mva[trafo_index])*(net.trafo.vkr_percent[trafo_index]/100+net.trafo.vk_percent[trafo_index]/100*1j)
        #z = 0
        zabc = np.array([[z,0,0],[0,z,0],[0,0,z]])
        zabc_trafos.append(zabc)

    for line_index in net.line.index:
        # Fortecue transformation matrix
        A = np.array([  [1,1,1],
                    [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j],
                    [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j]])
        invA = 1/3*np.array([[1,1,1], 
                        [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j],
                        [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j]])
        # Construccion of sequence impedance matrix from net.line
        z012 = np.array([[net.line.r0_ohm_per_km[line_index] + net.line.x0_ohm_per_km[line_index]*1j,0, 0],
        [0, net.line.r_ohm_per_km[line_index] + net.line.x_ohm_per_km[line_index]*1j, 0],
        [0,0,net.line.r_ohm_per_km[line_index] + net.line.x_ohm_per_km[line_index]*1j]])
        # Transformation to line impedance matrix
        zabc = A@z012@invA
        zabc_lines.append(zabc)

    return zabc_lines, zabc_trafos

def complex_voltages(net):
    v_comp = pd.DataFrame()
    v_comp['v_com_a'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_a_pu * ((net.res_bus_3ph.va_a_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_a_degree*pi/180).apply(sin)*1j)
    v_comp['v_com_b'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_b_pu * ((net.res_bus_3ph.va_b_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_b_degree*pi/180).apply(sin)*1j)
    v_comp['v_com_c'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_c_pu * ((net.res_bus_3ph.va_c_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_c_degree*pi/180).apply(sin)*1j)
    v_comp = v_comp
    return v_comp

def complex_currents(net, zabc_lines): #linea por linea
    v_com = complex_voltages(net)
    i_com = pd.DataFrame(columns=['i_com_a','i_com_b','i_com_c'])

    for i in net.line.index:
        zabc = zabc_lines[i]*net.line.length_km[i]
        
        v_from_p = v_com.loc[net.line.from_bus[i]]
        v_from = np.array([[v_from_p[0]],[v_from_p[1]],[v_from_p[2]]])
        
        v_to_p = v_com.loc[net.line.to_bus[i]]
        v_to = np.array([[v_to_p[0]],[v_to_p[1]],[v_to_p[2]]])
        
        iabc = np.linalg.inv(zabc)@(v_from-v_to)
        
        i_com.loc[i] = [iabc[0][0], iabc[1][0], iabc[2][0]]
    return i_com

def rlcp(net, zabc_lines, zabc_trafos): #todo el sistema
    i_com = complex_currents(net, zabc_lines)
    p_abc_lines = pd.DataFrame(columns=['losses_a', 'losses_b','losses_c'])
    p_abc_trafos = pd.DataFrame(columns=['losses_a', 'losses_b','losses_c'])

    for i in net.trafo.index:
        zabc = zabc_trafos[i]
        rabc = np.real(zabc)
        
        currents_lv_trafo = net.line[net.line.from_bus == net.trafo.lv_bus[i]]
        total_current_lv = i_com.loc[currents_lv_trafo.index].sum()
        iabc = np.array([[total_current_lv.i_com_a],[total_current_lv.i_com_b],[total_current_lv.i_com_c]])

        
        pabc = np.real(iabc*(rabc@np.conj(iabc)))
        p_abc_trafos.loc[i] = [pabc[0][0], pabc[1][0], pabc[2][0]]


    for i in net.line.index:
        zabc = zabc_lines[i]*net.line.length_km[i]
        rabc = np.real(zabc) 

        iabc = np.array([[i_com.i_com_a[i]],[i_com.i_com_b[i]],[i_com.i_com_c[i]]])

        pabc = np.real(iabc*(rabc@np.conj(iabc)))
        p_abc_lines.loc[i] = [pabc[0][0], pabc[1][0], pabc[2][0]]

    return p_abc_lines, p_abc_trafos

# ---- for BCDLA ---- #
def mat_incid(net):
    """
    Costruisce la matrice delle incidenze, considera i trasformatori come rami
    """
    rami_chiusi = net.line[net.line.in_service == True]
    rami_chiusi_trx = net.trafo[net.trafo.in_service == True]
    rami_aperti = net.line[net.line.in_service == False]

    rami_chiusi = rami_chiusi.reset_index(drop=True)
    rami_chiusi_trx = rami_chiusi_trx.reset_index(drop=True)
    net.line = pd.concat([rami_chiusi, rami_aperti], ignore_index= True)

    dim_L = len(rami_chiusi)+len(rami_chiusi_trx)
    mat_incid = np.zeros((dim_L,dim_L))
    cont = 0
    for i in rami_chiusi_trx.index:
        if rami_chiusi_trx['hv_bus'][i] != 0:
            mat_incid[cont][rami_chiusi_trx['hv_bus'][i]-1] = 1
        mat_incid[cont][rami_chiusi_trx['lv_bus'][i]-1] = -1
        cont += 1
    for i in rami_chiusi.index:   # para cada fila del df rami_chiusi
        if rami_chiusi['from_bus'][i] != 0:
            mat_incid[cont][rami_chiusi['from_bus'][i]-1] = 1
        mat_incid[cont][rami_chiusi['to_bus'][i]-1] = -1
        cont += 1
    return mat_incid

def shunt_currents(net, mat_gamma, i_long_b):
    """
    Computes shunt currents of nodes by phase and stores them in a df
    """
    
    i_lv_tr = pd.DataFrame(columns=['i_com_a','i_com_b','i_com_c'])
    for i in net.trafo.index:
        currents_lv_trafo = net.line[net.line.from_bus == net.trafo.lv_bus[i]]
        total_current_lv = i_long_b.loc[currents_lv_trafo.index].sum()
        iabc = np.array([[total_current_lv.i_com_a],[total_current_lv.i_com_b],[total_current_lv.i_com_c]])
        i_lv_tr.loc[i] = [iabc[0][0], iabc[1][0], iabc[2][0]]

    i_long = pd.concat([i_lv_tr, i_long_b], ignore_index= True)

    i_shunt = pd.DataFrame(columns=['i_sh_a','i_sh_b','i_sh_c'])
    g_tr = np.transpose(mat_gamma)
    aux = np.linalg.inv(g_tr)

    for i in range(3):
        i_long_ph = np.array(i_long.iloc[:, i])
        i_sh = np.dot(aux,i_long_ph)
        i_shunt.iloc[:,i] = i_sh
    
    return i_shunt

def bcdla_3ph(net, zabc_lines, zabc_trafos, i_long_b, i_shunt, gamma):
    # empty list for rxi product for each line
    rxi_all = []
    # empty dataframe with losses allocated to nodes
    losses_k = pd.DataFrame(columns=['losses_a', 'losses_b', 'losses_c'])

    for b in net.trafo.index: 
        # get matrix Rabc from Zabc of the trafo
        rabc = np.real(zabc_trafos[b])
        
        # lv current of the trafo
        currents_lv_trafo = net.line[net.line.from_bus == net.trafo.lv_bus[b]]
        total_current_lv = i_long_b.loc[currents_lv_trafo.index].sum()
        i_lv_tr = np.array([[total_current_lv.i_com_a],[total_current_lv.i_com_b],[total_current_lv.i_com_c]]).conjugate()
        
        # r x i product for the line
        rxi_b = rabc@i_lv_tr
        # add vector to list
        rxi_all.append(rxi_b)

    for b in net.line.index: 
        # get matrix Rabc from Zabc of the line
        rabc = np.real(zabc_lines[b]*net.line.length_km[b])
        # longitudinal current of the line
        i_long_br = np.array([[i_long_b.i_com_a[b]],[i_long_b.i_com_b[b]],[i_long_b.i_com_c[b]]]).conjugate()
        # r x i product for the line
        rxi_b = rabc@i_long_br
        # add vector to list
        rxi_all.append(rxi_b)

    for k in range(gamma.shape[0]):  
        # row of gamma matrix to act as a filter
        row_gamma = gamma[k,:]
        aux = []
        # for each node
        for i in range(len(row_gamma)):
            # filters rxi products keeping the ones that belong to the branches of the path from node k to the root
            aux.append(row_gamma[i]*rxi_all[i])
        # sum of rxi products in the path from node k to the root 
        aux = np.array(aux).sum(axis=0)

        # shunt current of the node k
        i_shunt_k = np.array([[i_shunt.i_sh_a[k]],[i_shunt.i_sh_b[k]],[i_shunt.i_sh_c[k]]])
        # product that get the allocated losses to node k by phase
        losses_k_arr = (i_shunt_k*aux).real
        # save allocated losses to node k in the BCDLA losses dataframe
        losses_k.loc[k] = [losses_k_arr[0][0], losses_k_arr[1][0], losses_k_arr[2][0]]
    return losses_k

# ---- 2 nodes network ---- # 
def basic_network_3ph():
        #  ---- Preliminari e sistemazione dati ----  #
    #Matrice per trasformata Fortescue
    A = np.array([  [1,1,1],
                    [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j],
                    [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j]])
    invA = 1/3*np.array([[1,1,1],
                        [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j],
                        [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j]])

    # Matrice di impedenze di fase -dopo Kron-  [ohm/miglio]
    Zabc = np.array([[.3375 + 1.0478j, .1535+.3849j, .1559+.5017j],[.1535+.3849j, .3414 +1.0348j, .1580+.4236j],[.1559+.5017j, .1580+.4236j, .3465+1.0179j]])
    Zabc = 3*Zabc #linea de 3 millas

    # Sistemo Zabc per considerare la linea trasposta
    Zabc[0,0] = Zabc[1,1] = Zabc[2,2] = (Zabc[0,0] + Zabc[1,1] + Zabc[2,2])/3
    Zabc[0,1] = Zabc[0,2] = Zabc[1,2] = (Zabc[0,1] + Zabc[0,2] + Zabc[1,2])/3
    Zabc[1,0] = Zabc[2,0] = Zabc[2,1] = Zabc[0,1]

    # Matrice delle impedenze in sequenza 012
    Z012 = invA@Zabc@A
    print(Z012)

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

    return net

net = basic_network_3ph()