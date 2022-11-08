import pandas as pd
# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)
from math import cos,sin,pi, acos, sqrt
import numpy as np
import pandapower as pp
import plotly.graph_objs as go
import time    

# ---- for both methods ---- # 
def get_zabc(net):
    zabc_lines = []
    zabc_trafos = []
    for trafo_index in net.trafo.index:
        z = (net.trafo.vn_lv_kv[trafo_index]**2/net.trafo.sn_mva[trafo_index])*(net.trafo.vkr_percent[trafo_index]/100+net.trafo.vk_percent[trafo_index]/100*1j)
        #z = 0
        zabc = np.array([[z,0,0],[0,z,0],[0,0,z]])
        zabc_trafos.append(zabc)
        
    # Fortecue transformation matrix
    A = np.array([  [1,1,1],
                [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j],
                [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j]])
    invA = 1/3*np.array([[1,1,1], 
                    [1, cos(2*pi/3)+ sin(2*pi/3)*1j, cos(4*pi/3)+sin(4*pi/3)*1j],
                    [1, cos(4*pi/3)+ sin(4*pi/3)*1j, cos(2*pi/3)+sin(2*pi/3)*1j]])
    
    for line_index in net.line.index:    
        # Construccion of sequence impedance matrix from net.line
        z012 = np.array([[net.line.r0_ohm_per_km[line_index] + net.line.x0_ohm_per_km[line_index]*1j,0, 0],
        [0, net.line.r_ohm_per_km[line_index] + net.line.x_ohm_per_km[line_index]*1j, 0],
        [0,0,net.line.r_ohm_per_km[line_index] + net.line.x_ohm_per_km[line_index]*1j]])
        # Transformation to line impedance matrix
        zabc = A@z012@invA
        zabc = zabc*net.line.length_km[line_index]
        zabc_lines.append(zabc)

    return zabc_lines, zabc_trafos

def complex_voltages(net):
    v_comp = pd.DataFrame()
    v_comp['v_com_a'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_a_pu * ((net.res_bus_3ph.va_a_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_a_degree*pi/180).apply(sin)*1j)
    v_comp['v_com_b'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_b_pu * ((net.res_bus_3ph.va_b_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_b_degree*pi/180).apply(sin)*1j)
    v_comp['v_com_c'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_c_pu * ((net.res_bus_3ph.va_c_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_c_degree*pi/180).apply(sin)*1j)
    v_comp = v_comp
    return v_comp

def complex_currents(net, yabc_lines): #linea por linea
    v_com = complex_voltages(net)
    num_lines = yabc_lines.shape[0]
    diff_v = v_com.loc[net.line.from_bus].reset_index(drop=True) - v_com.loc[net.line.to_bus].reset_index(drop=True)
    diff_v_array = diff_v.to_numpy().reshape((num_lines,3,1))
    iabc_array = yabc_lines@diff_v_array
    iabc_array = iabc_array.reshape((num_lines,3))
    i_com = pd.DataFrame(iabc_array, columns=['i_com_a','i_com_b','i_com_c'])
    return i_com

# ---- for RLCP ---- #
def rlcp(net, rabc_lines, rabc_trafos, i_long_b): #todo el sistema
    num_lines = i_long_b.shape[0]
    p_abc_trafos = pd.DataFrame(columns=['losses_a', 'losses_b','losses_c'])
    for i in net.trafo.index:
        rabc = rabc_trafos[i]
        
        currents_lv_trafo = net.line[net.line.from_bus == net.trafo.lv_bus[i]]
        total_current_lv = i_long_b.loc[currents_lv_trafo.index].sum()
        iabc = np.array([[total_current_lv.i_com_a],[total_current_lv.i_com_b],[total_current_lv.i_com_c]])
        
        pabc = np.real(iabc*(rabc@np.conj(iabc)))
        p_abc_trafos.loc[i] = [pabc[0][0], pabc[1][0], pabc[2][0]]

    i_long_b_array = i_long_b.to_numpy().reshape(num_lines,3,1)
    i_long_conj = np.conj(i_long_b_array).reshape(num_lines,3,1)

    p_abc_array = np.real(i_long_b_array*rabc_lines@i_long_conj).reshape((num_lines,3))
    p_abc_lines = pd.DataFrame(p_abc_array, columns=['losses_a', 'losses_b','losses_c'] )


    return p_abc_lines, p_abc_trafos

# ---- for BCDLA ---- #
def mat_incid(net):
    """
    Costruisce la matrice delle incidenze, solo considera trasformatori e linee come rami
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

    i_long_pha = np.array(i_long.iloc[:, 0])
    i_long_phb = np.array(i_long.iloc[:, 1])
    i_long_phc = np.array(i_long.iloc[:, 2])
    i_sha = np.dot(aux,i_long_pha)
    i_shb = np.dot(aux,i_long_phb)
    i_shc = np.dot(aux,i_long_phc)
    i_shunt.iloc[:,0] = i_sha
    i_shunt.iloc[:,1] = i_shb
    i_shunt.iloc[:,2] = i_shc
    
    return i_shunt

def bcdla_3ph(net, rabc_lines, rabc_trafos, i_long_b, i_shunt, gamma):
    num_lines = i_long_b.shape[0]
    num_nodes = gamma.shape[0]
    # zeros array for rxi product for each trafo
    rxi_trafos = np.zeros((len(net.trafo.index),3,1), dtype=complex)

    for b in net.trafo.index: 
        # get matrix Rabc from Zabc of the trafo
        rabc = rabc_trafos[b]
        # lv current of the trafo
        currents_lv_trafo = net.line[net.line.from_bus == net.trafo.lv_bus[b]]
        total_current_lv = i_long_b.loc[currents_lv_trafo.index].sum()
        i_lv_tr = np.array([[total_current_lv.i_com_a],[total_current_lv.i_com_b],[total_current_lv.i_com_c]]).conjugate()
        # r x i product for the trafo
        rxi_b = rabc@i_lv_tr
        # add vector to list
        rxi_trafos[b] = rxi_b

    i_long_br = i_long_b.to_numpy().conjugate().reshape(num_lines,3,1)
    rxi_lines = rabc_lines@i_long_br

    rxi_all = np.concatenate((rxi_trafos, rxi_lines), axis=0).reshape(num_nodes,3)
    aux = gamma@rxi_all
    aux = aux.reshape(num_nodes,3,1)

    i_shunt_array = i_shunt.to_numpy().reshape(num_nodes,3,1)
    losses_k_array = (i_shunt_array*aux).real.reshape(num_nodes,3)
    losses_k = pd.DataFrame(losses_k_array, columns=['losses_a', 'losses_b', 'losses_c'] )

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

def profilo_v(net, mat_gamma):
    aux = np.sum(mat_gamma, axis=0)
    nodi_terminal = np.where(aux == -1)
    percorsi_termiali = mat_gamma[nodi_terminal[0]]

    figurea = go.Figure()
    figureb = go.Figure()
    figurec = go.Figure()

    for i in range(percorsi_termiali.shape[0]):# 
        nodi_percorso = []
        distanze_percorso = [0]
        num_ramo = 0
        for j in percorsi_termiali[i][1:]:
            if j == -1:
                nodo_partenza = net.line.from_bus[num_ramo]
                lung_ramo = net.line.length_km[num_ramo]*1000
                nodi_percorso.append(nodo_partenza)
                distanze_percorso.append(lung_ramo + distanze_percorso[-1])
                num_ramo = num_ramo + 1
            else:
                num_ramo = num_ramo +1
        nodi_percorso.append(nodi_terminal[0][i]+1)
        
        profilo_a = pd.DataFrame()
        profilo_a['tensione'] = net.res_bus_3ph.loc[nodi_percorso, 'vm_a_pu']
        profilo_a['distanza'] = distanze_percorso
        profilo_a['nodo'] = nodi_percorso
        profilo_b = pd.DataFrame()
        profilo_b['tensione'] = net.res_bus_3ph.loc[nodi_percorso, 'vm_b_pu']
        profilo_b['distanza'] = distanze_percorso
        profilo_b['nodo'] = nodi_percorso
        profilo_c = pd.DataFrame()
        profilo_c['tensione'] = net.res_bus_3ph.loc[nodi_percorso, 'vm_c_pu']
        profilo_c['distanza'] = distanze_percorso
        profilo_c['nodo'] = nodi_percorso

        trace_a =  go.Scatter(x=profilo_a.distanza, y=profilo_a.tensione, mode='lines',
                name = f'Nodo terminale {nodi_terminal[0][i]+1}', text=profilo_a.nodo)
        trace_b =  go.Scatter(x=profilo_b.distanza, y=profilo_b.tensione, mode='lines',
                name = f'Nodo terminale {nodi_terminal[0][i]+1}', text=profilo_b.nodo)
        trace_c =  go.Scatter(x=profilo_c.distanza, y=profilo_c.tensione, mode='lines',
                name = f'Nodo terminale {nodi_terminal[0][i]+1}', text=profilo_c.nodo)
        
        figurea.add_trace(trace_a)
        figureb.add_trace(trace_b)
        figurec.add_trace(trace_c)

    figurea.update_traces(hovertemplate='Nodo: %{text} <br>Tensione pu: %{y} <br>Distanza [m]: %{x}')
    figurea.update_layout(title="Profilo di tensione fase A", xaxis_title="Distanza [m]", yaxis_title="Tensione pu", legend_title="Nodo terminale")
    figureb.update_traces(hovertemplate='Nodo: %{text} <br>Tensione pu: %{y} <br>Distanza [m]: %{x}')
    figureb.update_layout(title="Profilo di tensione fase B", xaxis_title="Distanza [m]", yaxis_title="Tensione pu", legend_title="Nodo terminale")
    figurec.update_traces(hovertemplate='Nodo: %{text} <br>Tensione pu: %{y} <br>Distanza [m]: %{x}')
    figurec.update_layout(title="Profilo di tensione fase C", xaxis_title="Distanza [m]", yaxis_title="Tensione pu", legend_title="Nodo terminale")

    figurea.show()
    figureb.show()
    figurec.show()

    return figurea

# ------- PV INCORPORATION ------ # 
def adjust_profile(profile_dir, time_space):
    # --- read file --- #
    profile_df = pd.read_excel(profile_dir, index_col = None)

    for column in profile_df.columns:
        name_check = column[0:8]
        if name_check == 'Unnamed:':
            profile_df = profile_df.drop([column], axis = 1)
    # ----------------- #
    
    # --- convertion of times to seconds --- #
    profile_df['secondo'] = 'NaN'
    for i in profile_df.index:
        profile_df['secondo'][i] = profile_df.ora[i].hour*3600 + profile_df.ora[i].minute*60 +profile_df.ora[i].second
    profile_df = profile_df.drop('ora', axis = 1)
    # -------------------------------------- #
    
    # --- completes original dataframe with 0 generation for night hours --- #
    num_data = len(profile_df)
    first_row = pd.DataFrame({'PAt[kW]': 0, 'secondo': 0}, index =[0])

    aux = (profile_df.secondo[0]//(time_space*60))*time_space*60
    second_row = pd.DataFrame({'PAt[kW]': 0, 'secondo': aux}, index =[0])

    aux = (profile_df.secondo[num_data-1]//(time_space*60)+1)*time_space*60
    second_tl_row = pd.DataFrame({'PAt[kW]': 0, 'secondo': aux}, index =[0])

    last_row = pd.DataFrame({'PAt[kW]': 0, 'secondo': 86400}, index =[0])

    profile_df = pd.concat([first_row, second_row, profile_df, second_tl_row,last_row]).reset_index(drop = True)
    # ----------------------------------------------------------------------- #

    # ---- interpolation ---- #
    adjusted_profile_df = pd.DataFrame(columns= profile_df.columns)
    adjusted_profile_df.loc[0] = [0,0]

    list_minutes = []
    list_power = []

    for i in range(1440//time_space+1):
        des_time_seconds = time_space*60*i
        list_minutes.append(des_time_seconds)
        start_range = time_space*60*i
        end_range = time_space*60*(i+1)

        closest = pd.DataFrame(np.abs(profile_df.secondo - des_time_seconds)).min(axis=1).idxmin()
        closest_time = profile_df.secondo[closest]

        time1 = 0
        time2 = 0
        power1 = 0
        power2 = 0
        dif_power = 0
        dif_time = 0

        if closest_time > des_time_seconds:
            time2 = closest_time
            time1 = profile_df.secondo[closest-1]
            power1 = profile_df['PAt[kW]'][closest-1]
            power2 = profile_df['PAt[kW]'][closest]
            dif_time = time2 - time1
            dif_power = power2 - power1

        if closest_time < des_time_seconds:
            time1 = closest_time
            time2 = profile_df.secondo[closest+1]
            power1 = profile_df['PAt[kW]'][closest]
            power2 = profile_df['PAt[kW]'][closest+1]
            dif_time = time2 - time1
            dif_power = power2 - power1

        rango = profile_df[(profile_df['secondo'] >= start_range) & (profile_df['secondo'] <= end_range)]

        if len(rango > 0):
            power_value = rango['PAt[kW]'].mean()
            list_power.append(power_value)
        else:
            power_value = dif_power/dif_time*(des_time_seconds - time1) + power1
            list_power.append(power_value)


        adjusted_profile_df.loc[i] =[power_value,des_time_seconds]
    # -------------------------------------------------------------- #

    # ---- vertical escalation ---- # 
    adjusted_profile_df['PAt[kW]'] = adjusted_profile_df['PAt[kW]']/adjusted_profile_df['PAt[kW]'].max()
    # ----------------------------- #

    # ---- recover hour as time indicator ---- # 
    adjusted_profile_df['ora'] = 'NaN'
    for i in adjusted_profile_df.index:
        adjusted_profile_df['ora'][i] = time.strftime("%H:%M:%S", time.gmtime(adjusted_profile_df.secondo[i]))    
    adjusted_profile_df = adjusted_profile_df.drop([0],axis=0)
    adjusted_profile_df.reset_index(drop=True, inplace=True)
    # ---------------------------------------- #
    
    
    return adjusted_profile_df
    