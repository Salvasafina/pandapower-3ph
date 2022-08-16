import pandas as pd
from math import cos,sin,pi, sqrt
import numpy as np

def get_zabc(net):
    zabc_lines = []
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
    return zabc_lines

def complex_voltages(net):
    v_comp = pd.DataFrame()
    v_comp['v_com_a'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_a_pu * ((net.res_bus_3ph.va_a_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_a_degree*pi/180).apply(sin)*1j)
    v_comp['v_com_b'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_b_pu * ((net.res_bus_3ph.va_b_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_b_degree*pi/180).apply(sin)*1j)
    v_comp['v_com_c'] = net.bus.vn_kv*1000/sqrt(3)*net.res_bus_3ph.vm_c_pu * ((net.res_bus_3ph.va_c_degree*pi/180).apply(cos) + (net.res_bus_3ph.va_c_degree*pi/180).apply(sin)*1j)
    v_comp = v_comp
    return v_comp

def complex_currents(net, impedances): #linea por linea 
    v_com = complex_voltages(net)
    i_com = pd.DataFrame(columns=['i_com_a','i_com_b','i_com_c'])

    for i in net.line.index:

        zabc = impedances[i]*net.line.length_km[i]

        v_from_p = v_com.loc[net.line.from_bus[i]]
        v_from = np.array([[v_from_p[0]],[v_from_p[1]],[v_from_p[2]]])

        v_to_p = v_com.loc[net.line.to_bus[i]]
        v_to = np.array([[v_to_p[0]],[v_to_p[1]],[v_to_p[2]]])

        iabc = np.linalg.inv(zabc)@(v_from-v_to)
        i_com.loc[i] = [iabc[0][0], iabc[1][0], iabc[2][0]]
    return i_com

def rlcp(net, impedances): #todo el sistema
    i_com = complex_currents(net, impedances)
    p_abc = pd.DataFrame(columns=['losses_a', 'losses_b','losses_c'])
    for i in net.line.index:
        zabc = impedances[i]*net.line.length_km[i]
        rabc = np.real(zabc) 
        iabc = np.array([[i_com.i_com_a[i]],[i_com.i_com_b[i]],[i_com.i_com_c[i]]])
        pabc = np.real(iabc*(rabc@np.conj(iabc)))
        p_abc.loc[i] = [pabc[0][0], pabc[1][0], pabc[2][0]]
    return p_abc
