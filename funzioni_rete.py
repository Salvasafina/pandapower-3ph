import scipy.io
import numpy as np
import pandas as pd
from math import sqrt, log , pi, e
from math import cos, sin, pi
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandapower as pp
import pandapower.plotting.plotly as pplotly

# -- funzioni create per scambio rami ---#
def crea_df_rami(file_mat):
    """
    Crea un DataFrame con l'informazione dei rami inclusi nel file .mat
    """
    info_rete = scipy.io.loadmat(file_mat)  # Si carica il file .mat
    mat_incid = info_rete['L']              # Salva la matrice delle incidenze è un numpy array
    mat_impedenza = info_rete['Zb_pu']      # Salva la matrice delle impedenze longitudinali
    nodo_partenza = []                      
    nodo_arrivo = []
    for ramo in mat_incid:                  # Per ogni ramo (ogni fila della matrice L)
        nodo_arrivo.append(np.where(ramo == -1)[0][0]+1) # Si ottiene il numero del nodo di arrivo
        try: 
            nodo_partenza.append(np.where(ramo == 1)[0][0]+1) # Si ottiene il numero del nodo di partenza
        except IndexError:      # Per evitare l'errore al cercare il nodo di partenza del primo ramo 
            nodo_partenza.append(0)

    resistenza_long = []
    reattanza_long =[]
    cont = 0 # Contatore per sempre stare sulla diagonale
    for ramo in mat_impedenza: # Per ogni ramo (ogni fila della matrice Zb)
        resistenza_long.append(ramo[cont].real) # Parte reale dell elemento della diagonale
        reattanza_long.append(ramo[cont].imag)  # Parte imaginaria dell elemento della diagonale
        cont +=1

    # I dati ricavati si salvano in un DataFrame di pandas
    rami_df = pd.DataFrame(info_rete['Ilim_Ampere'],columns=['I_limite'])
    rami_df['Nodo_Partenza'] = nodo_partenza
    rami_df['Nodo_Arrivo'] = nodo_arrivo
    rami_df['Resistenza_Long_pu'] = resistenza_long
    rami_df['Reattanza_Long_pu'] = reattanza_long
    rami_df['Condut_sh_pu'] = info_rete['Ysh_pu'].real
    rami_df['Sus_sh_pu'] = info_rete['Ysh_pu'].imag
    rami_df['In_servizio'] = [True for i in range(len(resistenza_long))]
    return rami_df

def crea_df_carichi(file_mat):
    """
    Crea un DataFrame con l'informazione dei carichi inclusi nel file .mat
    """
    info_rete = scipy.io.loadmat(file_mat)  # Si carica il file .mat

    # I dati ricavati si salvano in un DataFrame di pandas
    carichi_df = pd.DataFrame(range(1,info_rete['Scar_pu'].shape[0]+1), columns=['Nodo']) # Colonna del nodo a cui e connesso il carico
    carichi_df['P_carico_pu'] = info_rete['Scar_pu'].real
    carichi_df['Q_carico_pu'] = info_rete['Scar_pu'].imag
    carichi_df['Ammet_carico_pu'] = info_rete['Ycar_pu']
    return carichi_df

def dati_df_rete(file_mat):
    """
    Crea un DataFrame con le grandezze base e il numero di nodi della rete
    """
    info_rete = scipy.io.loadmat(file_mat)  # Si carica il file .mat

    rete_df = pd.DataFrame(info_rete['Sbase_kVA'], columns=['Sbase_kVA'])
    rete_df['Vbase_kV'] = info_rete['Vbase_kV']
    rete_df['Ibase_A'] = rete_df['Sbase_kVA'][0]/(sqrt(3)*rete_df['Vbase_kV'][0]).real
    rete_df['Zbase_ohm'] = (rete_df['Vbase_kV'][0]*1000)**2/(rete_df['Sbase_kVA'][0]*1000)
    rete_df['Num_nodi'] = info_rete['NumNodi']
    return rete_df

def aggiunge_ramo_df(rami_df, nuovi_rami):
    """
    Aggiunge dei rami nuovi, non inclusi nel file .mat, al DataFrame dei rami
    """
    for ramo in nuovi_rami:
        rami_df = rami_df.append(ramo, ignore_index = True)
    return rami_df

def rami_pp(rete, rami_df, rete_df):
    """
    Crea i rami della rete di PandaPower
    """
    for i in range(rami_df.shape[0]):
        from_bus = rami_df['Nodo_Partenza'][i]
        to_bus = rami_df['Nodo_Arrivo'][i]
        pp.create_line_from_parameters(rete,
        rami_df['Nodo_Partenza'][i], 
        rami_df['Nodo_Arrivo'][i],
        length_km = 1,
        r_ohm_per_km = rete_df['Zbase_ohm'][0]*rami_df['Resistenza_Long_pu'][i],
        x_ohm_per_km = rete_df['Zbase_ohm'][0]*rami_df['Reattanza_Long_pu'][i],
        max_i_ka = rami_df['I_limite'][i]/1000, 
        c_nf_per_km =  (rami_df['Sus_sh_pu'][i]/rete_df['Zbase_ohm'][0])/(2*pi*50)*10**9, 
        g_us_per_km= 0 , 
        in_service = rami_df['In_servizio'][i],
        name = str(from_bus) + '-' + str(to_bus),
        type= 'cs' )
    return rete.line

def crea_rete_pp (file_mat, nuovi_rami = None):
    """
    Crea una rete di PandaPower a partire del file. mat e i possibili rami adizionali

    Usa le funzioni: crea_rami(), crea_carichi(), dati_rete(), rami_pp(), aggiunge_ramo()
    """
    rami_df = crea_df_rami(file_mat)            # Chiamata alla funzione che crea il DataFrame dei rami
    carichi_df = crea_df_carichi(file_mat)      # Chiamata alla funzione che crea il DataFrame dei carichi
    rete_df = dati_df_rete(file_mat)            # Chiamata alla funzione che crea il DataFrame delle grandezze base

    #  ---- Costruzione della rete in Pandapower ---- #
    # Elemento rete vuoto
    rete = pp.create_empty_network()

    # Creazione dei nodi
    for i in range(rete_df['Num_nodi'][0]+1):
        pp.create_bus(rete, rete_df['Vbase_kV'],type='n', name = 'Bus ' + str(i))

    # Creazione rete esterna che attua come nodo slack
    pp.create_ext_grid(rete, bus= 0, vm_pu = 1, va_degree = 0)

    # Creazione dei rami originali e geodata dei nodi
    rete.line = rami_pp(rete, rami_df, rete_df)
    pplotly.create_generic_coordinates(rete)

    # Creazione dei rami inclusi quelli nuovi
    if nuovi_rami != None:
        rete.line = rete.line.drop(index = rete.line.index)
        rami_df = aggiunge_ramo_df(rami_df, nuovi_rami)
        rete.line = rami_pp(rete, rami_df, rete_df)

    # Creazione dei carichi
    for i in range(carichi_df.shape[0]):
        pp.create_load(rete, 
        carichi_df['Nodo'][i], 
        rete_df['Sbase_kVA'][0]*carichi_df['P_carico_pu'][i]/1000,
        rete_df['Sbase_kVA'][0]*carichi_df['Q_carico_pu'][i]/1000 )

    return rete

def line_coord(rete_pp):
    """
    Crea il DataFrame net.line_geodata, necessario per il corretto funzionamento della funzione schema(rete)
    """
    coords_rami = []
    for i in range(len(rete_pp.line)):
        x1 =  rete_pp.bus_geodata.loc[rete_pp.line.from_bus[i]]['x']
        y1 =  rete_pp.bus_geodata.loc[rete_pp.line.from_bus[i]]['y']
        x2 =  rete_pp.bus_geodata.loc[rete_pp.line.to_bus[i]]['x']
        y2 =  rete_pp.bus_geodata.loc[rete_pp.line.to_bus[i]]['y']
        coord = (np.array([[x1, y1] , [x2,y2]]))
        coords_rami.append(coord)
    rete_pp.line_geodata['coords'] = coords_rami
    return rete_pp.line_geodata

def schema_int_dopo_pf(rete):
    """
    Si definiscono le impostazioni per lo schema interattivo della rete
    Da usare dopo aver eseguito un power flow per vedere le tensioni in pu ai nodi,
    e il livello di carico dei rami
    """
    line_coord(rete)
    # Schema interattivo dei risultati, si crea un file chiamato 'temp-plot.html' che deve essere aperto con il navigatore
    open_lines = rete.line[rete.line.in_service == False].index
    closed_lines = rete.line[rete.line.in_service == True].index
    lc = pplotly.create_line_trace(rete, closed_lines, color="grey", width= 1.5,
     infofunc=pd.Series(index= rete.line.index, data= 'Ramo ' + rete.line.name + '<br>' + 'Liv carico: ' + rete.res_line.loading_percent.astype(str) + '%'))
    olc = pplotly.create_line_trace(rete, open_lines, color="red", 
     infofunc=pd.Series(index=rete.line.index, data= 'Ramo ' + rete.line.name + '<br>' + 'Ramo aperto'))
    
    bc = pplotly.create_bus_trace(rete, rete.bus.index, size=10, color="orange",
     infofunc=pd.Series(index= rete.bus.index, data = rete.bus.name + '<br>'+ 'Vn: '+rete.bus.vn_kv.astype(str) + 'kV'+ '<br>'+ 'Vpu: ' + rete.res_bus.vm_pu.astype(str)))
    egc = pplotly.create_bus_trace(rete, rete.ext_grid.index, size= 15, patch_type= 'square', color= 'orange', infofunc=pd.Series(index= rete.ext_grid.index, data = rete.bus.name + '<br>'+ 'Vn: '+rete.bus.vn_kv.astype(str) + 'kV'+ '<br>'+ 'Vpu: ' + rete.res_bus.vm_pu.astype(str)) )
    pplotly.draw_traces(bc + egc + lc + olc, auto_open= False)

def schema_int_base(rete):
    """
    Si definiscono le impostazioni per lo schema interattivo della rete, 
    mostra i nomi e tensione nominale dei nodi
    mostra il nome dei rami

    Da usare quando non si è eseguito ancora un power flow
    """
    line_coord(rete)
    # Schema interattivo dei risultati, si crea un file chiamato 'temp-plot.html' che deve essere aperto con il navigatore
    open_lines = rete.line[rete.line.in_service == False].index
    closed_lines = rete.line[rete.line.in_service == True].index
    lc = pplotly.create_line_trace(rete, closed_lines, color="grey", width= 1.5,
     infofunc=pd.Series(index= rete.line.index, data= 'Ramo ' + rete.line.name + '<br>'))
    olc = pplotly.create_line_trace(rete, open_lines, color="red", 
     infofunc=pd.Series(index=rete.line.index, data= 'Ramo ' + rete.line.name + '<br>' + 'Ramo aperto'))
    
    bc = pplotly.create_bus_trace(rete, rete.bus.index, size=10, color="orange",
     infofunc=pd.Series(index= rete.bus.index, data = rete.bus.name + '<br>'+ 'Vn: '+rete.bus.vn_kv.astype(str) + 'kV'+ '<br>'))
    egc = pplotly.create_bus_trace(rete, rete.ext_grid.index, size= 15, patch_type= 'square', color= 'orange', infofunc=pd.Series(index= rete.ext_grid.index, data = rete.bus.name + '<br>'+ 'Vn: '+rete.bus.vn_kv.astype(str) + 'kV'+ '<br>') )
    pplotly.draw_traces(bc + egc + lc + olc, auto_open= False)

def scambio_rami_pp(rete_pp):
    """
    Aggiorna la matrice delle incidenze e fa un scambio di rami
    """
    # aggionamento matrice delle incidenze
    rami_chiusi = rete_pp.line[rete_pp.line.in_service == True]
    rami_aperti = rete_pp.line[rete_pp.line.in_service == False]

    rami_chiusi = rami_chiusi.reset_index(drop=True)
    rete_pp.line = pd.concat([rami_chiusi, rami_aperti], ignore_index= True)

    dim_L = len(rami_chiusi)
    mat_incid = np.zeros((dim_L,dim_L))
    cont = 0
    for i in rami_chiusi.index:   # para cada fila del df rami_chiusi
        if rami_chiusi['from_bus'][i] != 0:
            mat_incid[cont][rami_chiusi['from_bus'][i]-1] = 1
        mat_incid[cont][rami_chiusi['to_bus'][i]-1] = -1
        cont += 1

    # matrice gamma a partire della matrice delle incidenze
    mat_gamma = np.linalg.inv(mat_incid)

    # scambio di rami
    rami_aperti = rete_pp.line[rete_pp.line.in_service == False]
    aperti = rami_aperti.index
    chiude = np.random.choice(aperti)
    ramo_chiudere = rete_pp.line[rete_pp.line.index == chiude]
    rete_pp.line['in_service'][ramo_chiudere.index[0]] = True

    sub_1 = np.where(mat_gamma[ramo_chiudere.from_bus-1] == -1)[1]          
    sub_2 = np.where(mat_gamma[ramo_chiudere.to_bus-1] == -1)[1]            
    no_comuni = list(set(sub_1) ^ set(sub_2))

    abre = np.random.choice(no_comuni)
    ramo_aprire = rete_pp.line[rete_pp.line.index ==  abre]
    rete_pp.line['in_service'][ramo_aprire.index[0]] = False

    return rete_pp.line

def funz_ob_pen_pp(rete_pp):
    """
    Fa il calcolo di una funzione obiettivo a partire dei risultati di un power flow
    """
    perdite = rete_pp.res_line.pl_mw.sum()
    viol_tens_min = rete_pp.res_bus[rete_pp.res_bus.vm_pu < 0.97]
    viol_tens_max = rete_pp.res_bus[rete_pp.res_bus.vm_pu > 1.03]
    viol_liv_car = rete_pp.res_line[rete_pp.res_line.loading_percent/100 > 1]

    f_obb = perdite*(1 + 5*sum(0.97 - viol_tens_min.vm_pu) + 5*sum(viol_tens_max.vm_pu - 1.03) + 5*sum(viol_liv_car.loading_percent/100 - 1))
  
    return f_obb

def scambio_con_verifica_pp(rete_pp, cm = 1):
    """
    Fa un scambio di rami e verifica se la nuova soluzione viene accetata 

    Usa le funzioni: funz_ob_pen(), scambio_rami().
    Quando si usa la funzione scambio_rami() si aggiorna la matrice delle incidenze.

    Risolve la rete originale e calcola la funzione obiettivo penalizzata, fa un scambio, risolve la rete nuova e calcola la 
    funzione obiettivo penalizzata.
    Se miglioramento: la configurazione della rete resta quella della rete nuova.
    Se peggioramento: si fa una prova di accetazione, se rifiutata, la configurazione della rete torna a quella di partenza,
                      se accetata, la configurazione della rete resta quella della rete nuova.
    """
    pp.runpp(rete_pp, algorithm='bfsw', init='flat')
    f_obb = funz_ob_pen_pp(rete_pp)
    config_rete = rete_pp.line.copy()
    scambio_rami_pp(rete_pp)

    pp.runpp(rete_pp, algorithm='bfsw', init='flat')
    f_obb2 = funz_ob_pen_pp(rete_pp)

    if f_obb2 < f_obb:
        print('Miglioramento')
    else:
        peggioramento = f_obb2 - f_obb
        controllo = e**(-1*peggioramento/cm)
        r = np.random.random()
        if controllo <= r:
            print('Configurazione con peggioramento rifiutata')
            rete_pp.line = config_rete
            pp.runpp(rete_pp, algorithm='bfsw', init='flat')
        else:
            print('Configurazione con peggioramento accettata')
    return rete_pp.line
    
def iniz_controllo_pp (rete_pp):
    """
    Calcola il valore iniziale del parametro di controllo Cm per un possibile Simulated Annealing

    Usa le funzioni: funz_ob_pen(), scambio_rami()
    Quando si usa la funzione scambio_rami() si aggiorna la matrice delle incidenze.
    
    Calcola la funzione obiettivo penalizzata della configurazione base, fa scambi e valutazioni della funzione
    obiettivo penalizzata fino raggiungere 10 peggioramenti e tornado alla rete base, calcola il peggioramento medio, 
    calcola il parametro C0 considerando una prob iniziale di 0.5
    """
    peggioramento = 0
    cont = 0
    config_rete = rete_pp.line.copy()
    pp.runpp(rete_pp, algorithm='bfsw', init='flat')
    f_obb = funz_ob_pen_pp(rete_pp)
    while cont < 10:
        scambio_rami_pp(rete_pp)
        pp.runpp(rete_pp, algorithm='bfsw', init='flat')
        f_obb2 = funz_ob_pen_pp(rete_pp)

        if f_obb2 > f_obb:
            peggioramento = peggioramento + f_obb2 - f_obb
            cont += 1
        rete_pp.line = config_rete
    
    pegg_med = peggioramento/cont
    c0 = pegg_med/(log(1/0.5))
    return c0

#-- funzioni create per uso di profili --#

def carichi_separati_df(file_mat):
    """
    Crea un dataframe con le potenze nominali (in P e Q) di ogni tipo di carico in ogni nodo
    """
    # Si carica il file .mat
    info_profili = scipy.io.loadmat(file_mat)

    # I dati ricavati si salvano in un DataFrame di pandas
    # Colonna del nodo a cui e connesso il carico
    carichi_df = pd.DataFrame(range(1,info_profili['NumNodi'][0][0]+1), columns=['Nodo'])
    # Colonne delle P per tipo ad ogni nodo
    carichi_df['P_res'] = info_profili['P_res']
    carichi_df['P_ind'] = info_profili['P_ind']
    carichi_df['P_agr'] = info_profili['P_agr']
    # Colonne delle Q per tipo ad ogni nodo
    carichi_df['Q_res'] = info_profili['Q_res']
    carichi_df['Q_ind'] = info_profili['Q_ind']
    carichi_df['Q_agr'] = info_profili['Q_agr']
    return carichi_df

def carica_profili_df(file_mat):
    """
    Crea un dataframe dei profili di carico e generazione PV
    """
    info_profili = scipy.io.loadmat(file_mat)  # Si carica il file .mat
    profili_df = pd.DataFrame(info_profili['res'],columns=['res'])
    profili_df['ind'] = info_profili['ind']
    profili_df['agr'] = info_profili['agr']
    profili_df['PV'] = info_profili['PV']
    return profili_df

def carichi_pp(rete, carichi_df):
    """
    Ad ogni nodo crea un carico di ogni tipo nella rete di pandapower
    """
    rete.load.drop(rete.load.index, inplace=True)

    for i in range(carichi_df.shape[0]):
        pp.create_load(rete, carichi_df['Nodo'][i], carichi_df['P_res'][i]/1000, carichi_df['Q_res'][i]/1000, 
        name= 'Res '+str(carichi_df['Nodo'][i]))
    for i in range(carichi_df.shape[0]):   
        pp.create_load(rete, carichi_df['Nodo'][i], carichi_df['P_ind'][i]/1000, carichi_df['Q_ind'][i]/1000, 
        name= 'Ind '+str(carichi_df['Nodo'][i]))
    for i in range(carichi_df.shape[0]):
        pp.create_load(rete, carichi_df['Nodo'][i], carichi_df['P_agr'][i]/1000, carichi_df['Q_agr'][i]/1000, 
        name= 'Agr '+str(carichi_df['Nodo'][i]))
        
def gen_pv_df(file_mat):
    """
    Crea un dataframe con le potenze nominali (in P e Q) di generazione PV ad ogni nodo
    """
    # Si carica il file .mat
    info_profili = scipy.io.loadmat(file_mat)
    gen_pv_df = pd.DataFrame(range(1,info_profili['NumNodi'][0][0]+1), columns=['Nodo'])
    gen_pv_df['P_PV'] = info_profili['P_PV']
    gen_pv_df['Q_PV'] = info_profili['Q_PV']

    return gen_pv_df

def gen_pv_pp(rete, gen_pv_df):
    """
    Ad ogni nodo crea un generatore statico nella rete di pandapower
    """
    rete.sgen.drop(rete.sgen.index, inplace=True)
    for i in range(gen_pv_df.shape[0]): 
        pp.create_sgen(rete, gen_pv_df['Nodo'][i], -1*gen_pv_df['P_PV'][i]/1000, gen_pv_df['Q_PV'][i]/1000, 
        name= 'PV '+str(gen_pv_df['Nodo'][i]))

# -- funzioni create per allocazione perdite --#

def mat_incid(rete_pp):
    """
    Costruisce la matrice delle incidenze, considera i trasformatori come rami
    """
    rami_chiusi = rete_pp.line[rete_pp.line.in_service == True]
    rami_chiusi_trx = rete_pp.trafo[rete_pp.trafo.in_service == True]
    rami_aperti = rete_pp.line[rete_pp.line.in_service == False]

    rami_chiusi = rami_chiusi.reset_index(drop=True)
    rami_chiusi_trx = rami_chiusi_trx.reset_index(drop=True)
    rete_pp.line = pd.concat([rami_chiusi, rami_aperti], ignore_index= True)

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

def correnti_rami(line_df, line_res_df, rete_df):
    """
    Calcola la corrente longitudinale nei rami in complessa e in pu
    """
    i_long= np.zeros((line_df.shape[0], 1),dtype= complex)
    for ramo in range(line_df.shape[0]):
        # impedenza longitudinale in ohm
        z_ramo = line_df.r_ohm_per_km[ramo]*line_df.length_km[ramo] + (line_df.x_ohm_per_km[ramo]*line_df.length_km[ramo])*1j  
        # tensione nel nodo di partenza (in complessa e in pu)
        v_from = line_res_df.vm_from_pu[ramo]*cos(line_res_df.va_from_degree[ramo]*pi/180) + line_res_df.vm_from_pu[ramo]*sin(line_res_df.va_from_degree[ramo]*pi/180)*1j
        # tensione nel nodo di arrivo (in complessa e in pu)
        v_to = line_res_df.vm_to_pu[ramo]*cos(line_res_df.va_to_degree[ramo]*pi/180) + line_res_df.vm_to_pu[ramo]*sin(line_res_df.va_to_degree[ramo]*pi/180)*1j
        # corrente longitudinale (in complessa e in pu)
        i_long[ramo] = (v_from - v_to)/(z_ramo/rete_df['Zbase_ohm'][0])
    return i_long

def resistenza_long(line_df,rete_df):
    """
    Ottiene un vettore con le resistenze longitudinali dei rami (in pu)
    """
    imp_real = np.zeros((line_df.shape[0],1))
    for i in range(line_df.shape[0]):
        imp_real[i]= line_df.r_ohm_per_km[i]*line_df.length_km[i]/rete_df['Zbase_ohm'][0]
    return imp_real

def correnti_nodi(mat_gamma,i_long):
    """
    Calcola le correnti shunt dei nodi (in complessa e in pu)
    """
    g_tr = np.transpose(mat_gamma)
    aux = np.linalg.inv(g_tr)
    i_nodi = np.dot(aux,i_long)
    return i_nodi