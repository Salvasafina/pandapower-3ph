import pandapower.networks as nw
import matplotlib.pyplot as plt
import pandapower as pp
import pandapower.plotting as pplot
from pandapower.pf.runpp_3ph import runpp_3ph

# load IEEE European LV network from Pandapower
net = nw.ieee_european_lv_asymmetric()

# fictional loads added just to be plotted
pp.create_loads(net= net, buses = net.asymmetric_load.bus, p_mw = 0, q_mvar = 0)

# number of elements by type present on the network
print(net)

# real asymmetric loads
print(net.asymmetric_load)

# base scenario power flow
runpp_3ph(net)

# power flow results 
print(net.res_bus_3ph)

# network diagram including loads
pplot.simple_plot(net, show_plot=False, bus_size= .4, plot_loads= True, load_size= .8)
print(net.res_bus_3ph)
plt.show()
