import matplotlib.pyplot as plt 
import matplotlib.cbook as cbook
import matplotlib
import numpy as np
# matplotlib.rcParams['text.usetex'] = True
# matplotlib.use('PS')
#2 matplotlib.rcParams['text.latex.unicode']=True
# from matplotlib.axes.Axes import inset_axes, indicate_inset_zoom, imshow
# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
# from mpl_toolkits.axes_grid1.inset_locator import mark_inset
plt.style.use('ggplot')
data_x=[]
data_y=[]
# fname2 = cbook.get_sample_data('data_x_x2_x3.csv', asfileobj=False)
# with cbook.get_sample_data('test_bisection.txt') as file:
array = np.loadtxt('test_bisection.txt')

# for i in range(100):
#     data_x.append(i/100)
#     data_y.append((i/100)**2)

data_x = array[:,0]
data_y = array[:,3]
fig1,ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(data_x,data_y)
ax1.set_xlim(0,0.003)
ax1.set_xlabel(r"$r/a_0$")
ax1.set_ylabel(r"$P(r)$")
# ax1.set_title("Wavefunction for state 3D of "+r"$^{16}_{8}O$"+" under potential "+r"$V_0=-20 MeV$")
ax1.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
 
axins = ax1.inset_axes([0.5,0.5,0.47,0.47])
# axins.plot(5.291772*10**(-11)*10**(15)*data_x,data_y)
axins.plot(data_x,data_y)
axins.vlines(x=5.7139*10**(-5),ymin=0.0,ymax=8.0)
axins.ticklabel_format(style='sci',axis='x',scilimits=(0,0))

x1,x2,y1,y2 = 0,0.00012,0,8
axins.set_xlim(x1,x2)
axins.set_ylim(y1,y2)
axins.xaxis.set_visible('False')
axins.yaxis.set_visible('False')
axins.annotate("Nuclear Radius\n      3.02fm ",(0.6*10**(-4),1))
ax1.indicate_inset_zoom(axins)
# ax1.indicate_inset(axins,1,1)
# print("hi")
# plt.show()
fig1.savefig('out.png',bbox_inches='tight')
