# Modeling the Evaporation of a Pure Water Film Using Finite Difference Methods
# Comparison with Experimental Results
# RH (March 29, 2025)

import xlrd
from scipy.special import erf
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)
from tqdm import tqdm
from iapws._iapws import _Liquid

# Data Importation
book = xlrd.open_workbook("Data_WaterFilm_1_RH0.xls")

print ("The number of worksheets is", book.nsheets)
print ("Worksheet name(s):", book.sheet_names())
sh = book.sheet_by_index(0)
print (sh.name, sh.nrows, sh.ncols)

NLines = sh.nrows
Temps = np.zeros([NLines])   #Time in s
Pesee = np.zeros([NLines])   #Mass in g

# Data
for i in range(0, NLines):
    Temps[i]=sh.cell_value(rowx=i, colx=0)
    Pesee[i]=sh.cell_value(rowx=i, colx=1)

TempH = Temps/60
Mpap = 0.0089        #paper mass in g

DM = (abs(Pesee-Mpap)/Mpap)*100
me = Pesee-Mpap      #Water Mass at t in g

def YAf(Z, xi):# Molar Fraction, L->infinite
    return (1+(yAinf/yA0)*erf(xi) - (1-yAinf/yA0)*erf(Z-xi))/(1+erf(xi))

def YAs(z, yA0, yAinf, s0):#Molar Fraction molaire, L->finite
    res = 1 - (1-yA0)*((1-yAinf)/(1-yA0))**((z-s0)/(L-s0))
    return res/yA0

def f(x, yA0, yAinf):
    return  (yA0/(1-yA0))*np.exp(-x**2)*(1-yAinf/yA0)/(1+erf(x))/np.pi**0.5-x

def Newton(f, x0, yA0, yAinf, eps = 1e-6, dx = 1e-5, N = 10000):
      for i in range(0, N):
          df = (f(x0+dx, yA0, yAinf)-f(x0, yA0, yAinf))/dx
          xi = x0 - f(x0, yA0, yAinf)/df
          if (abs(f(xi, yA0, yAinf)-f(x0, yA0, yAinf)) < eps):
              return (xi, i)
              break
          x0 = xi
      return (x0, i)

# System properties
D0, p_sat, p0 = 2.503e-5, 3141.8, 101325
a0, aL = 1.0, 0.05
yA0, yAinf = a0*p_sat/p0, aL*p_sat/p0
xi = Newton(f, 0, yA0, yAinf)[0]
T, p, R, Mw, Mas = 298.15, 101325, 8.31446, 18e-3, 29e-3
rho_w_L = _Liquid(T = T, P=0.101325)['rho']
beta_w_L = rho_w_L*R*T/(Mw*p)
Sigma = 8e-3*7e-3      # Paper area
print ("Paper area : Sigma = {0:1.3e} m^2".format(Sigma))
def YAInf(z, s0, t):
    #return (1-erf(z/(4*D0*t)**0.5-xi))/(1+erf(xi))
    return YAf((z-s0)/(4*D0*t)**0.5, xi)
L, s0 = 3.45e-3, (539.29e-6)/1
tau = L**2/D0
tau_0 = s0**2/D0
print("Characteristic time tau_0 = {0:1.3e} s".format(tau_0))
Fac_t = 100
tf = 60*60 # final time
trf = tf*tau # reduced final time
Nx, Nt = 20, int(200*Fac_t)
dx, dt = 1.0/Nx, tf/Nt
alpha = dt/dx**2.0
print("alpha = {0:1.3e}".format(alpha))
# Solutions vector
v = np.zeros([Nx+1, Nt+1]) # Vec. colonne solution chgt de variable
x = np.linspace(0, 1, Nx +1)
t = np.linspace(0, tf, Nt +1)
yA = np.zeros([Nx+1, Nt+1]) # Vec. colonne solution reelle (fraction molaire)
yAQS = np.zeros([Nx+1, Nt+1]) # Vec. colonne solution quasi-stationnaire (fraction molaire)
s = np.zeros([Nt+1]) # Front d'evaporation
s_rec = np.zeros([Nt+1]) # Front d'evaporation reconstruit
ds_dt = np.zeros([Nt+1]) # ds/dt
v_mol = np.zeros_like(s) # vitesse molaire moyenne dans la phase gazeuse
# drying rate pe
pe = -(p0*Mw*D0)*np.log((1-yAinf)/(1-yA0))/(R*T*rho_w_L*(L-s0))
# SA
s_AS = s0 + pe*t
s_AS[s_AS < 0] = 0
# QSA
Lc = L-s0
s_AQS = s0 - Lc*((1-2*pe*t/Lc)**0.5 - 1.0)
s_AQS[s_AQS < 0] = 0
# CLs and CIs
v[1:,0] = yAinf; v[0,:] = yA0 ; v[-1,:] = yAinf; s[0] = s0; s_rec[0] = s0
CL = np.zeros([Nx-1])

# some functions
def sn(v0n, v1n, snm1):# s[n]
    y1 = D0*dt*(v0n-v1n)/((beta_w_L*dx)*(1-yA0))
    y2 = (L+snm1)**2 + 4*(y1-L*snm1)
    return 0.5*(L+snm1-y2**0.5)

# Main loop
for n in tqdm(range(1, Nt+1), desc = 'Main Loop'):
    s[n] = sn(v[0,n-1], v[1, n-1], s[n-1])
    if s[n] < 0: # film is dried
        s[n] = 0
    dsSdt = (s[n] - s[n-1])/dt
    ds_dt[n] = dsSdt
    s_rec[n] = s_rec[n-1] + dsSdt*dt
    v_mol[n] = -beta_w_L*dsSdt
    # ABC matrix
    ABC = np.zeros([Nx-1,Nx-1])
    Bin = 2*alpha*D0+(L-s[n])**2
    xi_1_n = (1+beta_w_L-dx)*(L-s[n])*dsSdt
    A1n = -alpha*D0*(1-0.5*xi_1_n*dx/D0)
    CL[0] = A1n*yA0
    for i in range(0, Nx-1):
        xi_i_n = (1+beta_w_L-(i+1)*dx)*(L-s[n])*dsSdt
        Ain = -alpha*D0*(1-0.5*xi_i_n*dx/D0)
        Cin = -alpha*D0*(1+0.5*xi_i_n*dx/D0)
        if i == 0:#
            ABC[0,0] = Bin ; ABC[0, 1] = Cin
        elif i == Nx-2:#
            ABC[Nx-2, Nx-3] = Ain; ABC[Nx-2, Nx-2] = Bin
        else:#
            ABC[i, i-1] = Ain; ABC[i, i] = Bin; ABC[i, i+1] = Cin
    xi_Nxm1_n = (1+beta_w_L-(Nx-1)*dx)*(L-s[n])*dsSdt
    CNxm1n = -alpha*D0*(1+0.5*xi_Nxm1_n*dx/D0)
    CL[Nx-2] = CNxm1n*yAinf
    ABC_inv = LA.inv(ABC)
    v[1:-1, n] = ABC_inv@((L-s[n])**2*v[1:-1, n-1] - CL)
    zAs = np.linspace(s_AQS[n], L, Nx+1)
    yAQS[:,n] = yA0*YAs(zAs, yA0, yAinf, s_AQS[n])

# Graphical Figures
fig = plt.figure(1, figsize=(12,9))
n_fig = 3
plt.subplot(n_fig,1,1)
plt.plot(x, v[:, Nt], label = s)# v(x,t) Profile at t
plt.xlabel('$x$'); plt.ylabel('$v(x,t_f)$'); plt.grid();
# Real profile of yA(z,t) at t1
plt.subplot(n_fig,1,2)
nt1 = Nt
zt = x*(L-s[nt1]) + s[nt1]
z = x*L
plt.plot(zt*1e3, v[:, Nt], label = r'exact profile $y_w(z,t_f)$')# v(x,t) Profile at t
yw_lin = (yAinf - yA0)*(zt)/(L - s[nt1]) + yA0  # Approximate linear Profile at t
plt.plot(zt*1e3, yw_lin, 'k--', label = r'approximate linear profile $y^L_w(z,t)$')
zAs = np.linspace(s0, L, Nx+1)
plt.plot(zAs*1e3, yA0*YAInf(zAs, s0, tf), 'k*', label = r'$L_c\rightarrow\infty$  and $s(t)=s_0$')
# SA
plt.plot(zAs*1e3, yA0*YAs(zAs, yA0, yAinf, s0), 'k.', label = r'(SA), $L_c<\infty$ and $s(t)=s_0$')
# QSA
zAQs = np.linspace(s_AQS[nt1], L, Nx+1)
plt.plot(zAQs*1e3, yA0*YAs(zAQs, yA0, yAinf, s_AQS[nt1]), 'b--', label = r'(QSA), $L_c<\infty$ and $s(t)\neq\text{cste}$')
plt.xlabel('$z$ (in mm)'); plt.ylabel('$y_w(z,t_f)$'); plt.grid();
plt.legend()
# Film width s(t)
plt.subplot(n_fig,1,3)
mn, sec = 60, 1
ech_t = mn; str_ech_t = 'min'
plt.plot(t/ech_t, s/s0, label = 'exact resolution')#
plt.plot(t/ech_t, s_AS/s0, 'b-', label = 'SA')#
plt.plot(t/ech_t, s_AQS/s0, 'r--', label = 'QSA')#
nF, Dn = len(np.where(Temps<=tf)[0]), 2
# Exp. data
plt.plot((Temps[:nF:Dn]-0.0)/ech_t , me[:nF:Dn]/me[0],'k+',markerfacecolor='none',markeredgewidth=1.5,markersize=8,lw=1.25, label = 'Exp. results')
str_ax = '$t$ (in '+str_ech_t+')'; plt.xlabel(str_ax)
plt.ylabel('$s(t)/s_0$'); plt.grid();
plt.legend(); plt.tight_layout()
plt.show()
