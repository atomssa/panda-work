#!/usr/bin/python

from math import *
from matplotlib import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.integrate import quad

font = {'family' : 'serif', 'color'  : 'darkred', 'weight' : 'normal', 'size'   : 16}

#define some constants
#p_beam = 5.513
#e_beam = hypot(Mp, p_beam)
#s = 2*Mp*e_beam + 2*Mp*Mp
#sqrt_s = sqrt(s)

Mp = 0.938
Mp_sq = Mp*Mp
Mpi0 = 0.135
Mpi0_sq = Mpi0*Mpi0
Mjpsi = 3.097
Mjpsi2 = Mjpsi*Mjpsi
fpi = 0.093 # GeV
fpsi = 0.413 #GeV
fN = 0.005 #GeV
gpiNN = 13

alpha_s = 0.25
Mbar = 3.0 #GeV
C = pow(4.0*pi*alpha_s,3)*fN*fN*fpsi*10.0/fpi/81.0

M0 = 8776.88 #//GeV
#jpsi_pbarp_width = 1.2e-6;
#M0 = sqrt(jpsi_pbarp_width*243.0*pi*pow(Mbar,5)/pow(pi*alpha_s,6)/1280.0/pow(fpsi,2)/pow(fN,4))
#print "Mp = %4.3f" % Mp
#print "Mpi0 = %4.3f" % Mpi0
#print "Mjpsi = %4.3f" % Mjpsi
print "fpi = %4.3f" % fpi
print "fpsi = %4.3f" % fpsi
print "fN = %4.3f" % fN
print "gpiNN = %4.3f" % gpiNN
print "alpha_s = %4.3f" % alpha_s
print "Mbar = %4.3f" % Mbar
print "M0 = %4.3f" % M0
print "C = %4.3e" % C

# funciton definitions
def _lambda_sq(x, y, z):
    return (x*x)+(y*y)+(z*z)-(2*x*y)-(2*x*z)-(2*y*z)

def _xi_a(s):
    return pow(Mbar,2)/(2.0*s - pow(Mbar,2))

def _xi(Del_sq,s):
    return _xi_a(s);
    #return (pow(Mjpsi,2) - Del_sq - Mp_sq)/(2.0*s - pow(Mjpsi,2) + Del_sq - 3.0*Mp_sq)

def _DelT_sq(Del_sq,s):
    xi = _xi(Del_sq,s)
    return ((1-xi)/(1+xi)) * (Del_sq - 2.0*xi*( (Mp_sq/(1+xi)) - (Mpi0_sq/(1-xi)) ))

def _Del_sq_max(s):
    xi = _xi_a(s)
    return 2.0*xi*(Mp_sq*(xi-1) + Mpi0_sq*(xi+1))/(pow(xi,2)-1)

def _Del_sq_DelT_sq(DelT_sq,s):
    xi = _xi_a(s)
    return (DelT_sq*(1+xi)/(1-xi)) + 2*xi*((Mp_sq/(1+xi))-(Mpi0_sq/(1-xi)))

def _dsig_dDel_sq(Del_sq,s):
    pbarn = 1e12 #4.348e6 # 5.17377e6 #0.38941e9 # 4e6
    xi = _xi(Del_sq,s)
    DelT_sq = _DelT_sq(Del_sq,s)
    F1 = 1/16.0/pi/_lambda_sq(s,Mp_sq,Mp_sq)
    F2 = pow(C,2)*2.0*(1+xi)/4.0/xi/pow(Mbar, 6)
    I = fpi*gpiNN*Mp*(1-xi)*M0/(Del_sq - Mp_sq)/(1+xi)
    Iprim = fpi*gpiNN*Mp*M0/(Del_sq - Mp_sq)
    return pbarn * F1 * F2 * ( pow(I,2) - (DelT_sq*pow(Iprim,2)/Mp_sq) )

svals =    [12.2649, 15,   16.8745,  20,    24.346]
min_vals = [-0.45,   -0.5, -2.76,    -0.5,  -6.5]
max_vals = [0.61,    0.0,  0.46,     0.0,   0.32]
cols =     ['r',     'g',  'b',      'c',   'm']

#integration
def integrate(sel,min_val_idx):
    s = []
    xsect_tot = []
    xsect_tot_er = []
    for idx, _val in enumerate(svals):
        if sel==1 and (idx==1 or idx==3): continue
        if sel==2 and (idx==0 or idx==2 or idx==4): continue

        if (min_val_idx==1):
            min_vals[idx] = _Del_sq_DelT_sq(-0.5,_val)
            max_vals[idx] = _Del_sq_DelT_sq(0,_val)

        result = quad( _dsig_dDel_sq, min_vals[idx], max_vals[idx], args=( _val ) );
        s.append(_val)
        xsect_tot.append(result[0])
        xsect_tot_er.append(result[1])
        print "xsect(s= %4.2f, %4.2f < Delta^2 < %4.2f)= %4.2f" % (_val, min_vals[idx], max_vals[idx], result[0])
    fig, ax = plt.subplots()
    ax.plot(s, xsect_tot, 'b', label=r'$\sigma_{tot}$')
    plt.xlabel(r"$s[GeV^{2}]$", fontdict=font);
    plt.ylabel(r"$d\sigma_{tot}[pb]$",fontdict=font);

#plotting
def plot_del(ax, sel):
    #ax.set_yscale("log")
    npt = 100
    #svals = [12.2649, 15,  16.8745, 20,  24.346,    40,   50, 60,  70,  100, 150, 300, 600, 1000, 2000, 2500, 5000, 10000 ]
    #cols =  ['r',     'g', 'b',     'c', 'm', 'r',  'g', 'b', 'c', 'm', 'r', 'g', 'b', 'c', 'm',  'r',  'g',  'b',  'c',  ]
    for idx, _val in enumerate(svals):
        if sel==1 and (idx==1 or idx==3): continue
        if sel==2 and (idx==0 or idx==2 or idx==4): continue
        #if idx == 1 or idx == 3:
        #    Del_sq_min = _Del_sq_DelT_sq(-0.5,_val)
        #else:
        #    Del_sq_min = min_vals[idx]
        Del_sq_min = _Del_sq_DelT_sq(-0.5,_val)
        Del_sq_max = _Del_sq_max(_val)
        Del_sq=np.linspace(Del_sq_min,Del_sq_max,npt)
        #Del_sq=np.linspace(-0.6,0.5,npt)
        s = np.full(npt, _val)
        DelT_sq=map(_DelT_sq,Del_sq,s)
        xsect=map(_dsig_dDel_sq,Del_sq,s)
        ax.plot(DelT_sq, xsect, cols[idx], label=r'$s = %4.2f GeV^{2}$' % _val)
        #ax.plot(Del_sq, xsect, cols[idx], label=r'$s = %4.2f GeV^{2}$' % _val)
    plt.xlabel(r"$\Delta^{2}_{T}[GeV^{2}]$", fontdict=font);
    plt.ylabel(r"$d\sigma/d\Delta^{2}[pb/GeV^2]$",fontdict=font);
    #plt.axis([-0.51,0.01,50,250]);
    #plt.yticks(np.arange(min(xsect), max(xsect)+1, 10.0))
    legend = ax.legend(loc='best', shadow=True);
    frame = legend.get_frame();
    frame.set_facecolor('0.90')

def plot_s(ax):
    npt = 100
    a_s = np.linspace(10, 30, npt)
    Del_sq_max = map(_Del_sq_max,a_s)
    xsect = map(_dsig_dDel_sq, Del_sq_max, a_s)
    ax.plot(a_s, xsect, 'b', label=r'$\Delta^{2}_{T} = 0$')
    plt.xlabel(r"$s[GeV^{2}]$", fontdict=font);
    plt.ylabel(r"$d\sigma/d\Delta^{2}[pb/GeV^2]$",fontdict=font);
    legend = ax.legend(loc='best', shadow=True);
    frame = legend.get_frame();
    frame.set_facecolor('0.90')

def plot(idx=0, sel=0):
    #fig, (ax1, ax2) = plt.subplots(1,2)
    fig, ax = plt.subplots()
    if idx==0:
        plot_del(ax,sel)
    else:
        plot_s(ax)
    plt.show()

if __name__ == "__main__":
    plot()
