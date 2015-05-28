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

Mp = 0.94
Mp_sq = Mp*Mp
Mpi0 = 0.0 #0.135
#Mpi0 = 0.135
Mpi0_sq = Mpi0*Mpi0
Mjpsi = 3.097
Mjpsi2 = Mjpsi*Mjpsi
fpi = 0.093 # GeV
fpsi = 0.413 #GeV
#fpsi = 0.383 #GeV
fN = 0.005 #GeV^2
gpiNN = 13.21

alpha_s = 0.2520617559544493
Mbar = 3 #GeV
C = pow(4.0*pi*alpha_s,3)*fN*fN*fpsi*10.0/fpi/81.0
#C = 0.000435585

#(4 \[Pi]*\[Alpha]s)^3 (fpsi (fN)^2)/f\[Pi] 10/81;

M0 = 8776.88 #//GeV
#jpsi_pbarp_width = 0.42e-6;
#M0 = sqrt(jpsi_pbarp_width*243.0*pi*pow(Mbar,5)/pow(pi*alpha_s,6)/1280.0/pow(fpsi,2)/pow(fN,4))
print "MN = %4.3f" % Mp
print "Mpi0 = %4.3f" % Mpi0
#print "Mjpsi = %4.3f" % Mjpsi
print "fpi = %4.3f" % fpi
print "fpsi = %4.3f" % fpsi
print "fN = %4.3f" % fN
print "gpiNN = %4.3f" % gpiNN
print "alpha_s = %4.3f" % alpha_s
print "Mbar = %4.3f" % Mbar
print "M0 = %4.3f" % M0
print "C = %4.3e" % C
#print "CC = %4.3e" % CC

# funciton definitions
def _lambda_sq(x, y, z):
    return (x*x)+(y*y)+(z*z)-(2*x*y)-(2*x*z)-(2*y*z)

def _xi_a(s):
    return pow(Mbar,2)/(2.0*s - pow(Mbar,2))

def _xi(Del_sq,s):
    return (pow(Mjpsi,2) - Del_sq - Mp_sq)/(2.0*s - pow(Mjpsi,2) + Del_sq - 3.0*Mp_sq)

def _DelT_sq(Del_sq,s):
    xi = _xi_a(s)
    return ((1-xi)/(1+xi)) * (Del_sq - 2.0*xi*( (Mp_sq/(1+xi)) - (Mpi0_sq/(1-xi)) ))

def _Del_sq_max(s):
    xi = _xi_a(s)
    return 2.0*xi*(Mp_sq*(xi-1) + Mpi0_sq*(xi+1))/(pow(xi,2)-1)

def _Del_sq_DelT_sq(DelT_sq,s):
    xi = _xi_a(s)
    return (DelT_sq*(1+xi)/(1-xi)) + 2*xi*((Mp_sq/(1+xi))-(Mpi0_sq/(1-xi)))

def _dsig_dDel_sq(Del_sq,s):
    pbarn = 0.39e9 # 1gev^-2=3.9e9pbarn
    xi = _xi_a(s)
    DelT_sq = _DelT_sq(Del_sq,s)
    F1 = 1/16.0/pi/_lambda_sq(s,Mp_sq,Mp_sq)
    F2 = pow(C,2)*2.0*(1+xi)/4.0/xi/pow(Mbar, 8)
    I = fpi*gpiNN*Mp*(1-xi)*M0/(Del_sq - Mp_sq)/(1+xi)
    Iprim = fpi*gpiNN*Mp*M0/(Del_sq - Mp_sq)
    return pbarn * F1 * F2 * ( pow(I,2) - (DelT_sq*pow(Iprim,2)/Mp_sq) )

svals =    [12.2649, 15,   16.8745,  20,    24.346]
#min_vals = [-0.45,   -0.5, -2.76,    -0.5,  -6.5]
#min_vals = [-0.5,   -0.5,  -0.5,     -0.5,  -0.5]
min_vals = [-0.3,   -0.45,  -0.45,     -0.45,  -0.45]
max_vals = [0.6,    0.0,  0.46,     0.0,   0.32]
#max_vals = [0.0,    0.0,  0.0,     0.0,   0.0]
cols =     ['r',     'g',  'b',      'c',   'm']

def numevt():
    lumi = 2e3 # b-1
    br = 5.94e-2 # branching ratio
    f_t = [] # fraction within normalization window of |t| approx
    sig_normr = []
    sig_fullr = []
    nevt_normr = []
    nevt_fullr = []
    s_ = [12.2649, 16.8745, 24.346]
    # these are numbers derived from analytical formula for
    # values of t at cos_theta_cm = {-1, 0, 1}
    #{-1.50406, -0.443789, 0.616486}
    #{-5.96462, -2.75368, 0.457248}
    #{-13.2926, -6.48861, 0.31538}
    #normr = [ (-0.443789, 0.616486), (-1.0, 0.457248), (-1.0, 0.31538) ] # normalization winodw in Delta^2
    normr = [ (-0.443789, 0.616486), (-0.5, 0.457248), (-0.5, 0.31538) ] # normalization winodw in Delta^2
    fullr = [ (-0.443789, 0.616486), (-2.75368, 0.457248), (-6.48861, 0.31538)] # full range from 0 < cost < 1
    print s_
    for idx, s in enumerate(s_):
        print idx,s
        sig_normr.append(quad(_dsig_dDel_sq, normr[idx][0], normr[idx][1] , args=(s))[0])
        sig_fullr.append(quad(_dsig_dDel_sq, fullr[idx][0], fullr[idx][1] , args=(s))[0])
        f_t.append(sig_normr[idx]/sig_fullr[idx])
        nevt_normr.append(sig_normr[idx]*lumi*br*2)
        nevt_fullr.append(sig_fullr[idx]*lumi*br*2)
        print("s=%4.2f  xs(nr)=%4.2f nevt(nr)=%4.2f xs(nr)=%4.2f nevt(fr)=%4.2f"% (s, sig_normr[idx], nevt_normr[idx], sig_fullr[idx], nevt_fullr[idx]))

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
    plt.xlabel(r"$s[GeV^{2}]$", fontdict=font)
    plt.ylabel(r"$d\sigma_{tot}[pb]$",fontdict=font)

# test if pi0 mass hypothesis distorts Del^2 distributions too much
def plot_mpi0_test(iplab):
    global Mpi0
    global Mpi0_sq
    npt = 100
    _val = svals[iplab]
    s = np.full(npt, _val)

    fig, ax = plt.subplots()

    Del_sq_min = _Del_sq_DelT_sq(-.5,_val)
    Del_sq_max = _Del_sq_max(_val)
    Del_sq=np.linspace(Del_sq_min,Del_sq_max,npt)
    DelT_sq=map(_DelT_sq,Del_sq,s)
    xsect=map(_dsig_dDel_sq,Del_sq,s)
    ax.plot(DelT_sq, xsect, 'b', label=r'$s = %4.2f GeV^{2} M_{\pi^0}=%5.3f GeV/c^{2}$' %( _val, Mpi0))
    #ax.plot(Del_sq, xsect, cols[idx], label=r'$s = %4.2f GeV^{2}$' % _val)

    Mpi0 = 0.135
    Mpi0_sq = Mpi0*Mpi0
    Del_sq_min_135 = _Del_sq_DelT_sq(-.5,_val)
    Del_sq_max_135 = _Del_sq_max(_val)
    Del_sq_135 =np.linspace(Del_sq_min_135,Del_sq_max_135,npt)
    DelT_sq_135 =map(_DelT_sq,Del_sq_135,s)
    xsect_135 =map(_dsig_dDel_sq,Del_sq_135,s)

    ax.plot(DelT_sq_135, xsect_135, 'r', label=r'$s = %4.2f GeV^{2} M_{\pi^0}=%5.3f GeV/c^{2}$, ' % (_val, Mpi0))
    #ax.plot(Del_sq_135, xsect_135, cols[idx], label=r'$s = %4.2f GeV^{2}$' % _val)

    plt.xlabel(r"$\Delta^{2}_{T}[GeV^{2}]$", fontdict=font);
    plt.ylabel(r"$d\sigma/d\Delta^{2}[pb/GeV^2]$",fontdict=font);

    #plt.axis([-0.51,0.01,50,250]);
    #plt.yticks(np.arange(min(xsect), max(xsect)+1, 10.0))

    legend = ax.legend(loc='best', shadow=True);
    frame = legend.get_frame();
    frame.set_facecolor('0.90')

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
        Del_sq_min = _Del_sq_DelT_sq(-2.0,_val)
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

def plot_xsect_int():
    svals2 =    [12.2649, 16.8745,  24.346]
    nbins = 10;
    del2_bounds = [ [-0.3, 0.6], [-1.0, 0.46], [-1.0, 0.32] ]
    del2_ranges = []
    for b in del2_bounds:
        r = []
        for i in range(nbins):
            r.append([b[0]+i*(b[1]-b[0])/nbins, b[0]+((i+1)*(b[1]-b[0]))/nbins])
        del2_ranges.append(r)
    print del2_ranges;
    del2_mean = []
    xsect_int = []
    for iplab, blist in enumerate(del2_ranges):
        xs = []
        d2m = []
        for irange, bound in enumerate(blist):
            result = quad( _dsig_dDel_sq, bound[0], bound[1], args=(svals2[iplab]));
            print("plab= %4.2f, sig_tot(%4.2f<Delta^2<%4.2f) = %4.2f pm %4.2f" %(svals2[iplab], bound[0], bound[1], result[0], result[1]) )
            d2m.append((bound[0]+bound[1])/2.)
            xs.append(result[0])
        del2_mean.append(d2m)
        xsect_int.append(xs)

    fig, axes = plt.subplots(1,3,figsize=(14,6))
    for idx, ax in enumerate(axes):
        ax.plot(del2_mean[idx], xsect_int[idx], 'bo' if idx is 0 else 'rs' if idx is 1 else 'm^')

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
