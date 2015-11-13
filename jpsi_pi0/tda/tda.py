#!/usr/bin/python

from math import *
from matplotlib import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.integrate import quad

font = {'family' : 'serif', 'color'  : 'darkred', 'weight' : 'normal', 'size'   : 26}

#define some constants
#p_beam = 5.513
#e_beam = hypot(Mp, p_beam)
#s = 2*Mp*e_beam + 2*Mp*Mp
#sqrt_s = sqrt(s)

# used for mass simulations
# constants for prob. calculation
#_Mp = 0.938;
#_Mpi0 = 0.0; // This is required to reproduce x-sect dists. on PLB paper
#_Mjpsi2 = 9.59;
#pi = 3.14159;
#fpi = 0.093;// GeV
#gpiNN = 13;
#M0 = 8776.88; //GeV
#fphi = 0.413; //GeV
#fN = 0.005;//GeV
#alpha_s = 0.252062;
#M = 3.;//GeV
#C = pow(4.*pi*alpha_s,3)*(fN*fN*fphi/fpi)*(10./81.);

Mp = 0.938
Mp_sq = Mp*Mp
Mpi0 = 0.0
#Mpi0 = 0.135
Mpi0_sq = Mpi0*Mpi0
Mjpsi = 3.097
Mjpsi2 = Mjpsi*Mjpsi
fpi = 0.093 # GeV
fpsi = 0.413 #GeV
#fpsi = 0.383 #GeV
fN = 0.005 #GeV^2
gpiNN = 13.0 #.21
#alpha_s = 0.2520617559544493
alpha_s = 0.252062

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

def _lambda(x, y, z):
    return sqrt(_lambda_sq(x, y, z))

def _xi_a(s):
    return pow(Mbar,2)/(2.0*s - pow(Mbar,2))

def _xi_a_qdep(s, q):
    return pow(q,2)/(2.0*s - pow(q,2))

def _xi(Del_sq,s):
    return (pow(Mjpsi,2) - Del_sq - Mp_sq)/(2.0*s - pow(Mjpsi,2) + Del_sq - 3.0*Mp_sq)

def _xi_q(Del_sq,qsq,s):
    return (qsq - Del_sq - Mp_sq)/(2.0*s - qsq + Del_sq - 3.0*Mp_sq)

def _DelT_sq(Del_sq,s):
    xi = _xi_a(s)
    return ((1-xi)/(1+xi)) * (Del_sq - 2.0*xi*( (Mp_sq/(1+xi)) - (Mpi0_sq/(1-xi)) ))

def _Del_sq_max(s):
    xi = _xi_a(s)
    return 2.0*xi*(Mp_sq*(xi-1) + Mpi0_sq*(xi+1))/(pow(xi,2)-1)

def _Del_sq_DelT_sq(DelT_sq,s):
    xi = _xi_a(s)
    return (DelT_sq*(1+xi)/(1-xi)) + 2*xi*((Mp_sq/(1+xi))-(Mpi0_sq/(1-xi)))

def _t(qsq, costh, s):
    #t2=costh * sqrt(1-(4*Mp_sq/s)) * _lambda(s, qsq, Mpi0_sq)
    #print "%f" % (4*Mp_sq/s)
    f2=sqrt(1-(4*Mp_sq/s))
    f3=_lambda(s, qsq, Mpi0_sq)
    t2=costh * f2 * f3;
    return 0.5*(Mpi0_sq + t2 + 2*Mp_sq + qsq - s )

#def _t_max(qsq, s):
#    return _t(qsq, 1.0, s)
#
#def _t_min(qsq, s):
#    return _t(qsq, 0.5, s);

def _integrand(qsq, s, costh_min, costh_max):
    _tmin=_t(qsq,costh_min,s)
    _tmax=_t(qsq,costh_max,s)
    _ret = (_tmax-_tmin)/pow(qsq,5)
    #print "qsq=%f, qsqp5=%f, 1/qsqp5=%f, tmin=%f, tmax=%f, retval=%f"%(qsq,pow(qsq,5),1./pow(qsq,5),_tmin,_tmax,_ret)
    return _ret

def _cost(qsq, s, t):
    d1=sqrt(0.25*s-Mp_sq)
    d2=_lambda(qsq,s,Mpi0_sq)/sqrt(s)
    return -0.5*(2*Mp_sq+Mpi0_sq+qsq-s-2*t)/d1/d2

def sig_tot_epem(qsq_min, qsq_max, s):
    a_em = pow(4*pi/137,2)
    a_str = pow(4*pi*0.3,3)
    f_n_gstar = 5.3e-3
    den_gstar = pow(54,2)*pow(2*pi,3)*pow(fpi,2)
    const_K = a_em*a_str*pow(f_n_gstar,4)/den_gstar
    Isq=1.69e9
    K=const_K*Isq

    fbarn = 0.39e12 # 1gev^-2=3.9e9pbarn

    integral_d1=quad(_integrand, 3.0, 4.3, args=(5,0.5,1.0))[0]
    K = 1675*3*(5-4*Mp_sq)/8./integral_d1/fbarn

    re_integral_d1=quad(_integrand, 3.0, 4.3, args=(5,0.5,1.0))[0]
    xsect_d1 = fbarn*8*K*re_integral_d1/(5 - 4*Mp_sq)/3
    print "xsect_d1= %f"%xsect_d1

    integral_d2=quad(_integrand, 5.0, 9.0, args=(10,0.5,1.0))[0]
    xsect_d2 = fbarn*8*K*integral_d2/(10 - 4*Mp_sq)/3
    print "xsect_d2= %f"%xsect_d2

    normr = [ (-0.092, 0.59), (-1.3, 0.43), (-2.85, 0.3) ]
    s_ = [12.2649, 16.8745, 24.346]
    Br=0.0594
    _xs = [ 206800*Br, 309600*Br, 316300*Br]
    for idx, s in enumerate(s_):
        ct_min=_cost(Mjpsi2, s, normr[idx][0])
        ct_max=_cost(Mjpsi2, s, normr[idx][1])
        #integral=quad(_integrand, 2.96, 3.22, args=(s,ct_min,ct_max))[0]
        #integral=quad(_integrand, 2.94, 3.24, args=(s,ct_min,ct_max))[0]
        integral=quad(_integrand, pow(2.9,2), pow(3.3,2), args=(s,ct_min,ct_max))[0]
        xsect = fbarn*8*K*integral/(s - 4*Mp_sq)/3
        stob=_xs[idx]/xsect
        print "s=%5.2f dt=%5.3f, ct_min=%4.2f, ct_max=%4.2f xsect=%6.1f stob=%4.1f" % (s,normr[idx][1]-normr[idx][0],ct_min,ct_max,xsect,stob)

    integral=quad(_integrand, qsq_min, qsq_max, args=(s,0.5,1.0))[0]
    xsect = fbarn*8*K*integral/(s - 4*Mp_sq)/3
    return xsect

#def _dsig_dq_gstar(Del_sq, s, qsq):
#    pbarn = 0.39e9 # 1gev^-2=3.9e9pbarn
#    #xi = _xi_a_qdep(s,qsq)
#    #DelT_sq = _DelT_sq(Del_sq,s)
#    #I = fpi*gpiNN*Mp*(1-xi)*M0/(Del_sq - Mp_sq)/(1+xi)
#    #Iprim = fpi*gpiNN*Mp*M0/(Del_sq - Mp_sq)
#    K = 8.0 * const_K *  / (s - Mp_sq) / 3.0  #( pow(I,2) - (DelT_sq*pow(Iprim,2)/Mp_sq
#    var_K = 8*const_K/

#def _dsig_dq_gstar(Del_sq, s, qsq):
#    pbarn = 0.39e9 # 1gev^-2=3.9e9pbarn
#    xi = _xi_a_qdep(s,qsq)
#    DelT_sq = _DelT_sq(Del_sq,s)
#    I = fpi*gpiNN*Mp*(1-xi)*M0/(Del_sq - Mp_sq)/(1+xi)
#    Iprim = fpi*gpiNN*Mp*M0/(Del_sq - Mp_sq)
#    K = const_K * ( pow(I,2) - (DelT_sq*pow(Iprim,2)/Mp_sq) )
#    var_K = const_K/(s - Mp_sq)/pow(qsq,5)

def _dsig_dDel_sq(Del_sq,s):
    pbarn = 0.39e9 # 1gev^-2=3.9e9pbarn
    xi = _xi_a(s)
    DelT_sq = _DelT_sq(Del_sq,s)
    F1 = 1/16.0/pi/_lambda_sq(s,Mp_sq,Mp_sq)
    F2 = pow(C,2)*2.0*(1+xi)/4.0/xi/pow(Mbar, 8)
    I = fpi*gpiNN*Mp*(1-xi)*M0/(Del_sq - Mp_sq)/(1+xi)
    Iprim = fpi*gpiNN*Mp*M0/(Del_sq - Mp_sq)
    return pbarn * F1 * F2 * ( pow(I,2) - (DelT_sq*pow(Iprim,2)/Mp_sq) )

svals =    [12.25, 15,   16.8745,  20,    24.346]
#min_vals = [-0.45,   -0.5, -2.76,    -0.5,  -6.5]
#min_vals = [-0.092,   -0.5, -1.3,    -0.5,  -2.85]
min_vals = [-0.092,   -0.5, -1.0,    -0.5,  -1.0]
#min_vals = [-0.5,   -0.5,  -0.5,     -0.5,  -0.5]
#min_vals = [-0.3,   -0.45,  -0.45,     -0.45,  -0.45]
#max_vals = [0.6,    0.0,  0.46,     0.0,   0.32]
max_vals = [0.59,    0.0,  0.43,     0.0,   0.3]
#max_vals = [0.0,    0.0,  0.0,     0.0,   0.0]
cols =     ['r',     'g',  'b',      'c',   'k']

def plot_xi():

    costh_min=0.5
    costh_max=1.0
    npt = 10
    trange = []
    xi3 = []
    cols3 = []
    svals3 = []
    qvals3 = []

    rr=np.linspace(min_vals[0],max_vals[0],npt)
    trange.append(rr)
    xi3.append(map(lambda x: _xi_q(x, Mjpsi2, svals[0]), rr))
    cols3.append('k')
    svals3.append(svals[0])
    qvals3.append(Mjpsi2)

    rr=np.linspace(min_vals[2],max_vals[2],npt)
    trange.append(rr)
    xi3.append(map(lambda x: _xi_q(x, Mjpsi2, svals[2]), trange[1]))
    cols3.append('r')
    svals3.append(svals[2])
    qvals3.append(Mjpsi2)

    rr=np.linspace(min_vals[4],max_vals[4],npt)
    trange.append(rr)
    xi3.append(map(lambda x: _xi_q(x, Mjpsi2, svals[4]), trange[2]))
    cols3.append('b')
    svals3.append(svals[4])
    qvals3.append(Mjpsi2)

    nq=50
    for iq in range(0,nq+1):
        qval=3.0+(1.3*iq/nq)
        #print 'iq=%d, qval=%f'%(iq, qval)
        rr=np.linspace(_t(qval,costh_min,5.0),_t(qval,costh_max,5.0),npt)
        trange.append(rr)
        xi3.append(map(lambda x: _xi_q(x, qval, 5.0), rr))
        cols3.append('m')
        svals3.append(5.0)
        qvals3.append(qval)

    for iq in range(0,nq+1):
        qval=5.0+(4.0*iq/nq)
        #print 'iq=%d, qval=%f'%(iq, qval)
        rr=np.linspace(_t(qval,costh_min,10.0),_t(qval,costh_max,10.0),npt)
        trange.append(rr)
        xi3.append(map(lambda x: _xi_q(x, qval, 10.0), trange[3]))
        cols3.append('c')
        svals3.append(10.0)
        qvals3.append(qval)

    colors = []

    #cols3 =  ['r',     'g', 'b',     'c', 'm', 'k',  'y' ]
    #fig3, ax3 = plt.subplots(figsize=(12,8), dpi=100)
    fig3, ax3 = plt.subplots(figsize=(15,8))
    #fig3, ax3 = plt.subplots()
    for idx,_xi3 in enumerate(xi3):
        #print idx
        #print trange[idx]
        #print _xi3
        #ax3.plot(trange[idx], _xi3, cols3[idx], label=r'$s=%4.2g,Q^2=%4.2g$' %(svals3[idx], qvals3[idx]))
        if idx<3 or (idx-3)%(nq+1)==0:
            if idx < 3:
                line, = ax3.plot(trange[idx], _xi3, cols3[idx], label=r'$s \ = \ %4.2f \ GeV \ ^2 (\pi^{0}J/\psi)$' %(svals3[idx]), linewidth=3)
                colors.append(plt.getp(line,'color'))
            else:
                line, = ax3.plot(trange[idx], _xi3, cols3[idx], label=r'$s \ = \ %4.1f \ GeV \ ^2 (\pi^{0}\gamma^*)$' %(svals3[idx]), linewidth=3)
                colors.append(plt.getp(line,'color'))
            if idx==2 or idx==1:
                plt.text(trange[idx][0],_xi3[0]+0.01,r'$Q^2=%4.2g \ GeV \ ^2$' % (qvals3[idx]), fontsize=20)
            elif (idx-3)/(nq+1)==0:
                plt.text(trange[idx][npt-1]+0.04,_xi3[npt-1]-0.03,r'$Q^2=%4.2g \ GeV \ ^2$' % (qvals3[idx]), fontsize=20)
            else:
                plt.text(trange[idx][0]-0.5,_xi3[0],r'$Q^2=%4.2g \ GeV \ ^2$' % (qvals3[idx]), fontsize=20)
        else:
            if ((idx-3)%(nq+1))-(nq)==0:
                ax3.plot(trange[idx], _xi3, cols3[idx], linewidth=3)
                if (idx-3)/(nq+1)==0:
                    plt.text(trange[idx][npt-1]+0.05,_xi3[npt-1]-0.02,r'$Q^2=%4.2g \ GeV \ ^2$' % (qvals3[idx]), fontsize=20, backgroundcolor='w')
                else:
                    plt.text(trange[idx][0]-0.5,_xi3[0],r'$Q^2=%4.2g \ GeV \ ^2$' % (qvals3[idx]), fontsize=20)
            else:
                ax3.plot(trange[idx], _xi3, cols3[idx])


    plt.ylabel(r"$\xi$", fontdict=font);
    plt.xlabel(r"$t[GeV^2]$",fontdict=font);
    #plt.axis([-0.51,0.01,50,250]);
    #plt.yticks(np.arange(0.1, 1.0, 0.1))
    #plt.xticks(np.arange(-3.0, 1.5, 0.5))
    legend3 = ax3.legend(loc='best', shadow=True);

    frame = legend3.get_frame();
    frame.set_facecolor('0.96')

    for label in legend3.get_lines():
        label.set_linewidth(2.0)
    for tick in ax3.xaxis.get_major_ticks():
        tick.label.set_fontsize(22)
    for tick in ax3.yaxis.get_major_ticks():
        tick.label.set_fontsize(22)

    #style_axis(ax3)

    for color,text in zip(colors,legend3.get_texts()):
        text.set_color(color)
        text.set_fontsize(30);

    plt.show()
    fig3.savefig('xi_range.pdf')

def color_labels():
    x = np.arange(10)

    fig = plt.figure()
    ax = plt.subplot(111)

    colors = []
    for i in xrange(5):
        line, = ax.plot(x, i * x, label=r'$y = %ix\psi$' % i)
        colors.append(plt.getp(line,'color'))

        leg = ax.legend()

    for color,text in zip(colors,leg.get_texts()):
        text.set_color(color)

    plt.show()

def numevt(msv):
    lumi = [95.5,103.0,108.0] if msv else [2e3,2e3,2e3]
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
    #normr = [ (-0.443789, 0.616486), (-0.5, 0.457248), (-0.5, 0.31538) ] # normalization winodw in Delta^2
    #normr = [ (-0.092, 0.59), (-1.3, 0.43), (-2.85, 0.3) ] # normalization winodw in Delta^2
    normr = [ (-0.092, 0.59), (-1.0, 0.43), (-1.0, 0.3) ] # normalization winodw in Delta^2
    fullr = [ (-0.443789, 0.616486), (-2.75368, 0.457248), (-6.48861, 0.31538)] # full range from 0 < cost < 1
    print s_
    for idx, s in enumerate(s_):
        print idx,s
        sig_normr.append(quad(_dsig_dDel_sq, normr[idx][0], normr[idx][1] , args=(s))[0])
        sig_fullr.append(quad(_dsig_dDel_sq, fullr[idx][0], fullr[idx][1] , args=(s))[0])
        f_t.append(sig_normr[idx]/sig_fullr[idx])
        nevt_normr.append(sig_normr[idx]*lumi[idx]*br*2)
        nevt_fullr.append(sig_fullr[idx]*lumi[idx]*br*2)
        print("s=%4.2f  xs(nr)=%4.2f nevt(nr)=%4.2f xs(fr)=%4.2f nevt(fr)=%4.2f"% (s, 2*sig_normr[idx], nevt_normr[idx], 2*sig_fullr[idx], nevt_fullr[idx]))

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

def style_axis(_ax):
    _legend = _ax.legend(loc=2, shadow=True);
    _frame = _legend.get_frame();
    _frame.set_facecolor('0.90')
    for label in _legend.get_texts():
        label.set_fontsize(30);
    for label in _legend.get_lines():
        label.set_linewidth(2.0)
    for tick in _ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(22)
    for tick in _ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(22)

# test if pi0 mass hypothesis distorts Del^2 distributions too much
def plot_mpi0_test(iplab):
    global Mpi0
    global Mpi0_sq
    npt = 100
    _val = svals[iplab]
    s = np.full(npt, _val)

    fig, ax = plt.subplots(figsize=(12,8), dpi=1000)

    Mpi0_tmp = Mpi0

    Mpi0 = 0.0
    Mpi0_sq = Mpi0*Mpi0
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

    style_axis(ax)

    Mpi0 = Mpi0_tmp
    Mpi0_sq = Mpi0*Mpi0

#plotting
def plot_del(ax, sel):
    #ax.set_yscale("log")
    npt = 100
    #svals = [12.2649, 15,  16.8745, 20,  24.346,    40,   50, 60,  70,  100, 150, 300, 600, 1000, 2000, 2500, 5000, 10000 ]
    #cols =  ['r',     'g', 'b',     'c', 'm', 'r',  'g', 'b', 'c', 'm', 'r', 'g', 'b', 'c', 'm',  'r',  'g',  'b',  'c',  ]
    for idx, _val in enumerate(svals):
        if sel==1 and (idx==1 or idx==3): continue
        if sel==2 and (idx==0 or idx==2 or idx==4): continue
        if idx == 1 or idx == 3:
            Del_sq_min = _Del_sq_DelT_sq(-0.5,_val)
        else:
            Del_sq_min = min_vals[idx]
        #Del_sq_min = _Del_sq_DelT_sq(-2.0,_val)
        Del_sq_max = _Del_sq_max(_val)
        Del_sq=np.linspace(Del_sq_min,Del_sq_max,npt)
        #Del_sq=np.linspace(-0.6,0.5,npt)
        s = np.full(npt, _val)
        DelT_sq=map(_DelT_sq,Del_sq,s)
        xsect=map(_dsig_dDel_sq,Del_sq,s)
        ax.plot(DelT_sq, xsect, cols[idx], label=r'$s \ = \ %4.2f \ GeV \ ^{2}$' % _val, linewidth=3)
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
    #print del2_ranges;
    del2_mean = []
    xsect_int = []
    for iplab, blist in enumerate(del2_ranges):
        xs = []
        d2m = []
        for irange, bound in enumerate(blist):
            result = quad( _dsig_dDel_sq, bound[0], bound[1], args=(svals2[iplab]));
            #print("plab= %4.2f, sig_tot(%4.2f<Delta^2<%4.2f) = %4.2f pm %4.2f" %(svals2[iplab], bound[0], bound[1], result[0], result[1]) )
            d2m.append((bound[0]+bound[1])/2.)
            xs.append(result[0])
        del2_mean.append(d2m)
        xsect_int.append(xs)

    fig, axes = plt.subplots(1,3,figsize=(14,6))
    #style_axis(axes)
    for idx, ax in enumerate(axes):
        ax.plot(del2_mean[idx], xsect_int[idx], 'bo' if idx is 0 else 'rs' if idx is 1 else 'm^')

def plot(idx=0, sel=0):
    #fig, (ax1, ax2) = plt.subplots(1,2)
    fig, ax = plt.subplots(figsize=(12,8), dpi=1000)
    if idx==0:
        plot_del(ax,sel)
    else:
        plot_s(ax)

    legend = ax.legend(loc=2, shadow=True);
    frame = legend.get_frame();
    frame.set_facecolor('0.90')

    for label in legend.get_texts():
        label.set_fontsize(30);
    for label in legend.get_lines():
        label.set_linewidth(2.0)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(22)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(22)

    plt.show()
    #fig.savefig('tda_xsect_comp.pdf')

if __name__ == "__main__":
    plot()
