import numpy as np
from aux_functions import dlydlx
from scipy.optimize import fsolve
#from scipy.interpolate import interp1d
from aux_functions import k_b, mu, m_p, Grav, pi, sig_h2

def reconstruct_size_distribution(r,a,t,sig_g,sig_d,alpha,rho_s,T,M_star,v_f,a_0=1e-4,fix_pd=None):
    if fix_pd is not None: print('WARNING: fixing the inward diffusion slope')
    cs     = np.sqrt(k_b*T/mu/m_p)
    om     = np.sqrt(Grav*M_star/r**3)
    vk     = r*om
    gamma  = dlydlx(r,sig_g*np.sqrt(T)*om)
    p      = -dlydlx(r,sig_g)
    q      = -dlydlx(r,T)
    #
    # fragmentation size
    #
    b      = 3.*alpha*cs**2/v_f**2
    a_fr   = sig_g/(pi*rho_s)*(b-np.sqrt(b**2-4.))
    a_fr[np.isnan(a_fr)] = np.inf
    #
    # drift size
    #
    a_dr   = 0.55*2/pi*sig_d/rho_s*r**2.*(Grav*M_star/r**3)/(abs(gamma)*cs**2)
    #_a_drf = interp1d(np.log10(r),np.log10(a_dr))
    #a_drf  = lambda r: 10.**_a_drf(np.log10(r))
    #
    # transition to turblent velocities
    #
    Re = alpha*sig_g*sig_h2/(2*mu*m_p)
    a_bt = (8*sig_g/(pi*rho_s)*Re**-0.25*np.sqrt(mu*m_p/(3*pi*alpha))*(4*pi/3*rho_s)**-0.5)**0.4
    #
    # time dependent growth 
    #
    t_grow = sig_g/(om*sig_d)
    a_grow = a_0*np.exp(t/t_grow)
    #
    # the minimum of all of those
    #
    a_max  = np.minimum(np.minimum(a_fr,a_dr),a_grow)
    #
    # #################################
    # now fill in the size distribution
    # #################################
    #
    sig_1 = 1e-100*np.ones([len(a),len(r)]) # the fragment distribution
    sig_2 = 1e-100*np.ones(sig_1.shape)     # outward mixed fragments
    sig_3 = 1e-100*np.ones(sig_1.shape)     # drift distribution diffused
    #
    # find the radius where fragmentation stops
    #
    ir_f   = np.where(a_max==a_fr)[0]
    if len(ir_f)==0:
        ir_f = 0
    else:
        ir_f = ir_f[-1]
    r_f    = r[ir_f]
    #
    # FIRST: the fragment distribution
    #
    for ir in range(ir_f+1):
        ia_max = abs(a-a_max[ir]).argmin()
        i_bt = abs(a-a_bt[ir]).argmin()
        dist = (a/a[0])**1.5
        dist[i_bt:] = dist[i_bt]*(a[i_bt:]/a[i_bt])**0.375 # can be 1/4, or 1/2, or average 0.375, OO13 used 1/4
        dist[ia_max+1:] = 1e-100
        dist=dist/sum(dist)*sig_d[ir]
        sig_1[:,ir] = dist
    #
    # SECOND:
    # outward diffusing fragments
    #
    for ia in range(len(a)):
        #
        # outward diffusion from r_f
        #
        St = a[ia]*rho_s/sig_g*pi/2
        # the factor in front of alpha is experimental XXX
        i_mx = np.where(2*alpha/St/abs(gamma)>1-r_f/r)[0][-1]
        mask = np.arange(ir_f+1,i_mx+1)
        sig_2[ia,mask] = sig_1[ia,ir_f]*sig_g[mask]/sig_g[ir_f]*(r[mask]/r[ir_f])**-1.5
    #
    # THIRD:
    # the diffusive parts around the drift size
    #
    i_a_max= np.where(a>=a_max[ir_f])[0][0]
    for ia in range(i_a_max):
        #
        # find radius and surface density where a == a_max
        #
        _ir   = np.where((a_max<a[ia+1]) & (a_max>=a[ia]))[0]
        if len(_ir)==0:
            _ir = abs(a_max-a[ia]).argmin()
        else:
            _ir = _ir[-1]
        mask  = np.arange(max(0,_ir-2),min(len(r),_ir+3))
        f     = lambda _r: 10.**np.interp(np.log10(_r),np.log10(r[mask]),np.log10(a_max[mask]))-a[ia]
        _r    = fsolve(f,r[_ir])[0]
        _sigd = 10.**np.interp(np.log10(_r),np.log10(r),np.log10(sig_d+1e-100))
        _sigg = 10.**np.interp(np.log10(_r),np.log10(r),np.log10(sig_g+1e-100))
        #
        # inward diffusion from a_max
        #
        mask = np.arange(ir_f+1,_ir+1)
        #
        # analytical estimate (everything is a power-law):
        # v is approximated as v0 * (r/_r)**d
        #
        _ir2           = abs(r-0.5*_r).argmin()
        v0             = -a[ia]*rho_s*pi/(2*sig_g)*cs**2/vk*(p+(q+3.)/2.)
        d              = (p-q+0.5)
        v              = v0[_ir]*(r/_r)**d[_ir]
        pd_est         = 1./(2*alpha)*(v/cs*vk/cs-2*p*alpha+vk/cs*np.sqrt(   (v/cs)**2+4*(1+d-p)*v/vk*alpha+4*alpha*sig_d/sig_g  ))
        if fix_pd is not None:
            pd_est = fix_pd
        else:
            pd_est = pd_est[_ir2]
        sig_3[ia,mask] = _sigd*(r[mask]/_r)**pd_est
        #
        # outward diffusion
        #
        A=a[ia]*rho_s*pi*gamma[_ir]/(2*alpha[_ir]*_sigg*p[_ir])
        sig_3[ia,_ir+1:] = _sigd*np.exp(A*((r[_ir+1:]/_r)**p[_ir+1]-1.))
    #
    # normalize
    #
    for _ir in range(ir_f+1,len(r)):
        sig_3[:,_ir] = sig_3[:,_ir]/sum(sig_3[:,_ir])*(sig_d[_ir]-sum(sig_2[:,_ir]))
    #
    # add all and do a final normalization
    #
    sig_sol = sig_1+sig_2+sig_3
    #for _ir in range(len(r)):
    #    sig_sol[:,_ir] = sig_sol[:,_ir]/sum(sig_sol[:,_ir])*sig_d[_ir]
    return sig_sol,a_max,r_f,sig_1,sig_2,sig_3