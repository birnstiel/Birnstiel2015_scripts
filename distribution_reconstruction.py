# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from scipy.interpolate   import interp1d
print('New implementation')

def reconstruct_size_distribution(r,a,t,sig_g,sig_d,alpha,rho_s,T,M_star,v_f,a_0=1e-4,fix_pd=None):
    """
    Reconstructs the approximate size distribution based on the recipe of Birnstiel et al. 2015, ApJ.
    
    Arguments:
    ----------
    
    r : array
    :    radial grid [cm]
    
    a : array
    :    grain size grid, should have plenty of range and bins to work [cm]
    
    t : float
    :    time of the snapshot [s]
    
    sig_g : array
    :    gas surface densities on grid r [g cm^-2]
    
    sig_d : array
    :    dust surface densities on grid r [g cm^-2]
    
    alpha : array
    :    alpha parameter on grid r [-]
    
    rho_s : float
    :    material (bulk) density of the dust grains [g cm^-3]
    
    T : array
    :    temperature on grid r [K]
    
    M_star : float
    :    stellar mass [g]
    
    v_f : float
    :    fragmentation velocity [cm s^-1]
    
    Keywords:
    ---------
    
    a_0 : float
    :    initial particle size [cm]
    
    fix_pd : None | float
    :    float: set the inward diffusion slope to this values
         None:  calculate it
         
    Output:
    -------
    sig_sol,a_max,r_f,sig_1,sig_2,sig_3
    
    sig_sol : array
    :    2D grain size distribution on grid r x a
    
    a_max : array
    :    maximum particle size on grid r [cm]
    
    r_f : float
    :    critical fragmentation radius (see paper) [cm]
    
    sig_1, sig_2, sig_3 : array
    :    grain size distributions corresponding to the regions discussed in the paper
    """
    import warnings
    import numpy as np
    from aux_functions import dlydlx, k_b, mu, m_p, Grav, pi, sig_h2
    from uTILities import trace_line_though_grid
    from scipy.integrate import cumtrapz
    
    if fix_pd is not None: print('WARNING: fixing the inward diffusion slope')    
    warnings.warn('TODO: limit outward diffusion by drift')
    warnings.warn('TODO: apply the time limit to the maximum particle size')
    #
    # calculate derived quantities
    #
    alpha  = alpha*np.ones(len(r))
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
    #
    # time dependent growth 
    #
    t_grow = sig_g/(om*sig_d)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'overflow encountered in exp')
        a_grow = a_0*np.exp(t/t_grow)
    #
    # the minimum of all of those
    #
    a_max  = np.minimum(np.minimum(a_fr,a_dr),a_grow)
    if a_max.max()>a[-1]: raise ValueError('Maximum grain size larger than size grid. Increase upper end of size grid.')
    #
    # transition to turblent velocities
    #
    Re = alpha*sig_g*sig_h2/(2*mu*m_p)
    a_bt = (8*sig_g/(pi*rho_s)*Re**-0.25*np.sqrt(mu*m_p/(3*pi*alpha))*(4*pi/3*rho_s)**-0.5)**0.4
    #
    # apply the reconstruction recipes
    #
    sig_1     = 1e-100*np.ones([len(a),len(r)]) # the fragment distribution
    sig_2     = 1e-100*np.ones(sig_1.shape)     # outward mixed fragments
    sig_3     = 1e-100*np.ones(sig_1.shape)     # drift distribution diffused
    prev_frag = None
    #
    # Radial loop: to fill in all fragmentation regions
    #
    frag_mask = a_max==a_fr
    frag_idx  = np.where(frag_mask)[0]
    for ir in range(len(r)):
        ia_max = abs(a-a_max[ir]).argmin()
        if frag_mask[ir]:
            #
            # fragmentation case
            #
            i_bt = abs(a-a_bt[ir]).argmin()
            dist = (a/a[0])**1.5
            dist[i_bt:] = dist[i_bt]*(a[i_bt:]/a[i_bt])**0.375 # can be 1/4, or 1/2, or average 0.375, OO13 used 1/4
            dist[ia_max+1:] = 1e-100
            dist=dist/sum(dist)*sig_d[ir]
            sig_1[:,ir] = dist
            prev_frag = ir
        elif sig_g[ir]>1e-5:
            #
            # outward diffusion of fragments from index prev_frag
            #
            if prev_frag is not None:
                sig_2[:,ir] = sig_1[:,prev_frag]*sig_g[ir]/sig_g[prev_frag]*(r[ir]/r[prev_frag])**-1.5
                #
                # limit outward diffusion by drift
                #
                for ia in range(len(a)):
                    St    = a[ia]*rho_s/sig_g*pi/2.
                    vd    = 1./(St+1./St)*cs**2/vk*gamma
                    t_dri = r/abs(vd)
                    t_dif = (r-r[prev_frag])**2/(alpha*cs**2/om/(1.+St**2))
                    #
                    # the factor of 5 below changes how far out we take the
                    # diffusion, but doesn't affect the results too much
                    #
                    i_out = np.where(t_dif<5*t_dri)[0]
                    if len(i_out)!=0:
                        i_out=i_out[-1]
                        sig_2[ia,i_out:]=1e-100
                    else:
                        sig_2[ia,prev_frag+1:]=1e-100
    #
    # ---------------------
    # add up all the the radial approximations from all cells crossed by the drift limit
    # ---------------------
    #
    # find all intersected grid cells
    # we will assume that the interfaces are in the middle of the grid center
    # usually it's the other way around, but this doesn't really matter here
    #
    ri  = 0.5*(r[1:]+r[:-1])
    ai  = 0.5*(a[1:]+a[:-1])
    ri  = np.hstack((r[0]-(ri[0]-r[0]),ri,r[-1]+(r[-1]-ri[-1])))
    ai  = np.hstack((a[0]-(ai[0]-a[0]),ai,a[-1]+(a[-1]-ai[-1])))
    f   = interp1d(np.log10(np.hstack((0.5*ri[0],r))),np.log10(np.hstack((a_max[0],a_max))),bounds_error=False,fill_value=np.log10(a_max[-1]))
    res = trace_line_though_grid(np.log10(ri),np.log10(ai),f)
    
    for cell in res:
        ir,ia = cell
        #
        # if fragmentation limited: skip this cell
        #
        if frag_mask[ir]: continue
        #
        # if a_max increases outward in the inner disk, skip as well
        #
        if ir in [0,1,2,3] and a_max[ir+1]>a_max[ir]: continue
        #
        # find ra, our starting point, and dust and gas densities there
        #
        mask  = np.arange(max(ir-2,0),min(ir+2,len(r)-1)+1)
        mask  = mask[a_max[mask].argsort()] # needs to be sorted to work for interpolation?
        ra    = 10.**np.interp(np.log10(a[ia]),np.log10(a_max[mask]),np.log10(r[mask]))
        _sigd = 10.**np.interp(np.log10(ra),np.log10(r),np.log10(sig_d+1e-100))
        _sigg = 10.**np.interp(np.log10(ra),np.log10(r),np.log10(sig_g+1e-100))
        #
        # find previous and next fragmentation index
        #
        prev_frag = frag_idx[frag_idx<ir]
        if len(prev_frag)==0:
            prev_frag = 0
        else:
            prev_frag = prev_frag[-1]
        next_frag = frag_idx[frag_idx>ir]
        if len(next_frag)==0:
            next_frag = len(r)-1
        else:
            next_frag = next_frag[0]
        mask = np.arange(prev_frag+1,min(ir+1,len(r)))
        #
        # analytical estimate (everything is a power-law):
        # v is approximated as v0 * (r/ra)**d
        #
        v0     = -a[ia]*rho_s*pi/(2*sig_g)*cs**2/vk*(p+(q+3.)/2.)
        d      = (p-q+0.5)
        v      = v0[ir]*(r/ra)**d[ir]
        pd_est = 1./(2*alpha)*(v/cs*vk/cs-2*p*alpha+vk/cs*np.sqrt(   (v/cs)**2+4*(1+d-p)*v/vk*alpha+4*alpha*sig_d/sig_g  ))
        if fix_pd is not None:
            _pd_est = fix_pd
        else:
            _pd_est = pd_est[ir]        
        #sig_3[ia,mask] += _sigd*(r[mask]/ra)**_pd_est
        # integrate over slowly varying power-law
        sol = np.exp( np.minimum(709.,cumtrapz(pd_est[mask],x=np.log10(r[mask]),initial=0)))
        sol = sol/np.interp(ra,r[mask],sol)*_sigd
        sig_3[ia,mask] = np.maximum(sol,sig_3[ia,mask])
        #
        # outward diffusion
        #
        A = a[ia]*rho_s*pi*gamma[ir]/(2*alpha[ir]*_sigg*p[ir])
        sig_3[ia,ir+1:] = np.maximum(sig_3[ia,ir+1:] , _sigd*np.exp(A*((r[ir+1:]/ra)**p[ir+1]-1.)))
    #
    # normalize sig_3 (where it contributes).
    # as it can overlap with sig_2, but sig_2 density is normalized,
    # we will adapt the normalization to be sig_d-sig_2
    #
    sig_3_sum = sig_3.sum(0)
    mask      = sig_3_sum>1e-50
    sig_3[:,mask] = sig_3[:,mask]/sig_3_sum[mask]*(sig_d[mask]-sig_2.sum(0)[mask])
    #
    # add up all of them and normalize
    #
    sig_dr = sig_1+sig_2+sig_3
    sig_dr = sig_dr/sig_dr.sum(0)*sig_d
    
    return sig_dr,a_max,r[frag_idx[-1]],sig_1,sig_2,sig_3   
        
def test_reconstruction():
    from uTILities import trace_line_though_grid
    from constants import M_sun, AU, year, pi
    import numpy as np
    #
    # ================
    # set up the model
    # ================
    #
    nr = 100
    na = 200
    ri = np.logspace(-1,3,nr+1)*AU
    ai = np.logspace(-4,2,na+1)
    
    r = 0.5*(ri[1:]+ri[:-1])
    a = 0.5*(ai[1:]+ai[:-1])
    
    M_star = M_sun
    M_disk = 0.05 * M_star
    rc     = 60 * AU
    eps    = 0.01
    rho_s  = 1.2
    v_f    = 1e3
    alpha  = 1e-3
    
    sig_g  = r**-1*np.exp(-r/rc)
    sig_g  = M_disk*sig_g/np.trapz(2*pi*r*sig_g,x=r)
    sig_d  = sig_g*eps
    T      = 200*(r/AU)**-0.5
    #
    # call the reconstruction routine
    #
    sig_dr,a_max,_,_,_,_ = reconstruct_size_distribution(r,a,1e6*year,
        sig_g,sig_d,alpha,rho_s,T,M_star,v_f)
    #
    # ========
    # PLOTTING
    # ========
    #
    _,ax = plt.subplots()
    # 
    # plot the grid
    # 
    for _ai in ai: ax.plot(ri/AU,_ai*np.ones(len(ri)),'k')
    for _ri in ri: ax.plot(_ri/AU*np.ones(len(ai)),ai,'k')
    for _a in a: ax.plot(r/AU,_a*np.ones(len(r)),'kx')
    #
    # find all intersected grid cells
    #
    f   = interp1d(np.log10(r),np.log10(a_max),bounds_error=False,fill_value=np.log10(a_max[-1]))
    res = trace_line_though_grid(np.log10(ri),np.log10(ai),f)
    #
    # draw all intersected cells
    #
    for j,i in res:
        ax.plot(np.array([ri[j],ri[j+1],ri[j+1],ri[j],ri[j]])/AU,[ai[i],ai[i],ai[i+1],ai[i+1],ai[i]],'r')
    #
    # plot the distribution
    #
    mx = np.ceil(np.log10(sig_dr.max()))
    ax.contourf(r/AU,a,np.log10(sig_dr),np.arange(mx-10,mx+1),alpha=0.5)
    #
    # plot maximum particle size
    #
    ax.plot(r/AU,a_max,'r-+')
    
    ax.set_xlim(0.9*ri[0]/AU,1.1*ri[-1]/AU)
    ax.set_ylim(0.9*ai[0],1.1*ai[-1])
            
    ax.set_xscale('log')
    ax.set_yscale('log')