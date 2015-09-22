#!/usr/bin/env python
# coding: utf-8
from matplotlib.pyplot import gca,figure,tight_layout
from numpy import logspace,ones,abs,log10,where
from uTILities import better_plots

def mynote(x,y,dx,dy,text,c='k',ls='-',lw=1,align='tc',fontdict={},offset=0.0):
    gca().plot([x,x+dx],[y,y+dy],ls=ls,lw=lw,c=c)
    
    fontdict['color'] = c
    
    if align[0]=='t':
        va = 'top'
    elif align[0]=='b':
        va = 'bottom'
    elif align[0]=='c':
        va = 'center'
        
    if align[1]=='l':
        ha = 'left'
    elif align[1]=='r':
        ha= 'right'
    elif align[1]=='c':
        ha= 'center'
        
    if align=='tr':
        xoff = -abs(offset)
        yoff = 0
    if align=='tl':
        xoff = abs(offset)
        yoff = 0
    if align=='bl':
        xoff = 0
        yoff = abs(offset)
    if align=='br':
        xoff = 0
        yoff = -abs(offset)
        
        
    gca().text(x+dx+xoff,y+dy+yoff,text,fontdict={'color':c},verticalalignment=va,horizontalalignment=ha)
    
    gca().plot(x,y,marker='o',ms=8,mfc=c,mec='None')


def main():
    r       = logspace(-1,3,100)
    af      = (r/1.)**-.5
    rf_out  = 100
    ad      = 5.*(r/1.)**-(5./4.)
    rf      = r[abs(af-ad).argmin()]
    a_width = 20.
    r_width = 3.
    r_arrow = 5.
    ofs_fct = 15.
    
    better_plots(fs=20,sans=False)
    #cols = rcParams['axes.color_cycle']
    cols = ['#31a354','#e6550d','#08519c']
    f    = figure()
    ax   = gca()
    
    mask = r<=rf_out
    ax.loglog(r[mask],af[mask],lw=5,c=cols[1])
    ax.loglog(r,ad,lw=5,c=cols[0])
    
    ax.set_ylim(1e-4,1e2)
    
    ax.plot(rf*ones(2),ax.get_ylim(),'k--',lw=3)
    #
    # fragmentation area
    #
    mask = r<=rf
    ax.fill_between(r,af,y2=ax.get_ylim()[0]*ones(len(r)),        where=mask,alpha=0.7,color=cols[1])
    #
    # inward diffusion area
    #
    mask = r>=rf
    ax.fill_between(r,ad,y2=ad/a_width,                           where=mask,alpha=0.7,color=cols[0])
    #
    # outward diffusion area
    #
    mask = (r>=rf) & (r<=10.**((log10(rf)+3*log10(rf_out))/4.))
    ax.fill_between(r,ad/a_width,y2=ax.get_ylim()[0]*ones(len(r)),where=mask,alpha=0.7,color=cols[2])
    ax.annotate('II',xy=(10.**log10(r[where(mask)[0][[0,-1]]]).mean(),4e-4),xycoords='data',color='k',horizontalalignment='center')
    #
    # depleted region
    #
    mask = r>=10.**((log10(rf)+3*log10(rf_out))/4.)
    mask[where(mask==False)[0][-1]] = True # this is needed to avoid gaps
    ax.fill_between(r,ad/a_width,y2=ax.get_ylim()[0]*ones(len(r)),where=mask,alpha=0.7,color='0.5')
    ax.annotate('IV',xy=(2*r[where(mask)[0][0]],2e-4),xycoords='data',color='k',horizontalalignment='center')
    #
    # a_f annotation
    #
    i = abs(r-10.**(0.5*(log10(rf)+log10(rf_out)))).argmin()
    mynote(r[i],af[i],0,ax.get_ylim()[-1]/r_width,'$a_\mathrm{f}$',c='k',align='tr',fontdict={},offset=r[i]/ofs_fct) # bc
    #
    # r_f annotation
    #
    mynote(rf,1e1,-rf/1.5,0,'$r_\mathrm{f}$',c='k',align='bl',fontdict={},offset=rf/ofs_fct )
    #
    # a_d annotation
    #
    i = abs(r-10.**(0.5*(log10(rf)+log10(r[0])))).argmin()
    mynote(r[i],ad[i],0,ax.get_ylim()[-1]/r_width,'$a_\mathrm{d}$',c='k',align='tr',fontdict={},offset=r[i]/ofs_fct)
    #
    # fragmentation arrow down
    #
    opt = dict(color=cols[1],arrowstyle = 'simple,head_width=.75,head_length=.75',connectionstyle = 'arc3,rad=0')
    gca().annotate('',xy=(r[i],1e-3),xycoords='data',xytext =(r[i],af[i]),textcoords = 'data',arrowprops=opt,size=20)
    ax.annotate('I',xy=(r[i]*1.3,10.**(log10([af[i],1e-3]).mean())),xycoords='data',color='k',horizontalalignment='center')
    #
    # outward diffusion arrow
    #
    opt['color']=cols[2]
    gca().annotate('',xy=(r_arrow*rf,1e-3),xycoords='data',xytext =(rf,1e-3),textcoords = 'data',arrowprops=opt,size=20)
    #
    # inward diffusion arrow
    #
    i = abs(r-10.**((log10(rf)+2*log10(r[-1]))/3.)).argmin()
    opt['color']=cols[0]
    gca().annotate('',xy=(r[i]/r_arrow,ad[i]),xycoords='data',xytext =(r[i],ad[i]),textcoords = 'data',arrowprops=opt,size=20)
    ax.annotate('III',xy=(r[i]/r_arrow*2,ad[i]/2.5),xycoords='data',color='k')
    #
    # other styling
    #
    ax.minorticks_off()
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.set_xticklabels('')
    ax.set_yticklabels('')
    ax.set_xlabel('radius')
    ax.set_ylabel('grain size')
    #
    # save and show figure
    #
    tight_layout()
    f.savefig('fig1.pdf')
    return f
    
if __name__=='__main__':
    main()