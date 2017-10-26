#!/usr/bin/python2.7

import pdfpaper as pdf
import vid 
import matplotlib.pyplot    as plt
import matplotlib.gridspec  as grd
import numpy as np
import glob
from collections import defaultdict
from matplotlib.ticker import LogFormatter
from matplotlib.colors import LogNorm

plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':8})
plt.rc('text', usetex=True)
printmarg    = False
printjoint   = False
printmetrics = False
printvideo   = False
print1stvid  = False
printorient  = False

########## INPUTS ##########

filepath='./data/'
filt = 4   # default: 4  

pdftypes = ['dns','bvb5','bvb4','dirichlet','betadelta','sml']
switches = ['123'  ,'123' ,'123' ,'123','123' ,'123'       ]
suffixes = ['I1']
metrics  = ['l2norm', 'l1norm', 'yvar', 'covar-xy', 'covar-xz', 'covar-yz', 'xskew', 'yskew','DJS','DJSabs','success'] 
times    = ['0.19005E+00'] 

printvideo   = True

exact_integ = ['']

############################

# set up dictionary to store parameter guesses (only used for sml and bvb5) and labels for PDFtypes
guesses = defaultdict(dict)
for suffix in suffixes :
    for pdftype in pdftypes :
        guesses[suffix][pdftype] = [1,1,1,1,1,0,0]
names      = ['dns','bvb4','bvb5','dirichlet','betadelta','sml','sml3','sml4','smlc','bvb6']
labels     = ['DNS','CM-1' ,'BVB5-12','Dirichlet','Beta-2-Delta','SMLD','SMLD-3','SMLD-4a','SMLD-4b','BVB6']
namecolors = ['#000000','#00FFFF','#CC00FF','#0000FF','#FF0000','#00DD64','#0000FF','#FF0000','#00DD64','#DD0064']
names2labels = dict(zip(names,labels))
names2colors = dict(zip(names,namecolors))
colormap     = pdf.genColormap()
suffix2suffix = dict(zip( 
    ['I1','I2','I3','I4','I5','I6','P1' ,'P2' ,'P3' ,'P4','P5' ,'J1','J2','J3','J4'],
    ['I1','I3','I6','I4','I5','I2','PI1','PI2','PI3','P4','PL1','L3','L4','L1','L2']))

# find ties where data is present to run calculations
times_all = glob.glob(filepath + '/joint*')
for i in range(len(times_all)) : 
    times_all[i] = times_all[i][times_all[i].index('0.'):] 
times_all = sorted(times_all, key=float)
if times[0] == 'all' :
    times = times_all[times_all.index(times[1]):times[3]:times[2]]

# Prepare metrics file
if printmetrics :
    metrics = dict( zip( metrics, ([0]*len(pdftypes) for i in metrics) ) ) 
    line = ''.join('{:>14s}'.format(k) for k in (ty+'-'+sw for ty,sw in zip(pdftypes,switches)))
    line = '{:>14s}'.format('time') + line + '\n'
    for suffix in suffixes :
        for metric in metrics :
            filename = 'metrics/' + metric + '_' + suffix + '.txt'
            with open(filename,'w') as f:
                f.write(line)

# Loop over requested times, suffixes, and pdftypes
for time in times :
    # find the index of the time to get correct ensight file
    itime = times_all.index(time) + 1
    for jj in range(len(suffixes)) :
        suffix = suffixes[jj]
        print suffix, time
        varx='ZMIXA' + suffix
        vary='ZMIXB' + suffix

        # prepare for figures
        marginalx = 0
        marginaly = 0
        joint     = 0
        figsizejoint = (6.5,1.15)
        marginalx, axx = plt.subplots(1,3,sharey=True,figsize=(6.5,1.8))
        marginaly, axy = plt.subplots(1,3,sharey=True,figsize=(6.5,1.8))
        joint     = plt.figure('joint',figsize=figsizejoint)
        video     = plt.figure('video',figsize=( 12,5))
        orient    = plt.figure('orient',figsize=(6.5,1))
        nrows=2
        ncols=5

        dns = pdf.dnsdata([])

        fmean = 0; fvar = 0; covf = 0; zskew = 0

        # loop over PDF types, calculating PDF for the given condition and preparing the plots and metrics
        for ii in range(len(pdftypes)) :
            pdftype = pdftypes[ii]
            switch = switches[ii]
            #if pdftype == 'bvb6' and time == times[0]:
            #    guesses[suffix]['bvb6'][:] = 
            pdfobj = pdf.pdf(pdftype,filepath=filepath,time=time,varx=varx,vary=vary,DNSinit='true',filt=filt,switch=switch)
            pdfobj.fmean = fmean; pdfobj.fvar = fvar; pdfobj.covf = covf; pdfobj.zskew = zskew
            pdfobj.exact_integ = exact_integ
            pdfobj.PDFlabel = names2labels[pdfobj.PDFtype]
            if printorient  : guesses[suffix][pdftype] = [1,1,1,1,1,0,0]
            if printorient and pdftype=='bvb6' :
                guesses[suffix][pdftype] =  [ 121.6468291,124.35420632,3.27890227,  0, 3.70618777, 178.964164 , -237.21987634]#12
                #guesses[suffix][pdftype] = [4,125,122,-240,-1.67,178,0] #23
                #guesses[suffix][pdftype] = [ 121.6468291,3.27890227, 124.35420632,  0 ,-237.21987634, 178.964164,3.70618777 ] # 13

            pdfobj.calcPDFfull(guesses[suffix][pdftype])
            #guesses[suffix][pdftype][:] = pdfobj.params[:]
            if printmetrics :
                pdfobj.calcMetrics(dns)
                for metric in metrics: metrics[metric][ii] = pdfobj.metrics[metric]
            if printjoint :
                pdfobj.colormap = colormap
                pdfobj.plotPDF(joint,ii,'paper',ncols=ncols,nrows=nrows)
            if printorient :
                pdfobj.colormap = colormap
                pdfobj.plotPDF(orient,ii,'orient',ncols=ncols,nrows=nrows)
            if printvideo :
                pdfobj.colormap = colormap
                pdfobj.plotPDF(video,ii,'video',ncols=ncols,nrows=nrows)                    
            if printmarg :
                pdfobj.calcMarginal()
                pdfobj.plotPDFmarg(names2colors[pdftype],names2labels[pdftype],axx,axy)                    
            if pdfobj.PDFtype == 'dns':
                fmean = pdfobj.fmean ; fvar = pdfobj.fvar; covf = pdfobj.covf; zskew = pdfobj.zskew
            print pdftype
            #print pdfobj.params #, pdfobj.pdft, pdfobj.metrics['l1norm'], pdfobj.metrics['l2norm']
            print 'params', pdfobj.params
            #print pdfobj.pdfvals[:,0]
            #print pdfobj.pdfvals[0,:]

        # Make joint PDF contour plot with video frame from ensight
        if printvideo :
            if print1stvid :
                itime=0
                time='0.00000E+00'
            plt.figure('video')
            indextime=int(round(100*float(time)))
            vid.makevidframe(indextime,'video', 'bry',0,ncols=ncols,nrows=nrows,varx=varx,vary=vary)
            plt.figtext(0.125,0.885,'$t/\\tau_{eddy }= % 5.2f $' % (float(time)/0.2),size=14)
            plt.subplots_adjust(bottom=0.15, top=0.87, left=-0.01, right=0.91)
            plt.savefig('video/' + suffix + '_' + time + '.eps',format="eps")
            #vid.makevidframe(0, suffix, 'bry',0,ncols=ncols,nrows=nrows,varx=varx,vary=vary)
            #plt.savefig('joint/' + suffix + '_' + time + 'c.pdf',format="pdf")
            plt.clf()
        if printjoint :
            plt.figure('joint')
            time=int(round(100*float(time)))
            plt.subplots_adjust(bottom=0, top=0.76, left=0.2, right=1, wspace=0.08)
            vid.makevidframe(itime,'paper', 'bry',0,ncols=ncols,nrows=nrows,varx=varx,vary=vary)
            #plt.figtext(0.45,0.8, (suffix2suffix[suffix] + ': ' + '$t/\\tau_{eddy }= % 5.2f $' % (float(time)/0.2)) )
            plt.figtext(0.2,0.8, (suffix2suffix[suffix] + ': ' + '$t/\\tau_{eddy }= % 5.2f $' % (float(time)/0.2)) )
            plt.savefig('joint/' + suffix + '_' + time + '.eps',format="eps")
            plt.clf()
        if printorient :
            plt.figure('orient')
            plt.subplots_adjust(bottom=0, top=0.76, left=0.05, right=1, wspace=0.08)
            plt.figtext(0.5,0.8,'Components 1 and 2 Premixed',horizontalalignment='center')
            plt.savefig('orient/13.eps',format="eps")
            plt.clf()
            
        # Marginal Distribution plots
        if printmarg :
            plt.figure(marginalx.number)
            axx[0].set_ylabel('$P(Z_1)$')
            for kk in range(3) :
                plt.sca(axx[kk])
                plt.xlabel('$Z_1$')
                plt.ylim([0,10])
                plt.gca().set_yticks([0,2,4,6,8,10])
                #plt.gca().set_yticks([1,3,5,7,9],minor=True)
                plt.title( '$t/\\tau_{eddy }= % 5.2f $' % (float(time)/0.2) )
            plt.legend(frameon=False,fontsize=6)                
            plt.subplots_adjust(bottom=0.22, top=0.88, left=0.08, right=0.98)
            #plt.legend(frameon=False) #fontsize=14
            plt.savefig('marginal/x_' + suffix + '_' + time + '.eps',format="eps")
            plt.clf()
            
            plt.figure(marginaly.number)
            axy[0].set_ylabel('$P(Z_2)$')
            for kk in range(3) :
                plt.sca(axy[kk])
                plt.xlabel('$Z_2$')
                plt.ylim([0,10])
                plt.gca().set_yticks([0,2,4,6,8,10])
                #plt.gca().set_yticks([1,3,5,7],minor=True)
                plt.title( '$t/\\tau_{eddy }= % 5.2f $' % (float(time)/0.2) )
            plt.legend(frameon=False,fontsize=6)
            plt.subplots_adjust(bottom=0.22, top=0.88, left=0.08, right=0.98)            
            plt.savefig('marginal/y_' + suffix + '_' + time + '.eps',format="eps")
            plt.clf()

        # Print the metrics to file
        if printmetrics :
            for metric in metrics :
                line = ''.join('{:14.6g}'.format(k) for k in metrics[metric]) 
                line = '{:>14s}'.format(time) + line + '\n'
                filename = 'metrics/' + metric + '_' + suffix + '.txt'
                with open(filename,'a') as f:
                    f.write(line)

######################################################################################

















######################################################################################
                    
# Make header and footer for joint pdf figures
if printjoint :
    plt.figure('header',figsize=(figsizejoint[0],0.3))
    ii = 0
    for pdftype in pdftypes:
        ii = ii+1
        label = names2labels[pdftype]
        plt.subplot(1,6,ii)
        plt.gca().axis('off')
        plt.text(0.5,0.5,label,horizontalalignment='center',verticalalignment='center')
    plt.subplots_adjust(bottom=0.1, top=0.9, left=0.2, right=1, wspace=0.08)
    gs = grd.GridSpec(1,1)
    gs.update(left=0,right=0.15,bottom=0.1, top=0.9)
    plt.subplot(gs[0,0])
    plt.gca().axis('off')
    plt.text(0.5,0.5,'Scalar Field',horizontalalignment='center',verticalalignment='center')
    gs = grd.GridSpec(1,1)
    gs.update(left=0,right=1.0,bottom=0.95, top=1.0)
    plt.subplot(gs[0,0])
    plt.gca().axis('off')
    plt.plot([0,1],[0,0],'k')
    gs = grd.GridSpec(1,1)
    gs.update(left=0,right=1.0,bottom=0.0, top=0.05)
    plt.subplot(gs[0,0])
    plt.gca().axis('off')
    plt.plot([0,1],[0,0],'k')
    plt.savefig('joint/header.eps',format="eps")
    plt.clf()
    
    plt.figure('footer',figsize=(figsizejoint[0],0.82))
    for ii in range(6) :
        plt.subplot(1,6,ii+1)
        for spine in plt.gca().spines.values() : spine.set_visible(False)
        pad=-7
        plt.gca().annotate('$Z_1$', xy=(0.5, 0), xytext=(-pad,-plt.gca().xaxis.labelpad - pad ),
                    xycoords=plt.gca().xaxis.label, textcoords='offset points',
                    ha='right', va='center')
        plt.gca().yaxis.set_ticks([])
        plt.gca().xaxis.set_ticks([0,1])
        plt.gca().xaxis.set_ticks_position('none')
    plt.subplots_adjust(bottom=1.02, top=1.03, left=0.20, right=0.995, wspace=0.135)
    gs1 = grd.GridSpec(1,1)
    gs1.update(left=0,right=1.0,bottom=0.0, top=0.05)
    plt.subplot(gs1[0,0])
    plt.gca().axis('off')
    plt.plot([0,1],[0,0],'k')
    gs2 = grd.GridSpec(1,1)
    gs2.update(left=0.035,right=0.15,bottom=0.31, top=0.8)
    plt.subplot(gs2[0,0])
    vid.makelegend('bry',plt.gca())
    #plt.text( 1.2,-0.1,'$C_1$')
    #plt.text(-0.1, 1.12,'$C_2$')
    #plt.text(-0.55,-0.35,'$C_3$')

    axcb = plt.gcf().add_axes([0.45,0.45,0.3,0.12])
    minlog = np.floor(np.log10(pdfobj.minlog))
    maxlog = np.ceil(np.log10(pdfobj.maxlog))
    xx,yy = np.meshgrid([0,1],[0,1])
    zz = np.array([[1e-2,1e2],[1e-2,1e2]])
    grf1 = plt.contourf(xx,yy,zz,levels=pdfobj.levs, norm=LogNorm(vmin=pdfobj.minlog, vmax=pdfobj.maxlog),cmap=pdfobj.colormap)
    cbar = plt.colorbar(grf1,cax = axcb,cmap=pdfobj.colormap,
                        ticks=np.logspace(minlog,maxlog,maxlog-minlog+1), 
                        format=LogFormatter(10, labelOnlyBase=False),orientation='horizontal')
    cbar.set_label("\\vspace{-36pt} $P(Z_1,Z_2)$")
    plt.savefig('joint/footer.eps',format="eps")
    plt.clf()

if printorient :
    fig = plt.figure('headerO',figsize=(6.5,0.5))
    
    plt.subplot(1,9,1)
    plt.gca().axis('off')
    plt.text(0.5,0.5,'DNS',horizontalalignment='center',verticalalignment='center')
    
    plt.subplot(1,9,2)
    plt.gca().axis('off')
    plt.text(0.5,0.25,'-1',horizontalalignment='center',verticalalignment='center')
    plt.subplot(1,9,3)
    plt.gca().axis('off')
    plt.text(0.5,0.75,'CM',horizontalalignment='center',verticalalignment='center')
    plt.text(0.5,0.25,'-2',horizontalalignment='center',verticalalignment='center')
    plt.subplot(1,9,4)
    plt.gca().axis('off')
    plt.text(0.5,0.25,'-3',horizontalalignment='center',verticalalignment='center')
    
    plt.subplot(1,9,5)
    plt.gca().axis('off')
    plt.text(0.5,0.25,'-12',horizontalalignment='center',verticalalignment='center')
    plt.subplot(1,9,6)
    plt.gca().axis('off')
    plt.text(0.5,0.75,'BVB5',horizontalalignment='center',verticalalignment='center')
    plt.text(0.5,0.25,'-13',horizontalalignment='center',verticalalignment='center')
    plt.subplot(1,9,7)
    plt.gca().axis('off')
    plt.text(0.5,0.25,'-23',horizontalalignment='center',verticalalignment='center')
    
    plt.subplot(1,9,8)
    plt.gca().axis('off')
    plt.text(0.5,0.5,'BVB6',horizontalalignment='center',verticalalignment='center')
    plt.subplot(1,9,9)
    plt.gca().axis('off')
    plt.text(0.5,0.5,'SMLD',horizontalalignment='center',verticalalignment='center')
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=1, wspace=0.08)

    gs = grd.GridSpec(1,9)
    fig.add_subplot(gs[0,1:4])
    plt.gca().axis('off')
    plt.plot([0,1],[0,0],'k')
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=1, wspace=0.12)
    fig.add_subplot(gs[4:7])
    plt.gca().axis('off')
    plt.plot([0,1],[0,0],'k')
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=1, wspace=0.12)
    
    gs = grd.GridSpec(1,1)
    gs.update(left=0,right=1.0,bottom=0.95, top=1.0)
    plt.subplot(gs[0,0])
    plt.gca().axis('off')
    plt.plot([0,1],[0,0],'k')
    gs = grd.GridSpec(1,1)
    gs.update(left=0,right=1.0,bottom=0.0, top=0.05)
    plt.subplot(gs[0,0])
    plt.gca().axis('off')
    plt.plot([0,1],[0,0],'k')
    
    plt.savefig('orient/headerO.eps',format="eps")
    plt.clf()
    
    plt.figure('footer',figsize=(6.5,0.82))
    for ii in range(9) :
        plt.subplot(1,9,ii+1)
        for spine in plt.gca().spines.values() : spine.set_visible(False)
        pad=-7
        plt.gca().annotate('$Z_1$', xy=(0.5, 0), xytext=(-pad,-plt.gca().xaxis.labelpad - pad ),
                    xycoords=plt.gca().xaxis.label, textcoords='offset points',
                    ha='right', va='center')
        plt.gca().yaxis.set_ticks([])
        plt.gca().xaxis.set_ticks([0,1])
        plt.gca().xaxis.set_ticks_position('none')
    plt.subplots_adjust(bottom=1.02, top=1.03, left=0.05, right=0.995, wspace=0.135)
    gs1 = grd.GridSpec(1,1)
    gs1.update(left=0,right=1.0,bottom=0.0, top=0.05)
    plt.subplot(gs1[0,0])
    plt.gca().axis('off')
    plt.plot([0,1],[0,0],'k')
    axcb = plt.gcf().add_axes([0.35,0.45,0.3,0.12])
    minlog = np.floor(np.log10(pdfobj.minlog))
    maxlog = np.ceil(np.log10(pdfobj.maxlog))
    xx,yy = np.meshgrid([0,1],[0,1])
    zz = np.array([[1e-2,1e2],[1e-2,1e2]])
    grf1 = plt.contourf(xx,yy,zz,levels=pdfobj.levs, norm=LogNorm(vmin=pdfobj.minlog, vmax=pdfobj.maxlog),cmap=pdfobj.colormap)
    cbar = plt.colorbar(grf1,cax = axcb,cmap=pdfobj.colormap,
                        ticks=np.logspace(minlog,maxlog,maxlog-minlog+1), 
                        format=LogFormatter(10, labelOnlyBase=False),orientation='horizontal')
    cbar.set_label("\\vspace{-36pt} $P(Z_1,Z_2)$")
    plt.savefig('orient/footerO.eps',format="eps")
    plt.clf()
