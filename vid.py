#!/usr/bin/python2.7

import numpy                as np
import matplotlib.pyplot    as plt
import matplotlib.patches   as patches
import matplotlib.gridspec  as grd 

# make triangle legend
def makelegend(colorscheme,ax) :

    # create color scheme
    if colorscheme == 'wok' : # white orange black
        colorA=(0.0,0.0,0.0)
        colorB=(1.0,143.0/255.0,0.0)
        colorC=(1.0,1.0,1.0)
    elif colorscheme == 'special' : # mix to Princeton
        colorA=(1.0,1.0,0.0)
        colorB=(1.0,0.1686,0.0)
        colorC=(225.0/245.0,0.0,3.0*37.0/255.0)
    else :# Blue Red Yellow
        colorA=(0.0,0.0,1.0)
        colorB=(1.0,0.0,0.0)
        colorC=(1.0,1.0,0.0)

    #fig, ax = plt.subplots()
    ax.set_aspect('equal')
    #ax.set_xlabel('$Z_1$')
    #ax.set_ylabel('$Z_2$')
    pad=-10
    plt.gca().annotate('$Z_1$', xy=(0.5, 0), xytext=(5,-plt.gca().xaxis.labelpad - pad  ),
                       xycoords=plt.gca().xaxis.label, textcoords='offset points',
                       ha='right', va='center')
    pad=-12.5
    plt.gca().annotate('$Z_2$', xy=(0, 0.5), xytext=(-plt.gca().yaxis.labelpad - pad, 0),
                       xycoords=plt.gca().yaxis.label, textcoords='offset points',
                       ha='right', va='center')
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(top='off',right='off')
    n=50
    hstep = 1.0/n
    vstep = 1.0/n
    for i in range(n) :
        for j in range(n-i) :
            x_vert = [i*hstep, (i    )*hstep, (i+1)*hstep]
            y_vert = [j*vstep, (j+1.0)*vstep,  j   *vstep]
            ZA = np.mean(x_vert)
            ZB = np.mean(y_vert)
            ZC = 1.0 - ZA - ZB
            color  = [max(min(ZA*colorA[0]+ZB*colorB[0]+ZC*colorC[0],1),0),
                      max(min(ZA*colorA[1]+ZB*colorB[1]+ZC*colorC[1],1),0),
                      max(min(ZA*colorA[2]+ZB*colorB[2]+ZC*colorC[2],1),0),
                      1                                                     ]
            # x_vert += j*0.5*hstep*np.ones(3)
            ax.fill(x_vert,y_vert,facecolor=color,edgecolor=color)
        for j in range(n-i-1) :
            x_vert = [(i    )*hstep, (i+1.0)*hstep, (i+1.0)*hstep]
            y_vert = [(j+1.0)*vstep, (j    )*vstep, (j+1.0)*vstep]
            ZA = np.mean(x_vert)
            ZB = np.mean(y_vert)
            ZC = 1.0 - ZA - ZB
            color  = [max(min(ZA*colorA[0]+ZB*colorB[0]+ZC*colorC[0],1),0),
                      max(min(ZA*colorA[1]+ZB*colorB[1]+ZC*colorC[1],1),0),
                      max(min(ZA*colorA[2]+ZB*colorB[2]+ZC*colorC[2],1),0),
                      1                                                     ]
            # x_vert[:] += j*0.5*hstep*np.ones(3)
            ax.fill(x_vert,y_vert,facecolor=color,edgecolor=color)
    plt.plot([0,1],[1,0],'k')
    





# make time varying plots from ensight data
def makevidframe(idat, condition, colorscheme, ii, nrows, ncols, varx, vary, fig='') :
    #simulation parameters (nneded because data is stored in domains by processor)
    nx=256
    nprocx = 8
    nprocy = 8
    nprocz = 8
    nx_ = nx/nprocx
    ny_ = nx/nprocy
    nz_ = nx/nprocz

    # extract data from ensight files
    dataA = np.zeros((nx,nx))
    dataB = np.zeros((nx,nx))
    data = {varx : dataA, vary : dataB}
    datnumber = str(idat+1).zfill(6)
    for ifi in [varx, vary] :
        with open("./data/" + ifi + "." + datnumber,"rb") as f :
            f.read(80) # var name
            f.read(80) # part
            np.fromfile(f,dtype=np.int32,count=1) # int
            f.read(80) #ndarray
            for jp in range(nprocy) :
                for ip in range(nprocx) :
                    for i in range(nx_) :
                        data[ifi][jp*ny_:(jp+1)*ny_,ip*nx_+i] = np.fromfile(f,dtype=np.float32,count=ny_)[:] 
                    f.seek(4*nx_*ny_*(nz_-1)           ,1) # seek to end of processor data
                    f.seek(4*nx_*ny_* nz_   *(nprocz-1),1) # seek to next relevant processor
        f.close()


    # create color scheme
    if colorscheme == 'wok' : # white orange black
        colorA=(0.0,0.0,0.0)
        colorB=(1.0,143.0/255.0,0.0)
        colorC=(1.0,1.0,1.0)
    elif colorscheme == 'special' : # mix to Princeton
        colorA=(1.0,1.0,0.0)
        colorB=(1.0,0.1686,0.0)
        colorC=(225.0/245.0,0.0,3.0*37.0/255.0)
    else :# Blue Red Yellow
        colorA=(0.0,0.0,1.0)
        colorB=(1.0,0.0,0.0)
        colorC=(1.0,1.0,0.0)

    # convert data from mixture fractions into RGB colors
    color=np.zeros((nx,nx,4))
    xvec=range(nx)
    xx,yy = np.meshgrid(xvec,xvec)    
    step   = 1.0/255
    for i in range(nx) :
        for j in range(nx) :
            x_min  = float(i)*step
            y_min  = float(j)*step
            x_vert = [x_min, x_min     , x_min+step, x_min+step]
            y_vert = [y_min, y_min+step, y_min+step, y_min     ]
            ZA     = max(min(data[varx][j,i],1),0)
            ZB     = max(min(data[vary][j,i],1),0)
            ZC     = max(min(1 - ZA - ZB   ,1),0)
            color[j,i,:]  = [max(min(ZA*colorA[0]+ZB*colorB[0]+ZC*colorC[0],1),0),
                             max(min(ZA*colorA[1]+ZB*colorB[1]+ZC*colorC[1],1),0),
                             max(min(ZA*colorA[2]+ZB*colorB[2]+ZC*colorC[2],1),0),
                             1                                                     ]
    # make plot
    if condition == 'icspaper' :
        plt.subplot(nrows,ncols,ii)
    elif condition == 'video' :
        plt.subplot2grid((nrows,ncols),(0,0),rowspan=2,colspan=2)
    elif condition == 'paper' :
        gs = grd.GridSpec(1,1)
        gs.update(left=0,right=0.15,bottom=0.0, top=0.87)
        plt.subplot(gs[0,0]) #171
    image = plt.imshow(color,aspect='equal',extent=(0,0.5,0.5,1)) #.set_ticklabels([])
    image.axes.get_xaxis().set_ticklabels([])
    image.axes.get_yaxis().set_ticklabels([])

