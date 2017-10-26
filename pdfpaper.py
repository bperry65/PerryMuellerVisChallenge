import numpy as np
import math as math
import scipy.stats as stat
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import LogFormatter
from scipy.optimize import root
from scipy.special import beta
from scipy.special import betaln
from scipy.integrate import dblquad

class dnsdata(object) :
     def __init__(self,data) :
         self.data = data
     def redefdata(self,data) :
         self.data = data

########################################################
#                   PDF Object
########################################################

class pdf(object):

    def __init__(self, PDFtype, xmean=0.5, xvar=0, ymean=0.25, yvar=0, covar=0, nx=0, ny=0, DNSinit='false',filepath='',time='0.0',varx='ZMIXA',vary='ZMIXB',filt=1,switch='F'):
        # inputs
        self.PDFtype = PDFtype
        self.PDFlabel = ''
        self.filepath = filepath
        self.timestr = time
        self.time = float(time)
        self.varx = varx
        self.vary = vary
        self.filt = filt
        self.switch = switch
        self.error = 0
        if self.PDFtype == 'dns':
            self.initDNSparams(self.filepath,self.time,self.varx,self.vary)
            if (self.nx == 0) and (self.ny == 0) :
                print "__init__ : warning : nx and ny not found in input_box file, using values from function input"
                self.nx = nx
                self.ny = ny
        elif DNSinit == 'true' :
            self.initDNSparams(self.filepath,self.time,self.varx,self.vary)
            if (nx != 0) or (self.nx == 0):
                self.nx = nx
            if (ny != 0) or (self.ny == 0):
                self.ny = ny
        else:
            self.xmean   = xmean
            self.xvar    = xvar
            self.ymean   = ymean
            self.yvar    = yvar
            self.covar   = covar
            self.nx      = nx
            self.ny      = ny

        # other values
        if self.filt == 0 :
            self.dx      = 1.0/(self.nx)
            self.dy      = 1.0/(self.ny)
        else :
            self.dx      = 1.0/(self.nx - 1.0)
            self.dy      = 1.0/(self.ny - 1.0) 

        # data storage
        self.pdfvals = np.zeros((self.ny,self.nx))
        self.pdfarea = np.ones((self.ny,self.nx))
        self.pdfarea = self.pdfarea*self.dx*self.dy
        if self.filt == 0 :
            self.xvec    = np.linspace(0.0+0.5/self.nx,1.0-0.5/self.nx,self.nx)
            self.yvec    = np.linspace(0.0+0.5/self.ny,1.0-0.5/self.ny,self.ny)
        else :
            self.xvec    = np.linspace(0.0,1.0,self.nx)
            self.yvec    = np.linspace(0.0,1.0,self.ny)
            self.pdfarea[0,0] = self.pdfarea[0,0] / 4.0
            self.pdfarea[1:self.nx-1,0] = self.pdfarea[1:self.nx-1,0] / 2.0
            self.pdfarea[0,1:self.nx-1] = self.pdfarea[0,1:self.nx-1] / 2.0
            for i in range(1,self.nx-1) :
                self.pdfarea[i,self.nx-i-1] = self.pdfarea[i,self.nx-i-1] / 2.0
            self.pdfarea[0,self.nx-1] = self.pdfarea[0,self.nx-1] / 8.0
            self.pdfarea[self.nx-1,0] = self.pdfarea[self.nx-1,0] / 8.0
        self.xmarg   = np.zeros(self.nx)
        self.ymarg   = np.zeros(self.ny)
        self.xmargB  = np.zeros(self.nx)
        self.ymargB  = np.zeros(self.ny) 
        self.params  = np.zeros(7)

        # error checking
        if self.PDFtype not in ['dirichlet','bvb5','betadelta','sml','sml3','sml4','smlc','bvb4','bvb6','dns'] :
            raise RuntimeError("PDFtype must be one of: dirichlet, betadelta, bvb[45], sml[34c], dns \
                                \n %s was specified" % self.PDFtype)
        if (self.xmean > 1 or self.xmean < 0):
            raise RuntimeError("xmean must be between 0 and 1")
        if (self.ymean > 1 or self.ymean < 0):
            raise RuntimeError("ymean must be between 0 and 1")  

    ###################################################################

    # Initialize means and variances from values in the monitor/hit_scalar file
    def initDNSparams(self,filepath,time,variable1,variable2):
        # open file and get data
        hit_scalar = np.genfromtxt(filepath + "/hit_scalar",skip_header=1)
        openfile = open(filepath + "/hit_scalar",'r')
        colnames = openfile.readline().strip().split()
        openfile.close()
        timecol = 1
        
        # find data corresponding the requested time
        if hit_scalar[0,timecol] > time :
            raise RuntimeError("initDNSparams : time requested is before start of given file")
        if hit_scalar[0,timecol] == time :
            ilow = 0
            ihigh = 0
            alpha = 1
        for i in range(1,len(hit_scalar)) :
            if hit_scalar[i,timecol] >= time :
                ilow = i-1
                ihigh = i
                alpha = (time - hit_scalar[i-1,timecol])/(hit_scalar[i,timecol] - hit_scalar[i-1,timecol])
                break
        else :  # attached to for loop, occurs when loop reaches end without break
            raise RuntimeError("initDNSparams : time requested is after end of given file")    
        
        # get variables
        found = 0
        for i in range(len(colnames)) :
            if colnames[i] == 'avg_' + variable1 :
                found = found+1
                self.xmean = hit_scalar[ilow,i] + alpha * (hit_scalar[ihigh,i] - hit_scalar[ilow,i])
            if colnames[i] == 'avg_' + variable2 :
                found = found+10
                self.ymean = hit_scalar[ilow,i] + alpha * (hit_scalar[ihigh,i] - hit_scalar[ilow,i])
            if colnames[i] == 'var_' + variable1 :
                found = found+100
                self.xvar = hit_scalar[ilow,i] + alpha * (hit_scalar[ihigh,i] - hit_scalar[ilow,i])
            if colnames[i] == 'var_' + variable2 :
                found = found+1000
                self.yvar = hit_scalar[ilow,i] + alpha * (hit_scalar[ihigh,i] - hit_scalar[ilow,i])
            if colnames[i] == 'cov_' + variable1[4:].strip() + '_' + variable2[4:].strip()  or \
                     colnames[i] == 'cov_' + variable2[4:].strip() + '_' + variable1[4:].strip()  :
                found = found+10000
                self.covar = hit_scalar[ilow,i] + alpha * (hit_scalar[ihigh,i] - hit_scalar[ilow,i])
            if colnames[i] == 'skw_' + variable1 :
                found = found+100000
                self.xskew = hit_scalar[ilow,i] + alpha * (hit_scalar[ihigh,i] - hit_scalar[ilow,i])
            if colnames[i] == 'skw_' + variable2 :
                found = found+1000000
                self.yskew = hit_scalar[ilow,i] + alpha * (hit_scalar[ihigh,i] - hit_scalar[ilow,i])
        if found != 1111111 :
            raise RuntimeError("initDNSparams : could not find statistics for requested variables")

        # get nx and ny
        for line in open(filepath + "/input_box") :
            if 'Number of bins' in line :
                self.nx = self.ny = self.nbins = int(line.split(':')[-1])
                if self.filt != 0 :
                    self.nx = self.nx // self.filt + 1
                    self.ny = self.ny // self.filt + 1
                elif (self.nbins % self.filt) != 0 :
                    raise RuntimeError("initDNSparams : number of points not divisible by filter")
                break
        else :
            self.nx = self.ny = 0

        # if self.switch == 'x' :
        zmean = 1 - self.xmean - self.ymean
        zvar = self.xvar + self.yvar + 2 * self.covar
        covaryz = 0.5*(self.xvar - self.yvar - zvar)
        self.xmean = zmean
        self.xvar = zvar
        self.covar = covaryz
        # elif self.switch == 'y' :
        # zmean = 1 - self.xmean - self.ymean
        # zvar = self.xvar + self.yvar + 2 * self.covar
        # covarxz = 0.5*(self.yvar - self.xvar - zvar)
        # self.ymean = zmean
        # self.yvar = zvar
        # self.covar = covarxz
        # elif self.switch == 'z' :
        # tmp1 = self.xmean
        # tmp2 = self.xvar
        # self.xmean = self.ymean
        # self.xvar  = self.yvar
        # self.ymean = tmp1
        # self.yvar  = tmp2
        if self.switch == '312' :
            zmean = 1 - self.xmean - self.ymean
            zvar = self.xvar + self.yvar + 2 * self.covar
            covarxz = 0.5*(self.yvar - self.xvar - zvar)
            self.ymean = self.xmean
            self.yvar = self.xvar
            self.xmean = zmean
            self.xvar = zvar
            self.covar = covarxz
        elif self.switch == '231' :
            zmean = 1 - self.xmean - self.ymean
            zvar = self.xvar + self.yvar + 2 * self.covar
            covaryz = 0.5*(self.xvar - self.yvar - zvar)
            self.xmean = self.ymean
            self.xvar = self.yvar
            self.ymean = zmean
            self.yvar = zvar
            self.covar = covaryz

    ###################################################################

    # Calculate value of PDF at a given location in x,y space for given input parameters
    def calcPDFval(self,x,y,parameters):
        if self.PDFtype == 'dirichlet':
            a0 = parameters[0]; a1 = parameters[1]; a2 = parameters[2]
            logC = parameters[5]
            #value = pre * ( x**(a0-1) * y**(a1-1) * (1-x-y)**(a2-1) )  
            value = math.exp( min( logC + (a0-1)*math.log(x)   + (a1-1)*math.log(y) + (a2-1)*math.log(1-x-y),500) )
            return value
        if self.PDFtype == 'bvb5':
            a0 = parameters[0]; a1 = parameters[1]; a2 = parameters[2];
            a3 = parameters[3]; a4 = parameters[4]; logC = parameters[5]
            #value = C * ( x**(a0-1) * y**(a1-1) * (1-x-y)**(a2-1) * (1-x)**(a3-1) * (1-y)**(a4-1) )
            value = math.exp( min(logC + (a0-1)*math.log(x  ) + (a1-1)*math.log(y  ) + (a2-1)*math.log(1-x-y) \
                                  + (a3-1)*math.log(1-x) + (a4-1)*math.log(1-y) , 500 ))
            return value
        if self.PDFtype == 'bvb6':
            a0 = parameters[0]; a1 = parameters[1]; a2 = parameters[2];
            a3 = parameters[3]; a4 = parameters[4]; logC = parameters[5]; a5 = parameters[6]
            #value = C * ( x**(a0-1) * y**(a1-1) * (1-x-y)**(a2-1) * (1-x)**(a3-1) * (1-y)**(a4-1) )
            value = math.exp( min(logC + (a0-1)*math.log(x  ) + (a1-1)*math.log(y  ) + (a2-1)*math.log(1-x-y) \
                                   + (a3-1)*math.log(1-x) + (a4-1)*math.log(1-y) + (a5-1)*math.log(x+y), 500 ))
            return value
        if self.PDFtype == 'bvb4':
            a0 = parameters[0]; a1 = parameters[1]; a2 = parameters[2];
            a3 = parameters[3]; logC = parameters[5]
            #value = C * ( x**(a0-1) * y**(a1-1) * (1-x-y)**(a2-1) * (1-x)**(a3-1) )
            value = math.exp( min(logC + (a0-1)*math.log(x  ) + (a1-1)*math.log(y  ) + (a2-1)*math.log(1-x-y) \
                                   + (a3-1)*math.log(1-x), 500) )
            return value
        if self.PDFtype == 'sml':
            logC = parameters[0]; a1 = parameters[1]; a2 = parameters[2];
            a3 = parameters[3]; a4 = parameters[4]; a5 = parameters[5]
            value = math.exp(min(logC + a1*x + a2*y + a3*x**2 + a4*y**2 + a5*x*y,500))
            return value
        if self.PDFtype == 'sml4':
            logC = parameters[0]; a1 = parameters[1]; a2 = parameters[2];
            a3 = parameters[3]; a4 = parameters[4]
            value = math.exp(logC + a1*x + a2*y + a3*x**2 + a4*y**2)
            return value
        if self.PDFtype == 'sml3':
            logC = parameters[0]; a1 = parameters[1]; a2 = parameters[2];
            a3 = parameters[3]
            value = math.exp(logC + a1*x + a2*y + a3*x**2)
            return value
        if self.PDFtype == 'smlc':
            logC = parameters[0]; a1 = parameters[1]; a2 = parameters[2];
            a3 = parameters[3]; a4 = parameters[4]
            value = math.exp(logC + a1*x + a2*y + a3*x**2 + a4*x*y)
            return value
        raise RuntimeError("calcPDFval not defined for PDF type %s" % self.PDFtype)

    ##################################################################

    # Calculate PDF values on boundary by explicit integration of approximate functions
    def calcPDFvals_boundary(self, parameters):
        a0 = parameters[0]; a1 = parameters[1]; a2 = parameters[2];
        a3 = parameters[3]; a4 = parameters[4]; logC = parameters[5]
        if self.PDFtype in ['dirichlet','bvb4','bvb5'] :
            if self.PDFtype not in self.exact_integ : #remove bvb5
                 self.pdfvals[0,0] = math.exp( min(logC - math.log(abs(a0*a1)) + (a0+a1)*math.log(self.dx/2.0),600) ) / self.pdfarea[0,0]
                 self.pdfvals[self.ny-1,0] = math.exp( min( logC + (a0+a2+a4-1.0)*math.log(self.dx/2.0) + betaln(a2,a0+1.0) \
                                                            - math.log(abs(a0*a4)) ,600)) / self.pdfarea[self.ny-1,0]
                 self.pdfvals[0,self.nx-1] = math.exp( min(logC + (a1+a2+a3-1.0)*math.log(self.dx/2.0) + betaln(a2,a1+1.0) \
                                                           - math.log(abs(a1*a3)) ,600)) / self.pdfarea[0,self.nx-1]
            else :
                 self.pdfvals[0,0]        , err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                          0,self.dx/2  , lambda y: 0,lambda y:self.dx/2)
                 self.pdfvals[self.ny-1,0], err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                          1-self.dx/2,1, lambda y: 0,lambda y:1-y      )
                 self.pdfvals[0,self.nx-1], err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                          0,self.dx/2  , lambda y: 1-self.dx/2,lambda y:1-y)
                 self.pdfvals[0,0]         = self.pdfvals[0,0]         / self.pdfarea[0,0] 
                 self.pdfvals[self.ny-1,0] = self.pdfvals[self.ny-1,0] / self.pdfarea[self.ny-1,0]
                 self.pdfvals[0,self.nx-1] = self.pdfvals[0,self.nx-1] / self.pdfarea[0,self.nx-1]
            for i in range(1,self.nx-1) :
                self.pdfvals[0,i] = math.exp(min(logC + (a1+1.0)*math.log(self.dx) - math.log(abs(a1)*2.0**a1) \
                    + (a0-1.0)*math.log(self.xvec[i]) + (a2+a3-2.0)*math.log(1.0-self.xvec[i]),600) ) / self.pdfarea[0,i] 
                self.pdfvals[i,0] = math.exp(min(logC + (a0+1.0)*math.log(self.dx) - math.log(abs(a0)*2.0**a0) \
                    + (a1-1.0)*math.log(self.yvec[i]) + (a2+a4-2.0)*math.log(1.0-self.yvec[i]),600) ) / self.pdfarea[i,0]
                self.pdfvals[i,self.nx-i-1] = math.exp( min(logC + (a0-1.0)*math.log(self.xvec[self.nx-i-1]) \
                    + (a1-1.0)*math.log(self.yvec[i]) + (a3-1.0)*math.log(1.0-self.xvec[self.nx-i-1]) \
                    + (a4-1.0)*math.log(1.0-self.yvec[i]) + (a2+1.0)*math.log(self.dx) - math.log(abs(a2*a2+a2)),200) ) \
                    / self.pdfarea[i,self.nx-i-1] 
            # constrain values on diagonal near corners where approx can go horribly wrong
            #if self.PDFtype in ['dirichlet','bvb4'] :
            if self.PDFtype not in self.exact_integ :
                 maxcor = 20
                 for i in range(1,maxcor) :
                      big =    max(self.pdfvals[0,self.nx-1], self.pdfvals[maxcor,self.nx-maxcor-1])
                      little = min(self.pdfvals[0,self.nx-1], self.pdfvals[maxcor,self.nx-maxcor-1])
                      self.pdfvals[i,self.nx-i-1] = max(min(self.pdfvals[i,self.nx-i-1],big),little)
                      big =    max(self.pdfvals[self.nx-1,0], self.pdfvals[self.nx-maxcor-1,maxcor])
                      little = min(self.pdfvals[self.nx-1,0], self.pdfvals[self.nx-maxcor-1,maxcor])
                      self.pdfvals[self.nx-i-1,i] = max(min(self.pdfvals[self.nx-i-1,i],big),little)
                      big =    max(self.pdfvals[0,        0], self.pdfvals[maxcor,               0])
                      little = min(self.pdfvals[0,        0], self.pdfvals[maxcor,               0])
                      self.pdfvals[        i,  0] = max(min(self.pdfvals[        i,  0],big),little)
                      big =    max(self.pdfvals[self.nx-1,0], self.pdfvals[self.nx-maxcor-1,     0])
                      little = min(self.pdfvals[self.nx-1,0], self.pdfvals[self.nx-maxcor-1,     0])
                      self.pdfvals[self.nx-i-1,0] = max(min(self.pdfvals[self.nx-i-1,0],big),little)
                      big =    max(self.pdfvals[0,        0], self.pdfvals[0,               maxcor])
                      little = min(self.pdfvals[0,        0], self.pdfvals[0,               maxcor])
                      self.pdfvals[0,          i] = max(min(self.pdfvals[0,          i],big),little)
                      big =    max(self.pdfvals[0,self.nx-1], self.pdfvals[0,     self.nx-maxcor-1])
                      little = min(self.pdfvals[0,self.nx-1], self.pdfvals[0,     self.nx-maxcor-1])
                      self.pdfvals[0,self.nx-i-1] = max(min(self.pdfvals[0,self.nx-i-1],big),little)
            else :
                 dx05 = self.dx/2
                 error = 100*np.ones((6,1))
                 count =     np.ones((6,1))
                 while error[0] >= 0.1 and count[0] <= self.nx/2-1:
                      i = int(count[0])
                      count[0] = count[0] + 1
                      pdfval_0i, err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                       0,dx05  , lambda y: self.xvec[i]-dx05, lambda y: self.xvec[i]+dx05)
                      pdfval_0i = pdfval_0i / self.pdfarea[0,i]
                      error[0] = abs(self.pdfvals[0,i]-pdfval_0i)/(pdfval_0i + 1e-9)
                      self.pdfvals[0,i] = pdfval_0i
                 while error[1] >= 0.1 and count[1] <= self.nx/2-1:
                      i = int(self.nx - count[1] - 1)
                      count[1] = count[1] + 1
                      pdfval_0i, err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                       0,dx05  , lambda y: self.xvec[i]-dx05, lambda y: self.xvec[i]+dx05)
                      pdfval_0i = pdfval_0i / self.pdfarea[0,i]
                      error[1] = abs(self.pdfvals[0,i]-pdfval_0i)/(pdfval_0i + 1e-9)
                      self.pdfvals[0,i] = pdfval_0i
                 while error[2] >= 0.1 and count[2] <= self.nx/2-1:
                      i = int(count[2])
                      count[2] = count[2] + 1                 
                      pdfval_i0, err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                       self.yvec[i]-dx05,self.yvec[i]+dx05, lambda y: 0, lambda y:dx05)
                      pdfval_i0 = pdfval_i0 / self.pdfarea[i,0]
                      error[2] = abs(self.pdfvals[i,0]-pdfval_i0)/(pdfval_i0 + 1e-9)
                      self.pdfvals[i,0] = pdfval_i0
                 while error[3] >= 0.1 and count[3] <= self.nx/2-1:
                      i = int(self.nx - count[3] - 1)
                      count[3] = count[3] + 1                 
                      pdfval_i0, err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                       self.yvec[i]-dx05,self.yvec[i]+dx05, lambda y: 0, lambda y:dx05)
                      pdfval_i0 = pdfval_i0 / self.pdfarea[i,0]
                      error[3] = abs(self.pdfvals[i,0]-pdfval_i0)/(pdfval_i0 + 1e-9)
                      self.pdfvals[i,0] = pdfval_i0
                 while error[4] >= 0.1 and count[4] <= self.nx/2-1:
                      i = int(count[4])
                      j = int(self.nx -i -1)
                      count[4] = count[4] + 1                     
                      pdfval_ji, err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                       self.yvec[j]-dx05,self.yvec[j]+dx05  , lambda y: self.xvec[i]-dx05 ,lambda y:1-y)
                      pdfval_ji = pdfval_ji / self.pdfarea[j,i]
                      error[4] = abs(self.pdfvals[j,i]-pdfval_ji)/(pdfval_ji + 1e-9)
                      self.pdfvals[j,i] = pdfval_ji
                 while error[5] >= 0.1 and count[5] <= self.nx/2-1:
                      j = int(count[5])
                      i = int(self.nx -j -1)
                      count[5] = count[5] + 1                     
                      pdfval_ji, err = dblquad(lambda x,y: self.calcPDFval(x,y,parameters), \
                                                       self.yvec[j]-dx05,self.yvec[j]+dx05  , lambda y: self.xvec[i]-dx05 ,lambda y:1-y)
                      pdfval_ji = pdfval_ji / self.pdfarea[j,i]
                      error[5] = abs(self.pdfvals[j,i]-pdfval_ji)/(pdfval_ji + 1e-9)
                      self.pdfvals[j,i] = pdfval_ji
            #print count.transpose()
            #print error.transpose()
            #print self.xvec[10], self.yvec[self.nx-10-1]
            #print self.pdfvals[0,10]*self.pdfarea[0,10]/math.exp(logC), self.pdfvals[self.nx-10-1,0]*self.pdfarea[self.nx-10-1,0]/math.exp(logC), self.pdfvals[self.nx-10-1,10]*self.pdfarea[self.nx-10-1,10]/math.exp(logC)
        if 'sml' in self.PDFtype :
            self.pdfvals[0,0]         = self.calcPDFval(            0,             0,parameters)
            self.pdfvals[self.ny-1,0] = self.calcPDFval(            0,             1,parameters)
            self.pdfvals[0,self.nx-1] = self.calcPDFval(            1,             0,parameters)
            for i in range(1,self.nx-1):
                self.pdfvals[0,i]     = self.calcPDFval(self.xvec[i],            0,parameters)
                self.pdfvals[i,0]     = self.calcPDFval(            0,self.yvec[i],parameters)
                self.pdfvals[i,self.nx-i-1]     = self.calcPDFval(self.xvec[self.nx-i-1],self.yvec[i],parameters)
            
    ##################################################################

    # Calculate parameters for desired distribution that give the specified moments
    def calcPDFparams(self, guess=np.zeros(7)):
        if self.PDFtype == 'betadelta':
            A = self.xmean*(1-self.xmean)/self.xvar - 1
            self.params[0] = self.xmean * A
            self.params[1] = (1-self.xmean) * A
            self.params[2] = self.ymean/(1-self.xmean)
            self.yvar  = self.xvar * self.params[2]*self.params[2]
            self.covar = -self.xvar * self.params[2]
            self.xskew = 2*(self.params[1]-self.params[0])*(self.params[1]+self.params[0]+1)**0.5 \
                / ((self.params[1]+self.params[0]+2)*(self.params[1]*self.params[0])**0.5)
            self.yskew = - self.xskew
            zvar = 2.0*self.covar + self.xvar + self.yvar
            self.covxz = 0.5*(self.yvar - self.xvar - zvar)
            self.covyz = 0.5*(self.xvar - self.yvar - zvar)
        if self.PDFtype == 'dirichlet':
            A = self.xmean*(1-self.xmean)/self.xvar - 1
            self.params[0] = self.xmean * A
            self.params[1] = self.ymean * A
            self.params[2] = (1-self.xmean-self.ymean)*A
            self.params[3] = 1.0
            self.params[4] = 1.0
            self.params[5] = (math.lgamma(self.params[0]+self.params[1]+self.params[2]) \
                                   - (math.lgamma(self.params[0]) + math.lgamma(self.params[1]) \
                                           + math.lgamma(self.params[2])))
            self.yvar  = self.ymean*(1.0-self.ymean)/(self.xmean*(1.0-self.xmean))*self.xvar
            self.covar = - self.ymean/(1-self.xmean)*self.xvar
        if self.PDFtype == 'bvb4':
            K = 1 - self.xmean + self.xvar/(1-self.xmean)
            self.params[0] = self.xmean * (self.xmean*(1-self.xmean)/self.xvar - 1) #a1
            self.params[1] = (K - self.ymean - self.yvar/self.ymean) / \
                             ( (1-self.xmean)/self.ymean*(self.ymean+self.yvar/self.ymean) - K) #b1
            self.params[2] = self.params[1] * ( (1-self.xmean)/self.ymean -1) #a1
            self.params[3] = self.params[0] * (1-self.xmean) / self.xmean #b2
            ##self.params[1] = self.params[0] * (1-self.xmean) / self.xmean
            ##self.params[2] = (K - self.ymean - self.yvar/self.ymean) / \
            ##                ( (1-self.xmean)/self.ymean*(self.ymean+self.yvar/self.ymean) - K)
            ##self.params[3] = self.params[2] * ( (1-self.xmean)/self.ymean -1)
            self.params[4] = 1.0
            self.params[5] = (math.lgamma(self.params[0] + self.params[3]) - math.lgamma(self.params[0]) - \
                                   math.lgamma(self.params[3]) ) + \
                             (math.lgamma(self.params[1] + self.params[2]) - math.lgamma(self.params[1]) - \
                                   math.lgamma(self.params[2]) )  
            ##self.params[5] = math.exp(math.lgamma(self.params[0] + self.params[1]) - math.lgamma(self.params[0]) - \
            ##                          math.lgamma(self.params[1]) ) * \
            ##                 math.exp(math.lgamma(self.params[2] + self.params[3]) - math.lgamma(self.params[2]) -
            ##                          math.lgamma(self.params[3]) )  
            self.covar = -self.ymean/(1-self.xmean)*self.xvar
            self.params[3] = self.params[3] - self.params[1] -self.params[2] + 1.0
        if (self.PDFtype in ['bvb6','bvb5', 'sml', 'sml4', 'sml3','smlc']):
            roots = root(self.convolute, guess)
            for i in range(7):
                self.params[i] = roots.x[i]
                self.converged = float(roots.success)

    ##################################################################

    # Calculate difference between desired statistics and statistics that result from the guess parameters
    def convolute(self,guess):
        answer  = [0,0,0,0,0,0,0]
        xtilde  = 0.0; x2tilde = 0.0; ytilde  = 0.0; y2tilde = 0.0; xytilde = 0.0
        x3tilde = 0.0; y3tilde = 0.0; z3tilde = 0.0; fmean   = 0.0; fvar    = 0.0; covf    = 0.0
        pdft = 0.0
        covar = 0.0; covxz = 0.0; covyz=0.0
        if self.PDFtype in ['bvb5','bvb6'] or 'sml' in self.PDFtype:
             for i in range(1,self.nx-1):
                  for j in range(1,self.ny-i-1):
                       self.pdfvals[j,i] = self.calcPDFval(self.xvec[i],self.yvec[j],guess)
             self.calcPDFvals_boundary(guess)
        for i in range(0,self.nx):
            for j in range(0,self.ny-i):
                zval    = 1.0 - self.xvec[i] - self.yvec[j]
                pdfwght = self.pdfvals[j,i] * self.pdfarea[j,i] 
                pdft    = pdft    +                             pdfwght
                xtilde  = xtilde  + self.xvec[i]              * pdfwght
                x2tilde = x2tilde + self.xvec[i]*self.xvec[i] * pdfwght
                ytilde  = ytilde  + self.yvec[j]              * pdfwght
                y2tilde = y2tilde + self.yvec[j]*self.yvec[j] * pdfwght
                covar   = covar   + self.xvec[i]*self.yvec[j] * pdfwght
                covxz   = covxz   + self.xvec[i]*zval         * pdfwght
                covyz   = covyz   + zval        *self.yvec[j] * pdfwght
        #for i in range(0,self.nx-1):
         #   for j in range(0,self.ny-i):
                #fmean   = fmean+ self.yvec[j]/(1-self.xvec[i])  * self.pdfvals[j,i] * self.pdfarea[j,i]
                #fvar    = fvar+(self.yvec[j]/(1-self.xvec[i]))**2  * self.pdfvals[j,i] * self.pdfarea[j,i]
                #fmean   = fmean+ self.yvec[j]/(self.yvec[j]+self.xvec[i])  * self.pdfvals[j,i] * self.pdfarea[j,i]
                #fvar    = fvar+(self.yvec[j]/(self.yvec[j]+self.xvec[i]))**2  * self.pdfvals[j,i] * self.pdfarea[j,i]
                #covf    = covf+self.yvec[j]/(1-self.xvec[i])*self.xvec[i]  * self.pdfvals[j,i] * self.pdfarea[j,i]
        divpdft = 1.0/pdft        
        xtilde  = xtilde *divpdft
        ytilde  = ytilde *divpdft
        x2tilde = x2tilde*divpdft
        y2tilde = y2tilde*divpdft
        covar   = covar  *divpdft
        covxz   = covxz  *divpdft
        covyz   = covyz  *divpdft
        fmean   = fmean  *divpdft
        fvar    = fvar   *divpdft
        covf    = covf   *divpdft
        ztilde = 1 -xtilde - ytilde
        for i in range(0,self.nx):
            for j in range(0,self.ny-i):
                x3tilde = x3tilde + (self.xvec[i] - xtilde)**3 * self.pdfvals[j,i] * self.pdfarea[j,i]
                y3tilde = y3tilde + (self.yvec[j] - ytilde)**3 * self.pdfvals[j,i] * self.pdfarea[j,i]
                z3tilde = z3tilde + (1-self.xvec[i]-self.yvec[j] - ztilde) **3 * self.pdfvals[j,i] * self.pdfarea[j,i]
        xskew = (x3tilde) * divpdft / (x2tilde - xtilde**2)**1.5
        yskew = (y3tilde) * divpdft / (y2tilde - ytilde**2)**1.5
        zskew = (z3tilde) * divpdft
        if self.PDFtype in ['dns'] :
             self.fmean = fmean
             self.fvar  = fvar - fmean**2
             self.covf  = covf - fmean*xtilde
             self.zskew = zskew
        answer[0] = (xtilde                   )             - self.xmean
        answer[1] = (x2tilde   - xtilde**2    )             - self.xvar
        answer[2] = (ytilde                   )             - self.ymean
        answer[3] = (y2tilde   - ytilde**2    )             - self.yvar
        answer[4] = (covar - xtilde*ytilde)                 - self.covar
        #answer[2] = fmean          - self.fmean
        #answer[3] = fvar - fmean**2 - self.fvar
        #answer[4] = covf - fmean*xtilde - self.covf
        answer[5] = pdft - 1.0
        self.pdft = pdft
        self.xskew=xskew
        self.yskew=yskew
        self.covxz = covxz - xtilde * (1-xtilde-ytilde)
        self.covyz = covyz - ytilde * (1-xtilde-ytilde)
        #print 'params'
        #print guess
        #print self.xvar, self.yvar
        #print self.fmean, self.fvar, self.covf
        if self.PDFtype in ['sml4','sml3'] :
            self.covar = covar - xtilde*ytilde
            answer[4] = guess[5]
            if self.PDFtype == 'sml3' :
                self.yvar = y2tilde - ytilde**2
                answer[3] = guess[4]
        if self.PDFtype == 'smlc' :
            self.yvar = y2tilde - ytilde**2
            answer[3] = guess[5]            
            self.covar = covar - xtilde*ytilde
            answer[4] = guess[4]
        if self.PDFtype == 'bvb6' :
            answer[6] = zskew - self.zskew
        else :
            answer[6] = guess[6]
        self.error = answer
        print self.error
        return answer


    ##################################################################

    # Fill self.pdfvals with the appropriate pdf
    def calcPDFfull(self,guess = [1,1,1,1,1,0]):
        self.converged = 1.0
        if self.PDFtype == 'betadelta':
            self.calcPDFparams()
            for i in range(self.nx):
                self.pdfvals[0,i] = stat.beta.pdf(self.xvec[i],self.params[1],self.params[0]) # pdf values
                self.pdfvals[1,self.nx -i-1] = self.xvec[i] * self.params[2] # y values
            self.pdfvals[0,0] = 2*self.pdfvals[0,1]-self.pdfvals[0,2]
            self.pdfvals[0,self.nx-1] = 2*self.pdfvals[0,self.nx-2] - self.pdfvals[0,self.nx-3]
            self.pdft = 1.0
        elif (self.PDFtype == 'dirichlet' or  self.PDFtype == 'bvb4') :
            self.calcPDFparams()
            for i in range(1,self.nx-1):
                for j in range(1,self.ny-1-i):
                    self.pdfvals[j,i] = self.calcPDFval(self.xvec[i], self.yvec[j], self.params)
            self.calcPDFvals_boundary(self.params)
            # calculate higher order moments by convolution and rescale
            approx = self.convolute(self.params)
            self.pdfvals = self.pdfvals / self.pdft
            #print 'xmean' , self.xmean, approx[0]
            #print 'xvar'  , self.xvar,  approx[1]
            #print 'ymean' , self.ymean, approx[2]
            #print 'yvar'  , self.yvar,  approx[3]
            #print 'covar' , self.covar, approx[4]
            #print 'sumPDF', '1.0',      approx[5]
        elif self.PDFtype in ['bvb5','bvb6']:
            self.calcPDFparams(guess)
        elif 'sml' in self.PDFtype :
            self.calcPDFparams(guess)
        elif self.PDFtype == 'dns':
            # get data from NGA hit dns outputs
            filename = self.filepath.replace('monitor','scalarPDF') + '/joint-' + self.timestr
            openfile = open(filename,'r')
            colnames = openfile.readline().strip().split()
            openfile.close()
            for i in range(len(colnames)) :
                if colnames[i] == self.varx + '-' + self.vary  or \
                         colnames[i] == self.vary + '-' + self.varx   :
                    coljoint = i
                    break
            else :
                raise RuntimeError("calcPDFfull : could not find statistics for requested variables")
            dnspdf = np.genfromtxt(filename,skip_header=1,usecols=(coljoint))
            # Allow for filtering of DNS data onto coarser x-y grid
            if self.filt == 0 :
                for kk in range(self.nbins*self.nbins) :
                    self.pdfvals[kk%self.nbins,kk//self.nbins] = dnspdf[kk]
                self.pdfvals = self.pdfvals / np.sum(self.pdfvals) /self.pdfarea
            
            elif self.filt == 1 :
                for kk in range(self.nbins*self.nbins) :
                    xindex =  (kk//self.nbins) 
                    yindex =  (kk% self.nbins)
                    self.pdfvals[yindex  ,xindex  ] = self.pdfvals[yindex  ,xindex  ] + dnspdf[kk]
                    self.pdfvals[yindex  ,xindex+1] = self.pdfvals[yindex  ,xindex+1] + dnspdf[kk]
                    self.pdfvals[yindex+1,xindex  ] = self.pdfvals[yindex+1,xindex  ] + dnspdf[kk]
                    self.pdfvals[yindex+1,xindex+1] = self.pdfvals[yindex+1,xindex+1] + dnspdf[kk]
                self.pdfvals = self.pdfvals / np.sum(self.pdfvals) / self.pdfarea
            elif self.filt > 1 :
                for kk in range(self.nbins*self.nbins) :
                    xindex =  (kk//self.nbins + self.filt/2) // self.filt
                    yindex =  (kk% self.nbins + self.filt/2) // self.filt
                    self.pdfvals[yindex,xindex] = self.pdfvals[yindex,xindex] + dnspdf[kk]
                self.pdfvals = self.pdfvals / np.sum(self.pdfvals) /self.pdfarea
            else :
                raise RuntimeError('calcPDFfull : filt cannot take a negative value for dns data initialization')
            if self.switch == 'x' :
                pdfvalstemp = np.zeros((self.ny,self.nx))
                for j in range(self.ny) :
                    for i in range(self.ny - j) :
                        pdfvalstemp[j,i] = self.pdfvals[j,self.ny-j-i-1]
                self.pdfvals = pdfvalstemp
            elif self.switch == 'y' :
                pdfvalstemp = np.zeros((self.ny,self.nx))
                for i in range(self.nx) :
                    for j in range(self.nx - i) :
                        pdfvalstemp[j,i] = self.pdfvals[self.nx-j-i-1,i]
                self.pdfvals = pdfvalstemp
            elif self.switch == 'z' :
                pdfvalstemp = np.zeros((self.ny,self.nx))
                for i in range(self.nx) :
                    for j in range(self.nx - i) :
                        pdfvalstemp[j,i] = self.pdfvals[i,j]
                self.pdfvals = pdfvalstemp

            # Approximate moments using integration of binned data to see how close it is
            # save skewness because it gets overwritten, then replace the value
            xskew = self.xskew; yskew = self.yskew
            approx = self.convolute([0,0,0,0,0,0,0])
            print 'xmean' , self.xmean, approx[0]
            print 'xvar'  , self.xvar,  approx[1]
            print 'ymean' , self.ymean, approx[2]
            print 'yvar'  , self.yvar,  approx[3]
            print 'covar' , self.covar, approx[4]
            print 'xskew' ,      xskew, self.xskew - xskew
            print 'yskew' ,      yskew, self.yskew - yskew 
            print 'sumPDF', '1.0',      approx[5]
            self.xskew = xskew; self.yskew = yskew

        if self.PDFtype not in ['dns','betadelta'] :
            pdfvalstemp = np.zeros((self.ny,self.nx))
            # in rotated coordinate system, rotate back to original - need to rescale to sum to 1
            sumpdf = 0.0
            if self.switch == '231' :
                for j in range(self.ny) :
                    for i in range(self.ny - j) :
                        pdfvalstemp[j,i] = self.pdfvals[j,self.ny-j-i-1]
                for i in range(self.nx) :
                    for j in range(self.nx - i) :
                        self.pdfvals[j,i] = pdfvalstemp[self.nx-j-i-1,i]
                        sumpdf = sumpdf + self.pdfvals[j,i]*self.pdfarea[j,i]
                self.pdfvals = self.pdfvals/sumpdf
            elif self.switch == '312' :
                for i in range(self.nx) :
                    for j in range(self.nx - i) :
                        pdfvalstemp[j,i] = self.pdfvals[self.nx-j-i-1,i]
                for j in range(self.ny) :
                    for i in range(self.ny - j) :
                        self.pdfvals[j,i] = pdfvalstemp[j,self.ny-j-i-1]
                        sumpdf = sumpdf + self.pdfvals[j,i]*self.pdfarea[j,i]
                self.pdfvals = self.pdfvals/sumpdf

     ##################################################################

    # calculate L2 norm of difference between PDF and DNS PDF
    def calcMetrics(self,dns):
        L2     = 0.0
        L1     = 0.0
        DJS    = 0.0
        DJSabs = 0.0
        dogs   = 0.0
        if self.PDFtype == 'dns' :
            dns.redefdata(self.pdfvals)
        if self.PDFtype not in ['betadelta'] :
            for i in range(0,self.nx) :
                for j in range(0,self.ny-i) :
                    dogs = dogs + self.pdfvals[j,i] * self.pdfarea[j,i]
                    L2 = L2 + (self.pdfvals[j,i] - dns.data[j,i]) ** 2 * self.pdfarea[j,i]
                    L1 = L1 + abs(self.pdfvals[j,i] - dns.data[j,i]) * self.pdfarea[j,i]
                    P = max( dns.data[j,i]   ,1e-50)
                    Q = max(self.pdfvals[j,i],1e-50)
                    M = 0.5 * (P+Q)
                    logP = 0.5 * P * math.log(P/M)
                    logQ = 0.5 * Q * math.log(Q/M)
                    DJS    = DJS    + (    logP +      logQ ) * self.pdfarea[j,i]
                    DJSabs = DJSabs + (abs(logP) + abs(logQ)) * self.pdfarea[j,i]
        #x3ord = self.xskew * self.xvar **1.5
        #y3ord = self.yskew * self.yvar **1.5
        zvar = 2.0*self.covar + self.xvar + self.yvar
        self.covxz = 0.5*(self.yvar - self.xvar - zvar)
        self.covyz = 0.5*(self.xvar - self.yvar - zvar)
        print  'running', self.pdft, dogs
        
        metrics = ['l2norm', 'l1norm', 'yvar'   , 'covar-xy', 'covar-xz', 'covar-yz', 'xskew'   , 'yskew'   , 'DJS', 'DJSabs','success'] 
        values  = [L2      , L1      , self.yvar, self.covar, self.covxz, self.covyz, self.xskew, self.yskew,  DJS ,  DJSabs ,self.converged]
        self.metrics = dict(zip(metrics,values))
        #print 'L2 for pdf', self.PDFtype, self.switch, '=', L2

   
    ##################################################################

    # plot PDFs
    def calcMarginal(self):

        if self.PDFtype not in ['betadelta'] :
             for i in range(self.nx):
                 for j in range(0,self.ny-i):
                     self.xmarg[i] = self.xmarg[i] + self.pdfvals[j,i] * self.dy
                 self.xmarg[i] = self.xmarg[i] - (self.pdfvals[0,i] + self.pdfvals[self.ny-i-1,i])/2.0*self.dy
             for j in range(self.ny):
                 for i in range(0,self.nx-j):
                     self.ymarg[j] = self.ymarg[j] + self.pdfvals[j,i] * self.dx
                 self.ymarg[j] = self.ymarg[j] - (self.pdfvals[j,0] + self.pdfvals[self.nx-j-1,j])/2.0*self.dx
             B1x = (self.xmean*(1-self.xmean)/self.xvar - 1)*self.xmean
             B2x = (self.xmean*(1-self.xmean)/self.xvar - 1)*(1-self.xmean)
             B1y = (self.ymean*(1-self.ymean)/self.yvar - 1)*self.ymean
             B2y = (self.ymean*(1-self.ymean)/self.yvar - 1)*(1-self.ymean)
        else :
            for i in range(self.nx):
                self.xmarg[i] = self.xmargB[i]
                self.ymarg[i] = stat.beta.pdf(self.yvec[i]*(1-self.xmean)/self.ymean,self.params[1] \
                                              ,self.params[0])*(1-self.xmean)/self.ymean
        if self.PDFtype in ['dns'] :
             for i in range(self.nx):
                 self.xmargB[i] = stat.beta.pdf(self.xvec[i],B1x,B2x)
                 self.ymargB[i] = stat.beta.pdf(self.yvec[i],B1y,B2y)

    ################################################################

    # plot PDFs
    def plotPDF(self,fig,ii,frametype,nrows=3,ncols=3, shift=2):
        plt.figure(fig.number)
        if frametype == 'video' :
             ax = plt.subplot2grid((nrows, ncols), (ii/(ncols-shift),shift+ii%(ncols-shift))) #ii/(ncols-2)
             ax.tick_params(labelsize=14)
             plt.text(0.5,0.97,self.PDFlabel,horizontalalignment='center',size=14)
             #if self.PDFtype not in ['betadelta', 'dns'] :
             #plt.text(0.25,0.83,'$||\epsilon||_1 = % 6.3f $' % self.metrics['l1norm'],)
             if (ii) // (ncols-shift) == nrows-1 :
                  plt.xlabel('Z$_1$',size=14)
             else :
                  ax.xaxis.set_ticklabels([])
             if ii % (ncols-shift) == 0 :
                  plt.ylabel('Z$_2$',size=14)
                  # pad=5
                  # ax.annotate(rows[(ii-1)//ncols], xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                  #    xycoords=ax.yaxis.label, textcoords='offset points',
                  #    size='large', ha='right', va='center')
             else :
                  ax.yaxis.set_ticklabels([])
        elif frametype == 'paper' :
             ax = plt.subplot(1,6,ii+1)
             ax.xaxis.set_ticklabels([])
             if ii == 0 :
                  pad=-12.5
                  ax.annotate('$Z_2$', xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                              xycoords=ax.yaxis.label, textcoords='offset points',
                              size=8, ha='right', va='center')
             else :
                  ax.yaxis.set_ticklabels([])
        elif frametype == 'orient' :
             ax = plt.subplot(1,9,ii+1)
             ax.xaxis.set_ticklabels([])
             if ii == 0 :
                  pad=-12.5
                  ax.annotate('$Z_2$', xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                              xycoords=ax.yaxis.label, textcoords='offset points',
                              size=8, ha='right', va='center')
             else :
                  ax.yaxis.set_ticklabels([])       
        else :
             raise RuntimeError("Oh shoot! wrong frametype specified")

        xx,yy = np.meshgrid(self.xvec,self.yvec)

        ax.xaxis.tick_bottom()
        ax.xaxis.set_ticks(np.arange(0,1.5,1.0))
        ax.yaxis.tick_left()
        ax.yaxis.set_ticks(np.arange(0,1.5,1.0))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_aspect('equal')

        self.minlog = 1e-2; minlog =self.minlog
        self.maxlog = 1e2; maxlog = self.maxlog
        nlevels = 100
        base = (maxlog/minlog)**(1.0/nlevels)
        self.levs = [0] * (nlevels+1); levs = self.levs
        for i in range(nlevels+1):
             levs[i] = minlog*pow(base,float(i))
        grf1 = plt.contourf(xx, yy, np.minimum(np.maximum(self.pdfvals,minlog),maxlog), 
                            levels=levs, norm=LogNorm(vmin=minlog, vmax=maxlog),cmap=self.colormap)
        points = [[1,0], [0,1], [1,1]]
        tri = patch.Polygon(points, ec='w', fc='w')
        plt.gca().add_patch(tri)
        points = [[0,0], [1,0], [0,1]]
        #tri2 = patch.Polygon(points, ec='k', fc='white')
        #plt.gca().add_patch(tri2)
        #self.pdfvals = np.minimum( self.pdfvals,0)
        if self.PDFtype == 'betadelta':
                  points = [[0,0], [1,0], [0,1]]
                  tri2 = patch.Polygon(points, ec='k', fc='white')
                  if not (self.switch == 'beta-u') :
                       plt.gca().add_patch(tri2)
                       points = np.array([self.xvec,self.pdfvals[1,:]]).T.reshape(-1,1,2)
                       segments = np.concatenate([points[:-1], points[1:]], axis=1)
                       grf1 = LineCollection(segments, cmap=self.colormap , norm=LogNorm(minlog,maxlog))
                       # Rescale probability to match 2d plots 
                       linewidth=4
                       bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                       wideness = bbox.width
                       rescale_fact = 72.0*wideness/linewidth
                       grf1.set_array((self.pdfvals[0,::-1])*rescale_fact)
                       grf1.set_linewidth(linewidth)
                       #grf1.set_zorder(101)
                       plt.gca().add_collection(grf1)
                  else :
                       yvecs = [[0]*len(self.xvec), 1.0 - self.xvec]
                       weights = [1.0-self.params[2], self.params[2] ]
                       grfs = [0] * 2
                       for iline in range(2) :
                            points = np.array([self.xvec,yvecs[iline]]).T.reshape(-1,1,2)
                            segments = np.concatenate([points[:-1], points[1:]], axis=1)
                            grfs[iline] = LineCollection(segments, cmap=self.colormap , norm=LogNorm(minlog,maxlog))
                            # Rescale probability to match 2d plots 
                            linewidth=8
                            bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                            wideness = bbox.width
                            rescale_fact = 72.0*wideness/linewidth*weights[iline]*2
                            grfs[iline].set_array((self.pdfvals[0,::-1])*rescale_fact)
                            grfs[iline].set_linewidth(linewidth)
                            grfs[iline].set_zorder(1)
                            plt.gca().add_collection(grfs[iline])
                       #points = [[1,1], [1,0], [0,1]]
                       #tri3 = patch.Polygon(points, fc='white')
                       #plt.gca().add_patch(tri3)

        grf2 = plt.plot((1,0),(0,1),color='k')

        if frametype == 'video' :
             axcb = fig.add_axes([0.91,0.65,0.016,0.25])
             axcb.set_aspect('auto')
             minlog = np.floor(np.log10(minlog))
             maxlog = np.ceil(np.log10(maxlog))
             cbar = plt.colorbar(grf1,cax=axcb,cmap=self.colormap,
                            ticks=np.logspace(minlog,maxlog,maxlog-minlog+1), 
                                              format=LogFormatter(10, labelOnlyBase=False))
             cbar.set_label('$P(Z_1,Z_2)$',size=14)
             cbar.ax.tick_params(labelsize=14)
    ##################################################################

    # plot marginal PDFs
    def plotPDFmarg(self, color, label, axx, axy):
         for qq in range(3) :
              if self.PDFtype == 'dns':
                   plt.sca(axx[qq])
                   plt.plot(self.xvec,self.xmarg,'o',label=label, color=color,markersize=2.8)
                   plt.plot(self.xvec,self.xmargB,label='Beta', color='black',linewidth=2.0)
                   plt.sca(axy[qq])
                   plt.plot(self.yvec,self.ymarg,'o',label=label, color=color,markersize=2.8)
                   plt.plot(self.yvec,self.ymargB,label='Beta', color='black',linewidth=1.0)
              else :
                   plt.sca(axx[qq])
                   plt.plot(self.xvec,self.xmarg,label=label, color=color, linewidth=1.0)
                   plt.sca(axy[qq])
                   plt.plot(self.yvec,self.ymarg,label=label, color=color, linewidth=1.0)
    
##################################################################

# create colormap for use in PDF contour plots
def genColormap():
    colors = ((1.0,1.0,1.0),
              (0.0,0.0,1.0),
              (0.0,1.0,1.0),
              (0.0,1.0,0.0),
              (1.0,1.0,0.0),
              (1.0,0.0,0.0))
    ncolor = len(colors)
    step = 1.0 / float(len(colors)-1)
    cdict = { 'red'   : [0]*ncolor,
              'green' : [0]*ncolor,  
              'blue'  : [0]*ncolor }
    for ic in range(ncolor):
         cdict['red'  ][ic] = (ic*step, colors[ic][0],colors[ic][0])
         cdict['green'][ic] = (ic*step, colors[ic][1],colors[ic][1])
         cdict['blue' ][ic] = (ic*step, colors[ic][2],colors[ic][2])

    colormap = LinearSegmentedColormap('whitebluered',cdict)
    return colormap

                            
