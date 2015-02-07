from pylab import *
#import scipy.interpolate as spintp
import scipy.interpolate as intp
import scipy.optimize
import pdb
import copy, itertools


#### need to impliment
# pass both counts and monitor into initialization
# parameter should return either counts or counts/mon
# add some assert commands to improve robustness of program

"""
s is single xye class instance
S is list of xye objects
"""

class xye(object):
    """General data processing class

    Parameters:
    -----------
            xye.x :         ndarray
                            x data
            xye.y :         ndarray
                            y data in raw counts
            xye.err :       ndarray
                            y error in raw counts
            xye.mon :       ndarray or single value in numpy array
                            montor values either for individal points or for
                            the overall spectrum. Can be either a real counter 
                            or per second. Default = 1
            xye.label :     string
                            label associate with xye. 
                            This is used in the plot legend
            xye.plotmon :   True | False (default)
                            If true data is plotted as counts over monitor    
    """ 
    def __init__(self, x, y, err, mon=array([1]), label='', plotmon=False):
        self.x = self.__checkForArray(x)
        self.y = self.__checkForArray(y)
        self.err = self.__checkForArray(err)
        mon = self.__checkForArray(mon)
        if len(mon) == 1:
            mon = tile(mon, shape(self.x))
        self.mon = mon
        self.label = label
        self.plotmon = plotmon
        self.__orderdata()
            
    def __orderdata(self):
        """Put data in order of increasing x"""
        xorder = argsort(self.x)
        self.x = self.x[xorder]
        self.y = self.y[xorder]
        self.err = self.err[xorder]
        self.mon = self.mon[xorder]

    def __checkForArray(self, a):
        """Check that input is a numpy array and return array. 
        if a list in input convert to numpy array and then return it"""
        if type(a) == list:
            a = array(a)
        elif type(a) == int or type(a) == float:
            a = array([a+ 0.0])
        if self.__isnumeric(a):
            return a
        else:
            print ("Is x, y, or err is not numeric")

    def __isnumeric(self, obj):
        """Check that input array is numeric."""
        try:
            obj+obj, obj-obj, obj*obj
        except ZeroDivisionError:
            return True
        except Exception:
            return False
        else:
            return True

    def __add__(self, obj):
        """ Return new class with obj added to input
        
        obj can be either class or constant        
        
        this method returns new object with modified label
        This can be used for quick plotting without modifying the class
        using the operator +
        """
        Stemp = xye(x=self.x, y=self.y, err=self.err, mon=self.mon,
                    plotmon=self.plotmon)
        Stemp.add(obj)
        if type(obj) == xye:
            newlabel = self.label + " + " + obj.label
        else:
            newlabel = self.label + " + " + str(obj)
        Stemp.label = newlabel
        return Stemp

    def __sub__(self, obj):
        """Return new class with obj subtracted from input
        
        obj can be either class or constant        
        
        this method returns new object with modified label
        This can be used for quick plotting without modifying the class
        using the operator -
        """
        Stemp = xye(x=self.x, y=self.y, err=self.err, mon=self.mon)
        Stemp.sub(obj)
        if type(obj) == xye:
            newlabel = self.label + " - " + obj.label
        else:
            newlabel = self.label + " - " + str(obj)
        Stemp.label = newlabel
        return Stemp

    def __mul__(self, obj):
        """Return new class with obj multiplied by input
        
        obj can be either class or constant        
        
        this method returns new object with modified label
        This can be used for quick plotting without modifying the class
        using the operator +
        """
        Stemp = xye(x=self.x, y=self.y, err=self.err, mon=self.mon,
                    plotmon=self.plotmon)
        Stemp.mul(obj)
        if type(obj) == xye:
            newlabel = self.label + " * " + obj.label
        else:
            newlabel = self.label + " * " + str(obj)
        Stemp.label = newlabel
        return Stemp
 
    def __div__(self, obj):
        """ Return new class with obj divided by input
        
        obj can be either class or constant        
        
        this method returns a new object with modified label
        This can be used for quick plotting without modifying the class
        using the operator /
        """
        Stemp = xye(x=self.x, y=self.y, err=self.err, mon=self.mon,
                    plotmon=self.plotmon)
        Stemp.div(obj)
        if type(obj) == xye:
            newlabel = self.label + " / " + obj.label
        else:
            newlabel = self.label + " / " + str(obj)
        Stemp.label = newlabel
        return Stemp

    def __truediv__(self, obj):
        """ Return new class with obj divided by input
        
        obj can be either class or constant        
        
        this method returns a new object with modified label
        This can be used for quick plotting without modifying the class
        using the operator /
        """
        Stemp = xye(x=self.x, y=self.y, err=self.err, mon=self.mon,
                    plotmon=self.plotmon)
        Stemp.div(obj)
        if type(obj) == xye:
            newlabel = self.label + " / " + obj.label
        else:
            newlabel = self.label + " / " + str(obj)
        Stemp.label = newlabel
        return Stemp
        
    def __mod__(self, shift):
        """ Return new class with x shifted by input
        
        shiftx must be a float or int       
        
        this method returns a new object with modified label
        This can be used for quick plotting without modifying the class
        using the operator %
        """
        Stemp = xye(x=self.x, y=self.y, err=self.err, mon=self.mon, 
                    plotmon=self.plotmon)
        Stemp.shift_x(shift)
        newlabel = self.label + " % " + str(shift)
        Stemp.label = newlabel
        return Stemp
        
    def __and__(self, obj):
        """ Combine self and obj and return a new object with modified label
        This can be used for quick plotting without modifying the class using
        the operator +=
        """
        Stemp = xye(x=self.x, y=self.y, err=self.err, mon=self.mon, 
                    plotmon=self.plotmon)
        Stemp.combine(obj)
        newlabel = self.label + " += " + obj.label
        Stemp.label = newlabel
        return Stemp
    
    def gety(self, xvals, fill_value=NaN):
        """Use interpolation to return axis_to_get = y, err, mon at xvals"""
        func = intp.interp1d(self.x, self.y, kind='linear', 
                        bounds_error=False, fill_value=fill_value)
        return func(xvals)

    def geterr(self, xvals, fill_value=NaN):
        """Use interpolation to return axis_to_get = y, err, mon at xvals"""
        func = intp.interp1d(self.x, self.err, kind='linear', 
                        bounds_error=False, fill_value=fill_value)
        return func(xvals)

    def getmon(self, xvals, fill_value=NaN):
        """Use interpolation to return axis_to_get = y, err, mon at xvals"""
        func = intp.interp1d(self.x, self.mon, kind='linear', 
                        bounds_error=False, fill_value=fill_value)
        return func(xvals)

    def throwawayNaN(self):
        """Remove all Nan values from self"""
        inds = ~np.logical_or( np.logical_or(isnan(self.x), isnan(self.y)),
                              np.logical_or(isnan(self.err), isnan(self.mon)))
        self.x = self.x[inds]
        self.y = self.y[inds]
        self.err = self.err[inds]
        self.mon = self.mon[inds]

    def combine(self, obj):
        """Combine data in obj into self"""
        self.x = hstack((self.x, obj.x))
        self.y = hstack((self.y, obj.y))
        self.err = hstack((self.err, obj.err))
        self.mon = hstack((self.mon, obj.mon))
        self.__orderdata()
        
    def add(self, obj):
        """ Add obj onto xye class
        obj can be either class or constant
        Throws away data where spectrum
        to add does not overlap        
        """
        if type(obj) == xye:
            self.y = self.y + obj.gety(self.x)
            self.err = sqrt( self.err**2 + (obj.geterr(self.x))**2 )
            self.mon = self.mon + obj.getmon(self.x)
        elif self.__isnumeric(obj):
            self.y = self.y + obj
        else:
            raise Exception("argument to add neither xye nor number!")
        self.throwawayNaN()
        
    def sub(self, obj):
        """ Subtract obj from xye class
        obj can be either class or constant
        Does not alter the monitor values. Throws away data where spectrum
        to subtract does not overlap
        """
        if type(obj) == xye:                                     
            self.y = self.y - obj.gety(self.x)
            self.err = sqrt( self.err**2 + obj.gety(self.x)**2 )
            self.mon = self.mon + obj.getmon(self.x)
        elif self.__isnumeric(obj):
            self.y = self.y - obj
        else:
            raise Exception("argument to sub neither xye nor number!")

    def mul(self, obj):
        """ Multipy xye class by obj
        obj can be either class or constant
        Does not alter the monitor values. Throws away data where spectrum
        to multiply does not overlap
        """
        if type(obj) == xye:
            self.y = self.y * obj.gety(self.x)
            
            rel_err = sqrt( (self.err/self.y)**2 + 
            (obj.geterr(self.x)/obj.gety(self.x))**2 )
            self.err = self.y * rel_err
        elif self.__isnumeric(obj):
            self.y = self.y * obj
            self.err = self.err * obj
        else:
            raise Exception("argument to mul neither xye nor number!")
        self.throwawayNaN()

    def div(self, obj):
        """ Divide xye class by obj
        obj can be either class or constant
        Does not alter the monitor values. Throws away data where spectrum
        to divide does not overlap
        """
        if type(obj) == xye:
            self.y = self.y / obj.gety(self.x, fill_value=NaN)
            rel_err = sqrt( (self.err/self.y)**2 +
            (obj.geterr(self.x)/obj.gety(self.x))**2 )
            self.err = self.y * rel_err
            if any( obj.gety(self.x)  == 0):
                raise Exception("Attempt to divde by zero in self.div")
        elif self.__isnumeric(obj):
            if obj == 0:
                raise Exception("Attempt to divde by zero in self.div")
            self.y = self.y / obj
            self.err = self.err / obj
        else:
            raise Exception("argument to div neither xye nor number!")
        self.throwawayNaN()
        
    def get_intBG(self, xmin=-inf, xmax=+inf):
        """ Integrate up linear background based on trapezium rule."""
        if xmax <= xmin:
            raise Exception("wimax <= xmin")
        
        y_xmin = self.get_yval(xmin)
        y_xmax = self.get_yval(xmax)
        
        if y_xmin<0 or y_xmax<0:
            raise Exception("Negative y values")
        
        I_BG = (xmax - xmin) * (y_xmin + y_xmax) / 2
        
        return I_BG
    
    def get_intval(self, xmin=-inf, xmax=+inf):
        """ Return numerical integral of function between xmin and xmax.
        The integral is calculated using Trapezoidal rule method
        """
        self.__orderdata()
        Stemp = copy.copy(self)
        Stemp.cut(xmin,xmax)
        xcut = Stemp.x
        ycut = Stemp.y
        
        endind = len(xcut)
        dx = xcut[1:endind] - xcut[0:(endind-1)]
        dy = (ycut[0:(endind-1)] + ycut[1:endind])/2
        I = sum(dx*dy)
        return I
    
    def get_peakval(self, xmin=-inf, xmax=+inf):
        """ Get x, y, err values at the maximum Y value""" 
        self.__orderdata()
        Stemp = copy.copy(self)
        Stemp.cut(xmin,xmax)
        xcut = Stemp.x
        ycut = Stemp.y
        ecut = Stemp.e
        
        maxind = argmax(ycut)
        peakX = xcut[maxind]
        peakY = ycut[maxind]
        peakE = ecut[maxind]
        return peakX, peakY, peakE
    
    def get_overlap(self, S2, xmin=-inf, xmax=+inf):
        if xmax <= xmin:
            raise Exception("wimax <= xmin")
        self.__orderdata()
        S2.__orderdata()
        
        xin1 = self.x
        yin1 = self.y
        xin2 = S2.x
        yin2 = S2.y

        overlapmin = max(xmin, min(xin1), min(xin2) )
        overlapmax = min(xmax, max(xin1), max(xin2) )

        if overlapmax <= overlapmin:
            raise Exception("No overlap to perform ")
  
        # determine number step to oversample
        dx1 = mean(xin1[1:len(xin1)] - xin1[0:(len(xin1) - 1)])
        dx2 = mean(xin2[1:len(xin2)] - xin2[0:(len(xin2) - 1)])
        NoStep = 5* (overlapmax - overlapmin) / min(dx1, dx2)
        
        if NoStep <= 10:
            raise Exception("Insufficient number of steps for convolution")
        
        x = scipy.linspace(overlapmin, overlapmax, NoStep)
        dx = x[1]-x[0]
        
        fdatin1 = intp.interp1d(xin1, yin1, kind='linear',
                                  bounds_error = True)
        fdatin2 = intp.interp1d(xin2, yin2, kind='linear',
                                  bounds_error = True)                                                       
        y1 = fdatin1(x)
        y2 = fdatin2(x)
        
        if argmax(y1) == 0 or argmax(y2)== 0 or argmax(y1) == (len(y1)-1) or argmax(y2) == (len(y2)-1):
            raise Exception("Overlap region doesn't seem to contain a peak")
        
        # compute the cross-correlation between y1 and y2
        ycorr = scipy.correlate(y1, y2, mode='full')
        xcorr = scipy.linspace(0, len(ycorr)-1, num=len(ycorr))
        
        # define a gaussian fitting function where
        # p[0] = amplitude
        # p[1] = mean
        # p[2] = sigma
        fitfunc = lambda p, x: p[0]*scipy.exp(-(x-p[1])**2/(2.0*p[2]**2))
        errfunc = lambda p, x, y: fitfunc(p,x)-y
        
        # guess some fit parameters
        p0 = scipy.c_[max(ycorr), scipy.where(ycorr==max(ycorr))[0], 5]
        # fit a gaussian to the correlation function
        p1, success = scipy.optimize.leastsq(errfunc, p0.copy()[0],
                                             args=(xcorr,ycorr))
        
        # convert index to lag steps
        # the first point has index=0 but the largest (negative) lag
        # there is a simple mapping between index and lag
        nLags = p1[1]-(len(y1)-1)

        return -nLags*dx
        #x,  indexshift, dx, x, xcorr, ycorr, yselfcorr
    
    def get_binedges(self, binwidth):
        """Get x array covering spectrum separated by binwidth"""
        steps = floor((max(self.x) -  min(self.x)) / binwdith)
        return linspace(min(self.x), max(self.x), steps)
        
    def rebin(self, binedges):
        """Rebin the spectrum in bins set by edges"""
        xout = array([])
        yout = array([])
        errout = array([])
        monout = array([])
    
        for i in range(1, len(binedges)):
            inds = (self.x > binedges[i-1]) & (self.x <= binedges[i])
            ycurr = self.y[inds]
            errcurr = self.err[inds]
            moncurr = self.mon[inds]
        
            if len(ycurr>0):
                xout = hstack((xout, (binedges[i] + binedges[i-1])/2 ))
                yout = hstack((yout, np.sum(ycurr)))
                errout = hstack((errout, sqrt(np.sum(errcurr**2))))
                monout = hstack((monout, np.sum(moncurr)))

        self.x = xout
        self.y = yout
        self.err = errout
        self.mon = monout
        
    def get_peakfit(self, xmin=-inf, xmax=+inf):
        if xmax <= xmin:
            raise Exception("wimax <= xmin")
        #pdb.set_trace()
        Stemp = copy.copy(self)
        Stemp.cut(xmin=-inf, xmax=+inf)

        # gaussian function
        fitfunc = lambda p, x: p[0]* exp(-0.5 * (x - p[1])**2 / p[2]**2) + p[3] 
        # error function
        errfunc = lambda p, x, y: fitfunc(p, x) - y 
        
        # initial guess
        ind = argmax(Stemp.y)
        pic = Stemp.y[ind] - min(Stemp.y)
        cen = Stemp.x[ind]
        offset = min(Stemp.y)
        wid = abs(xmax - xmin)/4
        p0 = [pic, cen, wid, offset] # Initial guess for the parameters
        p1, success = scipy.optimize.leastsq(errfunc, p0[:],
                                             args=(Stemp.x, Stemp.y))
        if ((success == 1) or (success == 2) or (success == 3) or (success == 4)):
            print 'Fit was sucessful'
        else:
            raise Exception("Fit did not work") 
        return p1[1], p1[0]

    def shift_x(self, shift):
        """ Add the input value to all of self.x"""
        self.x = self.x + shift

    def shift_y(self, shift):
        """ Add the input value to all of self.y"""
        self.y = self.y + shift
        
    def normint(self, xmin=-inf, xmax=+inf):
        """ Normalize the spectrum by the integral between xmin and xmax"""
        if xmax <= xmin:
            raise Exception("wimax <= xmin")
        
        I = self.get_intval(xmin=xmin, xmax=xmax)
        self.y = self.y / I
        self.err = self.err / I
        
    def normint_subBG(self, xmin=-inf, xmax=+inf):
        """ Normalize the spectrum by the integral between xmin and xmax.
        A linear background is removed """
        if xmax <= xmin:
            raise Exception("wimax <= xmin")
        I = self.get_intval(xmin=xmin, xmax=xmax)
        I_BG = self.get_intBG(xmin=xmin, xmax=xmin)
        
        if I_BG>= I:
            raise Exception("Background bigger than peak!")
        
        self.y = self.y / (I - I_BG)
        self.err = self.err / (I - I_BG)
        
    def cut(self, xmin=-inf, xmax=+inf):
        """Throw away values set by xmin and xmax.
        If xmax < xmin central region is thrown away
        If xmax > xmin edges are thrown away"""
        if xmax < xmin:
            indices = (self.x > xmin) | (self.x < xmax)
        elif xmax > xmin:
            indices = (self.x > xmin) & (self.x < xmax)
        else:
            indices = isreal(self.x)
        
        self.x = self.x[indices]
        self.y = self.y[indices]
        self.err = self.err[indices]
        self.mon = self.mon[indices]
        
    def errorbar(self, *args, **kwargs):
        """ Method for errorbar plot of spectrum"""
        if self.plotmon:
            errorbar(x=self.x, y=self.y/self.mon, yerr=self.err/self.mon, label=self.label,
                 *args, **kwargs)
        else:
            errorbar(x=self.x, y=self.y, yerr=self.err, label=self.label,
                 *args, **kwargs)
        if self.label!=None:
            legend()

    def plot(self, *args, **kwargs):
        """ Method to plot spectrum """
        if self.plotmon:
            plot(self.x, self.y/self.mon, label=self.label, *args, **kwargs)
        else:
            plot(self.x, self.y, label=self.label, *args, **kwargs)
        if self.label!=None:
            legend()

##Some functions
        
def unpack(*obj):
    """convert spectra into list with no nesting"""
    if type(obj) == xye:
        S = [obj]
    else:
        def flatten(x):
            result = []
            for el in x:
                if hasattr(el, "__iter__") and not isinstance(el, basestring):
                    result.extend(flatten(el))
                else:
                    result.append(el)
            return result
        S = flatten(obj)
    return S
    
def isxye(obj):
    """ Print warning if input arguments aren't xye class instances. """
    S = unpack(obj)
    for s in S:
        if type(s) != xye:
            print "Some arguments are not xye class instances"

def addxye(*obj):
    """ sum spectra
    obj can be any list containing xye objects
    
    The resulting spectra will share the same axis as the first spectrum
    subsequent spectra are interpolated onto the first x-axis    
    
    """
    S = unpack(obj)
    isxye(S)
    
    if len(S) < 2:
        print "Only one spectrum was sent to addxye. The same spectrum will be returned"
    
    #Sout = xye(x=S[0].x, y=S[0].y, err=S[0].err, label=S[0].label)
    Sout = copy.copy(S[0])
    for s in S[1:]:
        Sout.add(s)
        Sout.label += ' + ' + s.label
    
    return Sout

def plotxye(obj, *args, **kwargs):
    """ plot spectra
    
    Arguments are:
    *obj spectra passed either as individual spectra or lists of spectra
    in any configuration
    additional arguments or keyword arguments are passed onto plot fucntion
    
    hold = True | False (Default True)
    If hold = False window is initially cleared. Subsequent plots are always
    kept regardless
    """

    if kwargs.has_key('hold'):
        if kwargs['hold'] == False:
            clf()
        del(kwargs['hold'])
    
    S = unpack(obj)
    isxye(S)
    
    # plot everything
    color=iter(cm.rainbow(np.linspace(0,1,len(S))))
    for s in S:
        c = next(color)
        if len(S) >7:
            s.plot(*args, c=c, **kwargs)
        else:
            s.plot(*args, **kwargs)

def errorbarxye(obj, *args, **kwargs):
    """ errorbar spectra
    
    Arguments are:
    *obj spectra passed either as individual spectrum or lists of spectra
    in any configuration
    additional arguments or keyword arguments are passed onto plot fucntion
    
    hold = True | False (Default True)
    If hold = False window is initially cleared. Subsequent plots are always
    kept regardless
    """
    #pdb.set_trace()
    if kwargs.has_key('hold'):
        if kwargs['hold'] == False:
            clf()
        del(kwargs['hold'])
    
    S = unpack(obj)
    isxye(S)

    # plot everything
    color=iter(cm.rainbow(np.linspace(0,1,len(S))))
    for s in S:
        if type(s) == xye:
            c=next(color)
            if len(S) >7:
                s.errorbar(*args, c=c, **kwargs)
            else:
                s.errorbar(*args, **kwargs)

def lineup(obj):
    """ Line up spectra
    input is list of xye class instances in any order
    output is list of xye class instances referenced to 1st class instance 
    """
    S = unpack(obj)
    for i in range(1,len(S)):
        shiftx = S[i-1].get_overlap(S[i], 0.8, 3.5)
        print "i = " + str(i) + "  shiftx = " + str(shiftx)
        S[i].shift_x(-shiftx)
    return S

def normintxye(obj, xmin=-inf, xmax=+inf):
    """ Line up spectra
    input is list of xye class instances in any order
    output is list of xye class instances referenced to 1st class instance 
    """
    S = unpack(obj)
    for i in range(len(S)):
        S[i].normint(xmin=xmin, xmax=xmax)
    return S