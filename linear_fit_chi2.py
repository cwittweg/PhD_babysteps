import numpy as np
from iminuit import Minuit
from iminuit.util import make_func_code

class Chi2Functor:
    def __init__(self,f,x,y,yerr):
        self.f = f
        self.x = x
        self.y = y
        self.yerr = yerr
        f_sig = describe(f)
        #this is how you fake function
        #signature dynamically
        self.func_code = make_func_code(f_sig[1:])#docking off independent variable
        self.func_defaults = None #this keeps np.vectorize happy
    def __call__(self,*arg):
        #notice that it accept variable length
        #positional arguments
        chi2 = sum(((y-self.f(x,*arg))**2)/(yerr**2) for x,y,yerr in zip(self.x,self.y,self.yerr))
        return chi2

x= ...
y= ...
yerr= ...
def linear(x,a,b):
    return np.add(a,np.multiply(b,x))
linear_chi2=Chi2Functor(linear,x,y,yerr)
m = Minuit(linear_chi2)
m.migrad()

ndf=len(x)-len(m.args)
chi2=m.fval
a_fit=m.fitarg['a'])
a_fit_error=m.fitarg['error_a']
