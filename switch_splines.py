from __future__ import division
from functools import partial
from numpy import *
from numpy import concatenate as conc
from scipy.interpolate import UnivariateSpline as spline
import sys



def mk_SwitchSpline(level1,switch1,delay1,level2,switch2=None,delay2=None,level3=None,maxT = 100, debug = False, smord = 3, smooth = 0):

	''' General two switch rate defined by spline interpolation of 5 points in the t-rate(t) plane given by: 
	    (0,level1),(switch1,level1),(switch1+delay1,level2),(switch2,level2),(switch2+delay2,level3)
		Caution with too short delays (boxes), induces spline oscillations at the borders!
		Returns a spline function.
	'''

	# only one switch
	if switch2 == None:
		# if delay is too small, then do first order spline
                if delay1 < 0.15:
                        smord = 1
                        smooth = 0.5
                        # raise Exception('delay smaller than 0.15!')

                if switch1 > maxT:
                        yvec = level1*ones(100)
                        xvec = linspace(0, maxT, 100)
                elif switch1 + delay1 > maxT:
                        n1 = 100

                        yvec = conc( [level1*ones(100),level2*ones(n1)] )
                        xvec = conc( [linspace(0,switch1,100),linspace(switch1+delay1, switch1+delay1+10,n1)] )
                        
                else:
                        n1 = (maxT-switch1)/(0.1*delay1)
                        n1 = n1 if n1>100 else 100

                        yvec = conc( [level1*ones(100),level2*ones(n1)] )
                        xvec = conc( [linspace(0,switch1,100),linspace(switch1+delay1,maxT,n1)] )

	else:
                if delay1 < 0.15 or delay2 < 0.15:
                        smord = 1
                        smooth = 0.5
                        # raise Exception('delay smaller than 0.15!')

		if maxT < switch2:
			print switch1,delay1,switch2
			raise Exception('maxT smaller than 2nd switching time!')
	
		if switch1 + delay1 > switch2:
                        print 'switch1', switch1, 'delay1', delay1, 'switch2', switch2
			raise Exception('First switch plus delay bigger than 2nd switch!')

		# adjust interpolation point by given delays, minimum 100 points!
		n1 = (switch2-switch1)/(0.1*delay1)
		n1 = n1 if n1>100 else 100

		n2 = (maxT-switch2)/(0.1*delay2)
		n2 = n2 if n2>100 else 100
		yvec = conc( [level1*ones(100),level2*ones(n1),level3*ones(n2)] )
		xvec = conc( [linspace(0,switch1,100),linspace(switch1+delay1,switch2,n1),linspace(switch2+delay2,maxT,n2)] )

	s = spline(xvec,yvec,k = smord, s = smooth)	

	if debug:
		return s,xvec,yvec  # to check the xvec and yvec generated from the parameters
	return s

