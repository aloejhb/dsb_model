def spline_pfunc(a):
    if a > 10:
        raise Exception('a should be less than 10!', ('spline_pfunc', a))

    return a ** 2

def mk_pfuncs(a):
    b = spline_pfunc(a)
    return b

def setAllPars(a):
    if a < 0:
        raise Exception('a should be more than 0!', a)

    c = None
    c = mk_pfuncs(a)
    
    return c

def objFunc(a):
    d = None
    try:
        d = setAllPars(a)
    except Exception as e:
        print 'Warning: value of a is invalid:'
        print str(e)
        print e.args[0]
        print e.args[1]
        return 100000

    return d
