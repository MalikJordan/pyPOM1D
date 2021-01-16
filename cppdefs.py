# THIS FILE IS INCLUDED IN ALL FILES
RELEASE = 'Created by BFM v. 5.1'

PATH_MAX = 255

stderr = 0
stdout = 6

# HANDY FOR WRITING
def STDOUT(text):
    print(text)

def STDERR(text):
    print(text)

# STANDARD OUTPUT FOR PARALLEL COMPUTATION
def LEVEL0():
    call0 = ''
    STDERR(call0)

def LEVEL1():
    call1 = '   '
    STDERR(call1)

def LEVEL2():
    call2 = '       '
    STDERR(call2)

def LEVEL3():
    call3 = '           '
    STDERR(call3)

def LEVEL4():
    call4 = '               '
    STDERR(call4)

def FATAL():
    callf = 'FATAL ERROR: '
    STDERR(callf)

def LINE():
    print('------------------------------------------------------------------------')


# SHAPE OF VARIABLES
POINT       = 0
Z_SHAPE     = 1
T_SHAPE     = 2
XY_SHAPE    = 3
XYT_SHAPE   = 4
XYZT_SHAPE  = 5
OCET_SHAPE  = 6
SURFT_SHAPE = 7
BOTT_SHAPE  = 8
G_SHAPE     = 9
XYZ_SHAPE   = 10

# CONSTANTS FOR AVERAGE COMPUTATIONS
INIT       = 0
MEAN       = 1
RESET      = 2
ACCUMULATE = 10

# TO AVOID DIVIDING BY ZERO
SMALL = 1E-08

# WHAT PRECISION WE WILL USE IN THIS COMPILATION
_ZERO_ = 0.0
_ONE_ = 1.0



