import math
import textwrap
#from sympy import *
import matplotlib.pyplot as plt
from numpy import genfromtxt
import csv

tcr_set = [round(x*0.1,1) for x in range(349,377,2)]
tsk_set = [round(x*0.1,1) for x in range(319,347,2)]
swr = [5,10,20,30,40,60,80,100]
bbfr = [12.6,10.08,7.56,5.04,2.52,1.26,.63,.315,.1575,.07875]
vd = [3.75, 7.5, 15, 30, 60, 90 , 120, 150]
vc = [1,0.8,0.6,0.4,0.2,0.1,0.05,0.025,0.0125,0.00625]
