
D = 122500
l = 0.5
space = 7.5
x_initial =  x
#laticeset = [cell()] * int(space/l)
laticeset = [0] * int(space/l) * n * m
t = 0

Tb = 300

tend = 2*Tb

C = []
C2 = []
C4 = []
T = []

px1 = float(sum(liststates))/float(len(liststates))
px0 = 1 - px1

py1 = float(sum(p4))/float(len(p4))
py0 = 1 - py1

if py1==0:
    py1=0.00000000001
    py0 =0.9999999999

py1x0 = float(sum(p2))/float(len(p2))
py0x0 = 1 - py1x0

if py1x0==0:
    py1x0=0.00000000001
    py0x0 =0.9999999999

#py1x1 = (px1*py1)/py1
py1x1 = py1+py1x0-(py1*py1x0)
py0x1 = 1 - py1x1

print px0,px1,py0,py1
print py1x0,py0x0,py1x1,py0x1
print py1*py1x0,py1*py0x0,py1*py1x1,py1*py0x1

from math import log

if (calc_capacity):
    ibits= px0*py0x0*(log(py0x0/(py0))/log(2))+px1*py0x1*(log(py0x1/(py0))/log(2))+px0*py1x0*(log(py1x0/(py1))/log(2))+px1*py1x1*(log(py1x1/(py1))/log(2))

    print ''
    print 'Informação Mútua:', ibits