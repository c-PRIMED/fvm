import sys
import math
from numpy import *

#***********************************************************************#
#operating pressure
pressure   = 83593
# beam mode frequency 
#frequency = 114415
frequency =759734
omega = 2*pi*frequency

# beam density
rho = 2330

# length (change to 1 for 2d)

length = 1e-4
#length = 1

# beam width
b = 35e-6

# beam height
h = 1e-6

#beam is given a velocity as U = A*cos(omega*t)
A = 0.1

## symmetry factor (Set to 1 when doing whole beam, 2 when doing half)
symmetryFactor = 2

#time step
#deltaT=5.0E-08;
deltaT = 1.0/(frequency*100.)

# how many cycles in force data are used to calculate damping factor
# the more cycles, the accurate the result is, but it takes longer to run

nCycles = 10

# write force to file option

forceOutput = True

# write forceExt to file option

forceExtOutput = True

# write forceIntegral to file option

forceIntegralOutput = True

#*************************************************************************#


#read pressure integral data from top.dat and bot.dat
#which give the force on both surfaces as a function of time
#
#fileBase = "/home/sm/prism-meshes/case1/unsteady-p=11145-prism-pIntegral-"
#fileBase = "/home/sm/prism-meshes/case1/fluent-results/P=83593-"
fileBase = "/home/sm/prism-meshes/3dbeam/mode2EV/new-beam-114k-p=45862-prism-pIntegral-"

fileSuffix = ".xy"
#fileSuffix = ".out"

fileName = fileBase + "bot" + fileSuffix

fileBot = open(fileName,"r")

lines = fileBot.readlines()

nLines = len(lines)

forceBot = zeros(nLines, float)

time = zeros(nLines, float)

count = 0

for line in lines:
    value = line.split()
    time[count] = float(value[0])
    forceBot[count] = float(value[1])
    count = count+1
#print time, forceBot

fileBot.close()


fileName = fileBase + "top" + fileSuffix

fileTop = open(fileName,"r")

lines = fileTop.readlines()

nLines = len(lines)

forceTop = zeros(nLines, float)

count = 0

for line in lines:
    value = line.split()
    #time[count] = float(value[0])
    forceTop[count] = float(value[1])
    count = count+1
#print time, forceTop

fileTop.close()

print ("read in pressure integral done!")

#*************************************************************************#
#calculate net force

force = zeros(nLines,float)

for i in range(0, nLines):
    force[i] = forceTop[i]+forceBot[i]
   
if (forceOutput):
    fileName = fileBase + "netForce.dat"
    file = open(fileName,"w")

    for i in range(0, nLines):
         file.write("%e\t" % time[i])
         file.write("%e\n" % force[i])

    file.close()

print ("calculate net force done!")
#*************************************************************************#
#extract a few full cycles in force data and extend to many cycles as needed#

checkCount = 0
checkPoint = zeros(50, int)
for i in range(1, nLines):
    if (force[i]*force[i-1] <= 0):
        checkPoint[checkCount] = i
        checkCount = checkCount+1

print 'check count = %s' % checkCount
startPoint = checkPoint[2]
endPoint = checkPoint[4]

lengthPerCycle = endPoint-startPoint

totalLength = lengthPerCycle*nCycles

timeExt = zeros(totalLength, float)
forceExt = zeros(totalLength, float)

#for i in range(0,startPoint):
#    timeExt[i]=time[i]
#    forceExt[i]=force[i]

for j in range(0,nCycles):
    for i in range(0, lengthPerCycle):
        timeExt[(j*lengthPerCycle+i)] = time[startPoint+i]+j*deltaT*lengthPerCycle
        forceExt[(j*lengthPerCycle+i)] = force [startPoint+i]

if (forceExtOutput):
    fileName = fileBase + "netForceExt.dat"
    file = open(fileName,"w")
    for i in range(0, totalLength):
        file.write("%e\t" % timeExt[i])
        file.write("%e\n" % forceExt[i])
    file.close()

print ("calculate force extension done!")
#*************************************************************************#
#calculate the time integral of force#

## we are essentially computing the Fourier coefficients, ie. assuming
## p(t) = A cos(omega*t) + B sin(omega*t)
##
## A and B can be found as
## A = 2*integral (p(t)*cos(omega*t))
## B = 2*integral (p(t)*sin(omega*t))
##
## with the integrals being over a large enough number of cycles


sumForce = zeros(totalLength, float)
sumForceI = zeros(totalLength, float)

for t in range(0, totalLength):
    for i in range(0, t):
        sumForce[t] = sumForce[t] + forceExt[i]*cos(omega*timeExt[i])*deltaT;
        sumForceI[t] = sumForceI[t] + forceExt[i]*sin(omega*timeExt[i])*deltaT;
    if(t > 0):
        sumForce[t] = sumForce[t]*2*symmetryFactor/(t*deltaT);
        sumForceI[t] = sumForceI[t]*2*symmetryFactor/(t*deltaT);

if (forceIntegralOutput):
    fileName = fileBase + "forceIntegral.dat"
    file = open(fileName,"w")
    for i in range(0, totalLength):
        file.write("%e\t" % timeExt[i])
        file.write("%e\n" % sumForce[i])
    file.close()
print ("calculate force integral done!")
#*************************************************************************#

#cut sumForce into N pieces, calculate the average, max and min value of last piece

N = 10
gap = int(totalLength/N)
average = 0.0
max = -1.0e10
min = 1.0e10

for i in range (totalLength-gap, totalLength):
    average = average + sumForce[i];
    if(sumForce[i] >= max):
        max = sumForce[i]
    if(sumForce[i] <= min):
        min = sumForce[i]

average = average/gap

## repeat for the imaginary part of pressure

averageI = 0.0
maxI = -1.0e10
minI = 1.0e10

for i in range (totalLength-gap, totalLength):
    averageI = averageI + sumForceI[i];
    if(sumForceI[i] >= maxI):
        maxI = sumForceI[i]
    if(sumForce[i] <= minI):
        minI = sumForceI[i]

averageI = averageI/gap


print ("calculate damping factor done!") 
#*************************************************************************#
# calculate damping factor


## computed from the code using the same method that computes A*p integral

## this is the value for mode 1
#AAdxdy = 2.283285e-15*2.


## this is the value for mode 2

AAdxdy = 1.693977e-13*2
#factor = 2*omega*rho*b*h*A*length
factor = 2*omega*rho*AAdxdy*h

factor -= averageI

print 'average %s' % average
print 'factor %s' % factor

zetaAvg = -average/factor

zetaMax = -max/factor

zetaMin = -min/factor

print zetaAvg, zetaMin, zetaMax

fileName = fileBase + "damping-ratio.dat"
file = open(fileName,"w")

file.write("pressure\t")
file.write("%f\n" % pressure)

file.write("frequency\t")
file.write("%f\n" % frequency)

file.write("beam density\t")
file.write("%f\n" % rho)

file.write("beam width\t")
file.write("%e\n" % b)

file.write("beam thickness\t")
file.write("%e\n" % h)

file.write("velocity magnitude\t")
file.write("%f\n" % A)

file.write("time step\t")
file.write("%e\n" % deltaT)

file.write("average damping ratio\t")
file.write("%e\n" % zetaAvg)

file.write("mininum damping ratio\t")
file.write("%e\n" % zetaMin)

file.write("maximum damping ratio\t")
file.write("%e\n" % zetaMax)

print ("results are written to damping-ratio.dat")
#*************************************************************************#
