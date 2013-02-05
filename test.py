import sys
import os
import re
import optparse
import scipy as S
import CijUtil
import espresso

cijdat = open("MgO.cijdat","r")

numStrainPatterns = (len(cijdat.readlines())-2)/4

finalCijs = S.zeros((21,1))
errors = S.zeros((21,1))

cijdat.seek(0)
latticeType,numsteps,TetrHigh,TrigHigh = cijdat.readline().split()
numsteps = int(numsteps)

latticeTypes = {0:"Unknown", 1:"Triclinic", 2:"Monoclinic", 3:"Orthorhombic", \
		4:"Tetragonal", 5:"Cubic", 6:"Trigonal-low", 7:"Trigonal-high/Hexagonal"}

symmetryType = latticeTypes[int(latticeType)]
magnitude = float(cijdat.readline())

for patt in range(1):
     for a in range(0,numsteps):
             pattern =cijdat.readline()
             line1 = cijdat.readline().split()
             line2 = cijdat.readline().split()
             line3 = cijdat.readline().split()
             if a == 0:
                     strain = S.array([float(line1[0]),float(line2[1]),float(line3[2]),2*float(line2[2]),2*float(line1[2]),2*float(line1[1])])
             else:
                     strain = S.row_stack(S.array((strain,[float(line1[0]),float(line2[1]),float(line3[2]),2*float(line2[2]),2*float(line1[2]),2*float(line1[1])])))
             (units, thisStress) = espresso.get_stress("MgO_cij__"+str(patt+1)+"__"+str(a+1))
             if a == 0:
                     stress = thisStress
             else:
                     stress = S.row_stack((stress,thisStress))

def __fit(index1,index2):
       from scipy import stats, sqrt, square
       print strain
       print stress
       (cijFitted,intercept,r,tt,stderr) = stats.linregress(strain[:,index2-1],stress[:,index1-1])
       if (S.__version__ < '0.7.0'):
           stderr = S.sqrt((numsteps * stderr**2)/(numsteps-2))
           error  = stderr/sqrt(sum(square(strain[:,index2-1])))
       else:
           fit_str = ((strain[index2-1,:] * cijFitted) + intercept)                
           error = sqrt((sum(square(stress[:,index1-1] - fit_str)) / \
                       (numsteps-2))/(sum(square(strain[:,index2-1]))))
       print 'Cij   ', cijFitted
       print 'Error   ', error
       print 'intercept   ', intercept
       return cijFitted, error

cij = S.zeros(21)
strainsUsed = S.zeros((6,1))
strainsUsed[0] = 1
strainsUsed[1]=0
strainsUsed[2]=0
strainsUsed[3]=1
strainsUsed[4]=0
strainsUsed[5]=0
finalCijs[0], errors[0] = __fit(1,1)
fit_21, fit_21_error = __fit(2,1)
fit_31, fit_31_error = __fit(3,1)
finalCijs[6] = (fit_21 + fit_31)/2
errors[6] = S.sqrt((fit_21_error**2)/4 + (fit_31_error**2)/4)
finalCijs[3], errors[3] = __fit(4,4)
c = S.matrix([[1, 7, 7, 0, 0, 0],
		[7, 1, 7, 0, 0, 0],
		[7, 7, 1, 0, 0, 0],
		[0, 0, 0, 4, 0, 0],
		[0, 0, 0, 0, 4, 0],
		[0, 0, 0, 0, 0, 4]])
finalCijMatrix = S.zeros((6,6))
for i in range(0,6):
     for j in range(0,6):
           index = int(c[i,j])
           if index > 0:
                 finalCijMatrix[i,j] = finalCijs[index-1]
           elif index < 0:
                 finalCijMatrix[i,j] = -finalCijs[-index-1]

print finalCijMatrix
