#This source code is licensed under the MIT license found in the
#MIT-LICENSE file in the root directory of this source tree.
import sys
import os
import csv
import math
import numpy as np
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
import random
#count the slices of a processed cif file difZ=12 A
#---------------------containers---------------------------------------------
class scattData(object):
    """docstring for scattData,get all the scattering data from each element."""

    def __init__(self, name, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5):
        self.name = name
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.a5 = a5
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        self.b4 = b4
        self.b5 = b5
    def showName(self):
        return self.name
    def showA1(self):
        return self.a1
    def showA2(self):
        return self.a2
    def showA3(self):
        return self.a3
    def showA4(self):
        return self.a4
    def showA5(self):
        return self.a5
    def showB1(self):
        return self.b1
    def showB2(self):
        return self.b2
    def showB3(self):
        return self.b3
    def showB4(self):
        return self.b4
    def showB5(self):
        return self.b5

class utilData(object):
    """docstring fo utilData in Angstroms."""

    def __init__(self, argElement, argX, argY, argZ):
        self.argElement = argElement
        self.argX = argX
        self.argY = argY
        self.argZ = argZ
    def showElement(self):
        return self.argElement
    def showX(self):
        return self.argX
    def showY(self):
        return self.argY
    def showZ(self):
        return self.argZ
#---------------------Start slices Library-------------------------------
def convertDataToColour(the2dReal,maxRealValue,theCut2):
    baseImage1=np.zeros((nextOne,nextOne),dtype = int)
    i=0
    while i<theCut2:
        j=0
        while j<theCut2:
            imagedata1=int((the2dReal[i,j]*100)/maxRealValue)
            baseImage1[i,j]=imagedata1
            j=j+1
        i=i+1
    return baseImage1
def powerOfTwo(target):
    if target > 1:
        for i in range(1, int(target)):
            if (2 ** i >= target):
                return 2 ** i
    else:
        return 1
def preCounterSizes(pathOrigin):
    allCounterX=[]
    allCounterY=[]
    allCounterZ=[]
    with open(pathOrigin) as f:
        for line in f:
            interest2=line.split()
            if interest2[0]=="ATOM" or interest2[0]=="HETATM":
                markX = float(interest2[10])
                markY = float(interest2[11])
                markZ = float(interest2[12])
                allCounterX.append(markX)
                allCounterY.append(markY)
                allCounterZ.append(markZ)
    minX=abs(min(allCounterX))
    maxX=abs(max(allCounterX))
    minY=abs(min(allCounterY))
    maxY=abs(max(allCounterY))
    allCounterZ.sort()
    theMinZ=allCounterZ[0]
    theSecondMaxZ=allCounterZ[-2]
    theMaxZ=allCounterZ[-1]
    biggest=0
    if minX > maxX:
        biggestX=minX
    else:
        biggestX=maxX
    if minY > maxY:
        biggestY=minY
    else:
        biggestY=maxY
    if biggestX > biggestY:
        biggest=biggestX
    else:
        biggest=biggestY
    outerLimits=[biggest,theMinZ,theSecondMaxZ,theMaxZ]
    return outerLimits
def qfunc(theX):
    subterm=0.5-(0.5*math.erf(theX/math.sqrt(2)))
    return subterm
def findLambda(volts):
    h = 6.62606957e-34
    emass = 9.10938291e-31
    squarelight = 8.987551787e+16
    elementalCharge = 1.60217657e-19
    subterm = (elementalCharge * volts)/(2 * emass * squarelight)
    term1 = h / math.sqrt(2 * emass * elementalCharge * volts)
    term2 = 1 / math.sqrt(1 + subterm)
    newlambda = term1 * term2
    return newlambda
def divideZero(number):
    if number==0:
        number=1
    return number
def convertToAngstroms(coordinate):
    if coordinate == 0:
        coordinate = 1
    theNewCoordinate=float(coordinate)*1e-10
    return theNewCoordinate
def elementsBrightResource(subArray3):
    allElements=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
    'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',
    'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru',
    'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr',
    'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W',
    'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
    'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf']
    res = []
    for atoms in subArray3:
        allGoods=atoms.showElement()
        res.append(allGoods)
    res2=list(set(res))
    return res2
def extractElements(elementals,pathLibrary):
    resources=[]
    with open(pathLibrary) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                if row[1] in elementals:
                    resources.append(scattData(row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11]))
                line_count += 1
    return resources
def getOneScatt(literal,arrayResources):
    for obj in arrayResources:
        if literal == obj.name:
            return obj
#---------------------Start electronic Potential Library-----------------
def twoDimensionArrayPotential(anotherOneLambda,subArray,elementsInUseData,big2,secondZ,lastZ):
    chemR1 = np.zeros((big2,big2),dtype=complex)
    planck2=4.390477714e-67
    massElectron=9.10938188e-31
    precharge=1-(2*math.pi*massElectron*anotherOneLambda)/planck2
    particleBand = (-planck2)/(2*math.pi*massElectron)
    squarePi=math.pow(math.pi,2)
    lastestz = convertToAngstroms(lastZ)
    secondLastestZ = convertToAngstroms(secondZ)
    for barrier in subArray:
        elemental=barrier.showElement()
        sca = getOneScatt(elemental,elementsInUseData)
        lastX2 = int(float(barrier.showX()))
        lastY2 = int(float(barrier.showY()))
        oX=convertToAngstroms(barrier.showX())
        oY=convertToAngstroms(barrier.showY())
        oZ=convertToAngstroms(barrier.showZ())
        u2 = pow(oX,2)+pow(oY,2)
        firstStage1 = float(sca.a1)/float(sca.b1)
        firstStage2 = float(sca.a2)/float(sca.b2)
        firstStage3 = float(sca.a3)/float(sca.b3)
        firstStage4 = float(sca.a4)/float(sca.b4)
        firstStage5 = float(sca.a5)/float(sca.b5)
        secondStage1a = math.cos(-squarePi*u2)/float(sca.b1)
        secondStage1b = math.sin(-squarePi*u2)/float(sca.b1)
        secondStage2a = math.cos(-squarePi*u2)/float(sca.b2)
        secondStage2b = math.sin(-squarePi*u2)/float(sca.b2)
        secondStage3a = math.cos(-squarePi*u2)/float(sca.b3)
        secondStage3b = math.sin(-squarePi*u2)/float(sca.b3)
        secondStage4a = math.cos(-squarePi*u2)/float(sca.b4)
        secondStage4b = math.sin(-squarePi*u2)/float(sca.b4)
        secondStage5a = math.cos(-squarePi*u2)/float(sca.b5)
        secondStage5b = math.sin(-squarePi*u2)/float(sca.b5)
        thirdInner1 = math.sqrt((2*squarePi)/float(sca.b1))*(lastestz-oZ)
        thirdInner2 = math.sqrt((2*squarePi)/float(sca.b2))*(lastestz-oZ)
        thirdInner3 = math.sqrt((2*squarePi)/float(sca.b3))*(lastestz-oZ)
        thirdInner4 = math.sqrt((2*squarePi)/float(sca.b4))*(lastestz-oZ)
        thirdInner5 = math.sqrt((2*squarePi)/float(sca.b5))*(lastestz-oZ)
        thirdQ1 = qfunc(thirdInner1)
        thirdQ2 = qfunc(thirdInner2)
        thirdQ3 = qfunc(thirdInner3)
        thirdQ4 = qfunc(thirdInner4)
        thirdQ5 = qfunc(thirdInner5)
        fourthInner1 = math.sqrt((2*squarePi)/float(sca.b1))*(secondLastestZ-oZ)
        fourthInner2 = math.sqrt((2*squarePi)/float(sca.b2))*(secondLastestZ-oZ)
        fourthInner3 = math.sqrt((2*squarePi)/float(sca.b3))*(secondLastestZ-oZ)
        fourthInner4 = math.sqrt((2*squarePi)/float(sca.b4))*(secondLastestZ-oZ)
        fourthInner5 = math.sqrt((2*squarePi)/float(sca.b5))*(secondLastestZ-oZ)
        fourthQ1 = qfunc(fourthInner1)
        fourthQ2 = qfunc(fourthInner2)
        fourthQ3 = qfunc(fourthInner3)
        fourthQ4 = qfunc(fourthInner4)
        fourthQ5 = qfunc(fourthInner5)
        sixthStage1 = thirdQ1 - fourthQ1
        sixthStage2 = thirdQ2 - fourthQ2
        sixthStage3 = thirdQ3 - fourthQ3
        sixthStage4 = thirdQ4 - fourthQ4
        sixthStage5 = thirdQ5 - fourthQ5
        pseudo1a = firstStage1 * secondStage1a * sixthStage1
        pseudo1b = firstStage1 * secondStage1b * sixthStage1
        pseudo2a = firstStage2 * secondStage2a * sixthStage2
        pseudo2b = firstStage2 * secondStage2b * sixthStage2
        pseudo3a = firstStage3 * secondStage3a * sixthStage3
        pseudo3b = firstStage3 * secondStage3b * sixthStage3
        pseudo4a = firstStage4 * secondStage4a * sixthStage4
        pseudo4b = firstStage4 * secondStage4b * sixthStage4
        pseudo5a = firstStage5 * secondStage5a * sixthStage5
        pseudo5b = firstStage5 * secondStage5b * sixthStage5
        pseudoPosA = pseudo1a + pseudo2a + pseudo3a + pseudo4a + pseudo5a
        pseudoPosB = pseudo1b + pseudo2b + pseudo3b + pseudo4b + pseudo5b
        electronicPotentialA = precharge * particleBand * pseudoPosA
        electronicPotentialB = precharge * particleBand * pseudoPosB
        electronicPotential = complex(electronicPotentialA,electronicPotentialB)
        chemR1[lastX2,lastY2]=electronicPotential
    return chemR1
#---------------------End electronic Potential Library-------------------
#---------------------Start Propagator Library---------------------------
def propagator(anotherLambda3,subArray2,big):
    preDone = np.zeros((big,big),dtype=complex)
    for wave in subArray2:
        lastX = int(float(wave.showX()))
        lastY = int(float(wave.showY()))
        oX=convertToAngstroms(wave.showX())
        oY=convertToAngstroms(wave.showY())
        oZ=convertToAngstroms(wave.showZ())
        u2 = (pow(oX,2)+pow(oY,2))/divideZero(oZ)
        divPrecharge=divideZero(oZ*anotherLambda3)
        precharge2 = 2*math.pi*oZ
        precharge2a = (-math.cos(precharge2))/(divPrecharge)
        precharge2b = (-math.sin(precharge2))/(divPrecharge)
        precharge3 = (math.pi*u2)/(divPrecharge)
        precharge4a = math.cos(precharge3)
        precharge4b = math.sin(precharge3)
        charge5a=precharge2a+precharge4a
        charge5b=precharge2b+precharge4b
        charge5=complex(charge5a,charge5b)
        preDone[lastX,lastY]=charge5
    return preDone
#---------------------End Propagator Library-----------------------------
def oneLayer(chemicalPotential,eTransmission,psiSignal,side):
    moc1=np.matmul(eTransmission,psiSignal)
    moc2=np.fft.fft2(moc1)
    moc3=np.fft.fftshift(moc2)
    moc4=np.matmul(moc3,chemicalPotential)
    moc5=np.fft.ifft2(moc4)
    return moc5
def generalExtractor(pathOrigin):
    all=[]
    with open(pathOrigin) as f:
        for line in f:
            interest2=line.split()
            if interest2[0]=="ATOM" or interest2[0]=="HETATM":
                elementData=interest2[2]
                oX = interest2[10]
                oY = interest2[11]
                oZ = interest2[12]
                all.append(utilData(elementData,oX,oY,oZ))
    return all
def generalLayerOneIterator(psiOrigin,lambBase,cifArray,elementsInUse4,scatterWindow,stageOneZ,stageTwoZ,secondMaximumZ,maximumZ):
    temporalLayer=[]
    for sliced in cifArray:
        term4=float(sliced.showZ())
        if stageOneZ <= term4 <= stageTwoZ:
            temporalLayer.append(sliced)
    #make layer chemicalPotential
    onePotential=twoDimensionArrayPotential(lambBase,temporalLayer,elementsInUse4,scatterWindow,secondMaximumZ,maximumZ)
    #make layer propagator
    onePropagator=propagator(lambBase,temporalLayer,scatterWindow)
    psiSlice=oneLayer(onePotential,onePropagator,psiOrigin,scatterWindow)
    return psiSlice
############################START#LAYERS###################################
pathTest=sys.argv[1]
pathLibrary=sys.argv[2]
voltageApplied=int(sys.argv[3])
theSizes=preCounterSizes(pathTest)
nextOne=powerOfTwo(theSizes[0])
blackSheep=findLambda(voltageApplied)
allText=generalExtractor(pathTest)
elementsInUse = elementsBrightResource(allText)
elementsInUseData = extractElements(elementsInUse,pathLibrary)
#start of the general generalLayerIterator
psi=np.identity(nextOne,dtype = complex)
theStart=float(theSizes[1])
secondBest=float(theSizes[2])
theEnd=float(theSizes[3])
while theStart <=theEnd:
    jump=theStart+12
    psi=generalLayerOneIterator(psi,blackSheep,allText,elementsInUseData,nextOne,theStart,jump,secondBest,theEnd)
    theStart=theStart+12
preImage1=np.fft.fft2(psi)
preImage2=abs(preImage1)
maxV=np.max(preImage2)
preImage3=convertDataToColour(preImage2,maxV,nextOne)
plt.imshow(preImage3)
plt.show()
