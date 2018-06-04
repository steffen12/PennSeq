#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict
import math, sys, os, re, pysam, time, copy
import numpy as np

###############################################################################
###  ARGUMENT SETTINGS
###############################################################################

# checking whether argument is valid or not
validArgList = ["-counts", "-readlen", "-precision", "-out"]
for argIndex in range(1,len(sys.argv)):
    if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList :
        print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
        sys.exit()
        
countFileExists = False
readLengthExists = False
precisionExists = False
outFileExists = False
argIndex = 1

while argIndex < len(sys.argv):
    if sys.argv[argIndex] == "-counts":  ## load in counts file
        argIndex += 1
        countFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        countTmp = sys.argv[argIndex].split("/")
        countFile = countFileAbsPath + "/" + countTmp[len(countTmp)-1]
        countFileExists = True
        argIndex += 1
    elif sys.argv[argIndex] == "-readlen":  ## load in read len
        argIndex += 1
        readLength = int(sys.argv[argIndex])
        readLengthExists = True
        argIndex += 1
    elif sys.argv[argIndex] == "-precision":  ## load in precision
        argIndex += 1
        diffMax = float(sys.argv[argIndex])
        precisionExists = True
        argIndex += 1
    elif sys.argv[argIndex] == "-out":  ## load in output file
        argIndex += 1
        outFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        outTmp = sys.argv[argIndex].split("/")
        outFile = outFileAbsPath + "/" + outTmp[len(outTmp)-1]
        outFileExists = True
        argIndex += 1                  

if (not countFileExists) or (not readLengthExists) or (not precisionExists) or (not outFileExists): ## lack enough arguments
    print("Please provide arguments:")
    print("-counts\tSegment count file")
    print("-readlen\tRead Length")
    print("-precision\tPrecision")
    print("-out\tOutput file")
    sys.exit()

countFileHandle = open(countFile, 'r')
OUT = open(outFile, 'w')


###########################################################################################################################
###  START TO ANALYZE DATA FOR EACH GENE ###
##########################################################################################################################

geneCount = 0

startTime = time.time()

countFileLine = countFileHandle.readline() #Analyze header
splitLine = countFileLine.strip().split("\t")

if len(splitLine) == 7:
    pairedEndMode = False
elif len(splitLine) == 10:
    pairedEndMode = True
else:
    print("Invalid Counts Input File")
    sys.exit()

OUT.write("GeneName\tIsoformName\tNumberOfReads\tRelativeAbundance\n") ## Header of Results

countFileLine = countFileHandle.readline()

while countFileLine != "":
    splitLine = countFileLine.strip().split("\t")
    currGene = splitLine[3] #Get Current Gene
    print(currGene)

    segmentIDs = []
    segmentCounts = []
    segmentLengths = []
    segmentIsoforms = dict()
    geneIsoforms = []

    while countFileLine != "" and splitLine[3] == currGene:
        segmentID = splitLine[0]
        segmentIDs.append(segmentID)
        segmentCounts.append(int(splitLine[1]))
        segmentLengths.append(int(splitLine[4]))

        segmentIsoforms[segmentID] = splitLine[6].split(",")
        for isoform in segmentIsoforms[segmentID]:
            if isoform not in geneIsoforms:
                geneIsoforms.append(isoform)

        countFileLine = countFileHandle.readline()
        splitLine = countFileLine.strip().split("\t")

    readCount = sum(segmentCounts)

    if readCount == 0:
        continue

    segmentIsoformIndicatorMatrix = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)), dtype="bool_")
    for i in range(len(segmentIDs)):
        segmentID = segmentIDs[i]
        segmentIsoformList = segmentIsoforms[segmentID]
        for isoform in segmentIsoformList:
            j = geneIsoforms.index(isoform)
            segmentIsoformIndicatorMatrix[i, j] = True

    effectiveSegmentLengths = np.zeros(shape=(len(segmentIDs)), dtype="int_")
    for i in range(len(segmentIDs)):
        effectiveSegmentLengths[i] = segmentLengths[i] - readLength + 1

    isoformCounts = np.zeros(shape=(len(geneIsoforms)), dtype="int_")
    isoformLengths = np.zeros(shape=(len(geneIsoforms)), dtype="int_")
    for j in range(len(geneIsoforms)):
        for i in range(len(segmentIDs)):
            if(segmentIsoformIndicatorMatrix[i,j]):
                isoformCounts[j] += segmentCounts[i]
                isoformLengths[j] += effectiveSegmentLengths[i]

    geneCount += 1
    
    ############################################################################################################################################
    ## Find H for each segment
    ############################################################################################################################################
    
    segmentH = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)))
    for i in range(len(segmentIDs)):
        for j in range(len(geneIsoforms)):
            if(segmentIsoformIndicatorMatrix[i,j]):
                segmentH[i, j] = float(segmentCounts[i])/(isoformCounts[j] * effectiveSegmentLengths[i])

    print("Time to Analyze Distribution", time.time() - startTime)

#     #####################################################################################################
#     ## EM algorithm
#     #####################################################################################################

    Thetas = np.full(shape=len(geneIsoforms), fill_value=1.0/len(geneIsoforms))
    Alpha = np.zeros(shape=len(geneIsoforms)) #This is Theta with ~ over it
    for j in range(len(geneIsoforms)):
        Alpha[j] = Thetas[j] * isoformLengths[j]
    sumAlpha = np.sum(Alpha)
    for j in range(len(geneIsoforms)):
        Alpha[j] /= sumAlpha
    oldAlpha = np.zeros(shape=len(geneIsoforms))

     #########################################################################################################
     ## iteration begins        
            
    diff = 1.0
    iterCount = 0
    Z = np.zeros(shape=(len(segmentIDs), len(geneIsoforms)))
    while diff > diffMax:
        #print(gene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+str(diff)+"\t"+str(tmpTime))

        #Expectation Step
        for i in range(len(segmentIDs)):
            probNumerators = np.zeros(shape=len(geneIsoforms))
            for j in range(len(geneIsoforms)):
                if segmentIsoformIndicatorMatrix[i, j]: 
                    probNumerators[j] = Alpha[j] * segmentH[i, j]
            probDenominator = sum(probNumerators)    
            for j in range(len(geneIsoforms)):
                Z[i, j] = segmentCounts[i] * probNumerators[j] / probDenominator

        #Maximization Step
        oldAlpha = np.copy(Alpha)
        Alpha = np.multiply(np.sum(Z, axis=0), 1/float(readCount))
        
        iterCount += 1

        diffArray = np.absolute(np.subtract(Alpha, oldAlpha))
        diff = diffArray.max()
                
    sumAlpha = np.sum(Alpha)
    if sumAlpha == 0: 
        continue

    Alpha = np.multiply(Alpha, 1.0 / sumAlpha)

    Thetas = np.zeros(shape=(len(geneIsoforms),))
    for j in range(len(Alpha)):
        Thetas[j] = Alpha[j] / isoformLengths[j]
    sumTheta = np.sum(Thetas)

    print(currGene+"\t"+str(iterCount)+" iterations\tDone!")

    
    for i in range(len(Alpha)):
        Thetas[i] /= sumTheta
        #print(currGene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+isoformNames[i]+"\t"+str(isoformRelativeAbundances[i])+"\t"+str(tmpTime))
        OUT.write(currGene+"\t"+geneIsoforms[i]+"\t"+str(isoformCounts[i])+"\t"+str(Thetas[i])+"\n") ## write results into specified file

    print("Time for EM", time.time() - startTime)

print(time.time() - startTime)
OUT.close()