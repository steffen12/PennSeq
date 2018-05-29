#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict
import math, sys, os, re, pysam, time, copy

# set up auto dictionary function, makes a dictionary with dictionary elements
def auto_dict():
    return defaultdict(auto_dict)

# load CIGAR information
def getCigarStringInformation(readCigar, readName, readNumber):
    cigarMatchRead = auto_dict() #Store CIGAR string characters
    cigarNumberRead = auto_dict() #Store CIGAR string numbers
    cigarMatchInfoCount = 0 #Number of elements in cigarMatchRead
    cigarNumberInfoCount = 0 #Number

    splitCigar = re.split("([A-Z])",readCigar[readName][readNumber]) #Split CIGAR String 
    for i in range(len(splitCigar)-1):
        if splitCigar[i].isalpha(): #If all characters are alphabetic
            cigarMatchRead[cigarMatchInfoCount] = splitCigar[i]
            cigarMatchInfoCount += 1
        else:
            cigarNumberRead[cigarMatchInfoCount] = int(splitCigar[i])
            cigarNumberInfoCount += 1

    return cigarMatchRead, cigarNumberRead, cigarMatchInfoCount, cigarNumberInfoCount

def getExonsReadMappedTo(cigarMatchInfoCount, cigarMatchRead, cigarNumberRead, exonStarts, base, numofExons, isoformNames, compatibleVector, exonIndicators):
    exonIndicatorRead = [0] * numofExons

    firstExonRead = -1
    lastExonRead = -1
    for i in range(cigarMatchInfoCount):
        
        if cigarMatchRead[i] == "M" or cigarMatchRead[i] == "I": ## matched CIGAR

            for j in range(1,cigarNumberRead[i]+1):
                tmpbase = base + j
                for k in range(len(exonStarts)):
                    if exonIndicatorRead[k] == 1: 
                        continue #Go to next loop iteration
                    if tmpbase >= exonStarts[k] and tmpbase <= exonEnds[k]:
                        if firstExonRead == -1:
                            firstExonRead = k 
                        lastExonRead = k
                        exonIndicatorRead[k] = 1 ## confirm that the read covers this exon

            base += cigarNumberRead[i] # jump to next match information

        if cigarMatchRead[i] == "N": ## skipping area
            base += cigarNumberRead[i] # jump to next match information directly

    for j in range(len(isoformNames)):
        for exonIndex in range(firstExonRead, lastExonRead+1):
            #print(exonIndicatorRead1[exonIndex], exonIndicators[isoformNames[j]][exonIndex])
            if exonIndicatorRead[exonIndex] != exonIndicators[isoformNames[j]][exonIndex]:
                compatibleVector[j] = 0
    tmpcount = sum(exonIndicatorRead)
    #print("Vectors")
    #print(compatibleVector)
    #print(compatibleVector2)

    return compatibleVector, tmpcount, firstExonRead, lastExonRead

# def getExonsReadMappedTo(cigarMatchInfoCount, cigarMatchRead, cigarNumberRead, exonStarts, base, numofExons, isoformNames, compatibleVector, exonIndicators):
#     exonIndicatorRead = [0] * numofExons

#     firstExonRead = -1
#     lastExonRead = -1
#     currentExon = 0
#     for i in range(cigarMatchInfoCount):
        
#         if cigarMatchRead[i] == "M" or cigarMatchRead[i] == "I": ## matched CIGAR

#             for j in range(1,cigarNumberRead[i]+1):
#                 tmpbase = base + j
#                 if tmpbase >= exonStarts[currentExon] and tmpbase <= exonEnds[currentExon]:
#                     if firstExonRead == -1:
#                         firstExonRead = currentExon 
#                     lastExonRead = currentExon
#                     exonIndicatorRead[currentExon] = 1 ## confirm that the read covers this exon
#                 while tmpbase > exonEnds[currentExon]:
#                     currentExon += 1
#             base += cigarNumberRead[i] # jump to next match information

#         if cigarMatchRead[i] == "N": ## skipping area
#             base += cigarNumberRead[i] # jump to next match information directly

#     for j in range(len(isoformNames)):
#         for exonIndex in range(firstExonRead, lastExonRead+1):
#             #print(exonIndicatorRead1[exonIndex], exonIndicators[isoformNames[j]][exonIndex])
#             if exonIndicatorRead[exonIndex] != exonIndicators[isoformNames[j]][exonIndex]:
#                 compatibleVector[j] = 0
#     tmpcount = sum(exonIndicatorRead)
#     #print("Vectors")
#     #print(compatibleVector)
#     #print(compatibleVector2)

#     return compatibleVector, tmpcount, firstExonRead, lastExonRead

###############################################################################
###  ARGUMENT SETTINGS
###############################################################################

# checking whether argument is valid or not
validArgList = ["-bam", "-ref", "-out"]
for argIndex in range(1,len(sys.argv)):
    if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList :
        print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
        sys.exit()
        

bamFileExists = False
refFileExists = False
outFileExists = False
argIndex = 1
while argIndex < len(sys.argv):
    if sys.argv[argIndex] == "-bam":  ## load in BAM file
        argIndex += 1
        bamFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        bamTmp = sys.argv[argIndex].split("/")
        bamFile = bamFileAbsPath + "/" + bamTmp[len(bamTmp)-1]
        bamFileExists = True
        argIndex += 1
    elif sys.argv[argIndex] == "-ref":  ## load in annotation file
        argIndex += 1
        refFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        refTmp = sys.argv[argIndex].split("/")
        refGeneFile = refFileAbsPath + "/" + refTmp[len(refTmp)-1]
        refFileExists = True
        argIndex += 1
    elif sys.argv[argIndex] == "-out":  ## load in output file
        argIndex += 1
        outFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        outTmp = sys.argv[argIndex].split("/")
        outFile = outFileAbsPath + "/" + outTmp[len(outTmp)-1]
        outFileExists = True
        argIndex += 1
        
                                                

if (not bamFileExists) or (not refFileExists) or (not outFileExists): ## lack enough arguments
    print("Please provide arguments:")
    print("-bam\tIndexed bam file")
    print("-ref\tGene annotation file")
    print("-out\tOutput file")
    sys.exit()


# load gene information
geneStructureInformation = auto_dict() #will make a dictionary object
geneLineCount = auto_dict()

with open(refGeneFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        tmpinf = line.split("\t")
        gene = tmpinf[0]
        
        if not bool(geneStructureInformation[gene]): #if first line of matrix
            geneLineCount[gene] = 0
            geneStructureInformation[gene][geneLineCount[gene]] = line
        else:
            geneLineCount[gene] += 1
            geneStructureInformation[gene][geneLineCount[gene]] = line


#####################################
## Using pysam to read in bam file !!
#####################################
bamFilePysam = pysam.Samfile(bamFile,"rb")


## RESULTS FILE
OUT = open(outFile, 'w')


###########################################################################################################################
###  START TO ANALYZE DATA FOR EACH GENE ###
##########################################################################################################################

geneCount = 0

startTime = time.time()

OUT.write("GeneName\tIsoformName\tNumberOfReads\tRelativeAbundance\n") ## Header of Results

for gene in geneStructureInformation: #This can be made parallel easily

    geneCount += 1
    tmpTime = (time.time() - startTime)/60.0
    
    
    sameReadCount = auto_dict() #Keep track of paired-end reads
    isPairedEnd = dict()
    readStart = auto_dict()
    readEnd = auto_dict()
    readCigar = auto_dict()

    numofExons = geneLineCount[gene]
    tmpgeneinf = geneStructureInformation[gene][0].split("\t") #Get first line of matrix
    geneChr = tmpgeneinf[1]
    geneStart = int(tmpgeneinf[3])
    geneEnd = int(tmpgeneinf[4])


    ## load all reads information which were mapped to the specific gene within this loop using pysam
    for read in bamFilePysam.fetch(geneChr, geneStart, geneEnd):
        line = str(read)
        tmpinf = line.split("\t")
        tmpReadName = tmpinf[0]
        tmpReadChr = geneChr
        tmpReadStart = int(tmpinf[3]) + 1 #Explain + 1 here
        tmpReadCigar = ""

        ## Adjust to different Pysam Version to get CIGAR string!! ##

        if ")]" in tmpinf[5]: ## vector format
            
            tmpinf[5] = tmpinf[5].rstrip(")]")
            tmpinf[5] = tmpinf[5].lstrip("[(")
            tmpinfcigar = tmpinf[5].split("), (")
            for cc in tmpinfcigar:
                ttcc = cc.split(", ")
                if ttcc[0] == "3":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "N"
                if ttcc[0] == "2":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "D"
                if ttcc[0] == "1":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "I"
                if ttcc[0] == "0":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "M"
                if not (ttcc[0] == "3" or ttcc[0] == "2" or ttcc[0] == "1" or ttcc[0] == "0"):
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "X"
        else:      ## 100M10N100M format
            tmpReadCigar = tmpinf[5]
                                    
        if not bool(sameReadCount[tmpReadName]): #if empty
            sameReadCount[tmpReadName] = 1
            isPairedEnd[tmpReadName] = False
        else:
            sameReadCount[tmpReadName] += 1
            isPairedEnd[tmpReadName] = True
                                        
        readStart[tmpReadName][sameReadCount[tmpReadName]] = tmpReadStart
        readCigar[tmpReadName][sameReadCount[tmpReadName]] = tmpReadCigar


    ## load structure information of the specific gene within this loop                    
                        
    tmpgeneinf[5] = tmpgeneinf[5].rstrip(",") #Column 5 contains transcript names
    isoformNames = tmpgeneinf[5].split(",")
    exonStarts = [None] * numofExons
    exonEnds = [None] * numofExons
    exonIndicators = auto_dict()
    
    #Iterate for each exon in gene
    for exonIndex in range(numofExons):
        tmpinf = geneStructureInformation[gene][exonIndex+1].split("\t") #+1 as first exon is in 2nd row
        exonStarts[exonIndex] = int(tmpinf[3])+1 #Explain +1
        exonEnds[exonIndex] = int(tmpinf[4])
        tmpinf[5] = tmpinf[5].rstrip(",") #Column 5 contains matrix values
        tmpExonIndicators = tmpinf[5].split(",")

        #For each isoform, make an array containing 1 if isoform has exon i, or 0 if it does not
        #Get column of matrix and store it
        for isoformIndex in range(len(isoformNames)):
            isoformName = isoformNames[isoformIndex]
            exonIndicators[isoformName][exonIndex] = int(tmpExonIndicators[isoformIndex])

    print("Time to Preprocess reads", time.time() - startTime)
    #########################################################################################################################################
    ## START TO ANALYZE EACH READ 
    ##################################################################################################################################################

    qualifiedRead = auto_dict()
    readCount = 0
    fragmentStart = auto_dict()
    fragmentEnd = auto_dict()
    fragmentFirstExon = dict()
    fragmentLastExon = dict()
    CompatibleMatrix = auto_dict()
    
    for readName in sameReadCount: #Iterate over all reads
        
        #Get CIGAR information
        cigarMatchRead1, cigarNumberRead1, cigarMatchInfoCount1, cigarNumberInfoCount1 = getCigarStringInformation(readCigar, readName, 1)
                                    
        #calculate read end positions by starting at readStart and adding "length" of CIGAR string,
        #which is the sum of all the numbers contained in it
        readEnd[readName][1] = readStart[readName][1]
        for i in range(cigarMatchInfoCount1):
            readEnd[readName][1] += cigarNumberRead1[i]
            
        #If paired end do for second
        if isPairedEnd[readName]:
            cigarMatchRead2, cigarNumberRead2, cigarMatchInfoCount2, cigarNumberInfoCount2 = getCigarStringInformation(readCigar, readName, 2)
            readEnd[readName][2] = readStart[readName][2]
            for i in range(cigarMatchInfoCount2):
                readEnd[readName][2] += cigarNumberRead2[i]

        # calculate fragment START and END positions
        if isPairedEnd[readName]: #paired end reads
            fragmentStart[readName] = min(readStart[readName][2], readStart[readName][1]) 
            fragmentEnd[readName] = max(readEnd[readName][1], readEnd[readName][2])
        else: #single end reads
            fragmentStart[readName] = readStart[readName][1]
            fragmentEnd[readName] = readEnd[readName][1]
            
        ##################################################################################################################################    
        ## Obtain compatible matrix of isoforms with respect to reads
        #################################################################################################################################
        
        #Check if read start is contained within gene, and for paired end reads check the pair
        if (readStart[readName][1] >= geneStart and readStart[readName][1] <= geneEnd) or (readStart[readName][2] >= geneStart and readStart[readName][2] <= geneEnd and isPairedEnd[readName]) :
            #Make sure that the CIGAR string was well formatted and for both reads (if there are two) this is true
            if cigarMatchInfoCount1 == cigarNumberInfoCount1 and ((not isPairedEnd[readName]) or cigarMatchInfoCount2 == cigarNumberInfoCount2):
                base1 = readStart[readName][1] - 1
                if isPairedEnd[readName]:
                    base2 = readStart[readName][2] - 1
                    exonIndicatorRead2 = [0] * numofExons
                compatibleVector = [1] * len(isoformNames)

                ##############################################################################################################################################
                ### SET UP COMPATIBLE INDICATOR VECTOR ###############
                ###############################################################################################################################################
                ## READ 1 ##
                # find exons where read 1 mapped to
                compatibleVector, tmpcount1, firstExonRead1, lastExonRead1 = getExonsReadMappedTo(cigarMatchInfoCount1, cigarMatchRead1, cigarNumberRead1, exonStarts, base1, numofExons, isoformNames, compatibleVector, exonIndicators)

                ## READ 2 ## SAME AS READ 1
                if isPairedEnd[readName]: 
                    compatibleVector, tmpcount2, firstExonRead2, lastExonRead2  = getExonsReadMappedTo(cigarMatchInfoCount2, cigarMatchRead2, cigarNumberRead2, exonStarts, base2, numofExons, isoformNames, compatibleVector, exonIndicators)
                    fragmentFirstExon[readName] = min(firstExonRead1, firstExonRead2)
                    fragmentLastExon[readName] = max(lastExonRead1, lastExonRead2)
                else:
                    fragmentFirstExon[readName] = firstExonRead1
                    fragmentLastExon[readName] = lastExonRead1

                ##################################################################################################################################################
                ## fill in compatible matrix ##
                if tmpcount1 > 0 or (isPairedEnd[readName] and tmpcount2 > 0):
                    readCount += 1
                    qualifiedRead[readName] = 1
                    for i in range(len(isoformNames)):
                        CompatibleMatrix[readName][isoformNames[i]] = compatibleVector[i]
                
    tmpCompatibleMatrix = copy.deepcopy(CompatibleMatrix)

    print("Time to Analyze reads", time.time() - startTime)

    ### COMPATIBLE MATRIX OBTAINED !!!
    ###############################################################################################################
    
    if readCount == 0: continue
    print(gene+"\t"+str(readCount)+" reads detected...")
    
    ##############################################################################################################
    ### ANALYZE EMPIRICAL READS DISTRIBUTION BASED ON NON-PARAMATRIC METHOD
    ##############################################################################################################

    positionsFragmentCovered = auto_dict()
    readsDistributionIsoform = auto_dict()
    readsDistributionIsoformKnown = auto_dict()
    denominatorKnown = auto_dict()
    denominator = auto_dict()
    
    for i in range(len(isoformNames)):
        for readName in qualifiedRead:
            ####################################################################################################
            ### CACULATE VALID FRAGMENT LOCATION ON THIS GENE ###################################
            firstExon = fragmentFirstExon[readName]
            lastExon = fragmentLastExon[readName]
            tmpStart = max(fragmentStart[readName], exonStarts[firstExon])
            tmpEnd = min(fragmentEnd[readName], exonEnds[lastExon])
            fragmentStart[readName] = tmpStart ## new starts and end position updated
            fragmentEnd[readName] = tmpEnd

        ############################################################################################################################################
        ## PILE UP AT EACH BASE PAIR POSITION
        ############################################################################################################################################
        

    for isoformIndex in range(len(isoformNames)):
        isoformName = isoformNames[isoformIndex]
        isoformOverlappingReadStarts = []
        isoformOverlappingReadEnds = []
        currentBaseDensity = 0
        denominator[isoformName] = 0
        startIndex = 0
        endIndex = 0
        numIsoformReads = 0
        for readName in qualifiedRead:
            if CompatibleMatrix[readName][isoformName] == 1: #If read could have come from isoform
                isoformOverlappingReadStarts.append(fragmentStart[readName])
                isoformOverlappingReadEnds.append(fragmentEnd[readName])
                numIsoformReads += 1
        isoformOverlappingReadStarts = sorted(isoformOverlappingReadStarts)
        isoformOverlappingReadEnds = sorted(isoformOverlappingReadEnds)
        for exonIndex in range(numofExons):
            if exonIndicators[isoformName][exonIndex] == 1: #If isoform has exon j
                #readExonStart = max(tmpStart, exonStarts[j])
                #readExonEnd = min(tmpEnd, exonEnds[j])
                for k in range(exonStarts[exonIndex], exonEnds[exonIndex]+1):
                    readsDistributionIsoform[isoformName][k] = 0
                    while(startIndex < numIsoformReads and isoformOverlappingReadStarts[startIndex] == k):
                        currentBaseDensity += 1
                        startIndex += 1
                    readsDistributionIsoform[isoformName][k] = currentBaseDensity
                    denominator[isoformName] += currentBaseDensity

                    while(endIndex < numIsoformReads and isoformOverlappingReadEnds[endIndex] == k):
                        currentBaseDensity -= 1
                        endIndex += 1

    print("Time to Analyze Distribution", time.time() - startTime)

    #####################################################################################################
    ## EM algorithm
    #####################################################################################################

    Alpha = [1.0] * len(isoformNames) #This is Theta with ~ over it
    oldAlpha = [None] * len(isoformNames)

    isoformLength = auto_dict()
    
    for i in range(len(isoformNames)):
        isoformLength[isoformNames[i]] = 0
        for j in range(len(exonStarts)):
            if exonIndicators[isoformNames[i]][j] == 1:
                isoformLength[isoformNames[i]] += exonEnds[j] - exonStarts[j] + 1  ## calculate isoform length

    print(isoformLength[isoformNames[0]])
    ########################################################################################################            
    ## UPDATE empirical probability of read mapped to the specific position

    Hfunction = auto_dict() #
    distributionRanges = auto_dict()
    numerator = auto_dict()
    rangeWidth = 10
    for i in range(len(isoformNames)):
        for exonIndex in range(numofExons):
            if exonIndicators[isoformNames[i]][exonIndex] == 1: #If isoform has exon j
                j = exonStarts[exonIndex]
                jStop = exonEnds[exonIndex] + 1
                while(j < jStop):
                    rangeTotal = 0
                    endJ = min(j + rangeWidth, exonEnds[exonIndex] + 1)
                    for k in range(j, endJ):
                        rangeTotal += readsDistributionIsoform[isoformNames[i]][k]
                    distributionRanges[isoformNames[i]][j] = rangeTotal
                    j += rangeWidth

    for i in range(len(isoformNames)):
        if not bool(denominator[isoformNames[i]]): 
            denominator[isoformNames[i]] = 0
        for readName in qualifiedRead:
            #print(isoformNames[i], readName)
            if CompatibleMatrix[readName][isoformNames[i]] == 1: #If read could have come from isoform
                numerator[readName][isoformNames[i]] = 0.0
                k = fragmentStart[readName]
                readLength = 0
                exonIndex = fragmentFirstExon[readName]
                firstExonStart = exonStarts[exonIndex]
                lastExonEnd = exonEnds[fragmentLastExon[readName]]
                while((k - firstExonStart) % rangeWidth != 0 and k <= exonEnds[exonIndex]):
                    #print(isoformNames[i], readName, numerator[readName][isoformNames[i]])
                    if bool(readsDistributionIsoform[isoformNames[i]][k]):
                        numerator[readName][isoformNames[i]] += float(readsDistributionIsoform[isoformNames[i]][k])
                    k += 1
                while k + rangeWidth<= fragmentEnd[readName]:
                    #print(isoformNames[i], readName, numerator[readName][isoformNames[i]])
                    if k > exonEnds[exonIndex]:
                        exonIndex += 1
                        k = exonStarts[exonIndex]
                    #if lociIndicators[isoformNames[i]][k] == 1 and bool(readsDistributionIsoform[isoformNames[i]][k]):
                    if bool(distributionRanges[isoformNames[i]][k]):
                        numerator[readName][isoformNames[i]] += float(distributionRanges[isoformNames[i]][k])
                    k += rangeWidth
                    readLength += 1
                while k <= fragmentEnd[readName]:
                    #print(isoformNames[i], readName, numerator[readName][isoformNames[i]])
                    #print(readsDistributionIsoform[isoformNames[i]][k+1]) 
                    #print(k)
                    #print(fragmentEnd[readName])
                    if bool(readsDistributionIsoform[isoformNames[i]][k]):
                        numerator[readName][isoformNames[i]] += float(readsDistributionIsoform[isoformNames[i]][k])
                    k += 1
                #print(isoformNames[i], readName, numerator[readName][isoformNames[i]])
                if denominator[isoformNames[i]] > 0 and numerator[readName][isoformNames[i]] >= 0: 
                    
                    Hfunction[readName][isoformNames[i]] = float(numerator[readName][isoformNames[i]]) / denominator[isoformNames[i]]
                else: 
                    Hfunction[readName][isoformNames[i]] = 0

    print("Time to Precompute", time.time() - startTime)
    #########################################################################################################
    ## iteration begins        
            
    diff = 1.0
    iterCount = 0
    d = 300
    diffMax = .0001
    while diff > diffMax:
        #print(gene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+str(diff)+"\t"+str(tmpTime))

        #Expectation Step
        if iterCount == 0:
            for readName in qualifiedRead:
                probNumerators = [0 for j in range(len(isoformNames))]
                for i in range(len(isoformNames)):
                    if CompatibleMatrix[readName][isoformNames[i]] == 1:
                        if (isoformLength[isoformNames[i]]-d+1) > 0: 
                            probNumerators[i] = Alpha[i] / (isoformLength[isoformNames[i]]-d+1)
                        else: 
                            probNumerators[i] = Alpha[i]
                probDenominator = sum(probNumerators)      
                if probDenominator != 0: 
                    for i in range(len(isoformNames)):
                        tmpCompatibleMatrix[readName][isoformNames[i]] = probNumerators[i] / probDenominator
        else:
            for readName in qualifiedRead:
                probNumerators = [0 for j in range(len(isoformNames))]
                for i in range(len(isoformNames)):
                    if CompatibleMatrix[readName][isoformNames[i]] == 1: 
                        probNumerators[i] = Alpha[i] * Hfunction[readName][isoformNames[i]]
                probDenominator = sum(probNumerators)    
                if probDenominator != 0: 
                    for i in range(len(isoformNames)):
                        tmpCompatibleMatrix[readName][isoformNames[i]] = probNumerators[i] / probDenominator

        #Maximization Step
        for i in range(len(isoformNames)):
            isoformReadCounts = 0.0
            for readName in qualifiedRead: 
                isoformReadCounts += tmpCompatibleMatrix[readName][isoformNames[i]]
            oldAlpha[i] = Alpha[i]
            Alpha[i] = isoformReadCounts / readCount #readCount is total number of reads per gene

        iterCount += 1

        diff = abs(Alpha[0] - oldAlpha[0])
        for t in range(1, len(Alpha)):
            if abs(Alpha[t] - oldAlpha[t]) > diff: 
                diff = abs(Alpha[t] - oldAlpha[t])

                
    sumAlpha = sum(Alpha)
    if sumAlpha == 0: continue

    for i in range(len(Alpha)):
        Alpha[i] = Alpha[i] / sumAlpha

    isoformRelativeAbundances = [None] * len(isoformNames)
    sumTheta = 0.0
    for i in range(len(Alpha)):
        sumTheta += Alpha[i] / isoformLength[isoformNames[i]]

    print(gene+"\t"+str(iterCount)+" iterations\tDone!")

    
    for i in range(len(Alpha)):
        isoformRelativeAbundances[i] = (Alpha[i] / isoformLength[isoformNames[i]]) / sumTheta
        #print(gene+"\t"+str(geneCount)+"\t"+str(iterCount)+"\t"+isoformNames[i]+"\t"+str(isoformRelativeAbundances[i])+"\t"+str(tmpTime))

        OUT.write(gene+"\t"+isoformNames[i]+"\t"+str(readCount)+"\t"+str(isoformRelativeAbundances[i])+"\n") ## write results into specified file

    print("Time for EM", time.time() - startTime)

print(time.time() - startTime)
OUT.close()
            
                            


                    
