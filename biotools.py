def calcMeltTemp(strand):
    # Count the number of different nucleotides
    numA = 0
    numT = 0
    numG = 0
    numC = 0
    for letter in strand:
        if (letter == "A"):
            numA += 1
        elif (letter == "T"):
            numT += 1
        elif (letter == "G"):
            numG += 1
        else:
            numC += 1
    numNuc = numA + numT + numG + numC
    # Calculate the melting temperature differently based on the size of the strand
    if (numNuc <= 13):
        meltTemp = (numA + numT) * 2 + (numG + numC) * 4
    else:
        meltTemp = 64.9 + 41 * (numG + numC - 16.4) / (numA + numT + numG + numC)

def calcGCContent(strand):
    # Count the number of G's and C's
    numG = 0
    numC = 0
    numNuc = 0
    for letter in strand:
        if (letter == "G"):
            numG += 1
        elif (letter == "C"):
            numC += 1
        numNuc += 1
    return (numG + numC) / numNuc * 100

    



