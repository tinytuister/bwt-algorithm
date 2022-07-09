# suquence alignment using BWT algorithm, bioinformatique project. EcoliGenome.fa as main genome and ReadSet5.1.fastq as our queries
# there are couple of capped ranges in reading the genome and amount of queries acquired just make a faster process
# modify the code in order to remove the genome and queries range caps
# examples of diffrent genome and readsets sizes follows :
# 500 lines of genome & 2500 queries = 54.8% success rate
# 1000 lines of genome & 5000 queries = 67.3% success rate
# 1000 lines of genome & 10000 queries = 67.0% success rate
# 5000 lines of genome & 50000 queries = 89.2% success rate
# space and time optimization is at its worst good luck :)

import math
# to keep the characters rank (Which appeared sooner)
charRank = []


def successRate(mycount, total):
    # Average success rate
    return '{:.1%}'.format(mycount/float(total))


def readGenome():
    # reads genome from .fa file and removes the first line and spaces
    with open('EcoliGenome.fa', 'r') as genome:
        first_line = genome.readline()
        genome = genome.read().replace(first_line, '')
        genome = genome.replace('\n', '')
    # [:500] is to cap genome sequence (using genome subsequence for faster calculation)
    return genome[:500]


def readQuery(firstColumn, lastColumn, lastToFirstColumn, mycount):
    # reads queries from .fastq file and searches for the query in genome sub sequence
    with open('ReadSet5.1.fastq', 'r') as query:
        # total queries
        total = 0
        # myrange to cap queries
        myrange = 10000
        # approximate number of queries
        print("Calculating Success Rate For ",
              math.floor(myrange / 4), " Queries . . .")
        for i, line in enumerate(query):
            # remove second condition to uncap
            if i % 4 == 1 and i < myrange:
                success = False
                top = len(firstColumn) - len(firstColumn)
                bot = len(firstColumn) - 1
                total = total+1
                # query[0,149] range is there for getting rid of spaces after sequence, every sequence is exactly 149 characters
                if patternMatching(success, firstColumn, lastColumn, lastToFirstColumn, line[0:149], top, bot) == True:
                    # counter for succeed queries
                    mycount = mycount + 1
    print(successRate(mycount, total))
    return


def bwTransform(myin):
    # making first and last column of bwt
    myin = myin + '$'
    bwtMatrix = [myin[i:] + myin[:i] for i in range(len(myin))]
    bwtMatrix = sorted(bwtMatrix)
    lastColumn = [row[-1:] for row in bwtMatrix]
    firstColumn = [row[:1] for row in bwtMatrix]
    bwt = ''.join(lastColumn)
    fColumn = ''.join(firstColumn)
    return bwt, fColumn


def makeLastToFirstColumn(lastColumn, firstColumn):
    # making last to first column for backward search
    # goes through all the values and ranks to find the exact match (I'm not sure if this is the right implementation)
    lastToFirstColumn = []
    for i in range(len(lastColumn)):
        currValue = lastColumn[i]["value"]
        currRank = lastColumn[i]["rank"]
        for ii in range(len(firstColumn)):
            if firstColumn[ii]["value"] == currValue and firstColumn[ii]["rank"] == currRank:
                lastToFirstColumn.append(ii)
    return lastToFirstColumn


def firstInColumn(queryColumn, query, top):
    # finds the first position where our query character was appeared and returns the position as i
    i = top
    for i in range(i, len(queryColumn)):
        if queryColumn[i]["value"] == query:
            return i


def lastInColumn(queryColumn, query, bot):
    # finds the last position where our query character was appeared and returns the position as i
    i = bot
    for i in reversed(range(len(queryColumn))):
        if queryColumn[i]["value"] == query and i <= bot:
            return i


def patternMatching(success, firstColumn, lastColumn, lastToFirstColumn, query, top, bot):
    # backward searching
    # using boolean variables to know if the algorithm failed to find the sequence
    # when top and bot are not available to find and are set to a nonetype, is considered a fail
    for i in reversed(range(len(query))):
        topDone = False
        botDone = False
        if firstInColumn(lastColumn, query[i], top) != None:
            top = lastToFirstColumn[firstInColumn(lastColumn, query[i], top)]
            topDone = True
        if lastInColumn(lastColumn, query[i], bot) != None:
            bot = lastToFirstColumn[lastInColumn(lastColumn, query[i], bot)]
            botDone = True
        if topDone == False or botDone == False:
            success = False
            break
        else:
            success = True
    return success


def mycounter(query, uniqueChar, i):
    # just a counter to find how many of the same character are appeared before the character (finding rank)
    myRank = 0
    queryChar = query[i]
    while i != -1:
        if uniqueChar[i]["value"] == queryChar:
            myRank = myRank + 1
        i = i - 1
    return myRank


def makeColumn(query):
    # makes two rows from first and last column of bwTransform(myin) with values and ranks
    mycolumn = {}
    uniqueChar = {}
    for i in range(len(query)):
        uniqueChar[i] = {
            "value": query[i]
        }

    for i in range(len(query)):
        mycolumn[i] = {
            "value": query[i],
            "rank": mycounter(query, uniqueChar, i)
        }
    return mycolumn


    # myin was used for assigning a desired short query
firstColumn = {}
lastColumn = {}
#myin = input("Enter a word : ")
# calls bwTransform on two empty rows of bwt(last column) and fColumn(first column)
bwt, fColumn = bwTransform(readGenome())
# makes columns for better representation and assining the ranks and values
lastColumn = makeColumn(bwt)
firstColumn = makeColumn(fColumn)
# makes last to first column for backward searching
lastToFirstColumn = makeLastToFirstColumn(lastColumn, firstColumn)
# just a simple representation of the bwt
for i in range(len(lastColumn)):
    print(firstColumn[i]["value"], firstColumn[i]
          ["rank"], "\t", lastColumn[i]["value"], lastColumn[i]["rank"])
print("\t \t")
# mycount is used to count success rate
mycount = 0
# calls readQuery() to read queries of .fastq and calling a backward search using them
readQuery(firstColumn, lastColumn, lastToFirstColumn, mycount)
