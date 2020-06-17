import array
import sys
import numpy
# import pysam
import argparse
import logging 
import time
logger = logging.getLogger()

"""See the comments below to see the code you need to complete.
"""

class MinimizerIndexer(object):
    """ Simple minimizer based substring-indexer. 
    
    Please read: https://doi.org/10.1093/bioinformatics/bth408
    
    Related to idea of min-hash index and other "sketch" methods.
    """
    def __init__(self, targetString, w, k, t):
        """ The target string is a string/array of form "[ACGT]*".
        
        Stores the lexicographically smallest k-mer in each window of length w, such that w >= k positions. This
        smallest k-mer is termed a minmer. 
        
        If a minmer occurs in the target sequence more than t times as a minmer then it is omitted from the index, i.e. if the given minmer (kmer) is a minmer
        in more than t different locations in the target string. Note, a minmer may be the minmer for more than t distinct windows
        and not be pruned, we remove minmers only if they have more than t distinct occurrences as minmers in the sequence.
        """
        
        self.targetString = targetString
        self.w = w
        self.k = k
        self.t = t # If a minmer occurs more than t times then its entry is removed from the index
        # This is a heuristic to remove repetitive minmers that would create many spurious alignments between
        # repeats
        
        # Hash of minmers to query locations, stored as a map whose keys
        # are minmers and whose values are lists of the start indexes of
        # occurrences of the corresponding minmer in the targetString, 
        # sorted in ascending order of index in the targetString.
        #
        # For example if k = 2 and w = 4 and targetString = "GATTACATTT"
        #
        # GATTACATTT
        # GATT (AT)
        #  ATTA (AT)
        #   TTAC (AC)
        #    TACA (AC)
        #     ACAT (AC)
        #      CATT (AT)
        #       ATTT (AT)
        #
        # then self.minimizerMap = { "AT":(1,6), "AC":(4,) }
        self.minimizerMap = {}
        # Code to complete to build index - you are free to define additional functions
        self.size = len(targetString)
        for nuc in range(self.size - self.w + 1):
            window = self.targetString[nuc: nuc + self.w]
            minimer = None # some holder for the minimizer

            # implement the window alogrithm from the reading.
            for k_position in range(self.w - k + 1):
                kmer = window[k_position: k_position + k]
                if minimer is None or kmer < minimer[0]: 
                    minimer = (kmer, nuc + k_position)
            self.minimizerMap.setdefault(minimer[0], set()).add(minimer[1])
            if len(self.minimizerMap[minimer[0]]) > t: del self.minimizerMap[minimer[0]]

    def getMatches(self, searchString):
        """ Iterates through search string finding minmers in searchString and
        yields their list of minmer occurrences in targetString, each as a pair of (x, (y,)*N),
        where x is the index in searchString and y is an occurrence in targetString.
        
        For example if k = 2 and w = 4 and targetString = "GATTACATTT" and searchString = "GATTTAC"
        then self.minimizerMap = { "AT":(1,6), "AC":(4,) }
        and getMatches will yield the following sequence:
        (1, (1,6)), (5, (4,))
        
        You will need to use the "yield" keyword
        """
        # Code to complete - you are free to define additional functions
        size = len(searchString)
        minimers = set()
        for nuc in range(size - self.w + 1):

            # builds the window size per iteration
            window = searchString[nuc : nuc + self.w]
            for k_position in range(self.w - self.k + 1):

                # looks at all k's in the window
                kmer = window[k_position:k_position + self.k]
                
                # if a new minimer appears, yield it
                if kmer in self.minimizerMap and kmer not in minimers:
                    minimers.add(kmer)
                    yield (nuc + k_position, self.minimizerMap[kmer])

class SeedCluster:
    """ Represents a set of seeds between two strings.
    """
    def __init__(self, seeds):
        """ Seeds is a list of pairs [ (x_1, y_1), (x_2, y_2), ..., ], each is an instance of a seed 
        (see static cluster seeds method below: static methods: https://realpython.com/blog/python/instance-class-and-static-methods-demystified/)
        """
        seeds = list(seeds)
        seeds.sort()
        self.seeds = seeds
        # Gather the minimum and maximum x and y coordinates
        self.minX = seeds[0][0]
        self.maxX = seeds[-1][0]
        ys = list(map(lambda x : x[1], seeds))
        self.minY = min(ys)
        self.maxY = max(ys)

    @staticmethod
    def clusterSeeds(seeds, l):
        """ Cluster seeds (k-mer instances) in two strings. This is a static constructor method that creates a set
        of SeedCluster instances.
        
        Here seeds is a list of tuples, each tuple has the form (x, (y_1, y_2, ... )), where x is the coordinate
        in the first string and y_1, y_2, ... are coordinates in the second string. Each pair of x and y_i
        is an occurence of a shared k-mer in both strings, termed a *seed*, such that the k-mer 
        occurrence starts at position x in the first string and starts at position y_i in the second string.
        The input seeds list contains no duplicates and is sorted in ascending order, 
        first by x coordinate (so each successive tuple will have a greater  
        x coordinate), and then each in tuple the y coordinates are sorted in ascending order.
        
        Two seeds (x_1, y_1), (x_2, y_2) are *close* if the absolute distances | x_2 - x_1 | and | y_2 - y_1 |
        are both less than or equal to l.   
        
        Consider a *seed graph* in which the nodes are the seeds, and there is an edge between two seeds if they
        are close. clusterSeeds returns the connected components of this graph
        (https://en.wikipedia.org/wiki/Connected_component_(graph_theory)).
        
        The return value is a Python set of SeedCluster object, each representing a connected component of seeds in the 
        seed graph.
        
        (QUESTION 1): The clustering of seeds is very simplistic. Can you suggest alternative strategies by
        which the seeds could be clustered, and what the potential benefits such alternative strategies could
        have? Consider the types of information you could use.

        You could turn all the (x (y_1, y_2, ..., y_n)) tuples in to (x, y) tuples and then from the first one see
        which satisfy the l distance requirement. That would make your first cluster, and then the remaining ones you
        do the process again relative to the first remaining tuple until all the tuples are clustered. After that see if there
        are common tuples in each cluster and then merge the clusters. This would remove the need to make edge class and graph class
        and just basic knowlege of unions and intersections (stats) to see if you can merge clusters (ex. intersection of two sets, if
        there is 1 or more intersecting tuples, merge.). The potential benefits of this is that the huge set up cost for the adjmatrix
        can be avoided.
        """ 
        
        # Code to complete - you are free to define other functions as you like
        
        # The function will work like this,
        # first change all the (x( y1, y2, .., yn)) tuples into (x, y) tuples.
        # then from that point make the adj matrix so that it starts with itself
        # 
        # Once the adjMatrix is set up the following will happen.
        # first a 2 seeds will compared.
        # if both are <= l distance away do this:
        #     check if they have different componenets:
        #         if they do merge the components then set
        #         all the seeds in that component to have the 
        #         same seed values: the sum of both components
        # 
        # I would like to cite James Casaletto for his help in office hours and
        # Alex Pearson and Kavya for thier helpful piazza posts on this probelm.

        # list for all the seeds
        adjMatrix = dict()
        Visited = dict()
        # first the tuple is turned into (x, y) tuples
        for x_value in range(len(seeds)):
            if seeds[x_value][1]:
                for y_value in seeds[x_value][1]:
                    adjMatrix.setdefault((seeds[x_value][0], y_value), []).append((seeds[x_value][0], y_value))
                    Visited[(seeds[x_value][0], y_value)] = False
        '''
        # now we cluster the seeds.
        for seed1 in range(len(list_of_seeds)):
            x_value, y_value = list_of_seeds[seed1]
            for seed2 in range(seed1, len(list_of_seeds)):
                x_value_2, y_value_2 = list_of_seeds[seed2]

                # if seeds are close get thier adj list
                if abs(x_value - x_value_2) <= l and abs(y_value - y_value_2) <= l:
                    adj_seeds_1 = adjMatrix[list_of_seeds[seed1]]
                    adj_seeds_2 = adjMatrix[list_of_seeds[seed2]]
                    
                    # if lists are different, concatenate.
                    if adj_seeds_1 != adj_seeds_2:
                        adj_seeds = adj_seeds_1 + adj_seeds_2
                       
                        # now set all seeds in adjMatrix equal to same components
                        for seed in adj_seeds_1:
                            adjMatrix[seed] = adj_seeds
                        for seed in adj_seeds_2:
                            adjMatrix[seed] = adj_seeds
        '''


        def Visit(value):
            ''' Visit the value and see whats near it. '''
            current = value
            Visited[current] = True
            for i in adjMatrix:
                if current == i : continue
                if abs(current[0] - i[0]) <= l and abs(current[1] - i[1]) <= l:
                    newList = list(set(adjMatrix[current] + adjMatrix[i]))
                    for item in (adjMatrix[i] + adjMatrix[current]):
                        adjMatrix[item] = newList
                    Visited[i] = True                   


        for seed in adjMatrix:
            if not Visited[seed]:
                Visit(seed)

        # turn the adj lists into seedCluster objects.
        set_of_clusters = set()

        for cluster in adjMatrix:
            set_of_clusters.add(SeedCluster(adjMatrix[cluster]))        
        return set_of_clusters 


class SmithWaterman(object):
    def __init__(self, string1, string2, gapScore=-2, matchScore=3, mismatchScore=-3):
        """ Finds an optimal local alignment of two strings.
        
        Implements the Smith-Waterman algorithm: 
        https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
        
        (QUESTION 2): The Smith-Waterman algorithm finds the globally optimal local alignment between to 
        strings, but requires O(|string1| * |string2|) time. Suggest alternative strategies you could implement
        to accelerate the finding of reasonable local alignments. What drawbacks might such alternatives have?

        The algorithm for Smith-Waterman has to deal with many tie breakers to find only one alignment.
        You could use seeds from the Minimizer class to align to each string. And use those positions to guide
        you the best alignment instead of backtracking in a matrix. The drawback could be that if the minimer
        occurs frequently (like near t), its placement might dampen your ability to get the best alignment. Like
        for example "A_long_long_time" and "this_is_a_long_long_time" the minimer "long" would not be helpful in 
        finding the alignment of interest since it occurs more than once. 
        """
        # set all inputs as attributes of class
        self.string1 = string1
        self.string2 = string2
        self.gap = gapScore
        self.matchScore = matchScore
        self.mismatchScore = mismatchScore

        # make a numpy array for the matrix
        self.matrix = numpy.zeros((len(string1) + 1, len(string2) + 1), dtype = int)

        # then set up the matrix and its values.
        for x_pos in range(1, len(string1) + 1):

            # set up the y axis
            for y_pos in range(1, len(string2) + 1):

                # check if we have a match or a mismatch
                if self.string1[x_pos - 1] == self.string2[y_pos - 1]:
                    SofIJ = self.matchScore
                else: SofIJ = self.mismatchScore

                # find the diagonal, horizontal, and veritcal max socres
                diag = self.matrix[x_pos - 1][y_pos - 1] + SofIJ
                vertical = max([self.matrix[x_pos][i] + (self.gap) * (y_pos - i) for i in range(y_pos)])
                horizontal = max([self.matrix[i][y_pos] + (self.gap) * (x_pos - i) for i in range(x_pos)])
                
                # make sure they are not < 0
                if horizontal < 0: horizontal = 0
                if vertical < 0: vertical = 0
                
                # set position equal to the max of the path before it
                self.matrix[x_pos][y_pos] = max([diag, horizontal, vertical])


    def getAlignment(self):
        """ Returns an optimal local alignment of two strings. Alignment
        is returned as an ordered list of aligned pairs.
        
        e.g. For the two strings GATTACA and CTACC an optimal local alignment
        is (GAT)TAC(A)
             (C)TAC(C)
        where the characters in brackets are unaligned. This alignment would be returned as
        [ (3, 1), (4, 2), (5, 3) ] 
        
        In cases where there is a tie between optimal sub-alignments use the following rule:
        Let (i, j) be a point in the edit matrix, if there is a tie between possible sub-alignments
        (e.g. you could chooose equally between different possibilities), choose the (i, j) to (i-1, j-1)
        (match) in preference, then the (i, j) to (i-1, j) (insert in string1) in preference and
        then (i, j) to (i, j-1) (insert in string2).
        """
        # Code to complete - generated by traceback through matrix to generate aligned pairs
        
        # find the position of the max_value
        max_value = self.getMaxAlignmentScore()
        max_pos = tuple(numpy.argwhere(self.matrix == max_value)[-1])
        x_pos = max_pos[0]; y_pos = max_pos[1]

        # array that holds the tuples
        path = list()

        # now find the path to the 0
        
        while self.matrix[x_pos][y_pos] != 0:
            
            # if diagonal is a match take that as priority
            if self.string1[x_pos - 1] == self.string2[y_pos - 1]:
                path.append((x_pos - 1, y_pos - 1))
                x_pos -=1; y_pos -= 1
                continue

            # finds the best horizontal alignment
            bestX = 0; bestY = y_pos - 1
            for i in range(x_pos - 1):
                if self.matrix[i][y_pos - 1] >= self.matrix[bestX][bestY]:
                    bestX = i
            
            # finds best vertical alignment
            bestX_vertical = x_pos - 1; bestY_vertical = 0
            for i in range(y_pos - 1):
                if self.matrix[x_pos - 1][i] >= self.matrix[bestX_vertical][bestY_vertical]:
                    bestY_vertical = i
            
            # if diagonal not satisfied, see which is better
            # the horizontal of vertical alignment.
            if self.matrix[bestX][bestY] < self.matrix[bestX_vertical][bestY_vertical]:
                path.append((bestX_vertical, bestY_vertical))
                x_pos = bestX_vertical; y_pos = bestY_vertical
            else:
                path.append((bestX, bestY))
                x_pos = bestX; y_pos = bestY

        return path[::-1] # reversed because we want origin to highest element.
    
    def getMaxAlignmentScore(self):
        """ Returns the maximum alignment score
        """
        # get max of each row
        # max_scores = [max(i) for i in self.matrix]

        # return the max of the max vaules
        return numpy.max(self.matrix)

def simpleMap(targetString, minimizerIndex, queryString, config):
    """ Function takes a target string with precomputed minimizer index and a query string
    and returns the best alignment it finds between target and query, using the given options specified in config.
    
    Maps the string in both its forward and reverse complement orientations.
    
    (QUESTION 3): The code below is functional, but very slow. Suggest ways you could potentially accelerate it, 
    and note any drawbacks this might have.

    With all 3 algorithms essential in quadratic speed, time will quickly add up for the simpleMap function. I think
    running SeedCluster.clusterSeeds() with in get matches would speed up the class since the seeds produced in getMatches()
    is the input for the clusterSeeds() function. Another way could be to also do forward and reverseMaping at same time
    to double the speed, the trade off would be the increased amount of space needed to accomodate.
    """
    bestAlignment = [None]
    
    def mapForwards(queryString):
        """ Maps the query string forwards
        """
        # Find seed matches, aka "aligned kmers"
        seeds = list(minimizerIndex.getMatches(queryString))
        
        # For each cluster of seeds
        for seedCluster in SeedCluster.clusterSeeds(list(seeds), l=config.l):
            
            # Get substring of query and target to align
            queryStringStart = max(0, seedCluster.minX - config.c) # Inclusive coordinate
            queryStringEnd = min(len(queryString), seedCluster.maxX + config.k + config.c) # Exclusive coordinate
            querySubstring = queryString[queryStringStart:queryStringEnd]
            
            targetStringStart = max(0, seedCluster.minY - config.c) # Inclusive coordinate
            targetStringEnd = min(len(targetString), seedCluster.maxY + config.k + config.c) # Exclusive coordinate
            targetSubstring = targetString[targetStringStart:targetStringEnd]
            
            print( "target_aligning", targetStringStart, targetStringEnd, targetSubstring )
            print( "query_aligning", queryStringStart, queryStringEnd, querySubstring )
            
            # Align the genome and read substring
            alignment = SmithWaterman(targetSubstring, querySubstring, 
                                      gapScore=config.gapScore, 
                                      matchScore=config.matchScore,
                                      mismatchScore=config.mismatchScore)
            
            # Update best alignment if needed
            if bestAlignment[0] == None or alignment.getMaxAlignmentScore() > bestAlignment[0].getMaxAlignmentScore():
                bestAlignment[0] = alignment
        
        return bestAlignment
    
    def reverseComplement(string):
        """Computes the reverse complement of a string
        """
        rMap = { "A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
        return "".join(rMap[i] for i in string[::-1])
                
    # Run mapping forwards and reverse
    mapForwards(queryString)
    mapForwards(reverseComplement(queryString))
    
    return bestAlignment[0]

class Config():
    """ Minimal configuration class for handing around parameters
    """
    def __init__(self):
        self.w = 30
        self.k = 20
        self.t = 10
        self.l = 30
        self.c = 100
        self.gapScore=-2
        self.matchScore=3
        self.mismatchScore=-3
        self.logLevel = "INFO"
        
def main():
    # Read parameters
    config = Config()
    
    #Parse the inputs args/options
    parser = argparse.ArgumentParser(usage="target_fasta query_fastq [options]")

    parser.add_argument("target_fasta", type=str,
                        help="The target genome fasta file.")
    parser.add_argument("query_fastq", type=str,
                        help="The query sequences.")
    
    parser.add_argument("--w", dest="w", help="Length of minimizer window. Default=%s" % config.w, default=config.w)
    parser.add_argument("--k", dest="k", help="Length of k-mer. Default=%s" % config.k, default=config.k)
    parser.add_argument("--t", dest="t", help="Discard minmers that occur more frequently " 
                                            "in the target than t. Default=%s" % config.w, default=config.w)
    parser.add_argument("--l", dest="l", help="Cluster two minmers into the same cluster if within l bases of"
                                            " each other in both target and query. Default=%s" % config.l, default=config.l)
    parser.add_argument("--c", dest="c", help="Add this many bases to the prefix and suffix of a seed cluster in the"
                                            " target and query sequence. Default=%s" % config.c, default=config.c)
    parser.add_argument("--gapScore", dest="gapScore", help="Smith-Waterman gap-score. Default=%s" % 
                      config.gapScore, default=config.gapScore)
    parser.add_argument("--matchScore", dest="matchScore", help="Smith-Waterman match-score. Default=%s" % 
                      config.gapScore, default=config.gapScore)
    parser.add_argument("--mismatchScore", dest="mismatchScore", help="Smith-Waterman mismatch-score. Default=%s" % 
                      config.mismatchScore, default=config.mismatchScore)
    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" % 
                      config.logLevel, default=config.logLevel)
    parser.add_argument("-v", "--version", action='version', version=' %(prog)s 0.1')
    
    options = parser.parse_args()
    
    # Parse the log level
    numeric_level = getattr(logging, options.logLevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.logLevel)
    
    # Setup a logger
    logger.setLevel(numeric_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.debug("Established logger")
    
    startTime = time.time()
    
    # Parse the target sequence and read the first sequence
    with pysam.FastaFile(options.target_fasta) as targetFasta:
        targetString = targetFasta.fetch(targetFasta.references[0])
    logger.info("Parsed target string. Length: %s" % len(targetString))
    
    # Build minimizer index
    minimizerIndex = MinimizerIndexer(targetString.upper(), w=options.w, k=options.k, t=options.t)
    minmerInstances = sum(map(len, minimizerIndex.minimizerMap.values()))
    logger.info("Built minimizer index in %s seconds. #minmers: %s, #minmer instances: %s" %
                 ((time.time()-startTime), len(minimizerIndex.minimizerMap), minmerInstances))
    
    # Open the query files
    alignmentScores = [] # Array storing the alignment scores found
    with pysam.FastqFile(options.query_fastq) as queryFastq:
        # For each query string build alignment
        for query, queryIndex in zip(queryFastq, range(sys.maxsize)):
            print (queryIndex)
            alignment = simpleMap(targetString, minimizerIndex, query.sequence.upper(), config)
            alignmentScore = 0 if alignment is None else alignment.getMaxAlignmentScore()
            alignmentScores.append(alignmentScore)
            logger.debug("Mapped query sequence #%i, length: %s alignment_found?: %s "
                         "max_alignment_score: %s" % 
                         (queryIndex, len(query.sequence), alignment is not None, alignmentScore)) 
            # Comment this out to test on a subset
            #if queryIndex > 100:
            #    break
    
    # Print some stats
    logger.critical("Finished alignments in %s total seconds, average alignment score: %s" % 
                    (time.time()-startTime, float(sum(alignmentScores))/len(alignmentScores)))
    
if __name__ == '__main__':
    main()
