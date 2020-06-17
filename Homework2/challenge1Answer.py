import sys
import numpy

"""The following uses Python to challenge you to create an algorithm for finding
matches between a set of aligned strings. Minimal familiarity with Python is 
necessary, notably list and Numpy array slicing. 
"""

"""Problem 1.

Let X be a list of M binary strings (over the alphabet { 0, 1 }) each of length 
N. 

For integer 0<=i<=N we define an ith prefix sort as a lexicographic sort 
(here 0 precedes 1) of the set of ith prefixes: { x[:i] | x in X }.
Similarly an ith reverse prefix sort is a lexicographic sort of the set of
ith prefixes after each prefix is reversed.

Let A be an Mx(N+1) matrix such that for all 0<=i<M, 0<=j<=N, A[i,j] is the 
index in X of the ith string ordered by jth reverse prefix. To break ties 
(equal prefixes) the ordering of the strings in X is used. 

Complete code for the following function that computes A for a given X.

Here X is a Python list of Python strings. 
To represent A we use a 2D Numpy integer array.

Example:

>>> X = getRandomX() #This is in the challenge1UnitTest.py file
>>> X
['110', '000', '001', '010', '100', '001', '100'] #Binary strings, M=7 and N=3
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> 

Hint:
Column j (0 < j <= N) of the matrix can be constructed from column j-1 and the 
symbol in each sequence at index j-1.  

Question 1: In terms of M and N what is the asymptotic cost of your algorithm?
            The cost of my alogrithm is O(NM).
"""
def constructReversePrefixSortMatrix(X):
    #Creates the Mx(N+1) matrix
    A = numpy.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0])+1 ], dtype=int)
    
    #Code to write - you're free to define extra functions 
    #(inline or outside of this function) if you like.

    # adds the inital positions of the values for column 0
    for i in range(len(A)): A[i][0] = i
    # N is set as number of columns and M is a list 
    # of all the indicies form range 0 to M - 1 (the
    # number of rows in A)
    N = len(X[0])
    M = [i for i in range(len(X))]

    # for loop then follows Algorithm 1 in PBWT reading
    for col in range(N):
        a = []; b = []
        for row in M:
            if X[row][col] == '0': a.append(row)
            else: b.append(row)
        M = a + b
        for i in range(len(A)): A[i][col + 1] = M[i]

    
    return A

"""Problem 2: 

Following on from the previous problem, let Y be the MxN matrix such that for 
all 0 <= i < M, 0 <= j < N, Y[i,j] = X[A[i,j]][j].

Complete the following to construct Y for X. 

Hint: You can either use your solution to constructReversePrefixSortMatrix() 
or adapt the code from that algorithm to create Y without using 
constructReversePrefixSortMatrix().

Question 2: In terms of M and N what is the asymptotic cost of your algorithm?
            The asymtotic cost of my algorithm is O(NM)
            which found by NM from constructReversePrefixSortMatrix()
            and NM for assigning the values.
"""
def constructYFromX(X):
    #Creates the MxN matrix
    Y = numpy.empty(shape=[len(X), 0 if len(X) == 0 else len(X[0]) ], dtype=int)
    
    #Code to write - you're free to define extra functions
    #(inline or outside of this function) if you like.

    # calls constructReversePrefixSortMatrix to get A
    A = constructReversePrefixSortMatrix(X)
    N = len(X[0])
    M = len(X)

    # the for loop chain then just assigns values to Y
    # using A[row][col] value as the row and just col from
    # 0 to N - 1 as column value.

    for col in range(N):
        for row in range(M):
            Y[row][col] = X[A[row][col]][col]

    return Y



"""Problem 3.

Y is a transformation of X. Complete the following to construct X from Y, 
returning X as a list of strings as defined in problem 1.
Hint: This is the inverse of X to Y, but the code may look very similar.

Question 3a: In terms of M and N what is the asymptotic cost of your algorithm?
        The asymtotic cost of my algorithm is O(NM).

Question 3b: What could you use the transformation of Y for? 
Hint: consider the BWT.
        Y is useful for finding the original sequence X, alot like the
        BWT where you find the F postion of the L position in the BWT.

Question 3c: Can you come up with a more efficient data structure for storing Y?
        Numpy arrays are already much more efficent vs python lists and dictionaries.
        I can't imagine a much more efficent data structure.
"""
def constructXFromY(Y):
    #Creates the MxN matrix
    X = numpy.empty(shape=[len(Y), 0 if len(Y) == 0 else len(Y[0]) ], dtype=int)
    
    # the idea of M is to have the indicies of the row from 0 to M - 1
    # the only way I could think of keeping track of which 0 or 1 belongs to
    # which string was to impliment the same way we made A but insteaad we
    # use the columns of A as a map to find the repective 0 or 1 and to know
    # which indicy in the string it is in and which string the number
    # belongs to.

    M = [i for i in range(len(Y))]
    N = len(Y[0])
    for col in range(N):
        a = []; b = []
        for row in range(len(M)):
            # inital M stores the values in respective column in X
            for i in range(len(M)): X[M[i]][col] = Y[i][col]

            # a and b are then created by iterating through
            # the respective column 
            # if a -> the M[row] value will go in to a
            # else -> the M[row] value will go in to b

            if Y[row][col] == 0:
                a.append(M[row])
            else:
                b.append(M[row])

        # M is then repalced by a + b now with the values 0 to M - 1 sorted.
        M = a + b
            
    return list(map(lambda i : "".join(map(str, i)), X))#Convert back to a list of strings

"""Problem 4.

Define the common suffix of two strings to be the maximum length suffix shared 
by both strings, e.g. for "10110" and "10010" the common suffix is "10" because 
both end with "10" but not both "110" or both "010". 

Let D be a Mx(N+1) Numpy integer array such that for all 1<=i<M, 1<=j<=N, 
D[i,j] is the length of the common suffix between the substrings X[A[i,j]][:j] 
and X[A[i-1,j]][:j].  

Complete code for the following function that computes D for a given A.

Example:

>>> X = getRandomX()
>>> X
['110', '000', '001', '010', '100', '001', '100']
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> D = constructCommonSuffixMatrix(A, X)
>>> D
array([[0, 0, 0, 0],
       [0, 1, 2, 2],
       [0, 1, 2, 3],
       [0, 1, 1, 1],
       [0, 0, 2, 2],
       [0, 1, 0, 0],
       [0, 1, 1, 3]])

Hints: 

As before, column j (0 < j <= N) of the matrix can be constructed from column j-1 
and thesymbol in each sequence at index j-1.

For an efficient algorithm consider that the length of the common suffix 
between X[A[i,j]][:j] and X[A[i-k,j]][:j], for all 0<k<=i is 
min(D[i-k+1,j], D[i-k+2,j], ..., D[i,j]).

Question 4: In terms of M and N what is the asymptotic cost of your algorithm?
            The cost of my algorithm is O(NM).
"""
def constructCommonSuffixMatrix(A, X):
    D = numpy.zeros(shape=A.shape, dtype=int) #Creates the Mx(N+1) D matrix 

    # set up N and M
    N = len(A[0])
    M = len(X)
    for col in range(1, N):
        c = []; e = []
        p = 0; q = 0
        for row in range(M):
            # Checks which is the min value
            # p or value before D[row][column] + 1
            # then same for q and D[row][column - 1] + 1
            if D[row][col - 1] + 1 < p:
                p = D[row][col - 1] + 1
            if D[row][col - 1] + 1 < q:
                q = D[row][col - 1] + 1
            if X[A[row][col - 1]][col - 1] == '0':
                c.append(p)
                p = float('inf')
            else:
                e.append(q)
                q = float('inf')
        d = c + e
        for i in range(len(D)): D[i][col] = d[i]

    
    return D

# X = ['110', '000', '001', '010', '100', '001', '100']
# A = constructReversePrefixSortMatrix(X)
# print(constructCommonSuffixMatrix(A, X))

"""Problem 5.
    
For a pair of strings X[x], X[y], a long match ending at j is a common substring
of X[x] and X[y] that ends at j (so that X[x][j] != X[y][j] or j == N) that is longer
than a threshold 'minLength'. E.g. for strings "0010100" and "1110111" and length
threshold 2 (or 3) there is a long match "101" ending at 5.
    
The following algorithm enumerates for all long matches between all substrings of
X, except for simplicity those long matches that are not terminated at
the end of the strings.
    
Question 5a: What is the asymptotic cost of the algorithm in terms of M, N and the
number of long matches?
        With # of long matches as L, the asymptotic cost is O(max(NM, L)).
    
Question 5b: Can you see any major time efficiencies that could be gained by
refactoring?
        Yes, if we can run run it at same time as we run find the A and D 
        array values, we could achieve O(M) space.

Question 5c: Can you see any major space efficiencies that could be gained by
refactoring?
        Since b and c only hold inidcies of which strings have a 0 or 1 at
        position j. I think a numpy array would be more space efficient.

Question 5d: Can you imagine alternative algorithms to compute such matches?,
if so, what would be the asymptotic cost and space usage?
        No, because I would use a suffix tree for finding matches. Which I
        learned in lecture is rather space heavy (~15 bytes per node), which
        cancels out its desirable O(m) speed (m being 2n , and n is number 
        of nodes).

"""
def getLongMatches(X, minLength):
    assert minLength > 0
    
    A = constructReversePrefixSortMatrix(X)
    D = constructCommonSuffixMatrix(A, X)
    
    #For each column, in ascending order of column index
    for j in range(1, 0 if len(X) == 0 else len(X[0])):
        #Working arrays used to store indices of strings containing long matches
        #b is an array of strings that have a '0' at position j
        #c is an array of strings that have a '1' at position j
        #When reporting long matches we'll report all pairs of indices in b X c,
        #as these are the long matches that end at j.
        b, c = [], []
        
        #Iterate over the aligned symbols in column j in reverse prefix order
        for i in range(len(X)):
            #For each string in the order check if there is a long match between
            #it and the previous string.
            #If there isn't a long match then this implies that there can
            #be no long matches ending at j between sequences indices in A[:i,j]
            #and sequence indices in A[i:,j], thus we report all long matches
            #found so far and empty the arrays storing long matches.
            if D[i,j] < minLength:
                for x in b:
                    for y in c:
                        #The yield keyword converts the function into a
                        #generator - alternatively we could just to append to
                        #a list and return the list
                        
                        #We return the match as tuple of two sequence
                        #indices (ordered by order in X) and coordinate at which
                        #the match ends
                        yield (x, y, j) if x < y else (y, x, j)
                b, c = [], []
            
            #Partition the sequences by if they have '0' or '1' at position j.
            if X[A[i,j]][j] == '0':
                b.append(A[i,j])
            else:
                c.append(A[i,j])
        
        #Report any leftover long matches for the column
        for x in b:
            for y in c:
                yield (x, y, j) if x < y else (y, x, j)
