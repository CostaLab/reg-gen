import sys
#seq = [(1,1),(4,1),(8,1),(51,1),(4,1),(5,1),(4,1),(6,1),(4,1),(1,1),(100,1),(2,1),(5,1),(3,1),(4,1),(1,1),(10,1),(2,1),(5,1),(3,1)]
S,p,q = [],[],[]
prefixW = []
prefixD = []

def setU_i(upperBound):
    n = len(seq)
    j, U = n-1, [-1]*n
    for i in range(n-1, -1, -1):
        while width(i, j) > upperBound:
            j -= 1
        if j >= i: #fuer universal
            U[i] = j
        
    return U

def setL_i(lowerBound):
    n = len(seq)
    j, L = 0, [-1]*n
    for i in range(0, n):
        while j < n and width(i, j) < lowerBound:
            j += 1
        if j == n: break
        L[i] = j
    return L

def width(x, y):
    #return sum( map( lambda x: x [1], seq[x:y+1]))
    if x >= len(seq): x = len(seq)-1
    if y >= len(seq): y = len(seq)-1
    if x == y:
        return float( seq[x][1] )
    else:
        if x == 0:
            tmp = 0
        else:
            tmp = prefixW[x-1]
        return float( prefixW[y] - tmp )

def density(x, y):
    if x >= len(seq): x = len(seq)-1
    if y >= len(seq): y = len(seq)-1
    if x>y: return 0
    #return float(sum( map( lambda x: x[0], seq[x:y+1]) )) / sum( map( lambda x: x[1], seq[x:y+1]) )
    if x == y:
        return float(seq[y][0] / width(x,y))
    else:
        if x == 0:
            tmp = 0
        else:
            tmp = prefixD[x-1]
        return float( (prefixD[y] - tmp) / width(x, y) )
    
def LMatchInitialize(x, y):
    global p, S
    l = u = b = p[b] = y
    for i in range(y, x, -1):
        p[i] = i
        while p[i] < y and density(i,p[i]) <= density(p[i]+1,p[p[i]+1]):
            p[i] = p[p[i]+1]
        if S[p[i]] == -1: S[p[i]] = []
        S[p[i]].insert(0, i)
    
    return l, u, b

def _getMinRemoveS(S_u, l, b):
    for el in S_u:
        if el >= l:
            erg = el
            break
    
    for k in S_u:
        if k > erg:
            S_u.remove(k)
    
    return erg
    
def LMatchFind(i, L, l, u, b, x):
    global S, p
    while l > 1 + max(x, L[i]):
        l -= 1
        if p[l] >= u:
            b = l
    if 1 + max(x, L[i]) >= len(seq): #extra
        l = len(seq)
    if l < 1 + max(x, L[i]): #extra
        l = 1 + max(x, L[i])
    while u >= l and density(i, b-1) > density(i, p[b]):
        u = b-1
        if u >= l:
            b = _getMinRemoveS(S[u], l, b)
    
    return u, l, b

def maxIndex(lowerBound):
    for j in reversed(range(len(seq))):
        if width(j, len(seq)-1) >= lowerBound:
            return j

def maxDensity(g):
    max, pos = -1, -1
    for i in range(len(g)):
        if density(i,g[i]) > max:
            max = density(i,g[i]) 
            pos = (i,g[i]+1)
    return max, pos
    
def MaximumDensitySegmentL(lowerBound, x, y):
    global S, p
    L = setL_i(lowerBound)
    l, u, b = LMatchInitialize(x,y)
    s = maxIndex(lowerBound)
    g = [] 
    for i in range(s,-1,-1):
        if L[i] == len(seq):
            g.append( (i,len(seq)) )
        else:
            goodPartner = LMatchFind(i, L, l, u, b, x)
            g.append( (i, goodPartner ) )
    return maxDensity(g)
    
            
def UMatchInitialize(x, y):
    global q
    u = y
    for i in range(x+1, y+1):
        q[i] = i
        while q[i] > x and q[q[i]-1] != -1 and density(q[q[i]-1], q[i]-1) <= density(q[i], i):
            q[i] = q[q[i]-1]
    return u

def UMatchFind(i, u, U, x):
    while u > U[i]:
        u -= 1
    while u > x and density(i, q[u]-1) > density(i, u):
        u = q[u]-1
    return u
    
def _getPosibleBlock(L, U): #fuer universal
    for i in range(len(L)):
        if L[i] <= U[i]:
            return L[i], U[i], i
        else:
            print("H")
    
def ConstructBlocks(lowerBound, upperBound):
    blocks = []
    L, U = setL_i(lowerBound), setU_i(upperBound)
    s, i = maxIndex(lowerBound), 0
    #fuer universal: x, y = L[0], U[0]
    x, y, i = _getPosibleBlock(L,U) #fuer universal
    blocks.append((x,y))
    while y < U[s]:
        x = y+1
        while (i < s and L[i] < x) or L[i] > U[i]:
            i += 1
        y = U[i]
        blocks.append((x,y)) 
    return blocks

def _getNextBlock(blocks, blocksPara):
    (x,y) = blocks.pop()
    l, u, b, v = blocksPara.pop()
    return x, y, l, u, b, v

def _LUNotPossible(L, U): #fuer universal
    notPossible = 1
    for i in range(len(L)):
        if L[i] == -1: continue
        if L[i] <= U[i]:
            notPossible = 0
    return notPossible
          
def MaximumDensitySegmentLU(lowerBound,upperBound):
    L, U = setL_i(lowerBound), setU_i(upperBound)
    if _LUNotPossible(L,U): #fuer universal
        return -1, (-1,-1)
    blocks = ConstructBlocks(lowerBound, upperBound)
    #print(len(blocks))
    blocksPara = []
    s = maxIndex(lowerBound)
    g = [-1]*(s+1)
    for (i,j) in blocks:
        l, u, b = LMatchInitialize(i,j)
        v = UMatchInitialize(i,j)
        blocksPara.append((l,u,b,v))
    #print('t')
    x, y, l, u, b, v = _getNextBlock(blocks, blocksPara)
    
    for i in range(s,-1,-1):
        #print(len(blocks),file=sys.stderr)
        if L[i] > U[i]: continue #fuer universal
        if x > L[i] or y < L[i]:
            xOld, vOld = x, v 
            x, y, l, u, b, v = _getNextBlock(blocks, blocksPara)
            
        g[i], l, b = LMatchFind(i, L, l, u, b, x)
        u = g[i]
        if U[i] > y:
            alt = UMatchFind(i, vOld, U, xOld)
            vOld = alt
            if density(i,alt) >= density(i,g[i]):
                g[i] = alt

    return maxDensity(g)

def testNaiv(seqTmp, l, u):
    global seq, prefixW, prefixD
    seq, maxF, pos, prefixW, prefixD = seqTmp, -1, (-1,-1), [], []
    getPrefixWD()
    for i in range(0,len(seq)):
        j = i
        #mindest groesse fuer j
        while width(i,j) < l and j < len(seq):
            j += 1
        if j >= len(seq): break
        #j vergroessern
        while width(i,j) <= u:
            if j >= len(seq): break
            if density(i,j) > maxF:
                maxF, pos = density(i,j), (i,j)
            
            j += 1
    
    return (maxF, pos)

def _naivLU(l): #fuer universal
    maxF, pos = -1, (-1,-1)
    for i in range(len(seq)):
        if width(i,i) == 1 and density(i,i) > maxF:
            maxF, pos = density(i, i), (i,i)
    return maxF, pos

def getPrefixWD():
    global prefixW, prefixD
    prefixW.append(seq[0][1])
    prefixD.append(seq[0][0])
    for i in range(1,len(seq)):
        prefixW.append(seq[i][1]+prefixW[i-1])
        prefixD.append(seq[i][0]+prefixD[i-1])


def AlgGoldwasser(seqTmp, lowerBound, upperBound):
    """Find maximum density segment. Give sequence [(value, width)], and
    upper/lower bound for window to search in.
    Return tuple (max. density, c) with c = (s1,s2) tuple with start- and
    end coordinates of window (s1 including, s2 excluding)."""
    global seq, S, p, q, prefixW, prefixD
    seq, prefixW, prefixD = seqTmp, [], []
    S, p, q = [-1]*len(seq), [-1]*len(seq), [-1]*len(seq)
    getPrefixWD()
    #print(prefixW, prefixD)
    
    if lowerBound > upperBound:
        raise Exception('lower bound > upper bound')
    elif seq == []:
        raise Exception('empty sequence')
    elif lowerBound > len(seq):
        raise Exception('wrong bounds')
    elif lowerBound == upperBound == 1:
        return _naivLU(lowerBound)
    return MaximumDensitySegmentLU(lowerBound,upperBound)
    

    
if __name__ == '__main__':
    #global seq, S, p, q, prefixW, prefixD
    #prefixW,prefixD = [], []
    seq = [(2,1),(5,1),(8,1),(4,1),(3,1),(6,1),(6,1), (0,1), (0,1), (0,1), (0,1), (0,1), (9,1),(3,1),(7,1),(1,1),(2,1),(4,1),(3,1)]

    
    #S, p, q = [-1]*len(seq),[-1]*len(seq), [-1]*len(seq)
    #getPrefixWD()
    #seq = [(2,1),(5,1),(8,1),(4,1),(3,1),(6,1),(6,1),(9,1),(3,1),(7,1),(1,1),(2,1),(4,1),(3,1)]
    lowerBound, upperBound = 5, 6
    
    #a, b , c = LMatchInitialize(0,13)
    #print('t')
    #U = setU_i(upperBound)
    #L = setL_i(lowerBound)
    
    #S, p, q = [-1]*len(seq),[-1]*len(seq), [-1]*len(seq)
    print(seq)
    b = AlgGoldwasser(seq, lowerBound, upperBound)
    print(b)
    
    #print(testNaiv(seq, lowerBound, upperBound))
    