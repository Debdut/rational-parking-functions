
def minAreaDyckPath(dim_n,dim_m):
    #This function computes an m,n Dyck path of minimal area.
    #Here the m is the number of columns and the n is the number
    #of rows, which equals the length of the path.
    
    #Set empty min path
    min_path = [1..dim_n]
    
    #Set slope
    slpe = float(dim_n/dim_m)
    
    for i in range(dim_n):
        min_path[i] = floor(float(i/slpe))
        
    return min_path
    
#------------------------------------------------------------

def isDyckPath(path,dim_m):
    #This function checks whether or not a given vector
    #path is a m,n-Dyck Path where n = len(path), m = dim_m
    
    #Compute the minAreaDyckPath(len(path),dim_m))
    min_path=minAreaDyckPath(len(path),dim_m)
    
    #Sets a boolean flag
    isPath=true
    
    for i in range(len(path)):
        if path[i]>min_path[i]:
            isPath=false
            
    return isPath
    
#------------------------------------------------------------   
    
def areaDyckPath(path,dim_m):
    #This function computes the area of an 
    #(m,n)-Dyck Path, where m = dim_m and
    #n=len(path)
    
    #Compute the area of the minAreaDyckPath
    minSumArea = sum(minAreaDyckPath(len(path),dim_m))
    
    #Take the difference of the sum(path) and minSumArea
    return(minSumArea-sum(path))
    
#------------------------------------------------------------     
    
def armOfSquare(x_coord, y_coord, path):
    #This function computes the arm of a
    #single square above an (m,n)-Dyck path, where
    #m dim_m and n = len(path)
    
    return path[y_coord-1]-x_coord
    
#------------------------------------------------------------ 

def legOfSquare(x_coord, y_coord, path):
    #This function computes the leg of a
    #single square above an (m,n)-Dyck path, where
    #m dim_m and n = len(path)
    
    #Find the last entry of 'path' that is less
    #then x-coord
    i=0
    while path[i] < x_coord:
        i = i + 1
    
    return y_coord - 1 - i
    
#------------------------------------------------------------ 
    
def dinvCharacteristic(x_coord, y_coord, path, dim_m):
    #This function determines whether a square above
    #the Dyck path should count toward computing the
    #dinv of the said path
    
    #Retrieve the arm and leg statistics of the square
    arm = armOfSquare(x_coord,y_coord,path)
    leg = legOfSquare(x_coord,y_coord,path)    
    
    return (len(path)*arm < dim_m*(leg+1) and dim_m*leg < len(path)*(arm+1))
    
#------------------------------------------------------------ 
    
def dinvDyckPath(path,dim_m):
    #This function computes the dinv of a
    #(m,n)-Dyck path.  
    
    dinv = 0
    for i in range(len(path)):
        for j in [1..path[i]]:
            if dinvCharacteristic(j,i+1,path,dim_m): 
                dinv = dinv + 1
    
    return dinv
    
#------------------------------------------------------------ 
    
def ranksOnDyckPath(path,dim_m):
    #This function computes the ranks of the
    #squares immediately to the right of an
    #(m,n)-Dyck path.
    
    #Construct a default list
    ranks = range(len(path))
    
    for i in range(len(path)):
        ranks[i] = (-1 - path[i])*len(path) + i*dim_m
    
    return ranks
    
#------------------------------------------------------------ 
    
def maxParkingFunction(path,dim_m):
    #This function computes the parking function that
    #has a maximal dinv.  Since the Dyck path is already
    #stored in the variable 'path,' the function only
    #returns the parking function vector.
    
    #Construct a default list
    maxPF = range(len(path))
    
    #Find the ranks to the right of the Dyck path and
    #sort them.
    ranks = ranksOnDyckPath(path,dim_m)
    sort_ranks = sorted(ranks)
    
    #Compute the maxParkingFunction
    for i in range(len(path)):
       maxPF[ranks.index( sort_ranks[i] )] = i+1
       
    return maxPF
    
#------------------------------------------------------------ 

def tdinvParkingFunction(path,PF,dim_m):
    #This function calculates the tdinv of a given
    #parking function input as a Dyck path, 'path,'
    #and a parking vector, 'PF.'
    
    ranks = ranksOnDyckPath(path,dim_m)
    tdinv = 0
    
    for i in [1..len(path)]:
        rank_i = ranks[ PF.index(i) ]
        for j in [i+1..len(path)]:
            rank_j = ranks[ PF.index(j) ]
            if rank_i < rank_j and rank_j < rank_i + dim_m:
                tdinv = tdinv + 1
                
    return tdinv
    
#------------------------------------------------------------ 
    
def dinvParkingFunction(path,PF,dim_m):
    #This function calculates the dinv of a given
    #parking function input as a Dyck path, 'path,'
    #and a parking vector, 'PF.'
    
    path_dinv = dinvDyckPath(path,dim_m)
    tdinv     = tdinvParkingFunction(path,PF,dim_m)
    maxPF     = maxParkingFunction(path,dim_m)
    maxdinv   = tdinvParkingFunction(path,maxPF,dim_m)
    #print(path_dinv)
    #print(tdinv)
    #print(maxdinv)
    
    return path_dinv + tdinv - maxdinv

#------------------------------------------------------------ 
    
def getPermutation(path,PF,dim_m):
    #This function computes the permutation of 
    #a given parking function input as a Dyck path, 'path,'
    #and a parking vector, 'PF.'
    
    #First compute the ranks of the cars in the parking function
    car_ranks = ranksOnDyckPath(path,dim_m)
    sort_ranks = list(reversed( sorted(car_ranks) ))
    #print(car_ranks, sort_ranks)
    
    #Beginning with the largest rank and working towards the 
    #smallest, write down the cars in the PF
    permutation = [1..len(path)]
    for i in range(len(path)):
        permutation[i] = PF[car_ranks.index( sort_ranks[i] )]
        
    return permutation

#------------------------------------------------------------ 
    
def pides(path,PF,dim_m,adjusted=True):
    #This function returns the inverse descent set, as a list,
    #of a given parking function input as a Dyck path, 'path,'
    #and a parking vector, 'PF.'
    
    invDes = Permutation(getPermutation(path,PF,dim_m)).inverse().descents()
    if adjusted:
        invDes = [x+1 for x in invDes]
    
    return invDes
    

#------------------------------------------------------------ 
    
def decrementPath(location,path,min_path, all_paths):
    #This function computes the next 
    #'larger' path of a given path and 
    #the minimal path min_path
    
    new_path = list(path)
    new_path[location] = new_path[location] - 1
    if location > 0:
        
        if new_path[location] < new_path[location-1]:
            new_path[location] = min_path[location]
            new_path = decrementPath(location-1,new_path,min_path, all_paths)
    
    return new_path

#------------------------------------------------------------ 

def allDyckPaths(dim_n,dim_m):
    #This function computes all the Dyck paths
    #in a m(cols) by n(rows) grid
   
    min_path  = minAreaDyckPath(dim_n,dim_m)
    max_path  = [0]*dim_n
    path      = min_path
    all_paths = list()
    all_paths.append(min_path)
    #print(path, max_path, min_path)
    #print("all_paths:",all_paths)
    
    while (path != max_path):
        #print("all_paths:",all_paths)

        new_path = decrementPath(len(path)-1, path, min_path, all_paths)
        #path = [0,0]
        #print("new_path:",new_path)
        #print("all_paths:",all_paths)
        
        all_paths.append(new_path)  
        #print("all_paths:",all_paths)
        path = new_path
                        
    return all_paths
    
#------------------------------------------------------------ 
    
def isPFCompatible(PF,DyckPath):
    #This function tests whether or not a given
    #parking function is compatible with the given
    #DyckPath
    
    des_set = Permutation(PF).descents()
    
    is_compatible = True
    i = 0
    while ((is_compatible) and (i < len(des_set))):
        if (DyckPath[ des_set[i] ] == DyckPath[ des_set[i]+1 ]):
            is_compatible = False
        i=i+1
        
    return is_compatible

#------------------------------------------------------------ 

def allPFs(DyckPath):
    #This function computes all the parking functions
    #given a DyckPath.  It only returns the permutation
    #as the DyckPath is already known.
    
    id = [1..len(DyckPath)]
    all_perms = Arrangements(id,len(DyckPath)).list()
    
    all_pfs = []
    for perm in all_perms:
        if (isPFCompatible(perm, DyckPath)):
            all_pfs.append(perm)
    
    return(all_pfs)
    
#------------------------------------------------------------ 
    
def allDyckPathswithPFs(dim_n,dim_m):
    #This function computes all the parking functions
    #in a m(cols) by n(rows) grid
    
    all_PF_pairs = []
    
    #Compute all the Dyck Paths
    all_Dyck_paths = allDyckPaths(dim_n,dim_m)
    
    #For each Dyck_path in all_Dyck_paths, compute
    #all the parking functions
    for Dyck_path in all_Dyck_paths:
        presentations = allPFs(Dyck_path)
        all_PF_pairs.append([Dyck_path,presentations])
        
    return(all_PF_pairs)

#------------------------------------------------------------ 
    
def printAllPFs(dim_n,dim_m):
    #This function prints all the parking functions
    #in a m(cols) by n(rows) grid
    
    all_PF_pairs = allDyckPathswithPFs(dim_n,dim_m)
    
    for i in all_PF_pairs:
        print "Dyck Path: ", i[0]
        print "Parking Functions: ", i[1]
        print " "

#------------------------------------------------------------ 

def printHikitaPoly(dim_n,dim_m):
    #This function prints the Hikita polynomial
    #by computing all parking functions in the
    #dim_n by dim_m grid.  
    
    all_pairs = allDyckPathswithPFs(dim_n,dim_m)
    
    terms = []
    
    for pair in all_pairs:
        areaPF = areaDyckPath(pair[0],dim_m)
        for pf in pair[1]:
             dinvPF = dinvParkingFunction(pair[0],pf,dim_m) 
             pidesPF = pides(pair[0],pf,dim_m)
             terms.append( "t^" + str(areaPF) + " q^" + str(dinvPF) + " F_" + str(pidesPF) )

    Hikita_poly = terms[0]
    for i in [1..len(terms)-1]:
        Hikita_poly = Hikita_poly + "  +  " + terms[i]
                     
    return Hikita_poly
    
#------------------------------------------------------------
    
def HikitaPoly(dim_m,dim_n, expand = False):
    #This function prints the Hikita polynomial
    #by computing all parking functions in the
    #dim_n by dim_m grid.  It uses a fraction field
    
    R = QQ['q','t'].fraction_field()
    (q,t) = R.gens()
    
    Ht = SymmetricFunctions(R).macdonald().Ht()
    Fqt = QuasiSymmetricFunctions(Ht.base_ring()).F()
    
    all_pairs = allDyckPathswithPFs(dim_n,dim_m)
    
    HP = 0
    
    for pair in all_pairs:
        areaPF = areaDyckPath(pair[0],dim_m)
        for pf in pair[1]:
             dinvPF = dinvParkingFunction(pair[0],pf,dim_m)
             pidesPF = pides(pair[0],pf,dim_m,False)
             compo = Compositions().from_descents(pidesPF, dim_n)
             HP = HP + t^areaPF * q^dinvPF * Fqt[compo]
    
    if (expand == True):
        HP = HP.expand(dim_n)

    return HP


def mydinvParkingFunction(path,PF,dim_m):
    #This function calculates the tdinv of a given
    #parking function input as a Dyck path, 'path,'
    #and a parking vector, 'PF.'
    
    ranks = ranksOnDyckPath(path,dim_m)
    dinv = 0
    dim_n = len(path)
    
    for i in [1..dim_n]:
        rank_i = ranks[ PF.index(i) ]
        for j in [i+1..len(path)]:
            rank_j = ranks[ PF.index(j) ]
            if rank_i < rank_j and rank_j < rank_i + dim_n - dim_m:
                dinv = dinv + 1
                
    return dinv


def myHikitaPoly(dim_m,dim_n, expand = False):
    #This function prints the Hikita polynomial
    #by computing all parking functions in the
    #dim_n by dim_m grid.  It uses a fraction field
    
    R = QQ['q','t'].fraction_field()
    (q,t) = R.gens()
    
    Ht = SymmetricFunctions(R).macdonald().Ht()
    Fqt = QuasiSymmetricFunctions(Ht.base_ring()).F()
    
    all_pairs = allDyckPathswithPFs(dim_n,dim_m)
    
    HP = 0
    
    for pair in all_pairs:
        areaPF = areaDyckPath(pair[0],dim_m)
        for pf in pair[1]:
             dinvPF = mydinvParkingFunction(pair[0],pf,dim_m) 
             pidesPF = pides(pair[0],pf,dim_m,False)
             compo = Compositions().from_descents(pidesPF, dim_n)
             HP = HP + t^areaPF * q^dinvPF * Fqt[compo]
    
    if (expand == True):
        HP = HP.expand(dim_n)

    return HP



# Tests

K = QQ['q','t'].fraction_field()
s = SymmetricFunctions(K).schur()

m = 3
path = [0,0,0,0,1,1,2]
PF = [2,3,5,7,4,6,1]

isPFCompatible(PF,path)
path_dinv = dinvDyckPath(path,m)
tdinv     = tdinvParkingFunction(path,PF,m)
maxPF     = maxParkingFunction(path,m)
maxdinv   = tdinvParkingFunction(path,maxPF,m)
dinv      = dinvParkingFunction(path,PF,m)









