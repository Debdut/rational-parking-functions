import multiprocessing

class RationalDyckPath:
    """ Stores an (m,n)-Dyck path, i.e. an integer lattice path in the first quadrant
        that can only make northward and eastward steps, starting at the origin and ending at the
        point (m,n), and staying strictly above the line y=n/m x.
    """
    def __init__(self, width = 3, height = 5, path = [0]):
        """ Initialization is (m, n, path) where path is a list where path[i] is the number of
            eastward steps taken before the ith northward step, i.e. the number of cells in row
            i (from the bottom) that lie directly west of the path.
        """
        self.width = width
        self.height = height
        self.path = path
        if not self.is_valid():
            self.path = [0]*self.height

    def is_valid(self):
        """ Determine whether or not a path is valid.  A path is invalid if it doesn't have length
            n or if the path does not stay strictly above the line y=n/m x.
        """
        if not self.height == len(self.path):
            return False
        for iter in range(self.height):
            if iter*self.width - self.path[iter]*self.height < 0:
                return False
        return True

    def minimal_dyck_path(self):
        """ Computes the (m,n)-Dyck path of minimal area. """
        slope = float(self.height/self.width)
        return RationalDyckPath(self.width, self.height, [floor(float(i/slope)) for i in range(self.height)])

    def max_area(self):
        """ Compute the maximal area of all Dyck paths. """
        min_path = self.minimal_dyck_path()
        max_a = 0
        for iter in range(self.height):
            max_a += min_path.path[iter]
        return max_a
    def area(self):
        """ Compute the area of the Dyck path. """
        if gcd(self.height,self.width) == 1:
            return ((self.height-1) * (self.width-1)/2 - sum(self.path))
        else:
            return self.max_area() - sum(self.path)
    def dinv(self):
        """ Compute the dinv of the Dyck path. """
        path_dinv = 0
        for iter in range(self.height):
            for jter in [1..self.path[iter]]:
                if self.dinv_characteristic(jter,iter+1):
                    path_dinv = path_dinv + 1
        return path_dinv

    def arm(self, x_coord, y_coord):
        """ Compute the arm of a square. """
        return self.path[y_coord-1]-x_coord
    def leg(self, x_coord, y_coord):
        """ Compute the arm of a square. """
        iter = 0
        while self.path[iter] < x_coord:
            iter += 1
        return y_coord - 1 - iter
    def dinv_characteristic(self, x_coord, y_coord):
        """ This function determines whether a square above
            the Dyck path should count toward computing the
            dinv of the path.
        """
        arm = self.arm(x_coord,y_coord)
        leg = self.leg(x_coord,y_coord)
        if gcd(self.width,self.height) > 1:
            return (self.height*arm < self.width*(leg+1) and self.width*leg <= self.height*(arm+1))
        else:
            return (self.height*arm < self.width*(leg+1) and self.width*leg < self.height*(arm+1))

    def ranks_on_path(self):
        """ Compute the ranks of the squares immediately to the
            right of the Dyck path, ordered from bottom to top.
        """
        divisor = gcd(self.width,self.height)
        k = self.width / divisor
        if divisor > 1:
            return [(-1 - self.path[i])*self.height + (1/divisor)*(1-(self.path[i]/k).floor()) + i*self.width for i in range(self.height)]
        else:
            return [(-1 - self.path[i])*self.height + i*self.width for i in range(self.height)]
    def rank_word(self):
        """ Compute the rank word of the Dyck path. """
        return list(sorted(self.ranks_on_path()))

    def prettyprint(self):
        """ Print the Dyck path in a nicely presented way. """
        for i in reversed(self.path):
            for j in range(i):
                print ".",
            for j in range(self.width-i):
                print "x",
            print "\n"
        print "\n"

    def to_affine(self):
        """ Return the m-increasing affine permutation associated to the
            Dyck path.  This only holds if m and n are coprime.  """
        if gcd(self.width,self.height) > 1:
            return -1
        rank_word = sorted(self.ranks_on_path())
        area = self.area()
        window = [rank + self.height - area + 1 for rank in rank_word]
        A = AffinePermutationGroup(['A',self.height-1,1])
        return A(window)

    def all_parking_functions(self, dinv=-1):
        """ Return a list of all parking functions compatible with the Dyck path.

            If dinv >= 0 the function will only return parking functions with the
            given dinv.
        """
        # Determine where the eastward steps in the Dyck path are at.
        right_steps = list()
        for i in range(self.height-1):
            if not self.path[i] == self.path[i+1]:
                right_steps.append(i)

        # Determine whether each permutation has a descent somewhere else.
        S_n = Arrangements([1..self.height], self.height).list()
        all_pf = list()
        for sigma in S_n:
            valid = True
            for des in Permutation(sigma).descents():
                if not des in right_steps:
                    valid = False
                    break
            pf = RationalParkingFunction(self.width,self.height,self.path,sigma)
            if valid and (dinv < 0 or dinv == pf.dinv()):
                all_pf.append(pf)
        return all_pf


# This class stores an (m,n)-Dyck path and a parking function, PF, where PF is a permutation
# that can only have descents where the Dyck path has rightward steps.

# The initialization input is (m, n, path, PF), but if a Dyck path is already defined the
# class can be initialized as (dyck='', PF='').
class RationalParkingFunction:
    """ Store an (m,n)-Parking Function, a permutation on an (m,n)-Dyck path that can only
        have descents where the Dyck path has eastward steps.

        The initialization input is (m, n, path, PF) where PF is a permutation.  If a Dyck
        path pi is already defined then the parking function can be initialized as
        RationalParkingFunction(dyck=pi, PF). """
    def __init__(self, width = 3, height = 5, path = [0], PF = [1], dyck = -1):
        if dyck == -1:
            self.dyck = RationalDyckPath(width, height, path)
        else:
            self.dyck = dyck
        self.PF = Permutation(PF)
        if not self.is_valid():
            PF = self.maximal_parking_function()

    def is_valid(self):
        """ Determine if the parking function only has decents where the Dyck path has
            eastward steps. """
        des_set = self.PF.descents()

        if not self.dyck.height == self.PF.size():
            return False
        for des in des_set:
            if self.dyck.path[ des ] == self.dyck.path[ des+1 ]:
                return False
        return True

    def maximal_parking_function(self):
        """ Compute the parking function with maximal dinv statistic.  This is the
            parking function that corresponds to the underlying Dyck path.
        """
        maximal_PF = list()
        ranks = self.dyck.ranks_on_path()
        for rk in ranks:
            maximal_PF.append(sorted(ranks).index(rk)+1)
        return maximal_PF

    def to_permutation(self):
        """ Compute the S_n permutation of the parking function.

            Each parking function can be realized as an ordering of the ranks immediately
            right of the underlying Dyck path.  The ordering of the ranks is the
            S_n permutation of the parking function. 
        """
        # First compute the ranks of the cars in the parking function
        car_ranks = self.dyck.ranks_on_path()

        # Beginning with the largest rank and working towards the
        # smallest, write down the cars in the PF
        sigma = []
        for car in list( reversed( sorted( car_ranks ))):
            sigma.append(self.PF[ car_ranks.index(car) ])
        return Permutation(sigma)

    def to_affine(self):
        """ Return the corresponding m-increasing affine permutation. """
        # First compute the ranks of the cars in the parking function
        car_ranks = self.dyck.ranks_on_path()

        # Then we reorder based on the parking function
        sigma = list()
        for car in range(self.dyck.height):
            sigma.append(car_ranks[ self.PF.index(car+1) ])

        # Finally we add the constant '+ height - area + 1' to each rank
        window = [car + self.dyck.height - self.dyck.area() + 1 for car in sigma]
        A = AffinePermutationGroup(['A',self.dyck.height-1,1])
        return A(window)

    def dinv(self):
        """ Compute the dinv of the parking function. """
        sigma = self.to_affine()

        # Compute the number of m-bounded inversions within the window of
        # the affine permutation.
        m_inversions = 0
        for i in sigma[:-1]:
            for j in sigma[sigma.index(i)+1:]:
                if i > j and j + self.dyck.width > i:
                    m_inversions = m_inversions + 1

        # Then subtract this number from the dinv of the underlying Dyck path.
        return self.dyck.dinv() - m_inversions
    def area(self):
        """ Compute the area of the parking function. """
        return self.dyck.area()

    def prettyprint(self):
        """ Print the parking function in a nice way. """
        for i in range(self.dyck.height):
            steps = list(reversed(self.dyck.path))[i]
            for j in range(steps):
                print ".",
            print "%s" % list(reversed(self.PF))[i],
            for j in range(self.dyck.width-steps-1):
                print "x",
            print "\n"
        print "\n"

def all_rational_dyck_paths(width, height, above_path = -1, dinv = -1, area = -1, skip = -1):
    """ Compute all of the (m,n)-Dyck paths and return them as a list.

        If parameters such as above_path, dinv, and area are specified then the function
        will only return the Dyck paths that satisfy those parameters.
        above_path:  This will only return Dyck paths that lie above the specified
                     path (as a list of cells west of the path from bottom to top)
        dinv:        This will only return Dyck paths that have the stated dinv statistic
        area:        This will only return Dyck paths that have the stated area statistic
        skip:        This will only return Dyck paths that have the stated skip statistic,
                     recall skip + area + dinv = (width-1)(height-1)/2
    """
    if above_path == -1 or not RationalDyckPath(width,height,above_path).is_valid():
        above_path = RationalDyckPath(width,height,[0]*height).minimal_dyck_path().path
    partition = list(reversed(above_path))
    max_area  = sum(above_path)
    parts     = list()
    all_paths = list()
    for i in range(max_area+1):
        parts = parts + Partitions(i, outer=list(reversed(above_path))).list()
    for p in parts:
        new_p = list(p)
        while len(new_p) < height:
            new_p.append(0)
        new_path = RationalDyckPath(width,height,list(reversed(new_p)))
        if (dinv < 0 or new_path.dinv() == dinv) and (area < 0 or new_path.area() == area) and (skip < 0 or new_path.dinv()+new_path.area()+skip == (width-1)(height-1)/2):
            all_paths.append(new_path)
    return all_paths

def min_row(width,height):
    for iter in [2..height]:
        if (iter-1)*width % height == 0:
            return iter
    return -1
def all_modified_rational_dyck_paths(width, height):
    g = gcd(width,height)
    if g == 1:
        return all_rational_dyck_paths(width,height)
    p = RationalDyckPath(width,height).minimal_dyck_path()
    p.path[min_row(width,height)] -= 1
    return all_rational_dyck_paths(width,height,p.path)

def all_m_increasing_affine_permutations(m, window_size, area=-1, dinv=-1):
    """ Compute all of the m-increasing affine permutations in ~S_n.
        Return them as a list of permutations in window format.  If the parameters have non-negative
        values then only the affine permutations that satisfy those conditions will be returned.
        This is only valid if m and window_size are relatively prime.
    """
    if gcd(m,window_size) > 1:
        return -1
    all_permutations = list()

    all_paths = all_rational_dyck_paths(m, window_size)
    for path in all_paths:
        if area < 0 or path.area() == area:
            all_permutations = all_permutations + [pf.to_affine() for pf in path.all_parking_functions(dinv)]
    return all_permutations

def rational_qt_catalan(m,n):
    """ Compute the (m,n)-rational q,t-Catalan polynomial over Q. """
    Cmn = 0
    R = QQ['q','t'].fraction_field()
    (q,t) = R.gens()

    paths = all_rational_dyck_paths(m,n)
    for current_path in paths:
        Cmn = Cmn + t^current_path.area() * q^current_path.dinv()
    return Cmn
def modified_rational_qt_catalan(m,n):
    """ Compute the modified (m,n)-rational q,t-Catalan polynomial over Q. """
    Cmn = 0
    R = QQ['q','t'].fraction_field()
    (q,t) = R.gens()

    paths = all_modified_rational_dyck_paths(m,n)
    for current_path in paths:
        Cmn = Cmn + t^current_path.area() * q^current_path.dinv()
    return Cmn/t
def c(m,n):
    return rational_qt_catalan(m,n)
def mc(m,n):
    return modified_rational_qt_catalan(m,n)

def hikita_polynomial(m, n, expand=-1, term=-1):
    """ Compute the Hikita polynomial over Q for m,n coprime.
        expand: By default the Hikita polynomial will be presented in terms of symmetric or quasisymmetric
                bases.  If expand is set, then the polynomial will be express in 'expand' indeterminates.
        term:   This returns only the coefficient of the term F_term.
    """
    # Set the underlying framework of symmetric and quasisymmetric functions
    R = QQ['q','t'].fraction_field()
    (q,t) = R.gens()
    Fqt = QuasiSymmetricFunctions(R).F()

    # Compute the statistics on all parking functions for every Dyck path
    HP = 0
    dyck_paths = all_rational_dyck_paths(m,n)
    for path in dyck_paths:
        all_pfs = path.all_parking_functions()
        for pf in all_pfs:
            composition = Compositions().from_descents( pf.to_permutation().idescents(), n )
            if term == composition:
                HP = HP + t^pf.area() * q^pf.dinv()
            if term == -1:
                HP = HP + t^pf.area() * q^pf.dinv() * Fqt[composition]
    if not expand == -1:
        HP = HP.expand(expand)
    return HP

R = QQ['q','t'].fraction_field()
(q,t) = R.gens()
F = QuasiSymmetricFunctions(R).F()
s = SymmetricFunctions(R).s()


print "starting..."
print
h = s(hikita_polynomial(4,3).to_symmetric_function())
g = s(hikita_polynomial(1,3).to_symmetric_function())
print h-g.nabla()
print
print s([3]).nabla()
print s([2,1]).nabla()
print s([1,1,1]).nabla()
print
print "... done!"
#print hikita_polynomial (10,4)










