def all_disjoint(sets):
    r'''
    Checks whether a collection of subsets is disjoint

    INPUT: set of sets
    OUTPUT: True/False

    '''
    union = set()
    for s in sets:
        for x in s:
            if x in union:
                return False
            union.add(x)
    return True



def noninv(sigma,k):
    r'''
    Computes the k-th noninversion number

    INPUT: permutation sigma and k
    OUTPUT: k-th noninvserion number of sigma

    '''
    N = sigma.noninversions(2)
    noninv = 0
    m = len(N)
    T = Subsets(m,k)
    for t in T.list():
        sets = [N[i-1] for i in t]
        if all_disjoint(sets):
            noninv += 1
    return noninv




def gamma(n,lamb,k):
    r'''
    Computes the eigenvalue of nu_k corresponding to lambda and k

    INPUT: an integer n, a vector lambda representing the parts of a partition, and another integer k.
    OUTPUT: an integer computed via Reiner-Saliola-Welker's eigenvalue formula

    EXAMPLES

    sage: gamma(6,[4,2],2)
    616


    '''
    lamb = Partition(lamb)
    rho = SymmetricGroupRepresentation(lamb)
    chi = rho.to_character()
    foo = 0
    for s in SymmetricGroup(n):
        sigma = Permutation(s)
        # print sigma, sigma.number_of_noninversions(k)
        foo += noninv(sigma,k) * chi(s)
    return foo
