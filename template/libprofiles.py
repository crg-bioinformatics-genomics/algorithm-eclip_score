import numpy as np
import timeit
from scipy.optimize import curve_fit 
import Bio.SeqIO
import timeit
import random

AAs = "ARNDCQEGHILKMFPSTWYVX"



"""
    Curve fitting
"""
def sigmoidfunc(x, a, b, c, d):
    return a+b*np.tanh(c*x + d)    

def fit_curve(data):
    try:
        xs = data
        offset = -(xs[0] + xs[-1])/2
        ys = np.arange(data.size)/float(data.size)
        popt, pcov = curve_fit(sigmoidfunc, xs, ys, p0=(1,1,1,offset),maxfev=100000)
        return popt
    except RuntimeError, re:
        # probably not found a solution
        print "Warning - unable to find a fit during 100000 iterations"
        return np.array((0,0,0,0))


def calc_profile(profile):
    res = np.apply_along_axis(fit_curve, 1, profile)
    return res


def generate_seq(length):
    """
        Generates random sequence with X not being included
    """
    fun = lambda x: random.choice(AAs[:-1])
    return "".join(map(fun, xrange(length)))


"""
    Parse the table.txt containing scales we use, store them into an array
"""

def parse_scales(filename, AAs):
    """
        
    """
    with open(filename, "r") as f:
        lines = f.readlines()
        # stuff that is empty after strip is not really important for us, 
        # most likely an empty line
        lines = filter(lambda x: x.strip(), lines) 
    
    # variable definition
    titles = []
    no_properties = len(lines)/2
    # TODO - why are we using the extra line? ---removed now, let's see
    # property_table = np.zeros((no_properties+1, len(AAs)))
    property_table = np.zeros((no_properties, len(AAs)))
    
    # now go through the file and assign values
    for i, line in enumerate(lines):
        if i % 2 == 0:
            # scale title line
            titles.append(line.strip())
        else:
            # line of numbers
            numbers = map(float, line.split())
            property_table[i/2] = numbers + [0]
            # print numbers
            # print len(numbers)
    
    return (titles, property_table)



def score_proteins(records, AAs, property_table):
    """
        Applies scales to proteins
        Also converts the profile to uppercase for matching
    """
    scored = []
    for name,seq in records:
        indices =  map(lambda s: AAs.find(s), seq.upper())
        # TODO - what to assign to unknown chars?
        scored.append(property_table[:, indices])
    return scored


def scoremap_proteins(records, property_table):
    
    def scoreme(record):
        (name, seq) = record
        indices =  map(lambda s: AAs.find(s), seq)
        return property_table[:, indices]
    
    return map(scoreme, records)

def random_from_scored(scored_set, AAs, property_table):
    """
        Creates random based on what information it gets from
        already existing scored moments
    """
    extract_length = lambda item: item.shape[1]
    desired_lengths = map(extract_length, scored_set)
    
    # get the records now:
    records = [("gen", generate_seq(length)) for length in desired_lengths]
    return score_proteins(records, AAs, property_table)



def _timescore():
    score(records, property_table)

def _timescore2():
    scoremap(records, property_table)

def _timescorefuncs():
    t = timeit.Timer(_timescore)    
    print t.timeit(number=1000)
    
    t = timeit.Timer(_timescore2)    
    print t.timeit(number=1000)
    return


def _timewholething(records):
    res = map(smooth, score(records))
    t = timeit.Timer(timewholething)    
    print t.timeit(number=10)


def _time_smoothings(scored_records):
    test_mine = scored[0][1]    
    window_size = 7
    half_window = 3
    weight_step = 25
    
    def fun3():
        # proof-of-concept, very inefficient for real use
        x = []
        for i in range(half_window, len(test_mine) - half_window):
            tmp = 0
            for j in range(-half_window, half_window+1):
                tmp += (100-abs(j * weight_step))/float(100) * test_mine[i+j]
            tmp = tmp/float(window_size)
            x.append(tmp)
        return x
    
    
    multipliers = [0.25,0.5, 0.75, 1, 0.75,0.5, 0.25]
    def fun():
        arr = np.empty((window_size, len(test_mine)))
        for i in range(0, window_size):
            arr[i] = np.roll(test_mine, i) * multipliers[i]
        x = np.mean(arr, 0)[window_size-1:]
        return x
    
    def fun2():
        arr = np.empty((window_size, len(test_mine)))
        for i in range(0, window_size):
            arr[i] = np.roll(test_mine, i)
        x = np.average(arr, 0, multipliers)[window_size-1:]
        return x
    
    
    t = timeit.Timer(fun)    
    t.timeit(number=50000)    
    t = timeit.Timer(fun2)    
    t.timeit(number=50000) 
    t = timeit.Timer(fun3)    
    t.timeit(number=50000)
    # results  - fun and fun2 are roughly equivallent, fun3 is 10 times slower
 

def smooth(scale_values, window_size=7):
    # TODO - make the multipliers change according to the window size
    # multipliers = [0.25,0.5, 0.75, 1, 0.75,0.5, 0.25]
    
    one_half = window_size/2 + 1
    step = 1.0/one_half
    one_side = np.arange(0, 1, step).tolist()
    
    multipliers = one_side[1:] + [1] + one_side[::-1][:-1]
    
    # print multipliers
    
    # TODO - define the window sizes better
    scales, length = scale_values.shape
    length = length -window_size + 1
    try:
        # res = np.empty((scales, length))
        def _proc_scale(scale):        
            """
                Inner function smoothing individual scales
            """
            l = len(scale)
            arr = np.zeros((window_size, l))
            for i in range(0, window_size):
                # arr[i] = np.roll(scale, i) #old and slower way
                arr[i, i:] = scale[0:l-i]
            # print scale_values.shape, length, window_size
            return np.average(arr, 0, multipliers)[window_size-1:]
        try:
            return np.apply_along_axis(_proc_scale, 1, scale_values) 
        except FloatingPointError, fpe:
            print scale_values
            # import code; code.interact(local=locals())            
            raise
        
    except ValueError:
        # TODO - make explicit length check
        # probably the protein is too small
        return scale_values
        # import code; code.interact(local=locals())

def timesmooth(scored_values):
    def _timesmooth():
        map(smooth, scored)

    t = timeit.Timer(timesmooth)    
    print t.timeit(number=10)


MINIMUM_PROTEIN_LENGTH = 10
MAXIMUM_PROTEIN_LENGTH = 1000


def read_file_unprocessed(protfilename):
    """
        Enforces sizes and scores + smooths proteins
    """
    
    with open(protfilename, "r") as handle:
        data = Bio.SeqIO.parse(handle, "fasta")
        records_unfiltered = [(record.id, str(record.seq)) for record in data]
        # filter by length requirements
        records = filter_records(records_unfiltered)        
        return records
    


def read_score_smooth(AAs, property_table, protfilename):
    """
        Enforces sizes and scores + smooths proteins
    """
    
    with open(protfilename, "r") as handle:
        data = Bio.SeqIO.parse(handle, "fasta")
        records_unfiltered = [(record.id, str(record.seq)) for record in data]

        # filter by length requirements
        records = filter_records(records_unfiltered)
        
        scored = score_proteins(records, AAs, property_table)
        return (scored, map(smooth, scored))
    


def filter_records(records_unfiltered):
    return [record for record in records_unfiltered 
        if len(record[1])>=MINIMUM_PROTEIN_LENGTH and len(record[1]) <= MAXIMUM_PROTEIN_LENGTH]
    



# TODO - solve unknown/unexpected characters in the sequence and also solve capital/lowercase letters in the sequences
