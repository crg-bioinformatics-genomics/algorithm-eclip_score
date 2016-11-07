import numpy as np
import operator as op
import itertools as iter
import forkmap
import sys
import os
import scipy
import scipy.stats.stats as stats
from pylab import mean  

"""
Functions for the properties analysis, only one setting here - how many threads to use
"""

# this ensures that maximum of four threads run but also at most one per processor
numthreads = 4
perproc = 1


def print_simple_result(percentages, property_list, signs, NO_OF_PROPERTIES, count=20):
    items = sorted(zip(percentages, property_list, signs, range(0, NO_OF_PROPERTIES)))[::-1]
    for item in items[0:count]:    
        outstr = "%s %s %0.3f" % (
                                 item[2], 
                                 item[1].strip() ,
                                 item[0]
                             )
        print outstr

# now let's do the scale combination!
def property_decision_moments(res_1, res_2, THRESHOLD):
    
    
    elements_in_A, NO_OF_PROPERTIES, NO_OF_MOMENTS = res_1.shape
    elements_in_B = res_2.shape[0]
    
    coverages = []
    signs = []
    property_indices = dict()
    for scale in range(NO_OF_PROPERTIES):
        vals = [
                res_1[:,scale]>THRESHOLD,       # dataset1 - positive selection on 1
                res_1[:,scale]<=(1-THRESHOLD),  # dataset1 - negative selection on 1
                res_2[:,scale]<=(1-THRESHOLD),  # dataset2 - positive selection on 1
                res_2[:,scale]>THRESHOLD        # dataset2 - negative selection on 1
        ]
        mprcs = np.array(map(lambda el: np.mean(el, axis=0), vals))
        combinations = []
        indices_A = np.zeros(elements_in_A, np.dtype('bool'))
        indices_B = np.zeros(elements_in_B, np.dtype('bool'))
        for moment_in_question in range(NO_OF_MOMENTS):
            # positive 1
            larger_1 = sum(mprcs[:,moment_in_question][[0,2]] * \
                    (elements_in_A, elements_in_B)) / float(elements_in_A + elements_in_B)
            # print "(+) Larger 1:", larger_1
        
            # positive 2
            larger_2 = sum(mprcs[:,moment_in_question][[1,3]] * \
                    (elements_in_A, elements_in_B)) / float(elements_in_A + elements_in_B)
            # print "(-) Larger 2", larger_2
            # print ""
            if larger_1>larger_2:
                indices_A = np.logical_or(res_1[:,scale, moment_in_question]>THRESHOLD, indices_A)
                indices_B = np.logical_or(res_2[:,scale, moment_in_question]<=(1-THRESHOLD), indices_B)                
                combinations.append("+")
            else:
                indices_A = np.logical_or(res_1[:,scale, moment_in_question]<=(1-THRESHOLD), indices_A)
                indices_B = np.logical_or(res_2[:,scale, moment_in_question]>THRESHOLD, indices_B)
                combinations.append("-")
        signs.append(combinations)
        coverages.append(
            (np.mean(indices_A) * elements_in_A +  np.mean(indices_B) *\
            elements_in_B)/float(elements_in_A + elements_in_B)
        )
        property_indices[scale] = [np.where(indices_A)[0], np.where(indices_B)[0]]
    
        # print scale
    return (coverages, signs, property_indices)
    
    # tops_new = []
    # top_new = []
    # for item in sorted(zip(coverages, property_list, signs, range(0, NO_OF_PROPERTIES)))[::-1][0:30]:    
    #     outstr = "%s %s %0.3f" % (
    #                              ("+" if item[2] else "-"), 
    #                              item[1].strip() ,
    #                              item[0]
    #                          )
    #     print outstr, "".join((60 - len(outstr)) * [" "]), " ".join(item[2])
    #     tops_new.append(outstr)
    #     top_new.append(item[3])
    #     # print top


def property_decision(res_1, res_2, THRESHOLD, using_both = True):
    
    _, NO_OF_PROPERTIES, __ = res_1.shape
    
    new_prc = list()
    new_signs = list()
    property_indices = dict()
    for i in range(0, NO_OF_PROPERTIES):
        vals = [
                    res_1[:,i, 0]>THRESHOLD,       # dataset1 - positive selection on 1
                    res_1[:,i, 0]<=(1-THRESHOLD),  # dataset1 - negative selection on 1
                    res_2[:,i, 0]<=(1-THRESHOLD),  # dataset2 - positive selection on 1
                    res_2[:,i, 0]>THRESHOLD        # dataset2 - negative selection on 1
                ]            
        # which of the vals is the best?
        mprcs = map(mean, vals)
        if using_both:
            dec_prcs = [(mprcs[0] + mprcs[2])/2.0, (mprcs[1] + mprcs[3])/2.0]
        else:
            dec_prcs = mprcs
        best = max(dec_prcs)
        offset = dec_prcs.index(best)%2
        indices_1 = np.where(vals[offset])[0]
        indices_2 = np.where(vals[offset+2])[0]    
        property_indices[i] = [indices_1, indices_2]   
        new_prc.append(best)    
        new_signs.append(dec_prcs.index(best) in [0,2]) 
    # la.create_matrix("")
    new_signs = map(lambda x: "+" if x else "-", new_signs)    
    return new_prc, new_signs, property_indices



def compare_datasets(momentsA, momentsB):
    # infer initialise lengths and parameters
    elements_in_A, NO_OF_PROPERTIES,  NO_MOMENTS = momentsA.shape
    elements_in_B = momentsB.shape[0]
    
    
    # dataset 1
    res_1 = np.zeros((elements_in_A, NO_OF_PROPERTIES, NO_MOMENTS)) 
    for i in range(elements_in_A):
        # np.average because we are essentially counting no. of ones    
        res_1[i] = np.array(
                        [
                            np.average(momentsA[i, :, 0] > momentsB[:, :, 0], axis=0),
                            np.average(momentsA[i, :, 1] > momentsB[:, :, 1], axis=0),
                            np.average(momentsA[i, :, 2] > momentsB[:, :, 2], axis=0),
                            np.average(momentsA[i, :, 3] > momentsB[:, :, 3], axis=0)
                        ]
                    ).transpose() # and change the shape as we like
    
    
    
    # dataset 2
    res_2 = np.zeros((elements_in_B, NO_OF_PROPERTIES, NO_MOMENTS)) 
    for i in range(elements_in_B):
        res_2[i] = np.array(
                        [
                            np.average(momentsB[i, :, 0] > momentsA[:, :, 0], axis=0),
                            np.average(momentsB[i, :, 1] > momentsA[:, :, 1], axis=0),
                            np.average(momentsB[i, :, 2] > momentsA[:, :, 2], axis=0),
                            np.average(momentsB[i, :, 3] > momentsA[:, :, 3], axis=0)                                                                        
                        ]
                    ).transpose()
    
    return res_1, res_2


def generate_moment(dataset, NO_OF_PROPERTIES, NO_MOMENTS):
    element_count = len(dataset)
    moments = np.zeros((element_count, NO_OF_PROPERTIES, NO_MOMENTS))    
    # TODO debugging here only
    for row in range(element_count):
        moments[row, :, :] = np.array([
            scipy.mean(dataset[row][0:NO_OF_PROPERTIES,:], axis=1),
            # scipy.mean(dataset[row][0:NO_OF_PROPERTIES,:], axis=1),
            # scipy.mean(dataset[row][0:NO_OF_PROPERTIES,:], axis=1),
            # scipy.mean(dataset[row][0:NO_OF_PROPERTIES,:], axis=1),
            scipy.std(dataset[row][0:NO_OF_PROPERTIES,:], axis=1),        
            stats.skew(dataset[row][0:NO_OF_PROPERTIES,:], axis=1),
            stats.kurtosis(dataset[row][0:NO_OF_PROPERTIES,:], axis=1)
        ]).transpose()
    return moments
    

def generate_moments(sortedA, sortedB, NO_OF_PROPERTIES, NO_MOMENTS):    
    # NOT DRY - TODO separate the function to only process a single dataset at a time
    # elements_in_A = len(sortedA)
    #     elements_in_B = len(sortedB)
    #     
    #     # process the dataset A
    #     momentsA = np.zeros((elements_in_A, NO_OF_PROPERTIES, NO_MOMENTS))    
    #     for row in range(elements_in_A):
    #         momentsA[row, :, :] = np.array([
    #             scipy.mean(sortedA[row][0:NO_OF_PROPERTIES,:], axis=1),
    #             scipy.std(sortedA[row][0:NO_OF_PROPERTIES,:], axis=1),        
    #             stats.skew(sortedA[row][0:NO_OF_PROPERTIES,:], axis=1),
    #             stats.kurtosis(sortedA[row][0:NO_OF_PROPERTIES,:], axis=1)
    #         ]).transpose()
    #     
    #     # process the second dataset
    #     momentsB = np.zeros((elements_in_B, NO_OF_PROPERTIES, NO_MOMENTS))    
    #     for row in range(elements_in_B):
    #         momentsB[row, :, :] = np.array([
    #             scipy.mean(sortedB[row][0:NO_OF_PROPERTIES,:], axis=1),
    #             scipy.std(sortedB[row][0:NO_OF_PROPERTIES,:], axis=1),
    #             stats.skew(sortedB[row][0:NO_OF_PROPERTIES,:], axis=1),
    #             stats.kurtosis(sortedB[row][0:NO_OF_PROPERTIES,:], axis=1)
    #         ]).transpose()        
    momentsA = generate_moment(sortedA, NO_OF_PROPERTIES, NO_MOMENTS)
    momentsB = generate_moment(sortedB, NO_OF_PROPERTIES, NO_MOMENTS)
    
    return momentsA, momentsB





def every_third(mylist): 
    # TODO - to be re-coded using standard list operations [::3]
    # nice as a historic insight - what a bunch of inefficiency :)
    return map(
        lambda i: mylist[i],
        filter(
            lambda i: i%3 == 0,
            range(
                len(mylist)
            )
        )
    )

def get_path():
    path = os.path.dirname(__file__)
    sep = "/" if path else ""
    return path+sep

# TODO change the location!!
def property_table():    
    
    is_ascii = lambda let: ord(let) < 128
    
    with open(get_path()+"./scales/scales_generated.txt", "rU") as propfile:
        # replace typographical dashes to other type
        properties = map(
            lambda x: x.strip().replace("\xe2\x80\x93", "-"), 
            every_third(propfile.readlines())
        )                    
        return map(lambda i: str(i)+" "+ filter(is_ascii, properties[i]), range(0, len(properties)))    

def no_properties():
    return len(property_table())

def get_property_list():
    """
        Fetches the chosen property list from the settings subdirectory of the script
    """
    return property_table()


property_list =  get_property_list()
def list_winner(indices, cmb, prc, q_inds, groups,single_prc, signs):
    for i in indices:
        print "Combination: %s yields %0.8f percentage, consistency is %f, QI: %f" % (str(cmb[i]), prc[i], (float(q_inds[i])*float(q_inds[i]))/prc[i], q_inds[i])
        for key, val in groups[i].iteritems():
            props, mean, std =  val
            if props:
                print "%s: %0.3f %0.3f" % (key, mean, std)
                for prop in list(props):
                    print "%s %s %0.3f" % (
                                             ("+" if signs[prop] else "-"), 
                                             property_list[prop].strip() ,
                                             single_prc[prop]
                                         )
        print "-" * 20   
        print 
        # print 

        # if isinstance(cmb[i], int):
        #     items = [cmb[i]]
        # else:
        #     items = cmb[i]
        # items = sorted(items, key=lambda prop: (signs[prop], single_prc[prop]), reverse=True)
        # for prop in items:
        #     print "%s %s %0.3f" % (
        #                             ("+" if signs[prop] else "-"), 
        #                             property_list[prop].strip() ,
        #                             single_prc[prop]
        #                         )


from heapq import nlargest
def which_largest(perc, cmb, one_prc, count):
    """
        This function analyses the percentages and selects top 10 results
        It also sortes recursively, if the coverage is the same, then it is sorted by the individual 
        property discriminations
    """
    def comparison(index):
        return (perc[index],sorted(map(lambda j: one_prc[j], cmb[index])))
    
    return nlargest(count, range(0, len(perc)), comparison)


def limited_comb_expand(prcs, combs, chosen, limit):
    """
        This function first selects top 'limit' of combinations
        and then proceeds to combine each of them with new property
        from the list of available ones (variable chosen)
    """
    # selects ITERATION_LIMIT best entries to improve performance
    items = sorted(zip(prcs, combs))[-limit:]
    return [l[1] + [i] for i in chosen for l in items if i not in l[1]] 


def union_properties(which_set, comb, property_indices):
    """
        For a given set (0 or 1) and combination of properties, this 
        function produces a set containin union of all the selected 
        proteins for all of the properties
        
        each of the unions is calculated on the fly, which could be
        in future sped up by storing the unions for each of the combs
    """
    return \
        reduce(
            lambda cur_set, item: \
                cur_set.union(property_indices[item][which_set]), 
            comb,   # applying accross all elements of the combination
            set()   # empty set as an initialisator
        )




def get_chosen():
    """
        Fetches the chosen property list from the settings subdirectory of the script
        The file is in old format, converted to zero-based indexing now
    """
    with open(get_path()+"settings/chosen.txt", "rU") as chosenfile:
        chosen = chosenfile.readlines()
        return map(lambda x: int(x)-2, chosen)



# TODO completely redundant - we can just query a list with the ::2 etc.
def extract_odd(mylist): 
    return map(
        lambda i: mylist[i],
        filter(
            lambda i: i%2 == 1,
            range(
                len(mylist)
            )
        )
    )



def all_indices(value, qlist):
    indices = []
    idx = -1
    while 1:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices

def average(values):
    try:
        return sum(values)/len(values)
    except ZeroDivisionError:
        return 0



def create_matrix(filename, NO_OF_PROPERTIES):
    # extract numbers (or more precisely cast them properly)    
    extract_numbers = lambda x: map(float, x.split()[0:NO_OF_PROPERTIES*8])

    with open(filename, "rU") as datasetfile:
        dataset = datasetfile.readlines()

    return  np.array(
                map(
                    extract_odd,
                    map(extract_numbers, dataset)
                )
            ).reshape(len(dataset),NO_OF_PROPERTIES,4)

def ppmap_calc_chosen_wrapper((u, info)):
    """
        PP has special way of handling function and simple calc_chosen needs to be pre-populated beforehand
    """
    # (calc_chosen, num_matrix1, num_matrix2, valid_entries1, valid_entries2, (THRESHOLD, DATA_THRESHOLD)) = info
    calc_chosen = info[0]
    partial = calc_chosen(*info[1:])
    return partial(u)



# TODO write a unit test for this one
# very old and obsolete code now!
def calc_chosen(num_matrix1, num_matrix2, valid_entries1, valid_entries2, (THRESHOLD, DATA_THRESHOLD)): 
    """
        
        Obsolete now
    
    """
    # using absolute values as measures of signal strength - how to do it otherwise?
    # it sure does not make sense to use comparison between positive/negative numbers
    biggerf = lambda (ai,y):    map(
                                    lambda (r,s): r > s, 
                                    zip(ai, y)
                                )[0]
    
    smallerequalf = lambda (ai,y):  map(
                                        lambda (r,s): r <= s, 
                                        zip(ai, y)
                                    )[0]
                                        
    # different signs   
    # biggerf = lambda (ai,y): ai > y and np.sign(ai) != np.sign(y)
    # smallerequalf = lambda (ai,y): ai <= y and np.sign(ai) != np.sign(y)
    
    # biggerf = lambda (ai,y): abs(ai) > abs(y)
    # smallerequalf = lambda (ai,y): abs(ai) <= abs(y)
    
     
    # TODO re-write to make more clear and more efficient!
    # essentially does each element from a with each element from b comparison
    compfun = lambda a,b,compf: map( 
                                lambda ai:  sum(
                                                map(
                                                    compf, 
                                                    zip([ai] * len(b),b)
                                                )
                                            )/float(len(b)),
                                a
                            ) 
    
    
    # the tag is somehow obsolete - forkmap did not perform very well
    @forkmap.parallelizable(numthreads, perproc)      
    def fnc(u):
        # selecting the first column of the dataset1 and 2 and make them into a list    
        # TODO change the 0: selector to only include valid ids in the future
        # TODO comment why is this u-2?
        l = num_matrix1[list(valid_entries1),u-2].tolist()
        m = num_matrix2[list(valid_entries2),u-2].tolist()                                              
        
        try:
            data_1 = compfun(l,m, biggerf)
            data_2 = compfun(m,l, smallerequalf)
        except ZeroDivisionError:
            # print "Skipping property: "+str(u)
            return False
        
        ### Positive selection on dataset one
        # select positive set
        # data_1_pos = all_indices(True, map(lambda x: x > THRESHOLD, data_1))
        
        
        
        # Positive selection on dataset two
        # data_1_indices = all_indices(True, map(lambda x: x > THRESHOLD, data_2))            
        
        # TODO rename the variables here
        # is the calculated average larget than our THRESHOLD?
        data_1_indices = all_indices(True, map(lambda x: x > THRESHOLD, data_1))
        # convert back to global coordinates
        
        if data_1_indices:
            data_1_indices = op.itemgetter(*data_1_indices)(list(valid_entries1))
            # map(lambda index: ids1[index], data_1_indices_conv)            
            if type(data_1_indices) == int:
                data_1_indices = [data_1_indices]
        
        # is the 
        data_2_indices = all_indices(True, map(lambda x: x<= 1-THRESHOLD, data_2))
        
        if data_2_indices:        
            # convert back to global coordinates        
            data_2_indices = op.itemgetter(*data_2_indices)(list(valid_entries2))
            
            # sometimes the itemgetter returns an int of there is only a single entry
            if type(data_2_indices) == int:
                  data_2_indices = [data_2_indices]
        
        ps = average(data_1) >= 0.5
        
        # .. and now the other way around
        # 
        data_3_indices = all_indices(True, map(lambda x: x<= 1-THRESHOLD, map(lambda x: 1-x,data_1)))        
        data_4_indices = all_indices(True, map(lambda x: x > THRESHOLD, map(lambda x: 1-x, data_2)))
        
        if data_3_indices:        
            # convert back to global coordinates        
            data_3_indices = op.itemgetter(*data_3_indices)(list(valid_entries1))
            
            # sometimes the itemgetter returns an int of there is only a single entry
            if type(data_3_indices) == int:
                  data_3_indices = [data_3_indices]
        
        if data_4_indices:        
            # convert back to global coordinates        
            data_4_indices = op.itemgetter(*data_4_indices)(list(valid_entries2))
            
            # sometimes the itemgetter returns an int of there is only a single entry
            if type(data_4_indices) == int:
                 data_4_indices = [data_4_indices]
        
        ps2 = not average(map(lambda x: 1-x, data_2)) >= 0.5
        
        
        # is there enough signal in the sets to justify selection?
        # TODO - this needs to be re-done when elements of l are composite vectors now
        # if np.mean(l) < DATA_THRESHOLD:
        #     data_1_indices = []
        #     data_4_indices = []        
        # if np.mean(m) < DATA_THRESHOLD:
        #     data_2_indices = []
        #     data_3_indices = []
        
        # 
        # discriminated_one = len(list(data_1_indices)) + len(list(data_2_indices))
        # discriminated_two = len(list(data_3_indices)) + len(list(data_4_indices))
        
        # also think - how do we determine the best discrimination.. maybe percentages?
        # 
        discriminated_one = max(
            len(list(data_1_indices))/float(len(valid_entries1)), 
            len(list(data_2_indices))/float(len(valid_entries2))
        )
        discriminated_two = max(
            len(list(data_3_indices))/float(len(valid_entries1)), 
            len(list(data_4_indices))/float(len(valid_entries2))
        )
        
        # prop.txt equivalent now contains four values:
        prop_data = (average(data_1), u, discriminated_one, ps, 0)
        prop_data2 = (average(map(lambda x: 1-x,data_1)), u, discriminated_two, ps2, 1)
        # properties.append(prop_data)
        # chosen_indices[u] = ((data_1_indices, data_2_indices))
        # property_dict[u] = ps
        return [(prop_data, ((data_1_indices, data_2_indices)), ps),(prop_data2, ((data_3_indices, data_4_indices)), ps2)]
    return fnc



def equalsigns(signs):
    return reduce(
                lambda (x, y), r: (r,x==r and y),
                signs, 
                (signs[0], True)
            )[1]

    
def filter_equallist(proplist, properties, history):
    """ 
        Looks for sings of properties in proplist and first
        checks if they are all equal, if not, it raises exception.
        If they are, it then proceeds and filters the properties list
        removing all the entries where the sign does not match
    """
    signs = [item[1] for item in history if item[0] in proplist]
    if signs:  
        # consistency check
        if len(signs) > 1 and not equalsigns(signs): 
            raise Exception('Invalid history! ('+str(proplist)+')')
        
        # actual filter function
        cur_filter = lambda (_, id, __, ps, ___):                             \
                False if (id in proplist and signs[0] is not ps)         \
                else True
        properties = filter(cur_filter, properties)
    return properties

    
    