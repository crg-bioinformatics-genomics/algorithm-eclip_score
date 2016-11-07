import libanalysis as la
import numpy as np
import yaml

# hardcoded init values
def property_group_init(NO_OF_PROPERTIES, order=False):
    
    # import code; code.interact(local=locals())
    
    group_filename = "./scales/scales_groups_ordered.txt" if order else "./scales/scales_groups.txt"
    group_file_path = os.path.join(get_path(), group_filename)
    with open(group_file_path, "r") as propfile:
        property_groups_loaded = yaml.load(propfile)
    
    # TODO IMPORTANT - the use of negatively-contributing properties is obsolete now
    property_groups = dict()
    for key, values in property_groups_loaded.iteritems():
        property_groups[key] = (
            key,
            values,
            []
        )
    
    
    # property_groups = dict()
    # property_groups["hydrophobicity"] = (
    #     "Hydrophobicity/polarity", 
    #     set(range(0, 16+1)),     #hydrophobicity scales, 16 is included, hence the +1   
    #     set([34,35])             #polarity scales
    # )
    # 
    # property_groups["alpha"] = (
    #     "Alpha proteins",
    #     set([21,22,23]),
    #     set()
    # )
    # 
    # property_groups["beta"] = (
    #     "Beta proteins - sheet",
    #     set([24,25,26]),
    #     set()
    # )
    # 
    # property_groups["beta-turn"] = (
    #     "Beta proteins - turn",
    #     set([27,28,29]),
    #     set()
    # )
    # 
    # property_groups["disorder"] = (
    #     "Disorder",
    #     set([17,19, 30, 31, 32, 20]),     # TODO - check if all has a good sign
    #     set()
    # )


    
    # property_colours = {
    #     "alpha+" : ("red", ""),
    #     "beta+" : ("green", ""),
    #     # "beta-turn+" : ("lime", ""),
    #     "membrane+" : ("orange", ""),
    #     "hphob+" : ("purple", ""),
    #     "alpha-" : ("red", "\\"),
    #     "beta-" : ("green", "\\"),
    #     # "beta-turn-" : ("lime", "\\"),
    #     "membrane-" : ("orange", "\\"),
    #     "hphob-" : ("purple", "\\")
    # }

    property_all = dict()
    for key, val in property_groups.iteritems():
        for v in list(val[1]):
            # positive stuff
            property_all[v] = key + "+"
        for v in list(val[2]):
            property_all[v] = key + "-"
        
    

    all_properties = set(range(0, NO_OF_PROPERTIES))
    in_groups = reduce(
        lambda x, y: x.union(y),                    # union of the elements
        map(
            lambda x: set(x[1]).union(set(x[2])),   # select the positives and negatives from the property tuple
            property_groups.values()                
        ), 
        set()                                       # starting with the empty set
    )



    ungrouped_properties = all_properties.difference(in_groups)

    disallowed_properties = set()
    
    return property_groups, ungrouped_properties


# print "Ungrouped properties:"
# print ungrouped_properties
# 17, 19, 20, 30, 31, 32 - disorder
# also add burial to the hydrophobicity


# getting the percentages and signs data (the signs data is a bit redundant now)



def print_initial_consistency(bool_signs, percentages, weighted_consistency = True):
    signs_plusminus = np.array(map(lambda x: 1 if x else -1, bool_signs))
    prc_arr = np.array(percentages)
    
    
    for key, val in property_groups.iteritems():
        print "-" * 20    
        print key
        title, pos, neg = val
        if weighted_consistency:
            values = np.append((signs_plusminus * prc_arr)[list(pos)], (-1*signs_plusminus * prc_arr)[list(neg)])
        else:
            values = np.append((signs_plusminus * 1)[list(pos)], (-1*signs_plusminus * 1)[list(neg)])        
        print "Consistency: %0.3f, consensus sign: %d" % (np.abs(np.mean(values)), np.sign(np.mean(values)))
        print "Standard deviation: " + str(np.std(values))
        
        if np.abs(np.mean(values)) < np.std(values):
            print "a bit fishy this one!"
            # disallowed_properties = disallowed_properties.union(set(pos))
            # disallowed_properties = disallowed_properties.union(set(neg))
        # print np.array(new_prc)[pos]
        # print np.array(new_signs)[neg]



# TODO comment on this one
increase_threshold = lambda THRESHOLD: (1-THRESHOLD)/2 + THRESHOLD

import os
def get_path():
    path = os.path.dirname(__file__)
    sep = "/" if path else ""
    return path+sep


# deprecated!
def handle_combinations(res_1, res_2, combinations_calculate=3, THRESHOLD=0.75, weighted_consistency=False, ITERATION_LIMIT=500, combined=False):        
        
    
    
    elements_in_A, NO_OF_PROPERTIES, NO_OF_MOMENTS = res_1.shape
    elements_in_B = res_2.shape[0]
    
    property_groups, ungrouped_properties = property_group_init(NO_OF_PROPERTIES)
        
    if combined:
        def decision_fun(x, y, z):
            a, b, c = la.property_decision_moments(x, y, z)
            b = map(lambda x: x[0] == "+", b)
            return a,b,c
    else:
        decision_fun = la.property_decision
        
    percentages, bool_signs, property_indices = decision_fun(res_1, res_2, THRESHOLD)
    
    
    prc_arr = np.array(percentages)
    signs_plusminus = np.array(map(lambda x: 1 if x else -1, bool_signs))
    # nested break implementation
    from contextlib import contextmanager
    @contextmanager
    def nested_break():
        class NestedBreakException(Exception):
            pass
        try:
            yield NestedBreakException
        except NestedBreakException:
            pass
    
    
    all_properties = range(0, NO_OF_PROPERTIES)
    
    # initialise data stores
    prcs = dict()
    combs = dict()
    winners = dict()
    quality_indices = dict()
    groups_stats = dict()
    groups_consistencies=dict()
    
    # initialise first round
    combs[1] = map(lambda x: [x], all_properties)
    
    # allow breaking from a nested loop
    with nested_break() as mylabel:
        
        for num_prop in range(1, combinations_calculate+1):
            
            # TODO redo
            # a very convoluted way to make combinations unique
            # sorts the values and converts them to a string
            # then only takes unique strings
            combs[num_prop] = sorted(
                map(
                    lambda el: map(int, el),
                    map(
                        lambda x: x.split(","), 
                        np.unique(
                            map(
                                lambda lst: ",".join(map(str, sorted(lst))), 
                                combs[num_prop]
                            )
                        )
                    )            
                )
            )
            
            calc_percentage = lambda (set1, set2): max(len(set1)/float(elements_in_A), len(set2)/float(elements_in_B))
            calc_index = lambda (set1, set2, consistency): sqrt(max(len(set1)/float(elements_in_A), len(set2)/float(elements_in_B)) * consistency)
            
            
            # now calculate consistencies for the various groups            
            groups_stats[num_prop] = []
            groups_consistencies[num_prop] = []
            
            # for each combination, determine how consistent it is
            for comb in combs[num_prop]:
                consistencies = []
                scomb = set(comb)
                result = dict()
                for key, val in property_groups.iteritems():
                    title, pos, neg = val                
                    cur_pos = pos.intersection(scomb)
                    cur_neg = neg.intersection(scomb)
                    
                    if cur_pos or cur_neg:   
                        # calculate the consistency
                        if weighted_consistency:
                            values = np.append(
                                (signs_plusminus * prc_arr)[list(cur_pos)], 
                                (-1*signs_plusminus * prc_arr)[list(cur_neg)]
                            )
                        else:
                            values = np.append(
                                 (signs_plusminus * 1)[list(cur_pos)], 
                                 (-1*signs_plusminus * 1)[list(cur_neg)]
                             )
                        # print "Values:"
                        # print values                
                        
                        result[key] = (cur_pos.union(cur_neg), np.mean(values), np.std(values))
                        
                        
                        #TODO - only consider groups with more than one element selected 
                        # for the consistency calculation
                        if len(cur_pos) + len(cur_neg) > 1:
                            consistencies += [np.abs(np.mean(values))] * len(values)
                        # no need to loop any further, the groupings are exclusive
                if consistencies:                
                    consistency = np.mean(consistencies)
                else:
                    # only ungrouped, pretty much shows a strength only
                    if weighted_consistency:                    
                        consistency = np.mean(np.abs(prc_arr[comb]))
                    else:
                        # not weighted consistency, so we use one for single properties
                        consistency = 1
                # now handle the ungrouped
                result["ungrouped"] = (ungrouped_properties.intersection(scomb), np.nan, np.nan)
                
                # and store the entry for the comb
                groups_consistencies[num_prop] = groups_consistencies[num_prop] + [consistency]
                groups_stats[num_prop] = groups_stats[num_prop] + [result]
            
            # 
            itercount = 0
            while True:
                prcs[num_prop] = map(    
                        calc_percentage,
                        zip(
                            map(lambda comb: la.union_properties(0, comb, property_indices), combs[num_prop]),
                            map(lambda comb: la.union_properties(1, comb, property_indices), combs[num_prop])
                        )
                )
                quality_indices[num_prop] = list(np.sqrt(np.array(prcs[num_prop]) * np.array(groups_consistencies[num_prop])))
                
                # how do the prcs look like; do we need to increase the threshold?
                if 1 - sorted(prcs[num_prop])[-1] < 0.01:
                    print "increasing"
                    THRESHOLD = increase_threshold(THRESHOLD)
                    percentages, bool_signs, property_indices = decision_fun(res_1, res_2, THRESHOLD)
                else:
                    break #out of the while loop    
                if itercount < 10:
                    itercount += 1
                else:
                    print "Unable to find suitable threshold for current combination"
                    # break out of both loops
                    raise mylabel
            
            
            # calculate combination for the next iteration
            combs[num_prop+1] = la.limited_comb_expand(prcs[num_prop], combs[num_prop], all_properties, ITERATION_LIMIT)
            print "finished"
            
            winners[num_prop] = la.which_largest(prcs[num_prop], combs[num_prop], percentages, 2)
            
    return winners, combs, prcs, quality_indices, groups_stats, percentages, bool_signs
            
            
            



