# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 20:20:45 2016

@author: i7
"""

from Bio import SeqIO
from Bio.Restriction import *
from Bio import Restriction
from Bio.Restriction.Restriction_Dictionary import rest_dict
#from Bio.Restriction.Restriction import RestrictionBatch
from Bio.Restriction.Restriction_Dictionary import typedict
import numpy
import math
import sys

def RADdigest(sequence_file):
    #seq = SeqIO.read(open(sequence_file,'rU'),'fasta')
    #the whole restriction enzyme dictionary
    "EcoRI" in rest_dict.keys()
    
    #Create a list of enzymes in your digest
    multi = (BamHI, EcoRI, PstI)
    
    #open file
    output = open('test_fragments_assmbl.fasta', "wt")
    
    for record in SeqIO.parse(sequence_file, "fasta"):
        for enz in range(len(multi)):
            #obtain the cut site coordinates for each enzyme
            cut_sites = multi[enz].search(record.seq)
            print(len(cut_sites))
            # get the distance between cut sites
            dist = [cut_sites[r] - cut_sites[r-1] for r in range(1,len(cut_sites))]
            print(dist)
            # generate some summary stats
            dist_a = numpy.array(dist)
            # mean fragment size
            print ("mean fragment size: " + str(numpy.mean(dist_a)))
            # 95 % Confidence invterval around mean
            print ("confidence interval around mean: " + str(1.96 * numpy.std(dist_a, ddof=1)/math.sqrt(len(dist_a))))
            # # of fragments btw. 100 and 500 bp using numpy indexing!
            less_than_five_bool    = dist_a < 500
            less_than_five         = dist_a[less_than_five_bool]
            greater_than_one_bool  = less_than_five > 100
            one_to_five_hundred_bp = less_than_five[greater_than_one_bool]
            print ("# of fragments btwn 100 and 500 bp: " + str(len(one_to_five_hundred_bp)))
            #get the coordinates as tuples to extract the sequences of fragments
            for site in range(1,len(cut_sites)):
                frag = record.seq[cut_sites[site-1]:cut_sites[site]]
                output.write(str(frag) + '\n')
                #print  "%s    %d    %s" % (record.id, cut_sites[site], multi[enz])
                
def RADdigest2(sequence_file):                
    #or try to make it more general
    def create_enzyme(name):
        e_types = [x for t, (x, y) in typedict.items() if name in y][0]
        enzyme_types = tuple(getattr(Restriction.Restriction, x) for x in e_types)
        return Restriction.Restriction.RestrictionType(name, enzyme_types, rest_dict[name])
    
    seq = SeqIO.read(open(sequence_file,'rU'),'fasta')    
    enzyme_list = ["BamHI","EcoRI","PstI","MstI"]
    rb = reduce(lambda x, y: x + y, map(create_enzyme, enzyme_list))
    
    ana = Restriction.Analysis(rb, seq.seq, linear=True)
    print(ana.with_sites())

def commons():
    test=seq[0:10000]
    # take a look at the common restriction enzymes
    Restriction.CommOnly
    # do a restriction batch analysis on our test sequence with ALL common
    # enzymes
    Ana = Restriction.Analysis(Restriction.CommOnly, test.seq, linear=True)
    # look that the enzymes that cut (BLUNT) the test sequence in more than #
    # 10 spots
    results = []
    for r in Ana.blunt():
        if len(Ana.blunt()[r]) > 10:
            results.append([r, len(Ana.blunt()[r])])
    # sort by the 2nd element of the list
    for r in sorted(results, key=lambda record: record[1]):
        print r[0], r[1]
    # or you can do this:
    from operator import itemgetter
    for r in sorted(results, key=itemgetter(1)):
        print r[0], r[1]
    
if __name__ == '__main__':
    fasta_file = sys.argv[1]
    #RADdigest('E.coli.fasta')
    RADdigest(fasta_file)
    #RADdigest1('E.coli.fasta')
    #RADdigest2('E.coli.fasta')



 ## Read enzyme name from input. Useful for a UI - but were going .sh scripts 
    ##enzyme_name = input("Enter enzyme name:\n") # E.g EcoRI
    ##print (type(enzyme_name)) # <type 'str'>
    #
    ## Get RestrictionType by name
    #batch = RestrictionBatch()
    #batch.add(enzyme_name)
    #enzyme = batch.get(enzyme_name)
    #print (type(enzyme)) # RestrictionType
    
    #rb = create_enzyme("EcoRI") + create_enzyme("MstI")
    # Or if you have a long list of restriction enzymes:
    # enzyme_list = ["EcoRI", "MstI"]
    # rb = reduce(lambda x, y: x + y, map(create_enzyme, enzyme_list))
#def RADdigest1(sequence_file):
#    seq = SeqIO.read(open(sequence_file,'rU'),'fasta')
#    print ("seq length: " + str(len(seq)))
#
#    Restriction.BamHI.search(seq.seq)
#    sites = Restriction.BamHI.search(seq.seq)
#    # get the number of fragments produced by cutting
#    print ("# of fragments: " + str(len(sites)))
#    # sort the sites by position
#    sites.sort()
#    # get the distance between sites
#    dist = [sites[r] - sites[r-1] for r in xrange(1,len(sites))]
#    print(dist)
#    # generate some summary stats
#    dist_a = numpy.array(dist)
#    # mean fragment size
#    print ("mean fragment size: " + str(numpy.mean(dist_a)))
#    # 95 % Confidence invterval around mean
#    print ("confidence interval around mean: " + str(1.96 * numpy.std(dist_a, ddof=1)/math.sqrt(len(dist_a))))
#    
#    # # of fragments btw. 100 and 500 bp using numpy indexing!
#    less_than_five_bool    = dist_a < 500
#    less_than_five         = dist_a[less_than_five_bool]
#    greater_than_one_bool  = less_than_five > 100
#    one_to_five_hundred_bp = less_than_five[greater_than_one_bool]
#    
#    print ("# of fragments btwn 100 and 500 bp: " + str(len(one_to_five_hundred_bp)))