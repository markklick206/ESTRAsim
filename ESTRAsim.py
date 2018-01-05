# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 20:20:45 2016

@author: mk

Usage: python ESTRAsim.py <reference genome fasta file> -e <enzyme name> | -s <cutsite sequence>
Output: <file> formatted like Radiqual "rad_tags.txt" with cutsite sequence in filename

"""

from Bio import SeqIO
from Bio.Restriction import *
from Bio import Restriction
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio.Restriction.Restriction_Dictionary import typedict
import numpy
import math
import sys
import os

def RADdigest(sequence_file):
    #seq = SeqIO.read(open(sequence_file,'rU'),'fasta')
    #the whole restriction enzyme dictionary
    "EcoRI" in rest_dict.keys()
    
    #Create a list of enzymes in your digest
    multi = (BamHI, EcoRI, PstI)
    
    #Length of the fragments from the digest as boundaries
    frag_length = (49,100)
    
    for record in SeqIO.parse(sequence_file, "fasta"):
        for enz in range(len(multi)):
            #print recognition site for the enzyme
            print(str(multi[enz]) + "'s recognition sequence: " + multi[enz].site)
            
            #obtain the cut site coordinates for each enzyme
            cut_sites = multi[enz].search(record.seq)
            
            #print(len(cut_sites))
            # get the distance between cut sites
            dist = [cut_sites[r] - cut_sites[r-1] for r in range(1,len(cut_sites))]
            #print(dist)
            
            # generate some summary stats
            dist_a = numpy.array(dist)
            # mean fragment size
            print ("mean fragment size: " + str(numpy.mean(dist_a)))
            
            #open file
            rad_tags = open(str(sequence_file)+"_"+record.id +"_"+str(multi[enz])+"_"+multi[enz].site +".fasta", "wt")
            
            #get the coordinates as tuples to extract the sequences of fragments
            for site in range(1,len(cut_sites)):
                frag = record.seq[cut_sites[site-1]:cut_sites[site]]
                #print len(frag)
                #if len(frag) == frag_length:
                if len(frag) >= frag_length[0] and len(frag) <= frag_length[1]:
                    rad_tags.write(str(frag) + '\n')
                #print  "%s    %d    %s" % (record.id, cut_sites[site], multi[enz])

def append_rad_tag_files(rad_tags, rad_tags1):
    concat_rad_tags = open("concat_fastas_" + str(rad_tags) + str(rad_tags1), "wt")
    for line in open(rad_tags):
        concat_rad_tags.write(line)
    for line in open(rad_tags1):
        concat_rad_tags.write(line)
            
if __name__ == '__main__':
    fasta_file = sys.argv[1]
    
    RADdigest(fasta_file)
    #RADdigest(fasta_file,enzyme_name)
    #RADdigest('E.coli.fasta')
    
    append_rad_tag_files("speciesB.reference.fa_chr0_BamHI_GGATCC.fasta","speciesB.reference.fa_chr1_BamHI_GGATCC.fasta")
    append_rad_tag_files("speciesB.reference.fa_chr0_EcoRI_GAATTC.fasta","speciesB.reference.fa_chr1_EcoRI_GAATTC.fasta")
    append_rad_tag_files("speciesB.reference.fa_chr0_PstI_CTGCAG.fasta","speciesB.reference.fa_chr1_PstI_CTGCAG.fasta")
    
    
    
#def RADdigest2(sequence_file,enzyme_name):
#    #the whole restriction enzyme dictionary
#    enzyme_name in rest_dict.keys()
#    
#    #Length of the fragments from the digest as boundaries
#    frag_length = (49,100)
#    
#    for record in SeqIO.parse(sequence_file, "fasta"):
#        #print recognition site for the enzyme
#        print(str(enzyme_name) + "'s recognition sequence: " + enzyme_name.site)
#        
#        #obtain the cut site coordinates for each enzyme
#        cut_sites = enzyme_name.search(record.seq)
#        
#        #print(len(cut_sites))
#        # get the distance between cut sites
#        dist = [cut_sites[r] - cut_sites[r-1] for r in range(1,len(cut_sites))]
#        #print(dist)
#        
#        # generate some summary stats
#        dist_a = numpy.array(dist)
#        # mean fragment size
#        print ("mean fragment size: " + str(numpy.mean(dist_a)))
#        
#        #open file
#        rad_tags = open(str(sequence_file)+"_"+record.id +"_"+str(enzyme_name)+"_"+enzyme_name.site +".fasta", "wt")
#        
#        #get the coordinates as tuples to extract the sequences of fragments
#        for site in range(1,len(cut_sites)):
#            frag = record.seq[cut_sites[site-1]:cut_sites[site]]
#            #print len(frag)
#            if len(frag) >= frag_length[0] and len(frag) <= frag_length[1]:
#                rad_tags.write(str(frag) + '\n')
#            #print  "%s    %d    %s" % (record.id, cut_sites[site], multi[enz])