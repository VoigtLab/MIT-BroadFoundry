#!/usr/bin/env python
'''
Composite Part Sequences
========================
    Script to generate the reference sequences for composite part pools.

    Usage:
    ------
    python composite_sequences.py -R -p POOL_DESIGN_FILENAME -f PARTS_FASTA_FILENAME
        -i INVERSE_PRIMER_SEQS -L SELF_CIRC_SEQ > OUTPUT_FILENAME

    POOL_DESIGN_FILENAME    - Part variants used at each position in composite part constructs
    PARTS_FASTA_FILENAME    - FASTA file with all parts used in pool design file
    INVERSE_PRIMER_SEQS     - Primer sequences used to amplify fragments after circularization (forward,reverse)
    SELF_CIRC_SEQ           - Linking sequence used for self-circularization ('GCGGCCGC' for NotI)
    OUTPUT_FILENAME         - Output filename to save results to.
'''

# To run: python composite_sequences.py -p AmeR_2NOR.csv -f AllParts.fa -i GTCTACCaggtTAGTTGCGT,CACCTTTACGCACCTGATCA -L GCGGCCGC
#    > outputfilename.circ.fa
# For uncircularized linear fragments run the program without the -i and -L arguments
# For regular expressions run with "-R" added to the options

# Composite Sequences
# Copyright (C) 2016 by
# MIT
# All rights reserved.
# OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.


#(C) MIT, https://opensource.org/licenses/NPOSL-3.0
__license__ = 'OSI Non-Profit OSL 3.0'
__author__  = 'D. Benjamin Gordon, MIT-Broad Foundry, MIT'
__version__ = '1.0'

import sys
import argparse
import re
from string import maketrans


def main():
    parser = argparse.ArgumentParser(description="Generate reference sequences for composite part pools")


    args = parser.add_argument_group("Arguments")


    args.add_argument('--partsfile',      '-p', action="store", default=None, dest="PartsFile",     help='CSV file containing pool part composition')
    args.add_argument('--partSeqs',       '-f', action="store", default=None, dest="PartSeqFile",   help='FASTA file containing all parts')
    args.add_argument('--invprimers',     '-i', action="store", default=None, dest="InvPrimers",    help='"ACAGT,TCACGA" primers for inverse PCR')
    args.add_argument('--linker',         '-L', action="store", default=None, dest="Linker",        help='linker sequence')
    args.add_argument('--regex',          '-R', action="store_true", default=False,dest="RegEx",    help='output regular expressions')
    options = parser.parse_args()


    partsfile  = options.PartsFile
    seqfile    = options.PartSeqFile
    invprimers = options.InvPrimers
    linker     = options.Linker
    regexflag  = options.RegEx


    if invprimers:
        invprimers = invprimers.upper().split(',')


    if not partsfile or not seqfile:
        print "USAGE: Dialout_Modules.py -p partsfile.csv -f seqfile.fa"
        exit(1)


    id = partsfile.split('/')[-1].replace('.csv','')
    filenametoks = [x.lower() for x in re.split('[/\._]',partsfile.split('/')[-1])]
    #print filenametoks                             #Uncomment to add the parsing of POOL_DESIGN_FILENAME to output file header

    seqD = FastaLoad(seqfile)
    partslists = []
    cdsSet = set()


    GUESS_CDS_FROM_FILENAME = True
    tokidxs = []
    counter = 1
    for line in open(partsfile).readlines():
        parts = [x for x in line.strip().split(',') if x != '']
        partslists.append(parts)

        if GUESS_CDS_FROM_FILENAME:
            if len(parts) > 1:                      #found variable part (for determining naming)
                tokidxs.append(counter)
            elif parts[0].lower() in filenametoks:  #found part matching filename (for determining naming)
                tokidxs.append(counter)
        else:
            tokidxs.append(counter)                 #Add them all (for naming)
        counter += 1

        #print tokidxs                              #Uncomment to add positions used for construct names to output file header

    list_of_partslists = compute_list_of_partslists(partslists)

    #seqnamepairs will be replaced at the end of each iteration with a new list
    #that contains all combinations so far

    ans = []
    for partslists in list_of_partslists:
        seqnamepairs = [([],['module'])]


        for partslist in partslists:
            newseqnamepairs = []
            for seqlist,namelist in seqnamepairs:
                for part in partslist:
                    newseqlist = seqlist + [seqD[part]]
                    newnamelist = namelist + [part]
                    newseqnamepairs.append((newseqlist, newnamelist))
            seqnamepairs = newseqnamepairs
        ans.extend(seqnamepairs)


    if not invprimers:
        j = 0
        for seqlist,namelist in ans:
            #seqid = name.replace('module','%s_%d'%(id,i+1))
            toks = namelist
            seqid = '_'.join(toks[i] for i in tokidxs)
            seq = ''.join(seqlist)
            print '>%s\n%s'%(seqid,seq)
            j += 1
    else:
        j = 0
        for seqlist,namelist in ans:
            toks = namelist
            seqid = '_'.join(toks[i] for i in tokidxs)
            preseq = ''.join(seqlist).upper()


            try:
                fwds = [m[0] for m in align(invprimers[0].upper(),preseq) if m[1] == 'f']
                fwdpos = fwds[0]
            except:
                print j, invprimers[0]
                print preseq
                fwds = [m[0] for m in align(invprimers[0].upper(),preseq) if m[1] == 'f']
                print fwds
                fwdpos = fwds[0]


            revs = [m[0] for m in align(invprimers[1].upper(),preseq) if m[1] == 'r']
            revpos = revs[0]


            seq = preseq[fwdpos:] + linker + preseq[:revpos]
            if regexflag:
                seq = string2regex(seq)
            print '>%s\n%s'%(seqid,seq)
            j += 1


import re
def string2regex(s):
    laststart = 0
    s = s.upper()
    Nres = list(re.finditer('N+',s))
    newstring = ''
    for Nre in Nres:
        newstring = newstring + s[laststart:Nre.start()]
        Nlen = Nre.end() - Nre.start()
        newstring = newstring + '([ACGT]{%d})'%(Nlen)
        laststart = Nre.end()
    newstring = newstring + s[laststart:]
    return newstring


def test_string2regex():
    s = 'actNNNNTTTNNNTATATACTACNNNNNNNNNNNNNNC'
    print s
    print string2regex(s)
    s = 'actNNNNTTTNNNTATATACTACNNNNNNNNNNNNNN'
    print s
    print string2regex(s)
    s = 'NNNNTTTNNNTATATACTACNNNNNNNNNNNNNNC'
    print s
    print string2regex(s)


def compute_list_of_partslists(partslists):
    #Generate new partslists, only allowing f w f and r w r
    #There are allowed:


    #   RegTC_5 Bbs1 1 <more stuff> 9 Bbs1R
    #   RegTC_5 Bbs1 r2 <more stuff> r1 Bbs1R


    #and these are disallowed:


    #   RegTC_5 Bbs1 1 <more stuff> r1 Bbs1R
    #   RegTC_5 Bbs1 r2 <more stuff> 2 Bbs1R
    lopl = []


    #only Forward:
    ans = []
    for partslist in partslists:
        if all([part.replace('r','').isdigit() for part in partslist]):  #Only matches linkers
            ans.append([part for part in partslist if part.find('r')<0]) #Keep only forw
        else:
            ans.append(partslist[:])
    lopl.append(ans)


    #only reverse:
    ans = []
    for partslist in partslists:
        if all([part.replace('r','').isdigit() for part in partslist]):  #Only matches linkers
            ans.append([part for part in partslist if part.find('r')>=0]) #Keep only forw
        else:
            ans.append(partslist[:])
    lopl.append(ans)


    #print count_combinations(lopl[0]), lopl[0]
    #print count_combinations(lopl[1]), lopl[1]
    return lopl


def count_combinations(lol):
    ans = 1
    for sublist in lol:
        ans = ans * len(sublist)
    return ans


def FastaLoad(filename):
    D = {}
    f = open(filename)
    chunks = [c.strip() for c in f.read().split('>') if len(c) > 1]
    for chunk in chunks:
        lines = chunk.split('\n')
        name = lines[0].split()[0]
        seq  = ''.join([line.strip() for line in lines[1:]])
        D[name] = seq
        #print ">>%s<<%d>>\n[%s]"%(name,len(seq),seq)
    return D

def align(primer,seq):
    rcprimer = revcomplement(primer)
    forw = [(m.start(),"f")             for m in re.finditer('(?=%s)'%primer, seq)]
    rev  = [(m.start()+len(primer),"r") for m in re.finditer('(?=%s)'%rcprimer, seq)]
    return sorted(forw+rev)

#Global vars for revcomplement
_revcomplement_memo = {'A':'T'}
_revcompTBL = maketrans("AGCTagctWSKMYRnN", "TCGAtcgaWSMKTYnN")
def revcomplement(seq):
    """
    revcomplement(seq)
    A quick reverse-complement routine that memo-izes queries, understands
    IUPAC ambiguity codes, and preserves case.
    """
    global _revcomplement_memo
    try:
        rc = _revcomplement_memo[seq]
    except KeyError:
        _t = list(seq.translate(_revcompTBL))
        _t.reverse()
        rc = ''.join(_t)
        _revcomplement_memo[seq] = rc
        _revcomplement_memo[rc]  = seq
    return(rc)


if __name__ == '__main__': main()
