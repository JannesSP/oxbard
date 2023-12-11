#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import mappy as mp # https://pypi.org/project/mappy/
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from Bio import SeqIO
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Mappy cigar
# CIGAR = {
# 0 : alignment match, (M)
# 1 : insertion, (I)
# 2 : deletion, (D)
# 3 : skipped, (N)
# 4 : softclip, (S)
# 5 : hardclip, (H)
# 6 : padding, (P)
# 7 : sequence match (=)
# 8 : sequence mismatch (X)
# }

def parse() -> Namespace:
    # TODO i probably need the positions of barcodes within the given reference sequence (start, len)
    # then check the ref position when looping through the cigar string and check matchrate
    
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter) 
    parser.add_argument('basecalls', metavar = 'FASTQ', help='basecalls to classify')
    parser.add_argument('modRef', metavar = 'FASTA', help='reference with barcode of modified read')
    parser.add_argument('canRef', metavar = 'FASTA', help='reference with barcode of canonical read')
    parser.add_argument('outfile', metavar = 'PATH', help = 'output file')
    parser.add_argument('-f', '--frontBCInterval', metavar = '[x,y)', nargs="+", help='1-based interval for the barcode region in the 5\' end')
    parser.add_argument('-b', '--backBCInterval', metavar = '[x,y)', nargs="+", help='1-based interval for the barcode region in the 3\' end')
    parser.add_argument('--mr', metavar='FLOAT', type=float, default=0.5, help='match rate between 0.0 and 1.0')
    parser.add_argument('--strand', default=1, choices=[-1, 1], type=int, help='Which strand to look for, forward: 1, reverse complement: -1.')
    return parser.parse_args()

def analyseCigar(cigarArray : list, fBCInterval : tuple, bBCInterval : tuple, refPos : int):
    '''
    Iterate over cigar string and check for barcode regions.

    Parameters
    ----------
    cigarArray : list
        cigar array from mappy [(#bases, cigar_character)]
    fBCInterval : tuple
        region of front barcode on reference sequence (incl., excl.)
    bBCInterval : tuple
        region of back barcode on reference sequence (incl., excl.)

    Returns
    -------
    (fBCMatches, bBCMatches) : tuple
        fBCMatches : int
            number of bases matching to the front barcode region
        bBCMatches : int
            number of bases matching to the back barcode region
    '''

    # init counters
    fBCMatches = bBCMatches = numMatched = 0
    # queryPos = 0
    wasInFrontBC = wasInBackBC = hasMatched = False

    for cigar in cigarArray: # https://samtools.github.io/hts-specs/SAMv1.pdf

        if hasMatched:
            # count matches for front barcode
            if refPos in range(*fBCInterval):
                fBCMatches += min(refPos - fBCInterval[0] + 1, numMatched)
                wasInFrontBC = True
            # add last matched bases when leaving barcode region
            elif wasInFrontBC:
                if refPos - hasMatched in range(*fBCInterval):
                    fBCMatches += fBCInterval[1] - (refPos - numMatched)
                wasInFrontBC = False

            # count matches for back barcode
            if refPos in range(*bBCInterval):
                bBCMatches += min(refPos - bBCInterval[0] + 1, numMatched)
                wasInBackBC = True
            # add last matched bases when elaving barcode region
            elif wasInBackBC:
                if refPos - hasMatched in range(*bBCInterval):
                    bBCMatches += bBCInterval[1] - (refPos - numMatched)
                wasInBackBC = False

        hasMatched = False
        if cigar[1] in [0, 7, 8]: # M, =, X consumes both
            hasMatched = cigar[1] != 8
            numMatched = cigar[0]
            refPos += cigar[0]
            # queryPos += cigar[0]
        elif cigar[1] in [2, 3]: # D, N consumes ref
            refPos += cigar[0]
        # elif cigar[1] in [1, 4]: # I, S consumes query
            # queryPos += cigar[0]
        elif cigar[1] in [1, 4, 5, 6]: # H, P consomes nothing
            continue
        else:
            print(f'ERROR! Unknown cigar character! {cigar}')
            exit(1)

    return fBCMatches, bBCMatches

def mapReads(basecalls : str, mod_ref : str, can_ref : str, frontBCInterval : tuple, backBCInterval : tuple, matchRate : float, outfile : str, strand : int):

    # over1900 = under1900 = 0
    df = pd.DataFrame(columns=['readid', 'reference', 'read length', 'read quality', 'mapped ref span', 'first reference position mapped', 'last reference position mapped', 'front matches', 'back matches'])
    can_clas = set()
    mod_clas = set()
    mod_coverage = np.zeros(len(readFasta(mod_ref)['seq']), dtype=int)
    can_coverage = np.zeros(len(readFasta(can_ref)['seq']), dtype=int)

    # make intervals 0-based
    for i in range(len(frontBCInterval)):
        frontBCInterval[i] -= 1
        backBCInterval[i] -= 1
        # revert indices for other strand positions
        if strand == -1:
            frontBCInterval[i] = len(mod_coverage) - frontBCInterval[i]
            backBCInterval[i] = len(mod_coverage) - backBCInterval[i]
    if strand == -1:
        frontBCInterval = frontBCInterval[::-1]
        backBCInterval = backBCInterval[::-1]

    w = open(outfile, 'w')

    mod_aligner = mp.Aligner(mod_ref, preset = 'splice', n_threads = 12, k = 9, best_n = 1)  # load or build index
    can_aligner = mp.Aligner(can_ref, preset = 'splice', n_threads = 12, k = 9, best_n = 1)
    if not mod_aligner: raise Exception("ERROR: failed to load/build index in mod_aligner")
    if not can_aligner: raise Exception("ERROR: failed to load/build index in can_aligner")
    w.write('readID,classification,frontClassified,backClassified,numFrontBCMatches,ratioFront,numBackBCMatches,ratioBack,mapLen,fistMappedRefPos,lastMappedRefPos\n')
    
    for i, (name, seq, qual) in enumerate(mp.fastx_read(basecalls)): # read a fasta/q sequence

        qual = np.mean([ord(letter) - 33 for letter in qual])

        # cut polyA + adapter
        # saw in unclassified_align.fa, that the mod barcode with many A, C and Us could align to the adapter(?) sequence after (3') the "poly A", within the 300 unclassified reads the AAAAA was always the "poly A" before (5') the adapter(?) sequence
        # better cut it to avoid false positives/classified reads
        # try:
        #     seq = seq[:seq.index('AAAAAAAA')] # set by experience of looking at fastq sequences
        # except ValueError:
        #     pass

        # if len(seq) >= 1900:
        #     over1900+=1
        # else:
        #     under1900+=1

        if (i + 1) % 1000 == 0:
            print('mapping read', i + 1, end = '\r')

        for mod_hit in mod_aligner.map(seq): # traverse alignments

            if mod_hit.strand != strand: # skip revcomp mapped reads
                continue

            # check barcodes
            fBCMatches, bBCMatches = analyseCigar(mod_hit.cigar, frontBCInterval, backBCInterval, mod_hit.r_st)
            frontClassified = fBCMatches >= matchRate * (frontBCInterval[1] - frontBCInterval[0])
            backClassified = bBCMatches >= matchRate * (backBCInterval[1] - backBCInterval[0])
            w.write(f'{name},mod,{frontClassified},{backClassified},{fBCMatches},{fBCMatches/(frontBCInterval[1] - frontBCInterval[0])},{bBCMatches},{bBCMatches/(backBCInterval[1] - backBCInterval[0])},{mod_hit.r_en-mod_hit.r_st},{mod_hit.r_st},{mod_hit.r_en}\n')
            if frontClassified or backClassified:
                mod_clas.add(name)
                mod_coverage[mod_hit.r_st : mod_hit.r_en + 1] += 1

            new_entry = pd.DataFrame({
                'readid' : [name],
                'reference' : ['mod'],
                'read length' : [len(seq)],
                'read quality' : [qual],
                'mapped ref span' : [mod_hit.r_en + 1 - mod_hit.r_st],
                'first reference position mapped' : [mod_hit.r_st],
                'last reference position mapped' : [mod_hit.r_en],
                'front matches' : [fBCMatches],
                'back matches' : [bBCMatches],
            })
            df = pd.concat((df, new_entry), ignore_index=True)

        for can_hit in can_aligner.map(seq): # traverse alignments

            if can_hit.strand != strand: # skip revcomp mapped reads
                continue

            fBCMatches, bBCMatches = analyseCigar(can_hit.cigar, frontBCInterval, backBCInterval, can_hit.r_st)
            frontClassified = fBCMatches >= matchRate * (frontBCInterval[1] - frontBCInterval[0])
            backClassified = bBCMatches >= matchRate * (backBCInterval[1] - backBCInterval[0])
            w.write(f'{name},can,{frontClassified},{backClassified},{fBCMatches},{fBCMatches/(frontBCInterval[1] - frontBCInterval[0])},{bBCMatches},{bBCMatches/(backBCInterval[1] - backBCInterval[0])},{can_hit.r_en-can_hit.r_st},{can_hit.r_st},{can_hit.r_en}\n')
            if frontClassified or backClassified:
                can_clas.add(name)
                can_coverage[can_hit.r_st : can_hit.r_en + 1] += 1

            new_entry = pd.DataFrame({
                'readid' : [name],
                'reference' : ['can'],
                'read length' : [len(seq)],
                'read quality' : [qual],
                'mapped ref span' : [can_hit.r_en + 1 - can_hit.r_st],
                'first reference position mapped' : [can_hit.r_st],
                'last reference position mapped' : [can_hit.r_en],
                'front matches' : [fBCMatches],
                'back matches' : [bBCMatches],
            })
            df = pd.concat((df, new_entry), ignore_index=True)

    w.close()

    # dir = os.path.dirname(outfile)
    modids = open(os.path.splitext(outfile)[0] + '_mod.ids', 'w')
    canids = open(os.path.splitext(outfile)[0] + '_can.ids', 'w')
    unclassified = open(os.path.splitext(outfile)[0] + '_unclassified.ids', 'w')

    for i, (name, seq, qual) in enumerate(mp.fastx_read(basecalls)): # read a fasta/q sequence
        if name in can_clas and name in mod_clas:
            unclassified.write(name + '\n')
        elif name in can_clas:
            canids.write(name + '\n')
        elif name in mod_clas:
            modids.write(name + '\n')
        else:
            unclassified.write(name + '\n')

    print(f'mapped {i + 1} reads')

    # print('Read Length >= 1900: ', over1900, ", < 1900: ", under1900)
    print(f'{os.path.splitext(os.path.basename(outfile))[0]}:\tClassified: mod: {len(mod_clas - can_clas)}, can: {len(can_clas - mod_clas)}, both: {len(can_clas & mod_clas)} total: {i+1}')
    w = open(os.path.splitext(outfile)[0] + '.log', 'w')
    w.write(f'Classified: mod: {len(mod_clas - can_clas)}, can: {len(can_clas - mod_clas)}, both: {len(can_clas & mod_clas)} total: {i+1}')
    w.close()
    with open(os.path.join(os.path.dirname(outfile), 'result.csv'), 'a') as a:
        a.write(f'{os.path.splitext(os.path.basename(outfile))[0]},mod,{len(mod_clas - can_clas)}\n')
        a.write(f'{os.path.splitext(os.path.basename(outfile))[0]},can,{len(can_clas - mod_clas)}\n')
        a.write(f'{os.path.splitext(os.path.basename(outfile))[0]},both,{len(can_clas & mod_clas)}\n')
    return mod_coverage, can_coverage, df

def plotCoverage(mod : np.ndarray, can : np.ndarray, file : str, frontBCInterval : list, backBCInterval : list) -> None:
    if not os.path.exists(os.path.dirname(file)):
        os.makedirs(os.path.dirname(file))

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 12), dpi = 300)
    ax[0].bar(
        x = np.arange(1, len(mod) + 1),
        height = mod)
    ax[0].grid(True, 'both', 'both', linewidth = 0.5, c = 'grey', alpha = 0.6)
    ax[0].set_title('mod')
    ax[1].bar(
        x = np.arange(1, len(can) + 1),
        height = can)
    ax[1].grid(True, 'both', 'both', linewidth = 0.5, c = 'grey', alpha = 0.6)
    ax[1].set_title('can')
    # for axis in ax:/
    for x in frontBCInterval + backBCInterval:
        ax[0].axvline(x, linewidth=1, linestyle='dashed', color='red', label='barcode')
        ax[1].axvline(x, linewidth=1, linestyle='dashed', color='red', label='barcode')
    fig.suptitle('Coverage')
    fig.supylabel('coverage')
    fig.supxlabel('position')
    # plt.grid(True, 'both', 'both', linewidth = 0.5, c = 'grey', alpha = 0.6)
    plt.tight_layout()
    plt.savefig(file)
    plt.close()

def readFasta(fasta : str) -> dict:
    '''
    Expects only on sequence in fasta
    
    Parameters
    ----------
    fasta : str
        path to fasta file

    Returns
    -------
    seq : dict
        keys:
            header, containes the sequence header
            seq, contains the nucleotide sequence
    '''
    record = next(SeqIO.parse(open(fasta),'fasta'))
    return {'header': record.id, 'seq': str(record.seq)}

def writeIds(ids : list, filepath : str) -> None:
    '''
    Write ids from list into file

    Parameters
    ----------
    ids : list
    filepath : str
    '''
    with open(filepath, 'w') as w:
        for id in ids:
            w.write(id + '\n')

def plotRest(data : pd.DataFrame, outfilebasename : str) -> None:
    # ['readid', 'reference', 'read length', 'read quality', 'mapped ref span', 'first reference position mapped', 'last reference position mapped', 'front matches', 'back matches']
    sns.set_theme()
    sns.kdeplot(data = data, x = 'read length', y = 'read quality', hue='reference', fill = True)
    plt.savefig(outfilebasename + '_readlength_vs_quality.pdf')
    plt.savefig(outfilebasename + '_readlength_vs_quality.svg')
    plt.cla()
    plt.close()

    sns.set_theme()
    sns.kdeplot(data = data, x = 'read length', y = 'mapped ref span', hue='reference', fill = True)
    plt.savefig(outfilebasename + '_readlength_vs_maplength.pdf')
    plt.savefig(outfilebasename + '_readlength_vs_maplength.svg')
    plt.cla()
    plt.close

    sns.set_theme()
    sns.kdeplot(data = data, x = 'read length', y = 'front matches', hue='reference', fill = True)
    plt.savefig(outfilebasename + '_readlength_vs_numfrontBCmatches.pdf')
    plt.savefig(outfilebasename + '_readlength_vs_numfrontBCmatches.svg')
    plt.cla()
    plt.close

    sns.set_theme()
    sns.kdeplot(data = data, x = 'read length', y = 'back matches', hue='reference', fill = True)
    plt.savefig(outfilebasename + '_readlength_vs_numbackBCmatches.pdf')
    plt.savefig(outfilebasename + '_readlength_vs_numbackBCmatches.svg')
    plt.cla()
    plt.close

    sns.set_theme()
    sns.histplot(data = data, x = 'first reference position mapped', hue='reference', log_scale=(False, True))
    plt.savefig(outfilebasename + '_firstmappedpos.pdf')
    plt.savefig(outfilebasename + '_firstmappedpos.svg')
    plt.cla()
    plt.close

    sns.set_theme()
    sns.histplot(data = data, x = 'last reference position mapped', hue='reference', log_scale=(False, True))
    plt.savefig(outfilebasename + '_lastmappedpos.pdf')
    plt.savefig(outfilebasename + '_lastmappedpos.svg')
    plt.cla()
    plt.close

def main() -> None:
    args = parse()
    outfile = args.outfile

    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    frontBCInterval = list(map(int, args.frontBCInterval))
    backBCInterval = list(map(int, args.backBCInterval))

    mod_coverage, can_coverage, data = mapReads(args.basecalls, args.modRef, args.canRef, frontBCInterval, backBCInterval, args.mr, outfile, args.strand)

    plotCoverage(mod_coverage, can_coverage, os.path.splitext(outfile)[0] + '_coverage.png', frontBCInterval, backBCInterval)
    plotRest(data, os.path.splitext(outfile)[0])

if __name__ == '__main__':
    main()
