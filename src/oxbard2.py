#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import numpy as np

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("forbarcode1")
    parser.add_argument("forbarcode2")
    parser.add_argument("revbarcode1")
    parser.add_argument("revbarcode2")
    parser.add_argument("reads")
    parser.add_argument("outpath")
    parser.add_argument("--threshold", type=int, default=8, help="Minimum number of matches for classification per barcode")
    parser.add_argument("--label1", type=str, default="barcode1", help="Name for barcode set 1")
    parser.add_argument("--label2", type=str, default="barcode2", help="Name for barcode set 2")
    return parser.parse_args()

def revcomp(seq : str) -> str:
    complement = {
        'A' : 'T',
        'C' : 'G',
        'G' : 'C',
        'T' : 'A'
    }
    seq = list(seq)
    for i in range(len(seq)):
        seq[i] = complement[seq[i]]
    seq = ''.join(seq)
    return seq[::-1]

def main() -> None:
    args = parse()
    forbc1 = args.forbarcode1.replace("U", "T")
    forbc2 = args.forbarcode2.replace("U", "T")
    revbc1 = args.revbarcode1.replace("U", "T")
    revbc2 = args.revbarcode2.replace("U", "T")
    reads = args.reads
    outpath = args.outpath
    # outcsv = open(os.path.join(outpath, "classification.csv"), "w")
    # outcsv.write("readid,strand,front_bc1_matches,back_bc1_matches,front_bc2_matches,back_bc2_matches\n")
    outbc1 = open(os.path.join(outpath, f"{args.label1}.ids"), "w")
    outbc2 = open(os.path.join(outpath, f"{args.label2}.ids"), "w")
    outboth = open(os.path.join(outpath, f"both.ids"), "w")
    reads_starting_with_polyA = open(os.path.join(outpath, f"polyAreads.ids"), "w")

    threshold = args.threshold
    
    aligner = PairwiseAligner()
    aligner.internal_gap_score = -1
    aligner.mismatch = -1

    counts = {args.label1 : 0, args.label2 : 0, "both" : 0, "sense" : 0, "antisense" : 0, "polyAreads" : 0}
    # [[forbc1_score, forbc1_matches], ...]
    senseScores = np.zeros((4, 2), dtype=int)
    antisenseScores = np.zeros((4, 2), dtype=int)
    # try:
    for record in SeqIO.parse(reads, "fastq"):
        # print(record.id)
        seq = str(record.seq).replace("U", "T")
        # first 30 nts
        front = seq[:50]
        # cut polyA
        try:
            if seq.index('AAAAAA') == 0:
                reads_starting_with_polyA.write(record.id + '\n')
                counts['polyAreads'] += 1
                continue
            seq = seq[:seq.index('AAAAAA')]
        except:
            pass
        # last 150 nts before the polyA (sometimes we see no poly A in the read but weird adapter sequences)
        back = seq[-150:]

        # FORWARD / SENSE Strand
        # align forward barcode at front of sequence
        alignments = aligner.align(front, forbc1)
        senseScores[0] = [alignments[0].score, alignments[0].counts()[1]]

        alignments = aligner.align(front, forbc2)
        senseScores[1] = [alignments[0].score, alignments[0].counts()[1]]

        # align rev barcode at back of sequence
        alignments = aligner.align(back, revbc1)
        senseScores[2] = [alignments[0].score, alignments[0].counts()[1]]

        alignments = aligner.align(back, revbc2)
        senseScores[3] = [alignments[0].score, alignments[0].counts()[1]]

        # REVERSE COMPLEMENT / ANTISENSE Strand
        alignments = aligner.align(front, revcomp(revbc1))
        antisenseScores[0] = [alignments[0].score, alignments[0].counts()[1]]

        alignments = aligner.align(front, revcomp(revbc2))
        antisenseScores[1] = [alignments[0].score, alignments[0].counts()[1]]

        alignments = aligner.align(back, revcomp(forbc1))
        antisenseScores[2] = [alignments[0].score, alignments[0].counts()[1]]

        alignments = aligner.align(back, revcomp(forbc2))
        antisenseScores[3] = [alignments[0].score, alignments[0].counts()[1]]

        if max(senseScores[:, 0]) > max(antisenseScores[:, 0]):
            # print("sense")
            sense = True
            strandScores = senseScores
        else:
            # print("antisense")
            sense = False
            strandScores = antisenseScores

        # print('barcode set 1 matches - front:', strandScores[0], 'back:', strandScores[2])
        # print('barcode set 2 matches - front:', strandScores[1], 'back:', strandScores[3])

        # number of matches per barcode for classification exceeds threshold
        # allow at max 3 gaps in barcode alignments
        forbc1_pass = strandScores[0, 1] >= threshold and strandScores[0, 0] + 3 >= strandScores[0, 1]
        forbc2_pass = strandScores[1, 1] >= threshold and strandScores[1, 0] + 3 >= strandScores[1, 1]
        revbc1_pass = strandScores[2, 1] >= threshold and strandScores[2, 0] + 3 >= strandScores[2, 1]
        revbc2_pass = strandScores[3, 1] >= threshold and strandScores[3, 0] + 3 >= strandScores[3, 1]

        # print('barcode set 1 - front: ', forbc1_pass, 'back: ', revbc1_pass)
        # print('barcode set 2 - front: ', forbc2_pass, 'back: ', revbc2_pass)

        bc1_pass = forbc1_pass or revbc1_pass
        bc2_pass = forbc2_pass or revbc2_pass

        if bc1_pass or bc2_pass:
            if sense:
                counts["sense"] += 1
            else:
                counts["antisense"] += 1

        # both barcode sets pass
        if bc1_pass and bc2_pass:
            counts["both"] += 1
            outboth.write(record.id + '\n')

        # barcode set 1 passes
        elif bc1_pass:
            counts[args.label1] += 1
            outbc1.write(record.id + '\n')

        # barcode set 2 passes
        elif bc2_pass:
            counts[args.label2] += 1
            outbc2.write(record.id + '\n')

    # except Exception as e:
    #     print(e)
    #     print(record.id)
    #     print(len(record.seq), record.seq)
    #     exit(1)

    print("Done")
    print(counts)

if __name__ == '__main__':
    main()