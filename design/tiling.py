"""
"""
import pandas as pd
from Bio import Seq
from Bio.SeqUtils import GC
from itertools import combinations

from utils import loadGenomeAsDict


def assembleGuide(targetSeq: str, strand: str, posOffset: int, sgRNALength: int):
    if strand == '+':
        rawSequence = targetSeq[posOffset:posOffset + 3 + sgRNALength].upper()
        sgSequence = 'G' + str(Seq.Seq(rawSequence[3:-1]).reverse_complement()).upper()
        return rawSequence, sgSequence

    elif strand == '-':
        rawSequence = targetSeq[posOffset + 1 - 3 - sgRNALength:posOffset + 1].upper()
        sgSequence = 'G' + rawSequence[1:-3].upper()
        return rawSequence, sgSequence


def findGuideInOffset(geneName: str, targetSeq: str,
                      posOffset: int, sgRNALength: int,
                      rangeLength: int, rangeStart: int):
    """Design and score guides
    [x] Find all NGG PAM guides in the selected regions
    [x] Filter if an N in the reference genome
    [x] 4-mer T/U repeats
    [x] More than 1 T/U in the last 4 nucleotides,
        which combine with the first few nucleotides of
        the gRNA scaffold to terminate Pol III transcription
    [x] Filter if <= 20% GC content and >= 90% GC content
    [ ] More than 40% total T/U content
    [ ] Filter on MIT specificity score > 50,
    [ ] 5-mer mononucleotide repeats
    [ ] Low-complexity sequences,
        defined as 10 nucleotides of 2-nt repeat,
        12 nucleotides of 3- or 4-nt repeats, or 18 nucleotides of 5- or 6-nt repeats

    Parameters
    ----------
    geneName : str
    targetSeq: str
    posOffset:
    sgRNALength : int
    rangeLength
    rangeStart
    """

    strand = None

    # search for guide
    if targetSeq[posOffset:posOffset + 1 + 1].upper() == 'CC' \
            and posOffset + 3 + sgRNALength < rangeLength \
            and 'N' not in \
            targetSeq[
            posOffset:posOffset + 3 + sgRNALength
            ].upper() \
            and 'AAAA' not in \
            targetSeq[
            posOffset:posOffset + 3 + sgRNALength
            ].upper() \
            and 'AA' not in \
            targetSeq[
            posOffset + sgRNALength:posOffset + 3 + sgRNALength
            ].upper() \
            and \
            20 < GC(targetSeq[posOffset:posOffset + 3 + sgRNALength]) < 90:

        foundGuide = True
        strand = '+'

    elif targetSeq[posOffset - 1:posOffset + 1].upper() == 'GG' \
            and posOffset - 3 - sgRNALength >= 0 \
            and 'N' not in \
            targetSeq[
            posOffset + 1 - 3 - sgRNALength:posOffset + 1
            ].upper() \
            and 'TTTT' not in \
            targetSeq[
            posOffset + 1 - 3 - sgRNALength:posOffset + 1
            ].upper() \
            and 'TT' not in \
            targetSeq[
            posOffset + 1 - 3:posOffset + 1
            ].upper() \
            and \
            20 < GC(targetSeq[posOffset + 1 - 3 - sgRNALength:posOffset + 1]) < 90:
        foundGuide = True
        strand = '-'
    else:
        foundGuide = False

    # assemble guide, if found
    if foundGuide:
        pamCoord = rangeStart + posOffset
        sgId = f'{geneName}_{strand}_{str(pamCoord)}'
        rawSequence, sgSequence = assembleGuide(targetSeq, strand, posOffset, sgRNALength)
        sgRNA = (sgId, geneName, sgSequence, rawSequence)
        sgRNAInfo = (sgId, geneName, sgRNALength + 3, pamCoord, 'not assigned', pamCoord, strand)
    else:
        sgRNA = None
        sgRNAInfo = None

    return sgRNA, sgRNAInfo


def findAllGuidesInRanges(geneName, chrom, rangeList, genomeDict, endBuffer=0, sgRNALength=20, PAM='NGG',
                          outFormat=None):
    """
    [x] Rearranges the columns and outputs a guidescan friendly table
    [ ] Rearranges the columns and outputs a BED file for viewing

    Parameters
    ----------
    geneName : str
    chrom : str
    rangeList : list
    genomeDict : dict
    endBuffer : int
    sgRNALength : int
    PAM: str
    outFormat: str
    """
    Library = []
    LibraryInfo = []

    for rangeStart, rangeEnd in rangeList:
        targetSeq = str(genomeDict[chrom].seq[rangeStart - endBuffer: rangeEnd + 1 + endBuffer])
        rangeLength = rangeEnd + 1 - rangeStart + (endBuffer * 2)

        for posOffset in range(rangeLength):
            sgRNA, sgRNAInfo = findGuideInOffset(geneName, targetSeq, posOffset, sgRNALength, rangeLength, rangeStart)
            if sgRNA and sgRNAInfo:
                Library.append(sgRNA)
                LibraryInfo.append(sgRNAInfo)

    LibraryTable = pd.DataFrame(Library, columns=[
        'sgId', 'geneName', 'sequence', 'genomic sequence'
    ]).set_index('sgId')
    LibraryInfoTable = pd.DataFrame(LibraryInfo, columns=[
        'sgId', 'geneName', 'length', 'pam coordinate', 'pass_score', 'position', 'strand'
    ]).set_index('sgId')

    if outFormat == "guidescan":
        # id, sequence, pam, chromosome, position, sense
        guidescanTable = pd.concat([LibraryInfoTable, LibraryTable[["sequence"]]], axis=1).reset_index()
        guidescanTable["chromosome"] = chrom
        guidescanTable["pam"] = PAM
        guidescanTable = guidescanTable[['sgId', 'sequence', 'pam', 'chromosome', 'position', 'strand']]
        guidescanTable.rename(columns={'sgId': 'id', "strand": 'sense'}, inplace=True)
        return guidescanTable
    else:
        return LibraryTable, LibraryInfoTable


def buildGuidePairs(LibraryInfoTable, footPrintRange=None):
    sgNameList = LibraryInfoTable.index.tolist()

    # Get all combinations of length 2
    pairs = [(sg1, sg2) for sg1, sg2 in combinations(sgNameList, 2)]

    df1 = LibraryInfoTable.loc[[sg1 for sg1, _ in pairs], ['position']].reset_index().add_prefix('sg1_')
    df2 = LibraryInfoTable.loc[[sg2 for _, sg2 in pairs], ['position']].reset_index().add_prefix('sg2_')

    out = pd.concat([df1, df2], axis=1)
    out['distance'] = out.sg2_position - out.sg1_position
    out.index = out.sg1_sgId + '::' + out.sg2_sgId

    print(f'All possible sgRNA pairs: {out.shape[0]}')

    # keep pairs if they present in given footPrintRange
    if footPrintRange:
        out = out[(out.distance.abs() >= footPrintRange[0]) & (out.distance.abs() <= footPrintRange[1])]
        print(f'Selected sgRNA pairs: {out.shape[0]}')

    return out
