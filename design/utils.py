"""Utility Functions
"""
import gzip
from Bio import SeqIO


def loadGenomeAsDict(genomeFasta):
    print('Loading genome file...')
    if '.gz' in genomeFasta:
        genomeDict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomeFasta, 'rt'), 'fasta'))
    else:
        genomeDict = SeqIO.to_dict(SeqIO.parse(genomeFasta, 'fasta'))
    print('Done\n')
    return genomeDict


def getPseudoIndices(table):
    return table.apply(lambda row: row.name[0][:6] == 'pseudo', axis=1)


def loadGencodeData(gencodeGTF, indexByENSG=True):
    pass
    # print('Loading annotation file...')
    # gencodeData = dict()
    # with open(gencodeGTF) as gencodeFile:
    #     for line in gencodeFile:
    #         if line[0] != '#':
    #             linesplit = line.strip().split('\t')
    #             attrsplit = linesplit[-1].strip('; ').split('; ')
    #             attrdict = {attr.split(' ')[0]: attr.split(' ')[1].strip('\"') for attr in attrsplit if
    #                         attr[:3] != 'tag'}
    #             attrdict['tags'] = [attr.split(' ')[1].strip('\"') for attr in attrsplit if attr[:3] == 'tag']
    #
    #             if indexByENSG:
    #                 dictKey = attrdict['gene_id'].split('.')[0]
    #             else:
    #                 dictKey = attrdict['gene_name']
    #
    #             # catch y-linked pseudoautosomal genes
    #             if 'PAR' in attrdict['tags'] and linesplit[0] == 'chrY':
    #                 continue
    #
    #             if linesplit[2] == 'gene':  # and attrdict['gene_type'] == 'protein_coding':
    #                 gencodeData[dictKey] = (
    #                     [linesplit[0], long(linesplit[3]), long(linesplit[4]), linesplit[6], attrdict], [])
    #             elif linesplit[2] == 'transcript':
    #                 gencodeData[dictKey][1].append(
    #                     [linesplit[0], long(linesplit[3]), long(linesplit[4]), linesplit[6], attrdict])
    #
    # print('Done\n')
    #
    # return gencodeData


def loadCageBedData(cageBedFile, matchList=['p1', 'p2']):
    pass
    # cageBedDict = {match: dict() for match in matchList}
    #
    # with open(cageBedFile) as infile:
    #     for line in infile:
    #         linesplit = line.strip().split('\t')
    #
    #         for name in linesplit[3].split(','):
    #             namesplit = name.split('@')
    #             if len(namesplit) == 2:
    #                 for match in matchList:
    #                     if namesplit[0] == match:
    #                         cageBedDict[match][namesplit[1]] = linesplit
    #
    # return cageBedDict


def matchPeakName(peakName, geneAliasList, promoterRank):
    pass
    # for peakString in peakName.split(','):
    #     peakSplit = peakString.split('@')
    #
    #     if len(peakSplit) == 2 \
    #             and peakSplit[0] == promoterRank \
    #             and peakSplit[1] in geneAliasList:
    #         return True
    #
    #     if len(peakSplit) > 2:
    #         print
    #         peakName
    #
    # return False


def generateAliasDict(hgncFile, gencodeData):
    pass
    # hgncTable = pd.read_csv(hgncFile, sep='\t', header=0).fillna('')
    #
    # geneToAliases = dict()
    # geneToENSG = dict()
    #
    # for i, row in hgncTable.iterrows():
    #     geneToAliases[row['Approved Symbol']] = [row['Approved Symbol']]
    #     geneToAliases[row['Approved Symbol']].extend(
    #         [] if len(row['Previous Symbols']) == 0 else [name.strip() for name in row['Previous Symbols'].split(',')])
    #     geneToAliases[row['Approved Symbol']].extend(
    #         [] if len(row['Synonyms']) == 0 else [name.strip() for name in row['Synonyms'].split(',')])
    #
    #     geneToENSG[row['Approved Symbol']] = row['Ensembl Gene ID']
    #
    # # for gene in gencodeData:
    # #   if gene not in geneToAliases:
    # #       geneToAliases[gene] = [gene]
    #
    # #   geneToAliases[gene].extend([tr[-1]['transcript_id'].split('.')[0] for tr in gencodeData[gene][1]])
    #
    # return geneToAliases, geneToENSG


def parseSgId(sgId):
    """Parse information from the sgRNA ID standard format
    """
    pass
    # parseDict = dict()
    #
    # # sublibrary
    # if len(sgId.split('=')) == 2:
    #     parseDict['Sublibrary'] = sgId.split('=')[0]
    #     remainingId = sgId.split('=')[1]
    # else:
    #     parseDict['Sublibrary'] = None
    #     remainingId = sgId
    #
    # # gene name and strand
    # underscoreSplit = remainingId.split('_')
    #
    # for i, item in enumerate(underscoreSplit):
    #     if item == '+':
    #         strand = '+'
    #         geneName = '_'.join(underscoreSplit[:i])
    #         remainingId = '_'.join(underscoreSplit[i + 1:])
    #         break
    #     elif item == '-':
    #         strand = '-'
    #         geneName = '_'.join(underscoreSplit[:i])
    #         remainingId = '_'.join(underscoreSplit[i + 1:])
    #         break
    #     else:
    #         continue
    #
    # parseDict['strand'] = strand
    # parseDict['gene_name'] = geneName
    #
    # # position
    # dotSplit = remainingId.split('.')
    # parseDict['position'] = int(dotSplit[0])
    # remainingId = '.'.join(dotSplit[1:])
    #
    # # length incl pam
    # dashSplit = remainingId.split('-')
    # parseDict['length'] = int(dashSplit[0])
    # remainingId = '-'.join(dashSplit[1:])
    #
    # # pass score
    # tildaSplit = remainingId.split('~')
    # parseDict['pass_score'] = tildaSplit[-1]
    # remainingId = '~'.join(tildaSplit[:-1])  # should always be length 1 anyway
    #
    # # transcripts
    # parseDict['transcript_list'] = remainingId.split(',')
    #
    # return parseDict


def parseAllSgIds(libraryTable):
    pass
    # sgInfoList = []
    # for sgId, row in libraryTable.iterrows():
    #     sgInfo = parseSgId(sgId)
    #
    #     # fix pam coordinates for -strand ??
    #     if sgInfo['strand'] == '-':
    #         sgInfo['pam coordinate'] = sgInfo['position']  # + 1
    #     else:
    #         sgInfo['pam coordinate'] = sgInfo['position']
    #
    #     sgInfoList.append(sgInfo)
    #
    # return pd.DataFrame(sgInfoList, index=libraryTable.index)
