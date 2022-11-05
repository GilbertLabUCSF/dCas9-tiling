import pandas as pd
from itertools import combinations


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
