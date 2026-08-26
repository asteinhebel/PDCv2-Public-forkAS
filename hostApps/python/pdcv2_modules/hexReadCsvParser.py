import numpy as np
import pandas as pd
from enum import IntEnum

# ------------------------------------------------------------------------
# --- Constants
# ------------------------------------------------------------------------
DSUM_CSV_SEP = ';'
DSUM_CSV_COL_NAMES = ["pdcIdx", "frameIdx", "dataIdx", "dsum"]
#DSUM_CSV_COL_TYPES = [np.uint8, np.uint32,  np.uint8,  np.uint16]
DSUM_CSV_COL_TYPES = [np.uint8, np.uint32,  np.uint32,  np.uint16]
DSUM_CSV_HEADER_LINE = DSUM_CSV_SEP.join(DSUM_CSV_COL_NAMES)
DSUM_CSV_COMMENT_CHAR = "#"

# ------------------------------------------------------------------------
# --- Functions
# ------------------------------------------------------------------------
def parseCsv(filename, sep=DSUM_CSV_SEP, comment=DSUM_CSV_COMMENT_CHAR):
    """
    Parse a csv file generated with hexRead.
    Set the proper dataType for each variable to reduce memory usage.
    """
    data = pd.read_csv(filename,
                       sep=sep,
                       comment=comment,
                       dtype=dict(zip(DSUM_CSV_COL_NAMES, DSUM_CSV_COL_TYPES)))
    if data["dataIdx"].max() < 256:
        # Default behavior is dataIdx from 0 to 127 (128 bins in PDC memory).
        # If using 'time' option, dataIdx is a timestamp (32 bits).
        # Optimize memory usage if possible.
        data["dataIdx"] = data["dataIdx"].astype(np.uint8)
    return data


def get_file_comments(filename,
                      comment=DSUM_CSV_COMMENT_CHAR,
                      keepCommentChar=False):
    """
    Associate each line with its comment.
    The comment associated with data is the last commented line before data lines
    """

    # first get the total number of lines in file
    with open(filename, "rb") as f:
        numLines = sum(1 for _ in f)

    # init an array to hold the line type of each data
    class LineId(IntEnum):
        empty = 0,
        header = 1,
        comment = 2,
        data = 3
    npLineId = np.zeros(numLines, dtype=np.uint8)

    # init list of comments with "N/A" for files without comments or data before first comment
    comments = ["N/A"]
    with open(filename, 'r', encoding='utf-8') as f:
        # skip header line
        for i, line in enumerate(f, 0): # 0-based indexing
            if line.startswith(DSUM_CSV_HEADER_LINE):
                # header line found
                npLineId[i] = LineId.header

            elif line.startswith(comment):
                # Line must start with a comment, because comment can be placed after data.
                # A comment after data is ignored
                npLineId[i] = LineId.comment
                if not keepCommentChar:
                    # replace only first occurence
                    line = line.replace(comment, '', 1).lstrip()
                comments.append(line.rstrip('\n'))
            elif not line.isspace():
                npLineId[i] = LineId.data

    commentIdx = (npLineId == LineId.comment).cumsum()[npLineId == LineId.data]

    return commentIdx, pd.Series(comments)

def comments_as_column(comments, commentsIdx):
    """
    From a list of comments found in the file, associated each row of data with its comment.
    """
    return comments[commentsIdx].values
