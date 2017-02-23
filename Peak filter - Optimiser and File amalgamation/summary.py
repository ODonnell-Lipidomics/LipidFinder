import pandas as pd
import numpy as np


def createSummaryDataFrame(lfd):
    """Create a summary DataFrame conatining only mean intensity for each sample and the
    M/Z, Time, Polarity and Column-type.

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has been fully processed and is ready for file amalgamation.

    """
    # Select out just mz and time from processed DataFrame - numbers will
    # always be the same
    temp = lfd.processedDataFrame.ix[
        :, lfd.firstRepOffset - 4:lfd.firstRepOffset - 2].copy()

    temp['Polarity'] = 'POS' if lfd.filePolarityMode == 'P' else 'NEG'
    temp['Column-type'] = 'Polar' if lfd.columnType == 'PO' else 'Non-polar'
    temp = pd.concat([temp, lfd.processedDataFrame.ix[
                     :, lfd.meanColList]], axis=1)
    # As we have taken out 2 and added 2 columns then lfd.firstRepOffset
    # reamins correct

    # reset index after drop
    temp.reset_index(inplace=True, drop=True)
    # Drop frames outside rt cutoffs
    temp.drop(np.where(((temp['Time'] < lfd.retentionTimeLowCutOff) | (
        temp['Time'] > lfd.retentionTimeHighCutOff)))[0], inplace=True)
    temp.sort_values(by=['MZ', 'Time'], inplace=True)
    lfd.summaryDataFrame = temp.reset_index(drop=True)
