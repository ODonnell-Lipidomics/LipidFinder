import pandas as pd
import numpy as np
import re


def performQCCalcs(lfd):
    """For the current 'processedDataFrame' the relative standard deviation (RSD) for the
    set of QC samples in each row is calculated and recorded. The number of sets where the
    RSD is lower than the lower QC RSD parameter (QCLowRSD) is recorded and the number of
    sets with an RSD lower than the higher QC RSD parameter (QCHighRSD) is also recorded.

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has NOT been previously feature clustered

    Returns:
        float: The ratio of the percentage of QC samples with an RSD less than QCLowRSD to the percentage of QC samples with an RSD less than QCHighRSD
    """
    # Column index of first QC sample
    indFirst = lfd.getFirstQCRepOffset()

    # Column index of last QC sample
    indLast = lfd.getLastQCRepOffset()

    # Insert mean and RSD columns for qc into dataFrame
    qcName = __colName__(lfd.processedDataFrame, indFirst)
    qcMeanName = qcName + '_mean'
    qcRSDName = qcName + '_RSD'

    # Get means and RSD of the qc samples
    means = lfd.processedDataFrame.ix[:, indFirst:indLast].mean(axis=1)

    # Pandas cython version of std is sample std (dof = 1) so using slower
    # numpy std as this is population (default ddof=0, degrees of freedom)
    RSDs = lfd.processedDataFrame.ix[
        :, indFirst:indLast].apply(np.std, axis=1) * 100 / means

    # Insert means and RSDs into processedDataFrame
    lfd.processedDataFrame.insert(indLast, qcMeanName, means)
    indLast += 1
    lfd.processedDataFrame.insert(indLast, qcRSDName, RSDs)

    lfd.qcData = lfd.processedDataFrame.copy()

    lfd.numberOfQCReps += 2

    # Calculate and report percentage lower RSD-QC Samples to Higher RSD-QC
    # Samples
    lowerRSDCount = (lfd.qcData[qcRSDName] < lfd.QCLowRSD).sum()
    upperRSDCount = (lfd.qcData[qcRSDName] < lfd.QCHighRSD).sum()
    propRSDLowerToUpper = round(lowerRSDCount / float(upperRSDCount) * 100, 1)
    return propRSDLowerToUpper


def __colName__(DataFrame, index):
    name = DataFrame.ix[:, index].name
    return re.sub('\d+$', "", name)
