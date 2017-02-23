import pandas as pd
import numpy as np


def reassignFrameMasses(lfd):
    """Assigns each mass in either a mass cluster or feature cluster to the mass of the
    row containing the highest sample mean intensity

    Args:
        lfd (LipidFinderData): A LipidFinderData object where the sample replicate means have been calculated.

    """
    if lfd.featureLevelMassAssignment:
        lfd.processedDataFrame = lfd.processedDataFrame.groupby(
            'Feature Cluster ID').apply(__maxGroupMass__, meanColList=lfd.meanColList)
    else:
        lfd.processedDataFrame = lfd.processedDataFrame.groupby(
            'mz Cluster ID').apply(__maxGroupMass__, meanColList=lfd.meanColList)
    lfd.reassignedFrameMasses = lfd.processedDataFrame.copy()


def __maxGroupMass__(groupMass, meanColList):
    maxIntInd = groupMass.ix[:, meanColList].max(axis=1).argmax()
    groupMass.ix[:, 'MZ'] = groupMass['MZ'][maxIntInd].copy()
    return groupMass
