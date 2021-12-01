"""
Geometry functions.
"""
from typing import List, Tuple

def calculateCentroid(
        pointCloud : List[Tuple[float, float, float]]
    ) -> Tuple[float, float, float]:
    """Calculate centroid of point cloud.

    Arguments
    --------------------------------------------------------------------------
    pointCloud (float 3-tuple list) -- list of xyz coordinates.

    Returns
    --------------------------------------------------------------------------
    centroid (float 3-tuple) -- centroid of points in point cloud.
    """
    numPoints = len(pointCloud)
    x, y, z = [], [], []
    for point in pointCloud:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    x, y, z = sum(x) / numPoints, sum(y) / numPoints, sum(z) / numPoints
    return x, y, z
