from __future__ import print_function
import math
class InvalidRangeException():
    pass

class Point:
    X = None
    Y = None

    def __init__(self, x, y):
        self.X = x
        self.Y = y


class Line:
    """An object representing a line made of 2 points."""

    Start = None
    End = None

    def __init__(self, start, end):
        """
        Construct the line from 2 points. Points must implement an attribute X and an attribute Y.
        :param start: The start point.
        :param end: The end point.
        """
        self.Start = start
        self.End = end

    def GetPointOnLine(self, ratio):
        """
        Returns the points at ratio-defined position.
        :param ratio: The ratio. Must be between 0 and 1.
        0 will return the starting point, 1 will return the endpoint.
        :return: A point
        """
        if(ratio < 0 or ratio > 1):
            raise InvalidRangeException()

        divider = 1 / ratio
        return Point((self.Start.X + self.End.X) / divider, (self.Start.Y + self.End.Y) / divider)

    def DivideLineIntoEqualsParts(self, totalParts):
        """Divide a line into an equal number of parts.
        :param: totalParts. The total number of part in which this line will be divided.
        :return: pointsSequence. An array of points. Points are returns as an array of shape [X,Y]
        """
        #http://stackoverflow.com/questions/3542402/partition-line-into-equal-parts
        pointsSequence = []
        if totalParts < 1:
            raise Exception("The number of parts must be greater than 1")

        step = 1 / totalParts
        totalSteps = 0
        while totalSteps <= 1:
            pX = self.Start.X * (1 - totalSteps) + self.End.X * totalSteps
            pY = self.Start.Y * (1 - totalSteps) + self.End.Y * totalSteps
            pointsSequence.append(Point(pX,pY))
            totalSteps = totalSteps + step

        return pointsSequence


    def DistanceSquared(self):
        """Returns the squared length of the line.
        :return: a float representing the squared length of the line.
        """
        return math.pow(self.Start.X - self.Start.Y, 2) + math.pow(self.End.X - self.End.Y, 2)

    def Distance(self):
        """Returns the length of the line"""
        return math.sqrt(self.DistanceSquared())

    def FindClosestPointOnLine(self, point):
        """
        Return the point on a line that is the closest to another point
        :param point:
        :return:
        """
        #http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment

        #If the line is an invalid geometry where both end points are the same, we can return any of the line point
        squaredDistance = self.DistanceSquared()

        if squaredDistance == 0:
            return [self.Start.X, self.Start.Y]

        ratio = ((point.X - self.Start.X) * (self.End.X - self.Start.X) + (point.Y - self.Start.Y) * (self.End.Y - self.Start.Y)) \
            / squaredDistance

        ratio = 1 if 1 < ratio else ratio
        ratio = ratio if ratio > 0 else 0

        return Point(
            self.Start.X + ratio * (self.End.X - self.Start.X),
            self.Start.Y + ratio * (self.End.Y - self.Start.Y)
        )



    def CurveSegment(self, point, otherLine, parts):

        curvePoints = [self.Start]
        #Get the list of intermediate points for segmentation
        points = self.DivideLineIntoEqualsParts(parts)

        #Iterate through points
        for p in points[1:-1]:
            #Get the closest point to the other line
            closestPointOnLine = otherLine.FindClosestPointOnLine(p)

            #Create 2 lines. connecting the intermediate points with the closest point and the closest line.
            pointToLine = Line(p, closestPointOnLine)
            pointToPoint = Line(point, p)

            #Get the distance-related information
            distanceToLine = pointToLine.Distance()
            distanceToPoint = pointToPoint.Distance()
            totalDistance = distanceToLine + distanceToPoint

            #Find the distance where the mid-point should be
            midDistance = totalDistance / 2

            #Set up logic to compute the line where the mid-distance point lies and where to find it.
            usePointToLine = True if midDistance > distanceToPoint else False
            locatedLine = pointToLine if usePointToLine == True else pointToPoint
            locatedDistance = distanceToPoint if usePointToLine == False else math.fabs(distanceToLine - distanceToPoint)
            pointRatio =  locatedDistance / locatedLine.Distance()

            #Get the point
            curvePoints.append(self.GetPointOnLine(pointRatio))

        curvePoints.append(self.End)
        return curvePoints