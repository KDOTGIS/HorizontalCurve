#-------------------------------------------------------------------------------
# Name:        KDOT_HorizontalCurveTool_Functions
# Purpose:     Principal functions for curve analyses
# Author:      URS Corporation
# Created:     2014
# Copyright:   (c) URS Corporation 2014

#-------------------------------------------------------------------------------
"""This module contains the functions and classes for this program"""
import os, sys, arcpy
import numpy as np
import pylab as plt
from numpy import linalg
from math import sin, cos, pi, exp, log, sqrt
import random, math

################################################################################
# Classes for Curve Fitting Program
class Line(object):
    """ Line class

    featID contains the original OID of the road segment
    name contains the Unique_ID of the road segment
    vertices contains every vertex in the road segment as a list of lists [[x1, y1, z1, m1], ... [xN, yN, zN, mN]]"""

    # method to initialize Line
    def __init__(self, featID, name, vertices):
        # OID of the road feature
        self.featID = featID
        # name contains the unique identifier string
        self.name = name
        # vertices contains [[x1, y1, z1, m1], ... [xN, yN, zN, mN]]
        self.vertices = vertices

class Curve(Line):
    """ A class to store the curve geometeries

    featID contains the original OID of the road segment
    name contains the Unique_ID of the road segment
    vertices contains every point in the road segment (even if it is not within a curve)
    curves constains a dictionary of curves along a road segment
    parameters contains a dictionary of curve parameters including in part the radius, center, and r-squared value """

    def __init__(self, featID, name, vertices, curves, parameters, curve_names, curve_directions):
        super(Curve, self).__init__(featID, name, vertices)
        # curves is a dictionary of keys and values where the keys are integers corresponding to curve number, and range of points (index values) within the curve
        # {1: [iStart1:iEnd1], 2: [iStart2:iEnd2], ..., N: [iStartN:iEndN]}
        self.curves = curves
        # params is a dictionary of keys and values where the keys are integers corresponding to the curve number, and values are empty initally
        # {1: [issolved1, radius1, center1, R-squared1], 2: [issolved2, radius2, center2, R-squared2], ..., N:  [issolvedN, radiusN, centerN, R-squaredN]}
        self.parameters = parameters
        #curve_names is a dictionary of keys and values where the keys are the same for the other members, and the values are the curve names:
        # {1: "curve1", 2: "curve2", ..., N: "curveN"}
        self.curve_names = curve_names
        #curve_directions is a dictionary of keys and values where the keys are the same for the other members, and the values are integers indications curve directions:
        # -1 is counter-clockwise, 1 is clock-wise. This is defined as a polar coordinate system where increasing values are counter-clockwise.
        # {1: direction1, 2: direction2, ..., N: directionN}
        self.curve_directions = curve_directions

################################################################################
def detailed_errors(thebool):
    # If TRUE, create detailed debug messages:
    global VERBOSE
    VERBOSE = thebool
    return VERBOSE

################################################################################
#Make a copy of a feature class with an EXTRA field containing the ORIGINAL object id
def make_copies_of_features(inputfeatures,  outputfeatures, thefield):
##    arcpy.AddMessage("funcs.make_copies_of_features")
    template = inputfeatures
    descin = arcpy.Describe(inputfeatures)
    inpath = descin.path
    inname = descin.name
##    arcpy.AddMessage('inpath: '+inpath)
##    arcpy.AddMessage('inname: '+inname)
    #Create the output feature class:
##    arcpy.AddMessage(outputfeatures)
    outpath, outname = os.path.split(outputfeatures)
##    arcpy.AddMessage('outpath: '+outpath)
##    arcpy.AddMessage('outname: '+outname)

    arcpy.CreateFeatureclass_management(outpath, outname, 'POLYLINE', template, 'ENABLED', 'DISABLED')
    arcpy.AddField_management(outputfeatures, thefield, "LONG", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    outfieldlist = arcpy.ListFields(outputfeatures)
    infieldlist = arcpy.ListFields(inputfeatures)

    endlist = []
    for fld in infieldlist:
        if(fld.type != "OID" and (fld.type != 'Geometry')):
            endlist.append(fld.name)
    #Create the inputfeatures field list. The first two members must be OID and SHAPE@
    sourcelist_names = ['OID@', 'SHAPE@'] + endlist
##    arcpy.AddMessage('sourcelist_names: '+ str(sourcelist_names))
    outlist_names = sourcelist_names + [thefield]

    ocursor = arcpy.da.InsertCursor(outputfeatures, outlist_names)
    with arcpy.da.SearchCursor(inputfeatures, sourcelist_names ) as icursor:
        for irow in icursor:
        #Create a new tuple by adding the value to it the value of the source OID@
##            arcpy.AddMessage(type(irow))
            #irow[0] is the OID value of the input data and is written to field "thefield"
            orow = irow + (irow[0],)
            newid = ocursor.insertRow(orow)
    return

################################################################################
def make_workspace_copy(inputfeatures,theworkspace,dotopologycheck,dosimplify,dosimplify_method, dosimplify_tolerance, thefield):
    """This function tests the input features for the topology error 'Must Be Single Part',
    and returns the Origin Feature's Object ID of the errant features to the calling module. Beware:
    the origing feature's object ID is that of the COPY of the input features. The object ID's of the copy
    may be different from the source "inputfeautures"!!!!. This is why the function passes back the name of the COPY so that the processing can
    continue on that feature class where the topologically errant features will be correctly identified
    by the values in the export topology errors geoprocessing tool."""

##    arcpy.AddMessage("funcs.make_workspace_copy")

    #Process the
    #roads with the simplify_line tool with the point_remove option at a tolerance of 0.001 meters so that redundant vertices on staight lines are removed.
    #If the user specifies their own parameters for simplify_line, THAT ARE NOT POINT_REMOVE AND THE TOLERANCE IS > 0.001 METERS, that is done additionally,
    #afterwards:

    #this section makes the feature class datasets, feature class names, and topology name:
    badfids =set()
    fdname = "KDOT_Topology_Check" #the feature dataset name for the topology check
    fdnamepath = theworkspace + "\\"+ fdname #the complete pathname of the feature dataset
    tpname = "CheckTopology" #the name of the topology
    topology_name = fdnamepath + "\\" + tpname #the complete pathname of the topology
##    arcpy.AddMessage("make_workspace_copy, fdnamepath: "+fdnamepath)
##    arcpy.AddMessage("make_workspace_copy, topology_name: "+topology_name)
    fcname = arcpy.ParseTableName(inputfeatures, theworkspace) #Split the inputfeatures to find the name from the path.
    namelist = fcname.split(", ") #the feature class name without the path. Used in creating a copy in the feature dataset.
##    arcpy.AddMessage('fcname = '+ namelist[2])
    topology_featureclass = fdnamepath +'\\' + namelist[2] + '_check' #the copy of inputfeatures used for the topology check
    topology_featureclass_errors = namelist[2] + '_errors' # the basename used for the export topology errors tool
##    arcpy.AddMessage(topology_featureclass)
    topology_featureclass_errors_line = fdnamepath +'\\' + namelist[2] + '_errors_line' #the output name of LINE errors from the export topology errors tool

    #Delete if the feature dataset currently exists:
    doesexistfd = arcpy.Exists(fdnamepath)
    try:
       if doesexistfd:
           arcpy.AddMessage('Previous topology check feature dataset exists. Now deleteing ')
           arcpy.Delete_management(fdnamepath)
    except arcpy.ExecuteError:
        print arcpy.GetMessages(2)
    except Exception as e:
        print e.args[0]

    #Re-create the topology feature dataset:
    arcpy.AddMessage('Generating the topology check scratch feature dataset')
    arcpy.CreateFeatureDataset_management(theworkspace, fdname, inputfeatures)

    #Make a copy of the input roads in the feature dataset that contains the topology:
    try:
        arcpy.AddMessage('Generating a copy of the input feature class in the scratch feature dataset')
        #This replaces the function "arcpy.CopyFeatures_management" so that we can retain the original FID:
##        make_copies_of_features(inputfeatures,  topology_featureclass, "Original_OID")
        make_copies_of_features(inputfeatures,  topology_featureclass, thefield)
##        arcpy.CopyFeatures_management(inputfeatures, topology_featureclass)
    except arcpy.ExecuteError:
        print arcpy.GetMessages(2)
    except Exception as e:
        print e.args[0]

    #Perform the topology check, if checked ON in input parameters:
##    if(VERBOSE): arcpy.AddMessage('make_workspace_copy, dotopology = ' + str(dotopologycheck))
##    if(dotopologycheck == True):
    if(str(dotopologycheck) == 'true'):
        arcpy.AddMessage('Creating the topology')
        arcpy.CreateTopology_management(fdnamepath, tpname)

        #Add the input roads to the topology
        arcpy.AddMessage('Adding the copy of the input features to the topology')
        arcpy.AddFeatureClassToTopology_management(topology_name, topology_featureclass, 1, 1)
        #Add a rule:
        arcpy.AddMessage('Adding rule "Must Be Single Part" to the topology')
        arcpy.AddRuleToTopology_management(topology_name,"Must Be Single Part (Line)", topology_featureclass)
        #Validate the topology:
        arcpy.AddMessage('Validating the topology')
        arcpy.ValidateTopology_management(topology_name)
        #Export the errant features to a feature class
        arcpy.AddMessage('Exporting the topologically-errant feature to feature class ' + topology_featureclass_errors)
        arcpy.ExportTopologyErrors_management(topology_name,fdnamepath,topology_featureclass_errors)
        arcpy.AddMessage("Completed exporting topology errors")

        #Extract the values from field "OriginObjectID". This is a field generated to identify the OID's of errant features. This is the ORIGINAL OID of the input features. Do NOT confuse this with "OriginalOID":
##        arcpy.AddMessage('Retrieving the object ID''s of the errant features')
        with arcpy.da.SearchCursor(topology_featureclass_errors_line,["OriginObjectID"]) as cursor:
            for row in cursor:
##                arcpy.AddMessage(str(row[0]))
                badfids.add(row[0])

    #Perform at the least, the default line simplification of 0.001 meters or 0.00328084 feet
    #SimplifyLine(mergedFeatures, simplifiedFeatures, dosimplify_method, dosimplify_tolerance, "RESOLVE_ERRORS", "KEEP_COLLAPSED_POINTS", "CHECK")
    simplified_featureclass = fdnamepath +'\\_simplified_roads'
    arcpy.SimplifyLine_cartography(topology_featureclass, simplified_featureclass, dosimplify_method, dosimplify_tolerance, False, False, False )

##    index_OID_Original, fieldlist =  ListFieldData(simplified_featureclass, thefield)
##    if(VERBOSE): arcpy.AddMessage("List of indecies and field names after copy, with original OID, and simplification of line:")
##    if(VERBOSE): arcpy.AddMessage("index for original OID field: "+thefield+" is "+str(index_OID_Original))

    arcpy.AddMessage('completed creating a workspace copy....')
##    arcpy.AddMessage('completed funcs.make_workspace_copy')
    return badfids, simplified_featureclass

################################################################################
def list_duplicates(seq):

    arcpy.AddMessage("funcs.list_duplicates")
    #http://stackoverflow.com/questions/9835762/find-and-list-duplicates-in-python-list
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set( x for x in seq if x in seen or seen_add(x) )
    # turn the set into a list (as requested)
    return list( seen_twice )

################################################################################
def break_linefeatures_into_Line_class_with_maxdistance(theFeatureClass, maxDistance, featurename_field, orig_id_field, badfids):
    """ Break line features into their point geometries and populate Line class

    theFeatureClass is and ESRI ArcGIS feature class """
    arcpy.AddMessage('Breaking lines into vertices.............')

    roadList = []

##    if(len(badfids) > 0):
##        arcpy.AddMessage("Features not used in analysis:")
##        arcpy.AddMessage("Original FID, road name")

    # Enter for loop for each feature
    # Remember that the SearchCursor function returns the row(columnindex) for the indicated fields;
    # here column 0 ,row(0) is OID, column 1, row(1) is SHAPE, column 2 (row(2) is the name for the unique road ID:
    #...................................................row(0), row(1),  row(2),            row(3)
    for row in arcpy.da.SearchCursor(theFeatureClass, ["OID@", "SHAPE@", featurename_field, orig_id_field]):
        # Print the current multipoint's ID
        #
        vList = []
##        arcpy.AddMessage("Feature {0}:".format(row[2]))
        featIDcopy = row[0]     #This is the OID of the copy, not the OID of the original source data.
        featID = row[3]     #This is the OID of the original source data.
        name = str(row[2])

        #Avoid processing the features that do not conform to the topology rule "must be singlepart", if chosen:
        if(featIDcopy in badfids):
            arcpy.AddMessage("Features not used in analysis (Original FID, road name): "+str(featID)+", "+name)
            continue

        #Proceed to processing the roads that are accepted:
        partnum = 0

        # Step through each part of the feature
        #
        for part in row[1]:
            # Print the part number
            #
##            arcpy.AddMessage("Part {0}:".format(partnum))

            # Step through each vertex in the feature
            #
            for pnt in part:
                if pnt:
                    # Print x,y coordinates of current point
                    #
##                    arcpy.AddMessage("{0}, {1}".format(pnt.X, pnt.Y, 0, pnt.M))
                    vList.append([pnt.X, pnt.Y, 0, pnt.M])
                else:
                    # If pnt is None, this represents an interior ring
                    #
##                    arcpy.AddMessage("Interior Ring:")
                    continue
            partnum += 1

        fList = []
        cumul_M = 0.0
        for i in xrange(len(vList) - 1):
            # X-coordinates
            x1 = vList[i][0]
            x2 = vList[i+1][0]

            # Y-coordinates
            y1 = vList[i][1]
            y2 = vList[i+1][1]

            # Calculate cumulative M value
            cumul_M = cumul_M + vList[i][3]

            # Calculate the distance between the two points
            d, dx, dy = calculate_distance(x1, y1, x2, y2)
            # if d <= maxDistance, then append current vList point to fList
##            arcpy.AddMessage("Distance: " + str(d))
            if d <= maxDistance:
                fList.append(vList[i])
            # if d > maxDistance, then calculate the insertion point(s), and then append the point(s) and the current vList point to fList
            else:
                fList.append(vList[i])
                # determine how many points to add between the two points where max distance is exceeded
                num_pts = int(d/maxDistance)
##                arcpy.AddMessage("Num point: " + str(num_pts))
                for i in range(num_pts):
                    maxD = maxDistance * (i + 1)  # add 1 here because otherwise the first i value will be zero, and the last i value will be (num_pts - 1)
##                    arcpy.AddMessage("MaxD: " + str(maxD))
                    new_pt = find_max_distance_point(x1, y1, x2, y2, maxD)
##                    arcpy.AddMessage("New Pt: " + str(new_pt))
                    fList.append([new_pt[0], new_pt[1], 0, cumul_M])
                    cumul_M = cumul_M + maxDistance
##        fList.append(vList[i+1])
##        arcpy.AddMessage(fList)
        # Test to see if the feature has more than 2 vertices... if it does, make the list of verticies (vList) into a Line instance
        if len(fList) > 2:
            # Create a Line class instance from the list of vertices for the current feature
            roadList.append(Line(featID, name, fList))
        else:
            arcpy.AddMessage(str(featID)+", "+name+", Line has only 2 vertices, skipped")
##            arcpy.AddMessage('Line has only 2 vertices, skipped')
        if(VERBOSE):
            arcpy.AddMessage("vlist output. Vertices before implemenatation of max_distance:")
            arcpy.AddMessage("x,y,z,m")
            for i, item in enumerate(vList):
              arcpy.AddMessage(str(item[0])+", "+str(item[1])+", "+str(item[2])+", "+str(item[3]))
            arcpy.AddMessage("flist output. Vertices before implemenatation of max_distance:")
            arcpy.AddMessage("x,y,z,m")
            for i, item in enumerate(fList):
              arcpy.AddMessage(str(item[0])+", "+str(item[1])+", "+str(item[2])+", "+str(item[3]))


    # return a list of Line instances
    return roadList

################################################################################
def populate_curves(featID, name, vertices, iStart, iEnd, geoList):
    """ Populate the curve information into the Curve class

    featID is the OID of the current Line instance from the original feature class
    name is the name of the current Line instance
    verticies are the vertices of the current Line instance as a list of lists [[x1, y1, z1, m1], ... [xN, yN, zN, mN]]
    iStart is the curve starting index in integer format
    iEnd is the curve ending index in integer format"""
##    arcpy.AddMessage("funcs.populate_curves")

    # zip the curve number and indices into a list of format (N, (iStart, iEnd))
    values = zip(iStart,iEnd)
    keys = range(1,(len(iStart)+1))

    # Create a dictionary from the list of the following format
    # {1: [iStart1, iEnd1], 2: [iStart2, iEnd2], ..., N: [iStartN, iEndN]}
    curves = dict(zip(keys, values))

    # Create an empty dictionary of the following format to be updated later
    # {1: [issolved1, radius1, center1, R-squared1], 2: [issolved2, radius2, center2, R-squared2], ..., N:  [issolvedN, radiusN, centerN, R-squaredN]}
    default_value = [None, None, None, None, None]
    parameters = dict.fromkeys(keys, default_value)
    curve_names = {}
    curve_directions = {}
##    for i, item in enumerate(zip(iStart, iEnd)):
    for i, item in enumerate(values):
##        cname = str(item[0]) + "_" + str(item[1]) + "_" + name
##        arcpy.AddMessage("i in calc curve_names and curve_directions: "+str(i))
        curve_names[i+1] = str(item[0]) + "_" + str(item[1]) + "_" + str(featID)
        curve_directions[i+1] = geoList[item[0]+1][3]

    # return the Curve instance
    return Curve(featID, name, vertices, curves, parameters, curve_names, curve_directions)

################################################################################
def calculate_distance(x1, y1, x2, y2):
    """ Calculate the distance between two points

    x1 is the first x-coordinate as float
    y1 is the first y-coordiante as float
    x2 is the second y-coordinate as float
    y2 is the second y-coordinate as float """
##    arcpy.AddMessage("funcs.calculate_distance")

    # Difference between x- and y- locations
    dx = x2 - x1
    dy = y2 - y1

    # Calculate the distance between the two points
    d = np.sqrt( (dx * dx) + (dy * dy) )

    # return the distance as a float
    return (d, dx, dy)

################################################################################
def find_max_distance_point(startX, startY, endX, endY, r):
    """ Calculate the end point when the maximum distance is exceeded for curve end

        startX, startY is the first point
        endX, endY is the second point
        r is the distance along a line in float format """

##    arcpy.AddMessage("funcs.find_max_distance_point")
    angle_md = polar_angle(startY, startX, endY, endX)
    x = math.cos(angle_md) * r + startX
    y = math.sin(angle_md) * r + startY

    # Convert the coordinates into a list called endpoint
    newpt = [x, y]

    # Return the endpoint as a list of length 1
    return newpt

################################################################################
def curve_flag(a,b):
    """ test to see if a curve has started """

    # set flag for in or out of curve.
    if a == b:
        # if cflag = 0, then a curve has not started (i.e. the length of the curveStart and curveEnd lists are equal)
        cflag = 0
    else:
        # if cflag = 1, then a curve HAS started (i.e. a curveStart point has been put in it's list, but the length of the curveEnd list has not changed, so the lists are of unequal length)
        cflag = 1

    return cflag

################################################################################
def detect_curves(road, geoList, threshold):
    """ Pinpoints curve start and end points

    road is the Line class instance
    geoList is a list of lists in format: [i, theta, change_angle, direction, pol_angle_a, pol_angle_b, d, dx, dy]
    threshold is a float value indicating the angle below which curvature in a Line instance will be ignored

    Rules:

        1. If a curve has not started, and a change in angle rises above the threshold, and the point at which this happens is not the last point on a line, then this is the start of a curve
        2. If a curve has started, and a change in  angle falls below the threshold, then this is the end of the curve
        3. If a curve has started, and the last change in angle on a line is > threshold (the change in angle does not fall back under the threshold from above), then this is the end of the curve   """
##    arcpy.AddMessage("funcs.detect_curves")

    curveStartX = []
    curveStartY = []
    curveEndX = []
    curveEndY = []
    kStart = []
    kEnd = []

    slength = len(curveStartX)
    elength = len(curveEndX)
    glength = len(geoList)-1

    for g in geoList:
        # index of the road segment
        i  = g[0]

        # angle change value for curve detection test
        change_angle = np.abs(g[2])

        # direction of the curve. This section needs to be explained.
        if i == 0:
            direction1 = 1
        else: direction1 = geoList[i-1][3]
        direction2 = g[3]

        # set flag for in or out of curve.
        cflag = curve_flag(slength, elength)

        # Find curve start
        # If a curve has not started, and a change in angle rises above the threshold, and the point at which this happens is not the first point on a line, then this is the start of a curve
        if cflag == 0 and i != 0 and (change_angle >= threshold):

            # append curve start x and y to the appropriate list
##            curveStartX.append(road.vertices[i][0])
##            curveStartY.append(road.vertices[i][1])
##            kStart.append(i)
##            slength = len(curveStartX)
##            elength = len(curveEndX)
##            # set flag for in or out of curve.
##            cflag = curve_flag(slength, elength)

            # append curve start x and y to the appropriate list
            curveStartX.append(road.vertices[i+1][0])
            curveStartY.append(road.vertices[i+1][1])
            kStart.append(i+1)
            slength = len(curveStartX)
            elength = len(curveEndX)
            # set flag for in or out of curve.
            cflag = curve_flag(slength, elength)

        # Find curve end
        # If a curve has started, and a change in angle falls below the threshold or the direction goes from clockwise to counter-clockwise or vice versa, then this is the end of the curve

        if cflag == 1 and i > 0 and ((change_angle < threshold) or (direction1 != direction2) or (i == glength)):
            # measure number of points in the potential curve.  if threePTID is 3 or less, it does not count as a curve
            threePtID = i - kStart[-1]
##            if threePtID >= 3:
            if threePtID >= 2:
##                arcpy.AddMessage("Non-3pt: " +str(threePtID))
                # append curve end x and y to the appropriate list
                curveEndX.append(road.vertices[i][0])
                curveEndY.append(road.vertices[i][1])
##                curveEndX.append(road.vertices[i+1][0])
##                curveEndY.append(road.vertices[i+1][1])
                kEnd.append(i)
                slength = len(curveStartX)
                elength = len(curveEndX)
                # set flag for in or out of curve.
                cflag = curve_flag(slength, elength)

            else:
##                arcpy.AddMessage("3 pts: " + str(kStart[-1]) + " " + str(i))
                curveStartX.pop()
                curveStartY.pop()
                kStart.pop()
                slength = len(curveStartX)
                elength = len(curveEndX)
                # set flag for in or out of curve.
                cflag = curve_flag(slength, elength)


##    arcpy.AddMessage(str(curveStartX) + "," + str(curveStartY))
##    arcpy.AddMessage(str(curveEndX) + "," + str(curveEndY))
##    arcpy.AddMessage("kstart: " + str(kStart))
##    arcpy.AddMessage("kend: " + str(kEnd))

    # return six lists
    # curveStartX is a list of curve start x-coordinates as floats
    # curveStartY is a list of curve start y-coordinates as floats
    # curveEndX is a list of curve end x-coordinates as floats
    # curveEndY is a list of curve end y-coordinates as floats
    # kStart is a list of the curve start point index values as integers
    # kEnd is a list of the curve end point index values as integers
    return (curveStartX, curveStartY, curveEndX, curveEndY, kStart, kEnd)

################################################################################
def convert_pointList_to_x_and_y_lists(ptList, flag):
    """ Convert a list of lists to a different format for specific functions

    ptList is a list of lists of every vertex in the current Line instance [[x1,y1,z1,m1], ... [xN,yN,zN,mN]] format
    flat is either the integer 1 or 2 and indicates the format the points should be exported as """

##    arcpy.AddMessage("funcs.convert_pointList_to_x_and_y_lists")
    easting = []
    northing = []

    # Read through the list of points, and append x and y values to separate lists
    for pt in ptList:
        easting.append(pt[0])
        northing.append(pt[1])

    # Return the x and y values in a specific format
    if flag == 1:
        # Return x and y values as two separate lists: [x0, x1, x2, ... xN] and [y0, y1, y2, ... yN]
        return (easting, northing)

    if flag == 2:
        # Return x and y values as a list of lists: [ [x0, y0], [x1, y1], [x2, y2], ... [xN, yN] ]
        return (zip(easting,northing))

    if flag != 1 and flag != 2:
        # If the flag is anything other than 1 or 2, the program generates an error because it does not have instructions on what to do
        return (arcpy.AddMessage("error in funcs.convert_pointList_to_x_and_y_lists.  check that flag value is either 1 or 2"))

################################################################################
def randomize_vertices(data, scale):
    """ add a small amount of random noise to curves that fit exactly to a perfect circle.  this is necessary because the hyperfit
    algorithm blows up if the fit between the observed and ideal curve is exact.

    data is a list of vertices as tuples [(2343119.0160308033, 1148169.4971354753), (2343159.3112258017, 1148167.0368385613)] """

##    arcpy.AddMessage("funcs.randomize_vertices")
    datarand = []
    for vtx in data:
        # Get the decimal portion of each coordinate value
        thefractionX = vtx[0]-int(vtx[0]); theintX = int(vtx[0])
        thefractionY = vtx[1]-int(vtx[1]); theintY = int(vtx[1])

        # Calculate a random perturbation of the decimal point using the predefined scale
        pt = random.uniform(thefractionX*-1, thefractionX)/scale + thefractionX + vtx[0], random.uniform(thefractionY*-1, thefractionY)/scale + thefractionY + vtx[1]

        # Add each perturbed point back to a list of data
        datarand.append(pt)

    # return the perturbed list of data for curve fitting by fit_hypersphere function
    return datarand

################################################################################
def hypersphere(c,rscale):
    """ For each curve in the Curve class instance parse out the x, and y values then pipe them into fit_hypersphere function

    c is the current Curve instance """

##    arcpy.AddMessage('funcs.hypersphere')

    # For each curve in the Curve class, read the index values for start and end of curve
    # 'n' is a counter variable to reference the key of the current curve object in the parameters dictionary
    # Suggestion: use iteritem() instead to get the key
    n = 1
##    arcpy.AddMessage(c.curves)
    for v in c.curves.values():        #Get the values (start and end of each curve in a line). c.curves is a dictionary.
        cStart =v[0]
        cEnd = v[1]+1

        # parse out the x, y values from c.vertices within the Curve class instance and store the actual coordinate points in the list "data"
        data = []
        for j in c.vertices[cStart:cEnd]:
            pt = j[0], j[1]
            data.append(pt)
##            arcpy.AddMessage(pt)
        #insert 2 additional vertices at the midpoints of 3-point curves:
        if(len(data) == 3):
##            arcpy.AddMessage("data before 3-point insertion")
##            arcpy.AddMessage(str(data))
            data = midpointinsertion(data)
##            arcpy.AddMessage("data after 3-point insertion")
##            arcpy.AddMessage(str(data))
        # Fit an ideal curve to the list of coordinates using the selected method(s) (default "Hyperaccuate method" aka "Hyper")
        method = 'Hyper'
        # The 'issolved' variable is a flag to indicate if the curve fitting algorithm was successful or not
        issolved = False
        try:
            radius, center = fit_hypersphere(data, method)
##                arcpy.AddMessage("Radius: " + str(radius) + " Center: " + str(center))
            issolved = True

        # if an exception is raised by the program, this indicates the observed curve points may be ideal as digitized
        # if this is the case, add a small measure of random noise to the coordinate points using the 'randomize_vertices' function
        except Exception as e:
            try:
                arcpy.AddMessage(e.args[0])
                datarand = randomize_vertices(data,rscale)
##                arcpy.AddMessage('using randomized values')
                radius, center = fit_hypersphere(datarand, method)
##                arcpy.AddMessage("Radius: " + str(radius) + " Center: " + str(center))
                issolved = True

           # if the randomize_vertices function is executed, and there is still a curve fitting problem, raise an exception
            except Exception as e:
                issolved = False
##            arcpy.AddMessage("Curve fit successful: " + str(issolved))

        # Calculate R-squared
        r2 = r_squared(radius, center, data)
##        arcpy.AddMessage('completed r_squared')

##        cname = str(cStart) + "_" + str(cEnd-1) + "_" + str(c.featID)
        # Update curve fitting parameters from the fit_hypersphere function into Curve class
        # {1: [issolved1, radius1, center1, R-squared1, CurveName1], 2: [issolved2, radius2, center2, R-squared2, CurveName2], ..., N:  [issolvedN, radiusN, centerN, R-squaredN, CurveNameN]}
##        params = [issolved, radius, center.tolist(), r2, cname]
        params = [issolved, radius, center.tolist(), r2, c.curve_names[n]]
        c.parameters[n] = params
##        arcpy.AddMessage(c.parameters)
        # Move to next curve
        n += 1

    # return to main function
    return

#################################################################################################################################################################
"""
        FitHypersphere.py

        fit_hypersphere(collection of tuples or lists of real numbers)
        will return a hypersphere of the same dimension as the tuples:
                (radius, (center))

        using the Hyper (hyperaccurate) algorithm of
        Ali Al-Sharadqah and Nikolai Chernov
        Error analysis for circle fitting algorithms
        Electronic Journal of Statistics
        Vol. 3 (2009) 886-911
        DOI: 10.1214/09-EJS419

        generalized to n dimensions

        Mon Apr 23 04:08:05 PDT 2012 Kevin Karplus

        Note: this version using SVD works with Hyper, Pratt, and Taubin methods.
        If you are not familiar with them, Hyper is probably your best choice.


        Creative Commons Attribution-ShareAlike 3.0 Unported License.
        http://creativecommons.org/licenses/by-sa/3.0/
"""
def fit_hypersphere(data, method="Hyper"):
    """returns a hypersphere of the same dimension as the
        collection of input tuples
                (radius, (center))

       Methods available for fitting are "algebraic" fitting methods
        Hyper   Al-Sharadqah and Chernov's Hyperfit algorithm
        Pratt   Vaughn Pratt's algorithm
        Taubin  G. Taubin's algorithm

       The following methods, though very similar, are not implemented yet,
          because the contraint matrix N would be singular,
          and so the N_inv computation is not doable.

        Kasa    Kasa's algorithm
    """
##    arcpy.AddMessage("funcs.FitHypersphere")
    num_points = len(data)
##    arcpy.AddMessage("DEBUG: num_points=" + str(num_points))

    if num_points==0:
        return (0,None)
    if num_points==1:
        return (0,data[0])
    dimen = len(data[0])        # dimensionality of hypersphere
##    arcpy.AddMessage(" DEBUG: dimen=" + str(dimen))

    if num_points<dimen+1:
        raise ValueError(\
            "Error: fit_hypersphere needs at least {} points to fit {}-dimensional sphere, but only given {}".format(dimen+1,dimen,num_points))

    # central dimen columns of matrix  (data - centroid)
    central = np.matrix(data, dtype=float)      # copy the data
    centroid = np.mean(central, axis=0)
    for row in central:
        row -= centroid
##    arcpy.AddMessage(" DEBUG: central=" + repr(central))

    # squared magnitude for each centered point, as a column vector
    square_mag= [sum(a*a for a in row.flat) for row in central]
    square_mag = np.matrix(square_mag).transpose()
##    arcpy.AddMessage(" DEBUG: square_mag=" + str(square_mag))

    if method=="Taubin":
        # matrix of normalized squared magnitudes, data
        mean_square = square_mag.mean()
        data_Z = np.bmat( [[(square_mag-mean_square)/(2*sqrt(mean_square)), central]])
##        arcpy.AddMessage(" DEBUG: data_Z=" + str(data_Z))
        u,s,v = linalg.svd(data_Z, full_matrices=False)
        param_vect= v[-1,:]
        params = [ x for x in np.asarray(param_vect)[0]]        # convert from (dimen+1) x 1 matrix to list
        params[0] /= 2*sqrt(mean_square)
        params.append(-mean_square*params[0])
        params=np.array(params)

    else:
        # matrix of squared magnitudes, data, 1s
        data_Z = np.bmat( [[square_mag, central, np.ones((num_points,1))]])
##        arcpy.AddMessage(" DEBUG: data_Z=" + str(data_Z))

        # SVD of data_Z
        # Note: numpy's linalg.svd returns data_Z = u * s * v
        #         not u*s*v.H as the Release 1.4.1 documentation claims.
        #         Newer documentation is correct.
        u,s,v = linalg.svd(data_Z, full_matrices=False)
##        arcpy.AddMessage(" DEBUG: u=" + repr(u))
##        arcpy.AddMessage(" DEBUG: s=" + repr(s))
##        arcpy.AddMessage(" DEBUG: v=" + repr(v))
##        arcpy.AddMessage(" DEBUG: v.I=" + repr(v.I))

        if s[-1]/s[0] < 1e-12:
            # singular case
            # param_vect as (dimen+2) x 1 matrix
            param_vect = v[-1,:]
            # Note: I get last ROW of v, while Chernov claims last COLUMN,
            # because of difference in definition of SVD for MATLAB and numpy
##            arcpy.AddMessage(" DEBUG: singular, param_vect=" + repr(param_vect))
##            arcpy.AddMessage(" DEBUG: data_Z*V=" + repr(data_Z*v))
##            arcpy.AddMessage(" DEBUG: data_Z*VI=" + repr(data_Z*v.I))
##            arcpy.AddMessage(" DEBUG: data_Z*A=" + repr(data_Z*v[:,-1]))
        else:
            Y = v.H*np.diag(s)*v
            Y_inv = v.H*np.diag([1./x for x in s])*v
##            arcpy.AddMessage(" DEBUG: Y=" + repr(Y))
##            arcpy.AddMessage(" DEBUG: Y.I=" + repr(Y.I), "\nY_inv=",repr(Y_inv))
            #Ninv is the inverse of the constraint matrix, after centroid has been removed
            Ninv = np.asmatrix(np.identity(dimen+2, dtype=float))
            if method=="Hyper":
                Ninv[0,0] = 0
                Ninv[0,-1]=0.5
                Ninv[-1,0]=0.5
                Ninv[-1,-1] = -2*square_mag.mean()
            elif method=="Pratt":
                Ninv[0,0] = 0
                Ninv[0,-1]=-0.5
                Ninv[-1,0]=-0.5
                Ninv[-1,-1]=0
            else:
                raise ValueError("Error: unknown method: {} should be 'Hyper', 'Pratt', or 'Taubin'")
##            arcpy.AddMessage(" DEBUG: Ninv=" + repr(Ninv))

            # get the eigenvector for the smallest positive eigenvalue
            matrix_for_eigen = Y*Ninv*Y
##            arcpy.AddMessage(" DEBUG: {} matrix_for_eigen=\n{}".format(method + repr(matrix_for_eigen)))
            eigen_vals,eigen_vects = linalg.eigh(matrix_for_eigen)
##            arcpy.AddMessage(" DEBUG: eigen_vals=" + repr(eigen_vals))
##            arcpy.AddMessage(" DEBUG: eigen_vects=" + repr(eigen_vects))
            positives = [x for x in eigen_vals if x>0]
            if len(positives)+1 != len(eigen_vals):
                # raise ValueError("Error: for method {} exactly one eigenvalue should be negative: {}".format(method,eigen_vals))
                arcpy.AddMessage(" Warning: for method {} exactly one eigenvalue should be negative: {}".format(method,eigen_vals))
            smallest_positive = min(positives)

            # chosen eigenvector as 1 x (dimen+2) matrix
            A_colvect =eigen_vects[:,list(eigen_vals).index(smallest_positive)]
##            arcpy.AddMessage(" DEBUG: A_colvect=" + repr(A_colvect))

            # now have to multiply by Y inverse
            param_vect = (Y_inv*A_colvect).transpose()
##            arcpy.AddMessage(" DEBUG: nonsingular, param_vect=" + repr(param_vect))
            params = np.asarray(param_vect)[0]  # convert from (dimen+2) x 1 matrix to array of (dimen+2)
##            arcpy.AddMessage(" DEBUG: params=" + repr(params))

    radius = 0.5* sqrt( sum(a*a for a in params[1:-1])- 4*params[0]*params[-1])/abs(params[0])
##    arcpy.AddMessage(radius)
    center = -0.5*params[1:-1]/params[0]
    center += np.asarray(centroid)[0]
##    arcpy.AddMessage(center)

    # return the radius as a float, and the center as a numpy array
    return (radius,center)

#################################################################################################################################################################
def polar_angle(startY, startX, endY, endX):

##    arcpy.AddMessage("funcs.polar_angle")
    theangle = math.atan2( endY - startY, endX - startX)
    if(theangle < 0): theangle = 2*pi + theangle
##    print(i+1, endXs[i], endYs[i], theangle, theangle/rads, names[i])
    return theangle
#################################################################################################################################################################
def generate_ideal_curve_with_rsquared(radius, c_x, c_y, startX, startY, endX, endY, deltaD, thedirection, obs_polyline, max_distance):
    """ Generate the points for an ideal curve with radius and center calculated by the hyperfit algorithm  """
##    arcpy.AddMessage("funcs.generate_ideal_curve")

    radsperdegree = math.pi/180.             #0.0174532925199
    degperrad = 180./math.pi                 # 57.29577951308232

    #Calculate the angle using the center point, and the start and end of the observed curve's start and end points:
    fullcircle = 2*math.pi                   #6.28318530718
    fullcircledeg = 360.

    deltarad = deltaD * radsperdegree

    angle_s = polar_angle(c_y, c_x, startY, startX)
    angle_e = polar_angle(c_y, c_x, endY, endX)

##    arcpy.AddMessage("angle_e before adjustment "+str(degperrad*angle_e))
    if(angle_e < angle_s and (thedirection == 1)):
            angle_e += fullcircle
    if(angle_e >= angle_s and (thedirection == -1)):
            angle_e -= fullcircle

    angle_range = abs(angle_e - angle_s)

    dec_nvertex = angle_range/deltarad
    nvertices = int(dec_nvertex)
    rmder = dec_nvertex - nvertices
    rmder_inrads = deltarad * rmder

##    arcpy.AddMessage("deltarad, deltaD, dec_nvertex, nvertices, rmder=, rmder as rads, rmder as degrees ")
##    arcpy.AddMessage(str(deltarad)+"," +str(deltaD)+"," +str(dec_nvertex)+"," +str(nvertices)+"," +str(rmder)+ "," + str(rmder_inrads)+ ", "+str(rmder_inrads*degperrad))
##    arcpy.AddMessage("deltaD, Radius: ")
##    arcpy.AddMessage(str(deltaD)+", "+str(radius))
##    arcpy.AddMessage("center X,Y: ")
##    arcpy.AddMessage(str(c_x)+", "+str(c_y))
##    arcpy.AddMessage('startX, startY, endX, endY = ')
##    arcpy.AddMessage(str(startX)+" "+str(startY) +" "+str(endX)+" "+str(endY))
##    arcpy.AddMessage("angle_s(rad), angle_s(deg), angle_e(rad), angle_e(deg), angle_range(deg)  ")
##    arcpy.AddMessage(str(angle_s)+", "+str(angle_s*degperrad)+", " + str(angle_e)+", "+ str(angle_e*degperrad)+", "+str(angle_range*degperrad))

    angles_list =[]
    p_angle = angle_s
    angles_list.append(p_angle)

    for i in xrange(0, nvertices):
##        arcpy.AddMessage("in angle calc "+str(i))
        p_angle += (deltarad * thedirection)
        angles_list.append(p_angle)
    p_angle += (rmder_inrads * thedirection)
    angles_list.append(p_angle)

##    arcpy.AddMessage("i, angles_list[i]")
##    for i in xrange(len(angles_list)):
##        arcpy.AddMessage(str(i)+", "+str(angles_list[i]))
    x = []
    y = []
    isect_point_list = []
##    arcpy.AddMessage("x_ctob, y_ctob")
    for i in xrange(len(angles_list)):
        # x and y point lists for center to best fit for intersection with observered curve, where
        # each x and y list has the x and y values from the center to the x and y values to best-fit, respectively (Center TO Bestfit).
        # these are intersected with the observed curve polyline (input parameter to function) to obtain a point for r-squared:
        xb = math.cos(angles_list[i]) * radius + c_x
        yb = math.sin(angles_list[i]) * radius + c_y
        x.append(xb)
        y.append(yb)
        #Create a line from the center to the current best-fit vertex, for intersection with polyline.
        x_ctob = math.cos(angles_list[i]) * max_distance + c_x
        y_ctob = math.sin(angles_list[i]) * max_distance + c_y
        # Put the arcpy.Point objects into an array, and then convert into a polyline
        ctob_array = arcpy.Array([arcpy.Point(c_x,c_y),arcpy.Point(x_ctob,y_ctob)])
        ctob_polyline = arcpy.Polyline(ctob_array)
        #Intersect the ctob_polyline (center of bestfit to max_distance) with the observed curve (obs_polyline):
        isect_points = ctob_polyline.intersect(obs_polyline, 1)
        isect_point = isect_points.getPart(0)
##        ismult = isect_point.isMultipart
##        arcpy.AddMessage(str(ismult))
        isect_point_list.append([isect_point.X, isect_point.Y])
##        arcpy.AddMessage(str(isect_point.X)+ ", "+str(isect_point.Y))


##        arcpy.AddMessage(str(i)+", "+str(angles_list[i])+", "+str(math.cos(angles_list[i])*radius)+", "+str(math.cos(angles_list[i]) * radius + center[0])+", "+str(math.sin(angles_list[i])*radius)+", "+str(math.sin(angles_list[i]) * radius + center[1]))

    return isect_point_list, (zip(x,y))

#################################################################################################################################################################
def r_squared(radius, center, data):
    """ Calculate the r-squared values between the observed and ideal curves

    radius is the radius of the ideal circle calculated by fit_hypersphere
    center is the center of the ideal circle calculated by fit_hypersphere
    data is a list of [x, y] points in list of list format  """

##    arcpy.AddMessage('funcs.r-squared')
    # Compute the distance between the center and each observed point
    dsum2 = 0
    rsum2 = 0
    npoints = len(data)
    for pt in data:
        # Calculate numerator
        ri, rix, riy = calculate_distance(center[0], center[1], pt[0], pt[1])
        di = (ri - radius) * (ri - radius)
        dsum2 += di
        # Calculate denominator
##        ri2 = (pt[0] - center[0]) * (pt[0] - center[0]) + (pt[1] - center[1]) * (pt[1] - center[1])
        mrad = max(ri, radius)
        ri2 = mrad**2
        rsum2 += ri2

    # Calculate correlation coefficient
    r2 = 1 - (dsum2 / rsum2)

    # return r-squared as a float
    return r2

################################################################################
def calculate_angle_geometry(pointList):
    """ Calculate the bearing between two line segments based on plane geometry

    pointList the r.vertices member of the current Line instance """

##    arcpy.AddMessage("funcs.calculate_geometric_bearing")
    line_geometry = []
    distance = []
    rads = math.pi/180.
    cnv = 180./math.pi
    if(VERBOSE):
        arcpy.AddMessage("length of array (from 0): " + str(len(pointList)))
        arcpy.AddMessage("Index, Xa, Xb, Xc, Ya, Yb, Yc:  ")

    for i in xrange(len(pointList) - 2):
        # X-coordinates
        Xa = pointList[i][0]
        Xb = pointList[i+1][0]
        Xc = pointList[i+2][0]

        # Y-coordinates
        Ya = pointList[i][1]
        Yb = pointList[i+1][1]
        Yc = pointList[i+2][1]
        if(VERBOSE): arcpy.AddMessage(str(i)+", "+str(Xa)+", "+str(Xb)+", "+str(Xc)+", "+str(Ya)+", "+str(Yb)+", "+str(Yc))

        # Calculate the bearing between the two points
        pol_angle_a = cnv*polar_angle(Ya, Xa, Yb, Xb)
        pol_angle_b = cnv*polar_angle(Yb, Xb, Yc, Xc)
        #Adjust angles if the comparison angles are in quadrant 4 and 1:
##        if((pol_angle_b-pol_angle_a) < -180. ):
##            pol_angle_b = pol_angle_b + 360.
##        elif((pol_angle_b-pol_angle_a) >= 180.):
##            pol_angle_a = pol_angle_a + 360.
        if((pol_angle_b-pol_angle_a) < -90. ):
            pol_angle_b = pol_angle_b + 360.
        elif((pol_angle_b-pol_angle_a) >= 90.):
            pol_angle_a = pol_angle_a + 360.

        angle_change = (pol_angle_b - pol_angle_a)
##        if(angle_change < 0.):
##            direction = 1.
##        elif(angle_change >= 0.):
##            direction = -1.
        if(angle_change < 0.):
            direction = -1.
        elif(angle_change >= 0.):
            direction = 1.

        d, dx, dy = calculate_distance(Xa, Ya, Xb, Yb)

        term = ((Xb - Xa)*(Xc - Xb) + (Yb - Ya)*(Yc - Yb))/(math.sqrt((Xb - Xa)**2 + (Yb - Ya)**2) * math.sqrt((Xc - Xb)**2 + (Yc - Yb)**2))
        # "term" should never be less than -1 or greater than 1. However, due to numerical precision problelms, occasionally this happens. The following will prevent this:
        if(term > 1.0):
            dotproduct = 1.0
        elif(term < -1.0):
            dotproduct = -1.0
        else:
            dotproduct = term

        theta = cnv * math.acos ( dotproduct )
        line_geometry.append([i, theta, angle_change, direction, pol_angle_a, pol_angle_b, d, dx, dy])

    if(VERBOSE):
        arcpy.AddMessage("i, theta, angle_change, direction, pol_angle_a, pol_angle_b, d, dx, dy")
        for i, item in enumerate(line_geometry):
            # i, theta, angle_change, direction, pol_angle_a, pol_angle_b, d, dx, dy
            arcpy.AddMessage(str(item[0])+", "+ str(item[1])+", "+ str(item[2])+", "+ str(item[3])+", "+ str(item[4])+", "+ str(item[5])+", "+ str(item[6])+", "+ str(item[7])+", "+ str(item[8]))

    # return thetaList (a list of bearings), and distance (a list of distances between points)
    return line_geometry

################################################################################
def midpointinsertion(ptlist):
    data = []
    if (len(ptlist) > 3):
        data = ptlist
        return ptlist
    else:
        x2 = ptlist[0][0] + (ptlist[1][0] - ptlist[0][0])/2.
        y2 = ptlist[0][1] + (ptlist[1][1] - ptlist[0][1])/2.
        x4 = ptlist[1][0] + (ptlist[2][0] - ptlist[1][0])/2.
        y4 = ptlist[1][1] + (ptlist[2][1] - ptlist[1][1])/2.
        data.append((ptlist[0][0], ptlist[0][1]))
        data.append((x2, y2))
        data.append((ptlist[1][0], ptlist[1][1]))
        data.append((x4, y4))
        data.append((ptlist[2][0], ptlist[2][1]))
    return data
