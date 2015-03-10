# Copyright:   (c) URS Corporation, 2014
# ArcGIS Version:   10.2
# Python Version:   2.7
#--------------------------------
import os
import sys
import arcpy
from arcpy import env
import KDOT_HorizontalCurveTool_Functions as funcs
# http://blogs.esri.com/esri/arcgis/2011/08/04/pythontemplate/
def do_analysis(*argv):
    """This function starts the analyses
    argv[0] = PolylineM roads feature class for curve detection
    argv[1] = workspace for output feature classes
    argv[2] = list of fields to concatenate for a unique identifier
    argv[3] = binary value to specify if topology check should be done ( 1 = yes, 0 = no)
    argv[4] = simplify tolerance
    argv[5] = threshold parameter for change in bearing to either start or stop a curve
    argv[6] = maximum distance parameter at which to truncate a curve
    argv[7] = Output feature class for lines feature class deliverables
    argv[8] = Output feature class for points feature class deliverables (disabled for now)
    argv[8] = Output feature class for hyperfit curves
    argv[9] = Angle resolution (in degrees) for the output feature class for the hyperfit curves """
    # Overwrite any existing feature classes of the same name
    arcpy.env.overwriteOutput = True
    try:
        #rscale defines the scale at which the decimal fraction of the vertices are randomized if the hyperfit fails. Hyperfit may fail if all the points are EXACTLY on a perfect curve.
        #A perfect curve may exist from digitization techniques. The randomization will only be done if the hyperfit fails.
        rscale = 100.

        mxd = arcpy.mapping.MapDocument("CURRENT")
        dataFrame = arcpy.mapping.ListDataFrames(mxd, "*")[0]

        VERBOSE = funcs.detailed_errors(False)
##        VERBOSE = funcs.detailed_errors(True)
##        arcpy.AddMessage("VERBOSE: "+str(VERBOSE))

#################################################################################################################################################################
        #THIS SECTION IDENTIFIES THE INPUT FEATURES, INPUT PARAMETERS, AND OUTPUT FEATURES, AND ALSO PREPARES THE OUTPUT FEATURES FOR WRITING:
        #
        #Obtain the inputroads features
        arcpy.AddMessage(" ")
        arcpy.AddMessage("Input parameters: ")
        inputroads = argv[0]

        # Check for indication of group layer
        arcpy.AddMessage("Input features: "+inputroads)
        if inputroads.find("\\"):
            arcpy.AddMessage('Input Feature Class in group layer')
            inputroads = inputroads.split("\\")[-1]
            arcpy.AddMessage("Input features: "+inputroads)

##        lyr_inputroads = "CRND Input Data"
##        #Add the output observed curves and best-fit curves to the dataframe:
##        arcpy.MakeFeatureLayer_management(inputroads, lyr_inputroads)
##        roadslayer = arcpy.mapping.Layer(lyr_inputroads)
##        arcpy.mapping.AddLayer(dataFrame, roadslayer, "AUTO_ARRANGE")

        #Obtain the spatial reference of the input roads:
        dsc = arcpy.Describe(inputroads)
##        dsc = arcpy.Describe(roadslayer)
        coord_sys = dsc.spatialReference

        #Obtain the workspace name
        theworkspace = argv[1]
        arcpy.AddMessage("Workspace: "+theworkspace)
        env.workspace = theworkspace

##        #The check to see if in an edit session is not yet functioning. Removed for now.
##        anedit = arcpy.da.Editor(theworkspace)
##        arcpy.AddMessage(str(anedit.isEditing))
##        if(str(anedit.isEditing) == "True"):
##            raise ineditsessionerror

        #Obtain the field for a unique identifier:
        featurename_field = argv[2]
        arcpy.AddMessage("Road name field identifier: "+featurename_field)

        #Deteremine if user wants to perform the topology check on the data:
        dotopologycheck = argv[3]
        arcpy.AddMessage('dotopologycheck '+str(dotopologycheck))

        #Perform simplify line? Now a required parameter.
##        dosimplify = argv[4]
        dosimplify = True
        #Simplify method: "point remove" only option. Method choice has been removed from parameter choices.
##        dosimplify_method = argv[5]
        dosimplify_method = "POINT_REMOVE"
        #Simplify tolerance (linear unit)
        dosimplify_tolerance = float(argv[4].split(" ")[0])
        arcpy.AddMessage("dosimplify_tolerance: "+str(dosimplify_tolerance))
        #Set the simplification tolerance to be at least 0.001 meters or 0.003280839 feet. This is the default xy tolerance in ArcGIS.
        if(dosimplify_tolerance < 0.003280839): dosimplify_tolerance = 0.003280839
##        arcpy.AddMessage("dosimplify_tolerance after test: "+str(dosimplify_tolerance))

        # Obtain the change in bearing threshold parameter
        threshold = round(float(argv[5]), 2)
        arcpy.AddMessage('Input threshold value: '+ str(threshold))

        # Obtain the maxDistance parameter. Defaults specified in toolbox parameters is 1040.03 # 1040.03 feet = 317 meters
        maxDistance = round(float(argv[6].split(" ")[0]),2)
        arcpy.AddMessage('Input maxDistance value: '+ str(maxDistance))

        #Obtain output features locations and names:
        #Output feature class lines for observed curves (spatial subsets of input roads)
        outputline_features = argv[7]
        arcpy.AddMessage('output line features of observed curves: '+outputline_features)
        #Output feature class points for observed curves (spatial subsets of input roads)
##        outputpoint_features = argv[8]
##        arcpy.AddMessage('output point features of observed curves: '+outputpoint_features)
        #Output feature class points for observed curves (spatial subsets of input roads)
        outputfit_features = argv[8]
        arcpy.AddMessage('output line features of fitted curves: '+outputfit_features)
        deltaD = float(argv[9])
        arcpy.AddMessage('angle resolution of fitted curves: '+str(deltaD))
        arcpy.AddMessage(" ")

#################################################################################################################################################################
        #Processing starts here:
        arcpy.AddMessage("*** BEGINNING PROCESSING ***")
        #Check the topology and return a list of bad id's and the location of a copy of the input roads. Also process the
        #roads with the simplify_line tool with the point_remove option at a tolerance of 0.001 meters so that redundant vertices on staight lines are removed.
        #If the user specifies their own parameters for simplify_line, THAT ARE NOT POINT_REMOVE AND THE TOLERANCE IS > 0.001 METERS, that is done additionally:
        orig_id_field = 'OriginalOID'
##        badfids, topology_featureclass = funcs.make_workspace_copy(roadslayer,theworkspace,dotopologycheck,dosimplify,dosimplify_method,dosimplify_tolerance, orig_id_field)
        badfids, topology_featureclass = funcs.make_workspace_copy(inputroads,theworkspace,dotopologycheck,dosimplify,dosimplify_method,dosimplify_tolerance, orig_id_field)

        outlinepath, outnameline = os.path.split(outputline_features)
        outpathfit, outnamefit = os.path.split(outputfit_features)
        outnameline_in_fdset = False
        outputfit_in_fdset = False

        addedfields = ["Original_OID", "Unique_ID", "CurveID"]
##        addedfields_ideal = ["Original_OID", "Unique_ID", "CurveID", "R_Squared"]
        addedfields_ideal = ["Original_OID", "Unique_ID", "CurveID", "R_Squared", "R_Squared_Dense", "X_Center", "Y_Center", "Radius", "X_Start", "Y_Start", "X_End", "Y_End"]

        # Lines: create a new line feature class of curves from the input vertices
        arcpy.CreateFeatureclass_management(outlinepath, outnameline, "POLYLINE")

        #Add new fields to the output feature line class:
        arcpy.AddField_management(outputline_features, addedfields[0], "LONG", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputline_features, addedfields[1], "TEXT", "", "", "250", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputline_features, addedfields[2], "TEXT", "", "", "250", "", "NULLABLE", "NON_REQUIRED", "")

        # Fitted line: create new feature class for ideal curve line
        arcpy.CreateFeatureclass_management(outpathfit, outnamefit, "POLYLINE")

##        arcpy.AddMessage(outlinepath)
##        arcpy.AddMessage(outpathfit)

        count_outlinepath = outlinepath.split(".gdb")[1].count("\\")
        count_outpathfit = outpathfit.split(".gdb")[1].count("\\")
        if(count_outlinepath == 1):
            outnameline_in_fdset = True
        else:
            outnameline_in_fdset = False

        if(count_outpathfit == 1):
            outputfit_in_fdset = True
        else:
            outputfit_in_fdset = False


##        arcpy.AddMessage(str(outnameline_in_fdset)+","+str(count_outlinepath))
##        arcpy.AddMessage(str(outputfit_in_fdset)+","+str(count_outpathfit))


        #Set the coordinate system ONLY if it is not in a FeatureDataset:
        if(outnameline_in_fdset == False): arcpy.DefineProjection_management(outnameline, coord_sys)
        if(outputfit_in_fdset == False): arcpy.DefineProjection_management(outnamefit, coord_sys)


        #Add new fields to the output feature line class:
        arcpy.AddField_management(outputfit_features, addedfields_ideal[0], "LONG", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[1], "TEXT", "", "", "250", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[2], "TEXT", "", "", "250", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[3], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[4], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[5], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[6], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[7], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[8], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[9], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[10], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(outputfit_features, addedfields_ideal[11], "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")

        #Prepare the insert cursor field types and names. The type for Fid and Shape must be "OID@" and "SHAPE@", not the name....
        cursor_line = arcpy.da.InsertCursor(outputline_features, ['SHAPE@'] + addedfields)
        cursor_fit = arcpy.da.InsertCursor(outputfit_features, ['SHAPE@'] + addedfields_ideal)

        ########################################################################
        #THIS SECTION EXTRACTS VERTICES FOR ALL INPUT ROADS, PUTS THEM AN LINE CLASS INSTANCES,  THEN IDENTIFIES THE CURVES WITHIN EACH, AND WRITES THEM OUT:
        # Break line features into xyzm list. This function will not include any roads that do not conform to the topology check, if topology check is selected:
        roadInstances = funcs.break_linefeatures_into_Line_class_with_maxdistance(topology_featureclass, maxDistance, featurename_field, orig_id_field, badfids)
##        arcpy.AddMessage('Completed roadInstances')
        curveInstances = []
        arcpy.AddMessage('CURVE DETECTION #########################################################')
        for road in roadInstances:
            if(VERBOSE):
                arcpy.AddMessage("Output of function calculate_angle_geometry")
                arcpy.AddMessage("Original FID, name: "+str(road.featID) +", "+str(road.name))
            # geoList = list of lists containing: [i, theta, angle_change, direction, pol_angle_a, pol_angle_b, d, dx, dy]
            geoList = funcs.calculate_angle_geometry(road.vertices)

            # Detect Curves
            # returns the curve starting x- and y- coordinates, the curve ending x- and y- coordinates, and the index of the curve start point and index of curve end point
##            arcpy.AddMessage(road.name)
            Xstart, Ystart, Xend, Yend, iStart, iEnd = funcs.detect_curves(road, geoList, threshold)
##            arcpy.AddMessage("Xstart, Ystart, Xend, Yend, iStart, iEnd: "+str(Xstart)+""+str(Ystart)+""+str(Xend)+""+str(Yend)+""+str(iStart)+""+str(iEnd))
            arcpy.AddMessage('completed curve detection for road name '+road.name)

            # Initialize curves into curve class
            # if there are curves on the line segment, then instantiate the Curve class
            # there is only one Curve instance for each line segment
            # the Curve instance contains a dictionary of the 1 or more curves on the line
            #
            # funcs.populate_curves and returns a Curve instance and PARTIALLY POPULATATES the curve instance. The parameters attribute of the class is populated
            # later within function hypersphere.
            curveInstance = funcs.populate_curves(road.featID, road.name, road.vertices, iStart, iEnd, geoList)
            curveInstances.append(curveInstance)
        arcpy.AddMessage('completed curve recognition......')

        #############################################################################################
        #Simultaneous cursors require an edit session:
        edit = arcpy.da.Editor(theworkspace)
        edit.startEditing(False, False)
        # Write the observed curves to a new points and lines feature classes
        # Compile points within all curves of this Line instance into a list for writing into a point and line feature class
        curve_points = []
        for crv in curveInstances:
            for key, value in crv.curves.iteritems():
                s_ref = value[0]
                e_ref = value[1]
##                arcpy.AddMessage("s_ref, e_ref "+str(s_ref)+", "+str(e_ref))
                cname = crv.curve_names[key]
                curve_points = crv.vertices[s_ref:e_ref+1]
                # Reformat point list of lists into separate x-list and y-list
                pts = funcs.convert_pointList_to_x_and_y_lists(curve_points ,2)
                # Create a new list of arcpy.Point objects from pts list
                a_pts = []
                [a_pts.append(arcpy.Point(k[0],k[1])) for k in pts]
##                for thepoint in a_pts:
##                    arcpy.AddMessage("x, y for observed curve in separte loop: " +str(a_pts))
                # Put the arcpy.Point objects into an array, and then convert into a polyline
                array = arcpy.Array(a_pts)
                polyline = arcpy.Polyline(array)
                # insert the new polyline geometery into the lines feature class
                cursor_line.insertRow([polyline, crv.featID, crv.name, cname])
        del cursor_line
        edit.stopEditing(True)
        arcpy.AddMessage('completed writing observed curves to the output feature class')

        #############################################################################################

        arcpy.AddMessage('HYPERFIT PROCESSING #####################################################')
        #THIS SECTION FITS AN IDEAL HORIZONTAL CURVE TO EACH DETECTED CURVE in the Curve Class:
        edit = arcpy.da.Editor(theworkspace)
        edit.startEditing(False, False)
        curve_points = []
        for roc in curveInstances:
            # Perform the "Hyperfit" algorithm for each curve.
            # THIS IS THE "CORE" OF THE PROGRAM WHERE WE CALCULATE THE BEST-FIT FOR EACH OBSERVED CURVE.
            # This function provides the interface to the hyperfit algorithm.
            # This function also inserts the parameters into the curve instances. View function hypersphere
            funcs.hypersphere(roc,rscale)
##            arcpy.AddMessage('completed funcs.hypersphere')

            #deltaD is the fractions of a degree with which is written the ideal fitted curves...
##            for i in range(len(roc.curves)):
            for key, value in roc.curves.iteritems():
                issolved = roc.parameters[key][0]
                radius = roc.parameters[key][1]
                center = roc.parameters[key][2]

                if(str(issolved) == "True"):
                    rfeatID = roc.featID
                    rname = roc.name
                    r_rsquared = roc.parameters[key][3]
                    cname = roc.parameters[key][4]
                    startpt = value[0]
                    endpt = value[1]
                    Xs = roc.vertices[startpt][0]
                    Ys = roc.vertices[startpt][1]
                    Xe = roc.vertices[endpt][0]
                    Ye = roc.vertices[endpt][1]

                    #Create the curve polyline to calculate R-squared from the perspective of the best-fit line:
                    curve_points = roc.vertices[startpt:endpt+1]
                    # Reformat point list of lists into separate x-list and y-list
                    pts = funcs.convert_pointList_to_x_and_y_lists(curve_points ,2)
                    # Create a new list of arcpy.Point objects from pts list
                    a_pts = []
                    [a_pts.append(arcpy.Point(k[0],k[1])) for k in pts]
                    # Put the arcpy.Point objects into an array, and then convert into a polyline
                    array = arcpy.Array(a_pts)
                    polyline = arcpy.Polyline(array)
                    #Calculate the maximum distance from the center to each point on the observed-curve. This is so that a maximum distance
                    #from the center to either an observed-curve or best-fit vertex will be used to obtain the intersection of a line from the center
                    #to the best-fit curve with the observed curve:
                    c_x = center[0]
                    c_y = center[1]
                    maxdist = []
                    for k in pts:
                        therad, dx, dy = funcs.calculate_distance(c_x, c_y, k[0], k[1])
##                        arcpy.AddMessage("c_x, c_y, k[0], k[1]: " + str(c_x) + ", " + str(c_y) + ", " + str(k[0]) + ", " + str(k[1]))
                        maxdist.append(therad)
                    themaxobs = max(maxdist)
                    usedradius = max(radius, themaxobs)
##                    arcpy.AddMessage("themaxobs, userradius = "+str(themaxobs)+", "+str(usedradius))

                    if r_rsquared == None:
                        r_rsquared = 0

##                    fitpoints =  funcs.generate_ideal_curve(radius, center, Xs, Ys, Xe, Ye, deltaD, roc.curve_directions[key])
                    isect_point_list, fitpoints =  funcs.generate_ideal_curve_with_rsquared(radius, c_x, c_y, Xs, Ys, Xe, Ye, deltaD, roc.curve_directions[key], polyline, usedradius)
                    ####################################################################
                    # Calculate R-squared using the vertices from the desified vertices on the observed using the radial location from the best-fit curve:
                    r2_dense = funcs.r_squared(radius, [c_x, c_y], isect_point_list)
##                    arcpy.AddMessage("curve, r_squared, r_squared_dense: "+cname+", "+str(r_rsquared)+", "+str(r2_dense))
                    # Write the fitted curves to a new lines feature classes
    ##                #IStart = a list of indices of the starting points of each curve. In otherwords, it's length is the number of curves...
    ##                arcpy.AddMessage(roc.featID[i])
                    if len(fitpoints) > 0:
                        pts = funcs.convert_pointList_to_x_and_y_lists(fitpoints ,2)
                         # Create a new list of arcpy.Point objects from pts list
                        a_pts = []
                        [a_pts.append(arcpy.Point(k[0],k[1])) for k in pts]
##                            arcpy.AddMessage("Point: " + str(k[0]) +", " + str(k[1]))
                        # Put the arcpy.Point objects into an array, and then convert into a polyline
                        array = arcpy.Array(a_pts)
                        polyline = arcpy.Polyline(array)
                        # insert the new polyline geometery into the lines feature class
                        #"Original_OID", "Unique_ID", "CurveID", "R_Squared", "R_Squared_Dense, "X_Center", "Y_Center", "Radius", "X_Start", "Y_Start", "X_End", "Y_End"
##                        cursor_fit.insertRow([polyline, rfeatID, rname, cname, r_rsquared])
                        cursor_fit.insertRow([polyline, rfeatID, rname, cname, r_rsquared, r2_dense, center[0], center[1], radius, Xs, Ys, Xe, Ye])
                    else: arcpy.AddMessage("Warning: no fitted points for: " + str(cname) + " on road " + str(roc.name))
                arcpy.AddMessage('completed hyperfit processing for road ' + str(roc.name)+ ", curve "+cname+", direction "+str(roc.curve_directions[key]))
        arcpy.AddMessage('completed hyperfit analysis......')
        arcpy.AddMessage('completed writing best-fit curves to the output feature class')
        arcpy.AddMessage('#########################################################')

        del cursor_fit
        edit.stopEditing(True)

        lyr_name_observed = str(outnameline)
        lyr_name_fit = str(outnamefit)
        #Add the output observed curves and best-fit curves to the dataframe:
        arcpy.MakeFeatureLayer_management(outputline_features, lyr_name_observed)
        layer = arcpy.mapping.Layer(lyr_name_observed)
        arcpy.mapping.AddLayer(dataFrame, layer, "AUTO_ARRANGE")
        arcpy.MakeFeatureLayer_management(outputfit_features, lyr_name_fit)
        layer2 = arcpy.mapping.Layer(lyr_name_fit)
        arcpy.mapping.AddLayer(dataFrame, layer2, "AUTO_ARRANGE")

    except arcpy.ExecuteError:
        arcpy.GetMessages(2)
    except Exception as e:
        arcpy.AddMessage(e.args[0])
##    except ineditsessionerror:
##        arcpy.AddMessage("You may not be in an editing session while this tool is active")
##        arcpy.AddMessage("Please discontinue the edit session before operating this tool")
# End do_analysis function

# This test allows the script to be used from the operating
# system command prompt (stand-alone), in a Python IDE,
# as a geoprocessing script tool, or as a module imported in
# another script
if __name__ == '__main__':
    # Arguments are optional
    argv = tuple(arcpy.GetParameterAsText(i) for i in range(arcpy.GetArgumentCount()))

if arcpy.CheckProduct("ArcEditor") == "AlreadyInitialized" or arcpy.CheckProduct("ArcInfo") == "AlreadyInitialized":
    arcpy.AddMessage("Current License Level: " + arcpy.ProductInfo())
    do_analysis(*argv)
else:
    arcpy.AddMessage("Current License Level: " + arcpy.ProductInfo())
    arcpy.AddMessage("Error: ArcEditor (ArcGIS for Desktop Standard) or ArcInfo (ArcGIS for Desktop Advanced) Required")
    arcpy.AddMessage("Program Terminated")







