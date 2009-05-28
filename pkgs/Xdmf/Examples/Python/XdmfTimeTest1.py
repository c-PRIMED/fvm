#!/bin/env python

from Xdmf import *

TimeTestTxt = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [] >

<Xdmf>
<Domain Name="Domain1" >
        <Grid Name="Parent One" GridType="Collection">
            <Time TimeType="List">
                <DataItem Format="XML" Dimensions="4" NumberType="Float">
                    0.0 0.1 0.2 0.3
                </DataItem>
            </Time>
            <Grid Name="One">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0    0.001 0 0   0.001 0.001 0
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Two">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.001    0.001 0 0.001   0.001 0.001 0.001
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Three">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.002    0.001 0 0.002   0.001 0.001 0.002
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Four">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.003    0.001 0 0.003   0.001 0.001 0.003
	        		</DataItem>
        		</Geometry>
            </Grid>
        </Grid>
        <Grid Name="Parent Two" GridType="Collection">
            <Time TimeType="Range">
                <DataItem Format="XML" Dimensions="2" NumberType="Float">
                    0.0 0.5
                </DataItem>
            </Time>
            <Grid Name="One">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0    0.001 0 0   0.001 0.001 0
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Two">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.001    0.001 0 0.001   0.001 0.001 0.001
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Three">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.002    0.001 0 0.002   0.001 0.001 0.002
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Four">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.003    0.001 0 0.003   0.001 0.001 0.003
	        		</DataItem>
        		</Geometry>
            </Grid>
        </Grid>
        <Grid Name="Parent Three" GridType="Collection">
            <Time TimeType="HyperSlab">
                <DataItem Format="XML" Dimensions="3" NumberType="Float">
                    0.0 0.05 4
                </DataItem>
            </Time>
            <Grid Name="One">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0    0.001 0 0   0.001 0.001 0
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Two">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.001    0.001 0 0.001   0.001 0.001 0.001
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Three">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.002    0.001 0 0.002   0.001 0.001 0.002
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Four">
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.003    0.001 0 0.003   0.001 0.001 0.003
	        		</DataItem>
        		</Geometry>
            </Grid>
        </Grid>
        <Grid Name="Parent Four" GridType="Collection">
            <Grid Name="One">
                <Time Value="0.0" />
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0    0.001 0 0   0.001 0.001 0
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Two">
                <Time Value="0.1" />
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.001    0.001 0 0.001   0.001 0.001 0.001
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Three">
                <Time Value="0.2" />
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.002    0.001 0 0.002   0.001 0.001 0.002
	        		</DataItem>
        		</Geometry>
            </Grid>
            <Grid Name="Four">
                <Time Value="0.3" />
		        <Topology Type="TRIANGLE"
        			NumberOfElements="1">
	        		<DataItem Format="XML" NumberType="Int"
		        	Dimensions="3">
        			0 1 2 
	        		</DataItem>
        		</Topology>
        		<Geometry Type="XYZ">
	        		<DataItem Format="XML" NumberType="Float"
		        	Dimensions="3 3">
        			0 0 0.003    0.001 0 0.003   0.001 0.001 0.003
	        		</DataItem>
        		</Geometry>
            </Grid>
        </Grid>
</Domain>
</Xdmf>
"""
######## Write out an Xdmf File ###############
fd = open('TimeTest.xmf', 'w')
fd.write(TimeTestTxt)
fd.close()
##############################################
dom = XdmfDOM()
dom.Parse('TimeTest.xmf')

### Loop thru the Four Collection Grids
for GridIndex in range(1, 5) :
    print 'Reading /Xdmf/Domain/Grid[%d]' % GridIndex
    ge = dom.FindElementByPath('/Xdmf/Domain/Grid[%d]' % GridIndex)
    g = XdmfGrid()
    g.SetDOM(dom)
    g.SetElement(ge)
    g.UpdateInformation()
    # Time is available after Light Data is Read
    t = g.GetTime()
    tt = t.GetTimeType()
    print 'Time Type = ', tt, ' : ', t.GetTimeTypeAsString()
    if(tt == XDMF_TIME_SINGLE) :
        print "Value = ", t.GetValue()
    elif (tt == XDMF_TIME_LIST) :
        print "List of Times = ", t.GetArray().GetValues()
    elif (tt == XDMF_TIME_RANGE) :
        print "Min/Max = ", t.GetArray().GetValues()
    elif (tt == XDMF_TIME_HYPERSLAB) :
        print "Start/Stride/Count = ", t.GetArray().GetValues()
    # Evaluate will populate the array with the indecies of the
    # XdmfGrids that are entirely in the range
    a = XdmfArray()
    if (t.Evaluate(g, a) == XDMF_SUCCESS) :
        print "Valid times = ", a.GetValues()
        print "MinMax = ", a.GetMinAsFloat64(), ',' , a.GetMaxAsFloat64()
    else :
        print 'No Valid Times ... checking all children'
        # The "Descend" Parameter will recursively check the children
        if (t.Evaluate(g, a, 1) == XDMF_SUCCESS) :
            print "Valid times = ", a.GetValues()
            print "MinMax = ", a.GetMinAsFloat64(), ',' , a.GetMaxAsFloat64()
    a = XdmfArray()
    # By default XdmfTime::Epsilon is 1e-7. Use XdmfTime::SetEpsilon()
    # to change it
    if(g.FindGridsInTimeRange(0.05, 0.1, a) == XDMF_SUCCESS) :
        print 'Grids in Range(0.05, 0.1) = ', a.GetValues()
    else :
        print "No times are in Range(0.05, 0.1)"
    # Read Heavy Data to See the DataItems
    g.Update()
    # print g.Serialize()
