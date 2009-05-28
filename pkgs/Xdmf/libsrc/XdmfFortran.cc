/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfFortran.cc,v 1.1 2009-05-20 20:15:35 kwleiter Exp $  */
/*  Date : $Date: 2009-05-20 20:15:35 $ */
/*  Version : $Revision: 1.1 $ */
/*                                                                 */
/*  Author:                                                        */
/*     Kenneth Leiter                                              */
/*     kenneth.leiter@arl.army.mil                                 */
/*     US Army Research Laboratory                                 */
/*     Aberdeen Proving Ground, MD                                 */
/*                                                                 */
/*     Copyright @ 2009 US Army Research Laboratory                */
/*     All Rights Reserved                                         */
/*     See Copyright.txt or http://www.arl.hpc.mil/ice for details */
/*                                                                 */
/*     This software is distributed WITHOUT ANY WARRANTY; without  */
/*     even the implied warranty of MERCHANTABILITY or FITNESS     */
/*     FOR A PARTICULAR PURPOSE.  See the above copyright notice   */
/*     for more information.                                       */
/*                                                                 */
/*******************************************************************/

#include <Xdmf.h>
#include <XdmfSet.h>

#include <sstream>
#include <map>
#include <stack>
#include <vector>

#include "XdmfFortran.h"

// This works with g77. Different compilers require different
// name mangling.
#define XdmfInit xdmfinit_
#define XdmfSetTime xdmfsettime_
#define XdmfAddCollection xdmfaddcollection_
#define XdmfCloseCollection xdmfclosecollection_
#define XdmfSetGridTopology xdmfsetgridtopology_
#define XdmfSetGridGeometry xdmfsetgridgeometry_
#define XdmfAddGridAttribute xdmfaddgridattribute_
#define XdmfAddArray xdmfaddarray_
#define XdmfWriteGrid xdmfwritegrid_
#define XdmfWriteToFile xdmfwritetofile_
#define XdmfSerialize xdmfserialize_

XdmfFortran::XdmfFortran()
{
	myDOM = new XdmfDOM();
	myRoot = new XdmfRoot();
	myDomain = new XdmfDomain();
	currentTime = -1;
	inCollection = false;
}

XdmfFortran::~XdmfFortran()
{

}

//
// C++ will mangle the name based on the argument list. This tells the
// compiler not to mangle the name so we can call it from 'C' (but
// really Fortran in this case)
//
extern "C" {

	long XdmfInit(char *outputName)
	{
		XdmfFortran * myPointer = new(XdmfFortran);
		myPointer->myRoot->SetDOM(myPointer->myDOM);
		myPointer->myRoot->Build();
		myPointer->myRoot->Insert(myPointer->myDomain);
		myPointer->myName = outputName;
		return (long)myPointer;
	}

	void WriteToXdmfArray(XdmfArray * array, XdmfPointer * data)
	{
		switch(array->GetNumberType()){
		case XDMF_INT8_TYPE :
			array->SetValues(0, (XdmfInt8*)data, array->GetNumberOfElements());
			return;
		case XDMF_INT16_TYPE :
			array->SetValues(0, (XdmfInt16*)data, array->GetNumberOfElements());
			return;
		case XDMF_INT32_TYPE :
			array->SetValues(0, (XdmfInt32*)data, array->GetNumberOfElements());
			return;
		case XDMF_INT64_TYPE :
			array->SetValues(0, (XdmfInt64*)data, array->GetNumberOfElements());
			return;
		case XDMF_FLOAT32_TYPE :
			array->SetValues(0, (XdmfFloat32*)data, array->GetNumberOfElements());
			return;
		case XDMF_FLOAT64_TYPE :
			array->SetValues(0, (XdmfFloat64*)data, array->GetNumberOfElements());
			return;
		case XDMF_UINT8_TYPE :
			array->SetValues(0, (XdmfUInt8*)data, array->GetNumberOfElements());
			return;
		case XDMF_UINT16_TYPE :
			array->SetValues(0, (XdmfUInt16*)data, array->GetNumberOfElements());
			return;
		case XDMF_UINT32_TYPE :
			array->SetValues(0, (XdmfUInt32*)data, array->GetNumberOfElements());
			return;
		default:
			array->SetValues(0, (XdmfFloat64*)data, array->GetNumberOfElements());
			return;
		}
	}

	void XdmfSetTime(long * pointer, double * t)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		myPointer->currentTime = *t;
	}

	void XdmfAddCollection(long * pointer, char * collectionType)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		XdmfGrid * currentCollection = new XdmfGrid();
		currentCollection->SetGridType(XDMF_GRID_COLLECTION);
		currentCollection->SetCollectionTypeFromString(collectionType);
		if (myPointer->inCollection)
		{
			myPointer->myCollections.top()->Insert(currentCollection);
		}
		else
		{
			myPointer->myDomain->Insert(currentCollection);
		}
		currentCollection->Build();
		myPointer->myCollections.push(currentCollection);
		myPointer->inCollection = true;
	}

	void XdmfCloseCollection(long * pointer)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		if(myPointer->inCollection)
		{
			myPointer->inCollection = false;
			myPointer->myCollections.pop();
		}
	}

	void XdmfSetGridTopology(long * pointer, char * topologyType, int * numberOfElements, XdmfInt32 * conns)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		myPointer->myTopology = new XdmfTopology();
		myPointer->myTopology->SetTopologyTypeFromString(topologyType);
		myPointer->myTopology->SetNumberOfElements(*numberOfElements);

		// Fortran is 1 based while c++ is 0 based so
		// Either subtract 1 from all connections or specify a BaseOffset
		//myPointer->myTopology->SetBaseOffset(1);

		// If you haven't assigned an XdmfArray, GetConnectivity() will create one.
		XdmfArray * myConnections = myPointer->myTopology->GetConnectivity();
		myConnections->SetNumberOfElements(*numberOfElements * myPointer->myTopology->GetNodesPerElement());
		myConnections->SetNumberType(XDMF_INT32_TYPE);
		myConnections->SetValues(0, conns, *numberOfElements * myPointer->myTopology->GetNodesPerElement());
	}

	void XdmfSetGridGeometry(long * pointer, char * geometryType, char * numberType, int * numberOfPoints, XdmfPointer * points)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		myPointer->myGeometry = new XdmfGeometry();
		myPointer->myGeometry->SetGeometryTypeFromString(geometryType);
		myPointer->myGeometry->SetNumberOfPoints(*numberOfPoints);

		XdmfArray * myPoints = myPointer->myGeometry->GetPoints();
		myPoints->SetNumberTypeFromString(numberType);

		switch(myPointer->myGeometry->GetGeometryType())
		{
			case XDMF_GEOMETRY_XYZ :
				myPoints->SetNumberOfElements(*numberOfPoints * 3);
				break;
			case XDMF_GEOMETRY_X_Y_Z :
				myPoints->SetNumberOfElements(*numberOfPoints * 3);
				break;
			case XDMF_GEOMETRY_XY :
				myPoints->SetNumberOfElements(*numberOfPoints * 2);
				break;
			case XDMF_GEOMETRY_X_Y :
				myPoints->SetNumberOfElements(*numberOfPoints * 2);
				break;
			case XDMF_GEOMETRY_VXVYVZ :
				//TODO: FIX THIS
				myPoints->SetNumberOfElements(*numberOfPoints * 3);
				break;
			case XDMF_GEOMETRY_ORIGIN_DXDYDZ :
				myPoints->SetNumberOfElements(6);
				break;
			default:
				myPoints->SetNumberOfElements(*numberOfPoints * 3);
				break;
		}
		WriteToXdmfArray(myPoints, points);
	}

	void XdmfAddGridAttribute(long * pointer, char * attributeName, char * numberType, char * attributeCenter, char * attributeType, int * numberOfPoints, XdmfPointer * data)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		XdmfAttribute * currAttribute = new XdmfAttribute();
		currAttribute->SetName(attributeName);
		currAttribute->SetAttributeCenterFromString(attributeCenter);
		currAttribute->SetAttributeTypeFromString(attributeType);

		XdmfArray * array = currAttribute->GetValues();
		array->SetNumberTypeFromString(numberType);
		array->SetNumberOfElements(*numberOfPoints);
		WriteToXdmfArray(array, data);
		myPointer->myAttributes.push_back(currAttribute);
	}

	void XdmfAddArray(long * pointer, char * name, char * numberType, int * numberOfValues, XdmfPointer * data)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;

		XdmfSet * currSet = new XdmfSet();
    	currSet->SetDOM(myPointer->myDOM);
       	currSet->SetSetType(XDMF_SET_TYPE_NODE);
       	currSet->SetName(name);

       	// Copy Elements from Set to XdmfArray
       	XdmfArray * array = currSet->GetIds();
       	array->SetNumberTypeFromString(numberType);
		array->SetNumberOfElements(*numberOfValues);
		std::stringstream heavyDataName;
		heavyDataName << myPointer->myName << ".h5:/" <<  name;;
		array->SetHeavyDataSetName(heavyDataName.str().c_str());
		WriteToXdmfArray(array, data);

		if (myPointer->inCollection)
		{
			myPointer->myCollections.top()->Insert(currSet);
		}
		else
		{
			XdmfGrid * myGrid = new XdmfGrid();
			myGrid->SetDOM(myPointer->myDOM);
	        myGrid->SetElement(myPointer->myDOM->FindElement("Domain"));
	        myGrid->Insert(currSet);
		}

        currSet->Build();
	}

	void XdmfWriteGrid(long * pointer, char * gridName)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		XdmfGrid * grid = new XdmfGrid();

		std::stringstream totalGridName;
		if(myPointer->myWrittenGrids.find(gridName) == myPointer->myWrittenGrids.end())
		{
			myPointer->myWrittenGrids[gridName] = 1;
			totalGridName << gridName;
		}
		else
		{
			myPointer->myWrittenGrids[gridName]++;
			totalGridName << gridName << "_" << myPointer->myWrittenGrids[gridName];
		}

		grid->SetName(totalGridName.str().c_str());

		//Modify HDF5 names so we aren't writing over top of our data!
		std::stringstream topologyDataName;
		topologyDataName << myPointer->myName << ".h5:/" <<  totalGridName.str() << "/Connections";
		myPointer->myTopology->GetConnectivity()->SetHeavyDataSetName(topologyDataName.str().c_str());
		grid->SetTopology(myPointer->myTopology);

		std::stringstream geometryDataName;
		geometryDataName << myPointer->myName << ".h5:/" <<  totalGridName.str() << "/XYZ";
		myPointer->myGeometry->GetPoints()->SetHeavyDataSetName(geometryDataName.str().c_str());
		grid->SetGeometry(myPointer->myGeometry);

		if (myPointer->inCollection)
		{
			myPointer->myCollections.top()->Insert(grid);
		}
		else
		{
			myPointer->myDomain->Insert(grid);
		}

		XdmfTime * t = new XdmfTime();
		if (myPointer->currentTime >= 0)
		{
			t->SetTimeType(XDMF_TIME_SINGLE);
			t->SetValue(myPointer->currentTime);
			grid->Insert(t);
		}

		while(myPointer->myAttributes.size() > 0)
		{
			XdmfAttribute * currAttribute = myPointer->myAttributes.back();

			std::stringstream attributeDataName;
			attributeDataName << myPointer->myName << ".h5:/" <<  totalGridName.str() << "/" << currAttribute->GetName();
			currAttribute->GetValues()->SetHeavyDataSetName(attributeDataName.str().c_str());

			grid->Insert(currAttribute);
			myPointer->myAttributes.pop_back();
		}

		grid->Build();
	}

	void XdmfWriteToFile(long * pointer)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		std::stringstream dataName;
		dataName << myPointer->myName << ".xmf";
		myPointer->myDOM->Write(dataName.str().c_str());
	}

	void XdmfSerialize(long * pointer)
	{
		XdmfFortran * myPointer = (XdmfFortran *)*pointer;
		cout << myPointer->myDOM->Serialize() << endl;
	}
}
