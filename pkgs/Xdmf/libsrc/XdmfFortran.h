/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfFortran.h,v 1.2 2009-05-20 20:19:45 kwleiter Exp $  */
/*  Date : $Date: 2009-05-20 20:19:45 $ */
/*  Version : $Revision: 1.2 $ */
/*                                                                 */
/*  Author:                                                        */
/*     Kenneth Leiter                                              */
/*     kenneth.leiter@arl.army.mil                                   */
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

#include <string>
#include <iostream>

#ifndef XDMFFORTRAN_H_
#define XDMFFORTRAN_H_

using std::cerr;
using std::cout;
using std::cin;
using std::endl;

class XdmfFortran{
public:
	XdmfFortran();
	~XdmfFortran();
	XdmfDOM * myDOM;
	XdmfRoot * myRoot;
	XdmfDomain * myDomain;
	XdmfTopology * myTopology;
	XdmfGeometry * myGeometry;
	std::stack<XdmfGrid*> myCollections;
	std::vector<XdmfAttribute*> myAttributes;
	std::map<char*, int> myWrittenGrids;
	char * myName;
	double currentTime;
	bool inCollection;
};

#endif /* XDMFFORTRAN_H_ */
