#include "XdmfSTLConverter.txx"

#include <XdmfDOM.h>
#include <XdmfDomain.h>
#include <XdmfRoot.h>

#include <vector>
#include <map>

int main( int argc, const char* argv[] )
{
	/*
	 * Create an STL vector and map and write out to XDMF
	 */

	XdmfSTLConverter * mySerializer = new XdmfSTLConverter();
	XdmfDOM * myDOM = new XdmfDOM();

	XdmfRoot * myRoot = new XdmfRoot();
	myRoot->SetDOM(myDOM);
	myRoot->Build();

	XdmfDomain * myDomain = new XdmfDomain();
	myRoot->Insert(myDomain);

	std::vector<int> intVector;
	int intTotal = 100 ;
	for (int i=1; i<6; i++)
	{
		intTotal += intTotal;
		intVector.push_back(intTotal);
	}

	std::map<int, double> myMap;
	myMap[12]= 13.785412486;
	myMap[60] = 45451.154987;
	myMap[-60] = 5231.554874;

	mySerializer->writeSetToXdmf(intVector, myDOM, "Vector");
	mySerializer->writeMapToXdmf(myMap, myDOM, "Map");

	cout << myDOM->Serialize() << endl;

	/*
	 * Now attempt to do the reverse.  Read XDMF information into STL containers:
	 */

	XdmfXmlNode domain = myDOM->FindElement("Domain");
	for (int i=0; i < myDOM->FindNumberOfElements("Set", domain); i++)
	{
	    XdmfSet * currSet = new XdmfSet();
	    currSet->SetDOM(myDOM);
	    currSet->SetElement(myDOM->FindElement("Set", i, domain));
	    currSet->Update();
	    if (currSet->GetNumberOfAttributes() == 0)
	    {
	    	// Must be Non-Associative
	    	std::vector<int> myData;
	    	mySerializer->getSetFromXdmf(myData, currSet);
	    	cout << "\nVECTOR CONTAINS:" << endl;
		    for (unsigned int j=0; j<myData.size(); j++)
		    {
		    	cout << myData[j] << endl;
		    }
	    }
	    else
	    {
	    	// Must be Associative
	    	std::map<int,double> myData;
	    	mySerializer->getMapFromXdmf(myData, currSet);
	    	cout << "\nMAP CONTAINS:" << endl;
	    	std::map<int,double>::iterator iter;
	    	for (iter = myData.begin(); iter != myData.end(); iter++)
	    	{
	    		cout << iter->first << '\t' << iter->second << endl;
	    	}
	    }
	}
}
