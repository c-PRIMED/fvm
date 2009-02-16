#include "Octree.h"
#include "NumType.h"

//#include <cstdio>
//#include <stdlib.h>
//#include <string.h>

typedef Octree::Point Point;
typedef Octree::Bounds Bounds;
typedef Octree::VectorT3Array VectorT3Array;
typedef Octree::VectorT3 VectorT3;
using namespace std;

inline double max ( double x, double y)
{
  return (x>y ? x:y);
}

inline double min ( double x, double y)
{
  return (x>y ? y:x);
}
// -----------------------------------------------------------------------------
// Construction -- "nullify" the class
// -----------------------------------------------------------------------------

Octree::Octree()
 {
  memset(_child, 0, sizeof(_child));
  _pointCount=0;
  _points=NULL;
  _center=0.0;
  _radius=0.0;
  _nodeType=0;
}

// -----------------------------------------------------------------------------
// Destruction -- free up memory
// -----------------------------------------------------------------------------
//???how to delete the entire octree//
Octree::~Octree()
{
  // delete[] _data;
}

// -----------------------------------------------------------------------------
// Build the octree
// -----------------------------------------------------------------------------
        // count: number of points
        // threshold:maximum number of points each node contains
        // maximumDepth: 
const   bool    Octree::build(Point *points,
			      const unsigned int count,
                              const unsigned int threshold,
                              const unsigned int maximumDepth,
                              const Bounds &bounds,
                              const unsigned int currentDepth
			      )
       
{
        //store the information for current node 
	_pointCount=count;
	_center=bounds.center;
	_radius=bounds.radius;
	_currentDepth=currentDepth;

        // The node is a leaf when...
        // 1. The number of points  <= the threshold
        // 2. We've recursed too deep into the tree
        //    (currentDepth >= maximumDepth)
        //
        //    NOTE: We specifically use ">=" for the depth comparison so that we
        //          can set the maximumDepth depth to 0 if we want a tree with
        //          no depth.

        if (count <= threshold || currentDepth >= maximumDepth)
        {
                // Just store the points in the node, making it a leaf

                _points = new Point [count];
		memcpy(_points, points, sizeof (Point)* count);

		_nodeType=1;
		    
	        return true;
        }
        // else, the current node is a intermediate node. 

	//set the data and node type
	_nodeType=0;
	_points=NULL;

	// then creat the child nodes
	

        unsigned int    childPointCounts[8]={0};	

        // Loop over all the points in current node, 
        // Classify each point to a child node

        for (unsigned int i = 0; i < count; i++)
        {
                // Current point

                Point   &p = points[i];

                // Center of this node

                const VectorT3 &c = bounds.center;

                // Here, we need to know which child each point belongs to. To
                // do this, we build an index into the _child[] array using the
                // relative position of the point to the center of the current
                // node

                p.code = 0;
		// bitwise OR x|y  each bit in x OR each bit in y
		// bitwise OR assignment x|=y assign x|y to x
                if (p.coordinate[0] > c[0]) p.code |= 1;
                if (p.coordinate[1] > c[1]) p.code |= 2;
                if (p.coordinate[2] > c[2]) p.code |= 4;

                //  keep track of how many points get assigned in each child 

                childPointCounts[p.code]++;
        }

	 // Recursively call build() for each of the 8 children

        for (int i = 0; i < 8; i++)
        {
                // Don't bother going any further if there aren't any points for
                // this child

                if (!childPointCounts[i]) continue;

                // Allocate the child

                _child[i] = new Octree;

                // Allocate a list of points that were coded JUST for this child
                // only

                Point   *newList = new Point [childPointCounts[i]];

                // Go through the input list of points and copy over the points
                // that were coded for this child

                Point   *ptr = newList;

                for (unsigned int j = 0; j < count; j++)
                {
                        if (points[j].code == i)
                        {
                                *ptr = points[j];
                                ptr++;
                        }
                }
                
                // At this point, we have a list of points that will belong to
                // this child node

		// assume newCount is the same as childPointCount
		int newCount=childPointCounts[i];

                // Generate a new bounding volume //
                // We use a table of offsets. These offsets determine where a
                // node is, relative to it's parent. So, for example, if want to
                // generate the bottom-left-rear (-x, -y, -z) child for a node,
                // we use (-1, -1, -1).
                // 
                // However, since the radius of a child is always half of its
                // parent's, we use a table of 0.5, rather than 1.0.
                // 
                // These values are stored the following table. 
		static VectorT3 boundsOffsetTable[8];
	
		boundsOffsetTable[0][0]=-0.5;
		boundsOffsetTable[0][1]=-0.5;
		boundsOffsetTable[0][2]=-0.5;
		
		boundsOffsetTable[1][0]=+0.5;
		boundsOffsetTable[1][1]=-0.5;
		boundsOffsetTable[1][2]=-0.5;
		
		boundsOffsetTable[2][0]=-0.5;
		boundsOffsetTable[2][1]=+0.5;
		boundsOffsetTable[2][2]=-0.5;

		boundsOffsetTable[3][0]=+0.5;
		boundsOffsetTable[3][1]=+0.5;
		boundsOffsetTable[3][2]=-0.5;

		boundsOffsetTable[4][0]=-0.5;
		boundsOffsetTable[4][1]=-0.5;
		boundsOffsetTable[4][2]=+0.5;

		boundsOffsetTable[5][0]=+0.5;
		boundsOffsetTable[5][1]=-0.5;
		boundsOffsetTable[5][2]=+0.5;
		
		boundsOffsetTable[6][0]=-0.5;
		boundsOffsetTable[6][1]=+0.5;
		boundsOffsetTable[6][2]=+0.5;

		boundsOffsetTable[7][0]=+0.5;
		boundsOffsetTable[7][1]=+0.5;
		boundsOffsetTable[7][2]=+0.5;


	
		/*	   {     {-0.5, -0.5, -0.5},
			{+0.5, -0.5, -0.5},
			{-0.5, +0.5, -0.5},
			{+0.5, +0.5, -0.5},
			{-0.5, -0.5, +0.5},
			{+0.5, -0.5, +0.5},
			{-0.5, +0.5, +0.5},
			{+0.5, +0.5, +0.5}
			};*/

                // Calculate our offset from the center of the parent's node to
                // the center of the child's node

                VectorT3   offset= boundsOffsetTable[i]*bounds.radius;
	
                // Create a new Bounds, with the center offset and half the
                // radius

                Bounds newBounds;
                newBounds.radius = bounds.radius * 0.5;
                newBounds.center = bounds.center + offset;

                // Recurse

		  _child[i]->build(newList, newCount, threshold, maximumDepth,
				   newBounds, currentDepth+1);

                // Clean up

                delete[] newList;
        }

	return true;
}

// -----------------------------------------------------------------------------
// Determine the [cubic] bounding volume for a set of points
// -----------------------------------------------------------------------------

const Octree::Bounds Octree::calcCubicBounds(const Point * points,
                                     const unsigned int count)
{
        // What will be returned to the caller

        Bounds  b;

        // Determine min/max of the given set of points

        VectorT3   min = points[0].coordinate;
        VectorT3   max = points[0].coordinate;

        for (unsigned int i = 0; i < count; i++)
        {
                const VectorT3 &p = points[i].coordinate;
                if (p[0] < min[0]) min[0] = p[0];
                if (p[1] < min[1]) min[1] = p[1];
                if (p[2] < min[2]) min[2] = p[2];
                if (p[0] > max[0]) max[0] = p[0];
		if (p[1] > max[1]) max[1] = p[1];
                if (p[2] > max[2]) max[2] = p[2];
		//if (p VectorT3::< min) min=p;
		//if (p VectorT3::> max) max=p;
		//min=setMin(min,p);
		//max=setMax(max,p);
        }

        // The radius of the volume (dimensions in each direction)

        VectorT3   radius;
	radius=(max-min)/2.;

        // Find the center of this space
	b.center=min+radius;
	

        // We want a CUBIC space. By this, I mean we want a bounding cube, not
        // just a bounding box. We already have the center, we just need a
        // radius that contains the entire volume. To do this, we find the
        // maxumum value of the radius' X/Y/Z components and use that

        b.radius = radius[0];
        if (b.radius < radius[1]) b.radius = radius[1];
        if (b.radius < radius[2]) b.radius = radius[2];

        // Done

        return b;
}

// -----------------------------------------------------------------------------
// Print out the octree leaf nodes
// -----------------------------------------------------------------------------

const bool Octree::report(FILE *fp)
{
  
  //if it is a leaf, then print out the data
  if (_nodeType==1){
    fprintf(fp,"currentDepth is  %i\n", _currentDepth);
    fprintf(fp,"node center is %f\t%f\t%f\n", _center[0],_center[1],_center[2]);
    fprintf(fp,"node radius is %f\n", _radius);
    fprintf(fp,"data in this node is ");
    for(unsigned int i=0; i<_pointCount; i++){
      fprintf(fp,"%i\t", _points[i].cellIndex);
    }
    fprintf(fp,"\n end of node\n");
  }
  //if it is a node, recursively traverse to child node
  if(_nodeType==0){
    for(int i=0; i<8;i++){
      //not guarantee a node has all 8 children
      if(!_child[i]) continue;
      else{
	_child[i]->report(fp);
      }
    }
  }
  return(true);
}
 


 /// <summary> A utility method to figure out the closest distance of a border
 /// to a point. If the point is inside the bounds, return 0.
 /// <returns> closest distance to the point. </returns>

const double Octree::borderDistance(VectorT3 coordinate){
  double x=coordinate[0];
  double y=coordinate[1];
  double z=coordinate[2];
  double nsdistance;
  double ewdistance;
  double fbdistance;
    
  double leftBound=_center[0]-_radius;
  double rightBound=_center[0]+_radius;
  double frontBound=_center[1]-_radius;
  double backBound=_center[1]+_radius;
  double bottomBound=_center[2]-_radius;
  double topBound=_center[2]+_radius;

  if(leftBound <= x && x <= rightBound)
    ewdistance = 0;
  else
    ewdistance = min(fabs(x-rightBound), fabs(x-leftBound));

  if(frontBound <= y && y <= backBound)
    fbdistance = 0;
  else
    fbdistance = min(fabs(y-backBound), fabs(y-frontBound));
 
  if(bottomBound <= z && z <= topBound)
    nsdistance = 0;
  else
    nsdistance = min(fabs(z-topBound), fabs(z-bottomBound));

  return (nsdistance * nsdistance +
                   ewdistance * ewdistance +
                   fbdistance * fbdistance);
  //note: it actually return the distance square
}

/// <summary> Get an object closest to a x/y/z. If there are branches at
/// this node, then the branches are searched. The branches are
/// checked first, to see if they are closer than the best distance
/// already found. If a closer object is found, bestDistance will
/// be updated with a new Double object that has the new distance.</summary>
/// <param name="x">left-right location in Octree Grid</param>
/// <returns> the object that matches the best distance, null if no closer objects were found.</returns>

const int Octree::getNode(VectorT3 coordinate,  double& shortestDistance)
{
  static int node;
  double distance;
  double closest;
  VectorT3 dR;

  //if it is a leaf, search all data in this leaf to find out the best distance
  if (_nodeType==1){
    for (unsigned int i=0;  i<_pointCount; i++){
      dR=_points[i].coordinate-coordinate;
      //dR=0;
      distance=mag2(dR);
      if(distance < shortestDistance*shortestDistance){
	shortestDistance=sqrt(distance);
        node=_points[i].cellIndex;
      }
    }
    return node;
  }
   //if it is a node, recursively traverse to child node
  else if(_nodeType==0){
     // Check the distance of the bounds of the branch,
     // versus the bestDistance. If there is a boundary that
     // is closer, then it is possible that another node has an
     // object that is closer.
    for (int i=0; i<8;i++){
      //not guarantee a node has all 8 children
      if(!_child[i]) continue;
      else{
	double childDistance=_child[i]->borderDistance(coordinate);
	if (childDistance < shortestDistance*shortestDistance){
	  int test=-1;
	  test=_child[i]->getNode(coordinate, shortestDistance);
	  if(test >= 0)
	    node=test;
	}
      }
    }
  }

  return(node);
}





// -----------------------------------------------------------------------------
// Get an node closest to a (x,y,z)
// <returns> the data in the node that matches the best distance,
//  null if no node were found.
// for this case, the data is an integer (cellIndex)
// -----------------------------------------------------------------------------    
	    
const int Octree::getNode(double x, double y, double z)
{
  int node;
  VectorT3 coordinate;
  double LargeNumber=1.0e20;
  coordinate[0]=x;
  coordinate[1]=y;
  coordinate[2]=z;
  node=Octree::getNode(coordinate, LargeNumber);
  return(node); 
}

const int Octree::getNode(VectorT3 coordinate)
{
  int node;
  double LargeNumber=1.0e20;
  node=Octree::getNode(coordinate, LargeNumber);
  return(node); 
}

/**********************************************************************/
/// <summary> Get all objects closest to a x/y/z within a radius. 
/// search mechanism similar to getNode(coordinate)
/**********************************************************************/

const vector<int> Octree::getNodes(VectorT3 coordinate,  double radius)
{
  static vector<int> node;
  double distance;
  VectorT3 dR;

  //if it is a leaf, search all data in this leaf to find out the best distance
  if (_nodeType==1){
    for (unsigned int i=0;  i<_pointCount; i++){
      dR=_points[i].coordinate-coordinate;
      distance=mag2(dR);
      if(distance < radius*radius){
	node.push_back(_points[i].cellIndex);
      }
    }
    return node;
  }
   //if it is a node, recursively traverse to child node
  else if(_nodeType==0){
     // Check the distance of the bounds of the branch,
     // versus the bestDistance. If there is a boundary that
     // is closer, then it is possible that another node has an
     // object that is closer.
    for (int i=0; i<8;i++){
      //not guarantee a node has all 8 children
      if(!_child[i]) continue;
      else{
	double childDistance=_child[i]->borderDistance(coordinate);
	if (childDistance < radius*radius){
	   _child[i]->getNodes(coordinate, radius);	 
	}
      }
    }
  }

  return(node);
}

/**********************************************************************/
/// <summary> Get all objects closest to a x/y/z within a radius. 
/// by simply looping over all the points
/**********************************************************************/

const int Octree::Naive_getNode(VectorT3 coordinate, int count, Point * points)
{
  
  //naive search by looping over all points
    double shortestDistance=100000;
    double shortestDistanceSqr=shortestDistance*shortestDistance;
    int cellIndex;

    for (int i=0; i<count; i++){
      VectorT3 dR=coordinate-points[i].coordinate;
      double distanceSqr=mag2(dR);
      if(distanceSqr < shortestDistanceSqr){
	shortestDistance=sqrt(distanceSqr);
	shortestDistanceSqr=shortestDistance*shortestDistance;
	cellIndex=points[i].cellIndex;
      }
    }
    return(cellIndex);
}


const vector<int> Octree::Naive_getNodes(VectorT3 coordinate, int count, Point * points, double radius)
{
  
  //naive search by looping over all points
   
    vector<int> cellIndexList;
    double radiusSqr=radius*radius;
    for (int i=0; i<count; i++){
      VectorT3 dR=coordinate-points[i].coordinate;
      double distanceSqr=mag2(dR);
      if(distanceSqr < radiusSqr){
	cellIndexList.push_back(points[i].cellIndex);
      }
    }
    return(cellIndexList);
}

const bool Octree::MPM_Points_Write(char *file)

{

  //here, we want to build up a solid circle 
  //the number of solid points and coordiantes are written to file
    FILE *fp;
   

    int nX=20, nY=20, nZ=1;
    double gapX=1.0/nX, gapY=1.0/nY, gapZ=1.0/nZ;
    double radius=0.2;
    VectorT3 center;
    center[0]=0.5;
    center[1]=0.5;
    center[2]=0.0;

    int count=0;
    VectorT3 temp;
    VectorT3 solidPoint[nX*nY*nZ];
    
    for(int i=0; i<nX; i++){
      for(int j=0; j<nY; j++){
	for(int k=0; k<nZ; k++){
	  temp[0]=i*gapX;
	  temp[1]=j*gapY;
	  temp[2]=k*gapZ;
	  VectorT3 ds=temp-center;
	  if(mag2(ds) <= radius*radius){
	    solidPoint[count][0]=temp[0];
	    solidPoint[count][1]=temp[1];
	    solidPoint[count][2]=temp[2];	  
	    count+=1;
	  }
	}
      }
    }

    fp=fopen(file,"w");
    fprintf(fp,"%i\n",count);
    for(int p=0; p<count; p++){
      fprintf(fp, "%lf\t%lf\t%lf\n", solidPoint[p][0],solidPoint[p][1],solidPoint[p][2]);
    } 
    fclose(fp);
    return (true);
}
    
const shared_ptr<VectorT3Array> Octree::MPM_Points_Read(char *file)

{
    FILE *fp;
    int nMPM;
    double x=0, y=0, z=0;

    fp=fopen(file,"r");
    fscanf(fp,"%i\n",&nMPM);
    //VectorT3 MPM_Points[nMPM];
    shared_ptr<VectorT3Array> MPM_Points ( new VectorT3Array(nMPM));
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &x, &y, &z);
      (*MPM_Points)[i][0]=x;
      (*MPM_Points)[i][1]=y;
      (*MPM_Points)[i][2]=z;
    }
    fclose(fp);
   
    return (MPM_Points);
}
