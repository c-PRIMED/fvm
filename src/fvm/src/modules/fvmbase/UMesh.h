// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _UMESH_H_
#define _UMESH_H_

#include "Mesh.h"
#include "Array.h"
#include "Vector.h"
#include "StorageSite.h"

class CRConnectivity;

#include "base_public.h"

class  BASE_PUBLIC UMesh : public Mesh
{
public:

  UMesh(const Args& args);

  virtual ~UMesh();

  DECLARE_HT("UMesh");

  const CRConnectivity& getConnectivity(const StorageSite& from,
                                        const StorageSite& to) const;

  const CRConnectivity& getAllFaceCells() const;
  const CRConnectivity& getFaceCells(const StorageSite& site) const;

  const Array<int>& getCellTypes() const;
  const Array<int>& getCellTypeCount() const;

  DECLARE_METHOD(maskCellsOfType);
  DECLARE_METHOD(getCellTypeCount);
  DECLARE_METHOD(getCellTypes);
  
private:
  const CRConnectivity& _allFaceCells;
  map<const StorageSite*,const CRConnectivity*> _faceCellsMap;
  mutable Array<int> *_cellTypes;
  mutable Array<int> *_cellTypeCount;
};

#endif
