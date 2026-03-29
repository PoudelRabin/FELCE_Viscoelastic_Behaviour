#pragma once
#include <FECore/FEPlotData.h>
#include <FECore/FEElement.h>
#include <FECore/units.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal LCE angle
//!
class FEPlotNodeLCEAngle : public FEPlotNodeData
{
public:
	FEPlotNodeLCEAngle(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) { }
	bool Save(FEMesh& m, FEDataStream& a);
};
//-----------------------------------------------------------------------------
//! Cosine of Nodal director angle
//!
class FEPlotCosineNodeLCEAngle : public FEPlotNodeData
{
public:
	FEPlotCosineNodeLCEAngle(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) { }
	bool Save(FEMesh& m, FEDataStream& a);
};
    
//=============================================================================
//							D O M A I N   D A T A
//=============================================================================
//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementLCEPK1Stress : public FEPlotDomainData
{
public:
	FEPlotElementLCEPK1Stress(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3F, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
	bool Save(FEDomain& dom, FEDataStream& a);
};
//-----------------------------------------------------------------------------
//! Element pressure
class FEPlotElementLCEPressure : public FEPlotDomainData
{
public:
	FEPlotElementLCEPressure(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) { SetUnits(UNIT_PRESSURE); }
	bool Save(FEDomain& dom, FEDataStream& a);
};
