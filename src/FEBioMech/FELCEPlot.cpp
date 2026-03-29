#include "stdafx.h"
#include "FELCEPlot.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FEReactiveFatigue.h"
#include "FEReactivePlasticity.h"
#include "FEReactivePlasticDamage.h"
#include "FERemodelingElasticMaterial.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FEElasticShellDomainOld.h"
#include "FEElasticEASShellDomain.h"
#include "FEElasticANSShellDomain.h"
#include "FEElasticMixture.h"
#include "FEElasticMultigeneration.h"
#include "FEUT4Domain.h"
#include "FEContactSurface.h"
#include "FERigidBody.h"
#include <FECore/FESPRProjection.h>
#include "FEUncoupledElasticMixture.h"
#include "FERigidMaterial.h"
#include "FEVolumeConstraint.h"
#include "FEFacet2FacetSliding.h"
#include "FEMortarSlidingContact.h"
#include "FEMechModel.h"
#include "FEPreStrainElastic.h"
#include <FECore/writeplot.h>
#include <FECore/FEDomainParameter.h>
#include <FECore/FEModel.h>
#include "FEDiscreteElasticMaterial.h"
#include "FEDiscreteElasticDomain.h"
#include "FEContinuousElasticDamage.h"
#include <FECore/FEMeshAdaptor.h> // for projectToNodes
#include "FESlidingInterface.h"
#include "FESlidingElasticInterface.h"
#include "FETiedContactSurface.h"
#include "FEReactiveVEMaterialPoint.h"
#include <FECore/FESurface.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FETrussDomain.h>
#include <FECore/FEElement.h>
#include "FELCEDomain.h"

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotNodeLCEAngle::Save(FEMesh& m, FEDataStream& a)
{
	FEModel* fem = GetFEModel();
	const int dof_LCEAngle = fem->GetDOFIndex("vart");

	writeNodalValues<double>(m, a, [=](const FENode& node) {
		return node.get(dof_LCEAngle);
		});
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotCosineNodeLCEAngle::Save(FEMesh& m, FEDataStream& a)
{
	FEModel* fem = GetFEModel();
	const int dof_LCEAngle = fem->GetDOFIndex("vart");

	writeNodalValues<double>(m, a, [=](const FENode& node) {
		return cos(node.get(dof_LCEAngle));
       
		});
	return true;
}


//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Store the average PK1 stress for each element. 
bool FEPlotElementLCEPK1Stress::Save(FEDomain& dom, FEDataStream& a)
{
    // get the domain
	FELCEDomain& Ldom = static_cast<FELCEDomain&>(dom);
	
    for (int i = 0; i < dom.Elements(); ++i)
	{
        FESolidElement& el = Ldom.Element(i); 
		mat3d AvgPK1 = Ldom.LCEAveragePK1Stress(el);
        
        a << AvgPK1;
    }
	
    
    return true;
}
//-----------------------------------------------------------------------------
//! Store the average PK1 stress for each element. 
bool FEPlotElementLCEPressure::Save(FEDomain& dom, FEDataStream& a)
{
	// get the domain
	FELCEDomain& Ldom = static_cast<FELCEDomain&>(dom);

	for (int i = 0; i < dom.Elements(); ++i)
	{
		FESolidElement& el = Ldom.Element(i);
		double Pressure = Ldom.PressureMix(el);

		a << Pressure;
	}

	return true;
}