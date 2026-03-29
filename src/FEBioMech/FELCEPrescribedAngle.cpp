#include "stdafx.h"
#include "FELCEPrescribedAngle.h"

//=======================================================================================
// NOTE(from original febio): I'm setting FEBoundaryCondition is the base class since I don't want to pull
//       in the parameters of FEPrescribedDOF. 
BEGIN_FECORE_CLASS(FELCEPrescribedAngle, FENodalBC)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list:director angle)");
	ADD_PARAMETER(m_scale, "value");
END_FECORE_CLASS();

FELCEPrescribedAngle::FELCEPrescribedAngle(FEModel* fem) : FEPrescribedDOF(fem)
{
}