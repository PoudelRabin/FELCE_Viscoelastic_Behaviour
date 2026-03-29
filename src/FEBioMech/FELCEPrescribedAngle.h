#pragma once
#include <FECore/FEPrescribedDOF.h>

class FELCEPrescribedAngle : public FEPrescribedDOF
{
public:
	FELCEPrescribedAngle(FEModel* fem);

private:
	DECLARE_FECORE_CLASS();
};
