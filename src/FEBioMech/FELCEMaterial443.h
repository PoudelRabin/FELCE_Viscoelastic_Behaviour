#pragma once
#include "FEElasticMaterial.h"
#include <FECore/FEMesh.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEModelParam.h>
#include <FECore/FENode.h>



//-----------------------------------------------------------------------------
//! LCE material with Neo-classical + frank energy

class FEBIOMECH_API FELCEMaterial443 : public FEElasticMaterial
{


public:
	FELCEMaterial443(FEModel* pfem);

public:
	FEParamDouble		m_mu;	//!< Shear modulus
	FEParamDouble		m_kfrank;	//!< frank constant
    FEParamDouble		m_a;        //!< stretch parameter
    FEParamDouble		m_bulk;     //!< bulk modulus
    
public:
    //! calculate stress at material point
    virtual mat3ds Stress(FEMaterialPoint& pt) override;

    //! calculate tangent stiffness at material point
    virtual tens4ds Tangent(FEMaterialPoint& pt) override;

    //! calculate the stresses and its derivatives at the material point
    void MaterialRoutine(FESolidElement el,int ngauss, matrix Bt, mat3d F, matrix t_e,
                double pressmix, matrix &PK1_vgt, matrix &Pi0,
               matrix &gt, matrix &dPdF_vgt, matrix &dPdt_vgt,
               matrix &dgdF_vgt, mat3d &dPidGradt, matrix &dgdt,
               matrix &dJdF_vgt);

    //Change 3by3 tensor to 9by1 voigt matrix
    void tens33_9(mat3d A, matrix &B);
               
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};    