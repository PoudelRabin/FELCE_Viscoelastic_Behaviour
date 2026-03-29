#pragma once
#include "FEElasticMaterial.h"
#include <FECore/FEMesh.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEModelParam.h>
#include <FECore/FENode.h>
#include "FELCEMaterialPoint.h"




//-----------------------------------------------------------------------------
//! LCE material with Neo-classical + frank energy

class FEBIOMECH_API FELCEMaterial444 : public FEElasticMaterial
{

//default constructer
public:
	FELCEMaterial444(FEModel* pfem);

//material parameters
public:
	FEParamDouble		m_mu_nc;	//!< Shear modulus of neoclassical 
    FEParamDouble		m_mu_nh;	//!< Shear modulus of neohokkean
    FEParamDouble		m_mu_neq_nc;	//!< Shear modulus of non-equilibrium network
	FEParamDouble		m_eta_u;	//!< viscocity of non-equilibrium network
    FEParamDouble		m_mu_neq_nh;	//!< Shear modulus of non-equilibrium rotation
    FEParamDouble		m_eta_t;	//!< viscocity of non-equilibrium rotation
	FEParamDouble		m_kfrank;	//!< frank constant
    FEParamDouble		m_a;        //!< stretch parameter


public:

    //! get the elastic base material
    FEElasticMaterial* GetBaseMaterial() { return m_Base; }

    //! Set the base material
    void SetBaseMaterial(FEElasticMaterial* pbase) {m_Base = pbase;}

    //! create material point data for this material
    FEMaterialPointData* CreateMaterialPointData() override ;
    
public:
    //! calculate stress at material point
    virtual mat3ds Stress(FEMaterialPoint& pt) override;

    //! calculate tangent stiffness at material point
    virtual tens4ds Tangent(FEMaterialPoint& pt) override;
    
    //! Calculate LCE stresses PK1_vgt, Pi0, gt and returns mat3d PK1
    mat3d StressLCE(FESolidElement el, int ngauss,  
        mat3d F, matrix Gradt, double pressmix, matrix &PK1_vgt, matrix &Pi0,
        matrix &gt);
            
    //! Calculate all the tangents needed for the LCE stiffness calculation
    void TangentLCE(FESolidElement el, int ngauss, 
          mat3d F, matrix Gradt, double pressmix, matrix &dPdF_vgt, 
          matrix &dPdt_vgt, matrix &dgdF_vgt, mat3d &dPidGradt, matrix &dgdt,
          matrix &dJdF_vgt);

    //Change 3by3 tensor to 9by1 voigt matrix
    void tens33_9(mat3d A, matrix &B);

    //change tens4d to 9by9 voigt matrix
    void tens3333_99(tens4d A, matrix& B);

private:
    FEElasticMaterial* m_Base;	//!< pointer to elastic solid material
 
public:
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};    