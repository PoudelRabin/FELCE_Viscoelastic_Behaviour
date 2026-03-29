#pragma once
#include <FECore/FESolidDomain.h>
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"
#include <FECore/FEDofList.h>
#include "FEElasticMaterialPoint.h"
#include "FELCEMaterial443.h"




//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! volume element of LCE 
class FEBIOMECH_API FE3FieldLCEDomain : public FESolidDomain, public FEElasticDomain
{

public:
    //! constructor
	FE3FieldLCEDomain(FEModel* pfem);

public:
    
    //! assignment operator
	FE3FieldLCEDomain& operator = (FE3FieldLCEDomain& d);
    
	//! activate
	void Activate() override;

	//! initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! Set flag for update for dynamic quantities
	void SetDynamicUpdateFlag(bool b);

	//! serialization
	void Serialize(DumpStream& ar) override;

	// get the total dof list
	const FEDofList& GetDOFList() const override;

public: // overrides from FEDomain

	//! get the material
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pm) override;
    
public: // overridden from FEElasticDomain

	// update stresses
	void Update(const FETimeInfo& tp) override;
    
    // update the stress of an element
	void UpdateElementStress(int iel, const FETimeInfo& tp) ;

	//! intertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, vector<double>& F) override;

	//! body forces
	void BodyForce(FEGlobalVector& R, FEBodyForce& BF) override;

	//! calculates inertial stiffness
	void MassMatrix(FELinearSystem& LS, double scale) override;

	//! body force stiffness
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;

	// calculate Global stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS) override;
    
    //! calculate the internal forces
	void InternalForces(FEGlobalVector& R) override;

protected:
	//! calculates the element internal residual forces
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

	//! Calculates element stiffness element matrix
	void ElementStiffness(FESolidElement& el, matrix& ke);

	//! Calculates the element stiffness and residual matrices
    void ElementRoutine(FESolidElement& el, double bulkmod, matrix &Kelem, matrix &Relem);
    
    //! Caluclates Bt matrix that interpolate Grad of theta
    void Btmatrix(FESolidElement el, int ngauss, matrix &Bt);
    
    //! Calcualtes Bu matrix that interpolate Grad of uflow
    void Bumatrix(matrix Bt, matrix &Bu);
    
    //! Calculate unassembled stiffness matrices
    void stiffnessmatrices(FESolidElement el, int ngauss,   
        matrix Bu, matrix Bt, double bulkmod,
        matrix dPdF_vgt, matrix dPdt_vgt,matrix dgdF_vgt,
        mat3d dPidGradt,matrix dgdt, matrix dJdF_vgt,  
        matrix &K11,double &K22, double &K23,double &K32, 
        matrix &K13, matrix &K31, matrix &K_ut, 
        matrix &K_tu,matrix &K_tt);
        
    //! Calculate unassembled residual matrices
    void residualmatrices(FESolidElement el, int ngauss, 
        matrix Bu, matrix Bt, matrix PK1_vgt, matrix Pi0,
        matrix gt, matrix &ru, matrix &rt);

	

protected:
    double              m_alphaf;
    double              m_alpham;
    double              m_beta;
	bool				m_update_dynamic;	//!< flag for updating quantities only used in dynamic analysis

	bool	m_secant_stress;	//!< use secant approximation to stress
	bool	m_secant_tangent;   //!< flag for using secant tangent

protected:
	FEDofList	m_dofU;		// displacement dofs
    FEDofList	m_dofDA;    // director angle
	FEDofList	m_dof;		// total dof list


private:
	FELCEMaterial443* m_pMat;

}; 


class FEStandard3FieldLCEDomain : public FE3FieldLCEDomain
{
public:
	FEStandard3FieldLCEDomain(FEModel* fem);

private:
	std::string		m_elemType;

	DECLARE_FECORE_CLASS();
};