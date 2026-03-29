#pragma once
#include <FECore/FESolidDomain.h>
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"
#include <FECore/FEDofList.h>
#include "FEElasticMaterialPoint.h"
#include "FELCEMaterial444.h"
#include "FELCEMaterialPoint.h"





//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! volume element of LCE 
class FEBIOMECH_API FELCEDomain : public FESolidDomain, public FEElasticDomain
{

public:
    //! constructor
	FELCEDomain(FEModel* pfem);

public:
    
    //! assignment operator
	FELCEDomain& operator = (FELCEDomain& d);
    
	//! activate
	void Activate() override;

	//! initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

	//! initialize material points in the domain
	void InitMaterialPoints() override;

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

	// update LCE quantities i.e current director angle at integration point
	void UpdateLCEquantities(int iel, const FETimeInfo& tp);

	// Calculate and update internal variable current and previous
	void UpdateInternalVariable_Av(int iel, const FETimeInfo& tp);

	// Calculate and update internal variable current and previous
	void UpdateInternalVariable_Fv(int iel, const FETimeInfo& tp);

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
public:
	//!calcualte PK1 stress at element centriod or say average of all gp (for visualization)
	mat3d LCEAveragePK1Stress(FESolidElement el);

	//! Calculates Pressure of mix formulation of an element el
	double PressureMix(FESolidElement& el);

protected:
	//! calculates the element internal residual forces
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

	//! Calculates element stiffness element matrix
	void ElementStiffness(FESolidElement& el, matrix& ke);

	//! Get the current dof values in nodes of element el
	void CurrentDofValues(FESolidElement& el, matrix& u_elem, matrix& t_elem);

	//! Calcualtes deformation Gradients
	mat3d DeformationGradientF(matrix u_elem, matrix Bt);

	//!Calculates Shape function 
	void NshapeFunction(FESolidElement& el, int ngauss, matrix& N_t);
    
    //! Caluclates Bt matrix that interpolate Grad of theta
    void Btmatrix(FESolidElement el, int ngauss, matrix &Bt);
    
    //! Calcualtes Bu matrix that interpolate Grad of uflow
    void Bumatrix(matrix Bt, matrix &Bu);
    
    //! Calculate unassembled stiffness matrices
    void stiffnessmatrices(FESolidElement el, int ngauss,
		matrix N_t, matrix Bu, matrix Bt, double bulkmod,
		matrix dPdF_vgt, matrix dPdt_vgt, matrix dgdF_vgt,
		mat3d dPidGradt, matrix dgdt, matrix dJdF_vgt,
		matrix& K11, double& K22, double& K23, double& K32,
		matrix& K13, matrix& K31, matrix& K_ut, matrix& K_tu,
		matrix& K_tt);
        
    //! Calculate unassembled residual matrices
    void residualmatrices(FESolidElement el, int ngauss,
		matrix N_t, matrix Bu, matrix Bt, matrix PK1_vgt,
		matrix Pi0, matrix gt, matrix& ru, matrix& rt);

	void tens33_9(mat3d A, matrix& B);
	void tens3333_99(tens4d A, matrix& B);





	

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
	FELCEMaterial444* m_pMat; //pointer to the LCEMaterial

}; 


class FEStandardLCEDomain : public FELCEDomain
{
public:
	FEStandardLCEDomain(FEModel* fem);

private:
	std::string		m_elemType;

	DECLARE_FECORE_CLASS();
};