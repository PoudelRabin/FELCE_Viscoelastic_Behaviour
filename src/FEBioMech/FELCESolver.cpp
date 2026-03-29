/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FELCESolver.h"
#include "FERigidConnector.h"
#include "FESlidingElasticInterface.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FE3FieldElasticShellDomain.h"
#include "FEElasticEASShellDomain.h"
#include "FEElasticANSShellDomain.h"
#include "FEBodyForce.h"
#include "FEResidualVector.h"
#include "FEUncoupledMaterial.h"
#include "FEContactInterface.h"
#include "FESSIShellDomain.h"
#include <FECore/log.h>
#include <FECore/DOFS.h>
#include <FECore/sys.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/vector.h>
#include "FESolidLinearSystem.h"
#include "FEBioMech.h"
#include "FESolidAnalysis.h"
#include "FETrussMaterial.h"
#include "FELinearTrussDomain.h"
#include "FELCEMaterial443.h"
#include "FE3FieldLCEDomain.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FELCESolver, FENewtonSolver)
	BEGIN_PARAM_GROUP("Nonlinear solver");	// make sure this matches FENewtonSolver. 
		ADD_PARAMETER(m_Dtol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "dtol"        );
		ADD_PARAMETER(m_Etol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "etol");
		ADD_PARAMETER(m_Rtol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol");
		ADD_PARAMETER(m_rhoi      , "rhoi"        );
		ADD_PARAMETER(m_alpha     , "alpha"       );
		ADD_PARAMETER(m_beta      , "beta"        );
		ADD_PARAMETER(m_gamma     , "gamma"       );
		ADD_PARAMETER(m_logSolve  , "logSolve"    );
		ADD_PARAMETER(m_arcLength , "arc_length"  );
		ADD_PARAMETER(m_al_scale  , "arc_length_scale");
	END_PARAM_GROUP();
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FELCESolver Construction
//
//Added by Rabin on 25.05.2023 (m_dofDA(pfem) on the below line)
FELCESolver::FELCESolver(FEModel* pfem) : FENewtonSolver(pfem), m_rigidSolver(pfem),\
m_dofU(pfem), m_dofDA(pfem)
{
	// default values
	m_Rtol = 0.0;	// deactivate residual convergence for value = 0.0
	m_Dtol = 0.001; // was 0.001
	m_Etol = 0.01;
	m_Rmin = 1.0e-20;
	m_Rmax = 0;	// not used if zero

	m_niter = 0;
	m_nreq = 0;

	m_logSolve = false;

	// default Newmark parameters (trapezoidal rule)
    m_rhoi = -2;
    m_alpha = m_alphaf = 1.0;
    m_alpham = 1.0;
	m_beta  = 0.25;
	m_gamma = 0.5;

	// arc-length parameters
	m_arcLength = ARC_LENGTH_METHOD::NONE; // no arc-length
	m_al_scale = 0.0;
	m_al_lam = 0.0;
	m_al_inc = 0.0;
	m_al_ds = 0.0;
	//Added by Rabin on 29.05.2023 to initialize m_al_gamma
	m_al_gamma = 0.0;

	//Added by Rabin on 29.05.2023 (copied from fechem PI solve)
	// we'll need a non-symmetric stiffness matrix
	m_msymm = REAL_UNSYMMETRIC;
	
    // get the DOF indices
	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT)); 

        //Added by Rabin on 25.05.2023
        m_dofDA.AddVariable(FEBioMech::GetVariableName(FEBioMech::DIRECTOR_ANGLE));
	}
}

//-----------------------------------------------------------------------------
FELCESolver::~FELCESolver()
{

}


//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FELCESolver
//
bool FELCESolver::Init()
{
	// initialize base class
	if (FENewtonSolver::Init() == false) return false;
	
	FEModel& fem = *GetFEModel();
	
    if (m_rhoi == -1) {
        // Euler integration
        m_alpha = m_alphaf = m_alpham = 1.0;
        m_beta = pow(1 + m_alpham - m_alphaf,2)/4;
        m_gamma = 0.5 + m_alpham - m_alphaf;
    }
    else if ((m_rhoi >= 0) && (m_rhoi <= 1)) {
        // Generalized-alpha integration (2nd order system)
        m_alpha = m_alphaf = 1.0/(1+m_rhoi);
        m_alpham = (2-m_rhoi)/(1+m_rhoi);
        m_beta = pow(1 + m_alpham - m_alphaf,2)/4;
        m_gamma = 0.5 + m_alpham - m_alphaf;
    }
    else {
        // for any other value of rhoi, use the user-defined alpha, beta, gamma parameters
        m_alphaf = m_alpham = m_alpha;
    }
    
	// allocate vectors
//	m_Fn.assign(m_neq, 0);
	m_Fr.assign(m_neq, 0);
	m_Ui.assign(m_neq, 0);
	m_Ut.assign(m_neq, 0);
	m_Uip.assign(m_neq, 0);

	// we need to fill the total displacement vector m_Ut
	FEMesh& mesh = fem.GetMesh();
	gather(m_Ut, mesh, m_dofU[0]);
	gather(m_Ut, mesh, m_dofU[1]);
	gather(m_Ut, mesh, m_dofU[2]);
	//Added by Rabin 29.05.2023
	gather(m_Ut, mesh, m_dofDA[0]);

	//for (int nr = 0; nr < 20; ++nr) {
	//    for (int nc = 0; nc < 1; ++nc) {
	//        feLog("%f\t", m_Ut[nr]);
	//    } feLog("\n");
	//} feLog("\n*****************************************************************\n");
	
	
	if (m_arcLength > 0)
	{
		m_Fint.assign(m_neq, 0.0);
		m_Fext.assign(m_neq, 0.0);
	}
    
	// set the dynamic update flag only if we are running a dynamic analysis
	//bool b = (fem.GetCurrentStep()->m_nanalysis == FESolidAnalysis::DYNAMIC ? true : false);
	//for (int i = 0; i < mesh.Domains(); ++i)
	//{
	//	FEElasticSolidDomain* d = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(i));
    //    FEElasticShellDomain* s = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(i));
    //    FEElasticEASShellDomain* seas = dynamic_cast<FEElasticEASShellDomain*>(&mesh.Domain(i));
    //    FEElasticANSShellDomain* sans = dynamic_cast<FEElasticANSShellDomain*>(&mesh.Domain(i));
	//	if (d) d->SetDynamicUpdateFlag(b);
    //    if (s) s->SetDynamicUpdateFlag(b);
    //    if (seas) seas->SetDynamicUpdateFlag(b);
    //    if (sans) sans->SetDynamicUpdateFlag(b);
	//}
	//
	return true;
	
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FELCESolver::Serialize(DumpStream& ar)
{
	// Serialize parameters
	FENewtonSolver::Serialize(ar);
	
	ar & m_nrhs;
	ar & m_niter;
	ar & m_nref & m_ntotref;
	ar & m_naug;
	ar & m_nreq;

	ar & m_alphaf;
	ar & m_alpham;

	ar & m_Ut & m_Ui;

	ar & m_arcLength;
	ar & m_al_scale;
	ar & m_al_lam & m_al_inc;// &m_al_ds;

	if (ar.IsLoading())
	{
//		m_Fn.assign(m_neq, 0);
		m_Fr.assign(m_neq, 0);
//		m_Ui.assign(m_neq, 0);
	}

	// serialize rigid solver
	//m_rigidSolver.Serialize(ar);
	
}

//-----------------------------------------------------------------------------
bool FELCESolver::InitEquations()
{
	// First call the base class.
	// This will initialize all equation numbers, except the rigid body equation numbers
	if (FENewtonSolver::InitEquations() == false) return false;
	
	// store the number of equations we currently have
	m_nreq = m_neq;

	// Next, we assign equation numbers to the rigid body degrees of freedom
	//int neq = m_rigidSolver.InitEquations(m_neq);
	//if (neq == -1) return false; 
	//else m_neq = neq;

	// Next, we add any Lagrange Multipliers
	//FEModel& fem = *GetFEModel();
	//for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	//{
	//	FENLConstraint* lmc = fem.NonlinearConstraint(i);
	//	if (lmc->IsActive())
	//	{
	//		m_neq += lmc->InitEquations(m_neq);
	//	}
	//}
	//for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	//{
	//	FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
	//	if (spc->IsActive())
	//	{
	//		m_neq += spc->InitEquations(m_neq);
	//	}
	//}
	
	// All initialization is done
	return true;
}


//-----------------------------------------------------------------------------
bool FELCESolver::InitEquations2()
{
	// First call the base class.
	// This will initialize all equation numbers, except the rigid body equation numbers
	
	if (FENewtonSolver::InitEquations2() == false) return false;
	
	// store the number of equations we currently have
	m_nreq = m_neq;

	// Next, we assign equation numbers to the rigid body degrees of freedom
	int neq = m_rigidSolver.InitEquations(m_neq);
	if (neq == -1) return false;
	else m_neq = neq;
	
	// Next, we add any Lagrange Multipliers
	//FEModel& fem = *GetFEModel();
	//for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	//{
	//	FENLConstraint* lmc = fem.NonlinearConstraint(i);
	//	if (lmc->IsActive())
	//	{
	//		m_neq += lmc->InitEquations(m_neq);
	//	}
	//}
	//for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	//{
	//	FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
	//	if (spc->IsActive())
	//	{
	//		m_neq += spc->InitEquations(m_neq);
	//	}
	//}
	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FELCESolver::UpdateKinematics(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// update rigid bodies
	//m_rigidSolver.UpdateRigidBodies(m_Ui, ui);

	// total displacements
	vector<double> U(m_Ut.size());
	for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];

	// update flexible nodes
	// translational dofs
	scatter(U, mesh, m_dofU[0]);
	scatter(U, mesh, m_dofU[1]);
	scatter(U, mesh, m_dofU[2]);
	
	//Added by Rabin on 29.05.2023
	//Director Angle dof
	scatter(U, mesh, m_dofDA[0]);


	// make sure the boundary conditions are fullfilled
	int nbcs = fem.BoundaryConditions();
	for (int i = 0; i<nbcs; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.Update();
	}


	// Update the spatial nodal positions
	// Don't update rigid nodes since they are already updated
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
        if (node.m_rid == -1) {
			node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
        }
		//node.m_dt = node.m_d0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
       // - node.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
	}
	
}

//-----------------------------------------------------------------------------
//! Update DOF increments
void FELCESolver::UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap)
{
	FEModel& fem = *GetFEModel();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();
    
    
	// update flexible nodes
	int n;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
        
		// displacement dofs
		// current position = initial + total at prev conv step + total increment so far + current increment
		if ((n = node.m_ID[m_dofU[0]]) >= 0) Ui[n] += ui[n];
		if ((n = node.m_ID[m_dofU[1]]) >= 0) Ui[n] += ui[n];
		if ((n = node.m_ID[m_dofU[2]]) >= 0) Ui[n] += ui[n];

		//Added by Rabin on 29.05.2023
		//Director Angle dofs
		if ((n = node.m_ID[m_dofDA[0]]) >= 0)
		{
			Ui[n] += ui[n];
			//Ui[n] = fmod(Ui[n], 2.0 * acos(-1));
			//if (Ui[n] > acos(-1)) Ui[n] -= 2.0 * acos(-1);
			//if (Ui[n] < -acos(-1)) Ui[n] += 2.0 * acos(-1);

		}



		
	}
	
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FELCESolver::Update(vector<double>& ui)
{
    FEModel& fem = *GetFEModel();
    FETimeInfo& tp = fem.GetTime();
    tp.currentIteration = m_niter;
    

	// update kinematics
	UpdateKinematics(ui);

	// update model state
	UpdateModel();
}
//-----------------------------------------------------------------------------
//! Updates the current state of the model
//! NOTE: The ui vector also contains prescribed displacement increments. Also note that this
//!       only works for a limited set of FEBio features (no rigid bodies!).
void FELCESolver::Update2(const vector<double>& ui)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// total displacements
	vector<double> U(m_Ut.size());
	for (size_t i = 0; i < m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];

	// update free nodes
	scatter3(U, mesh, m_dofU[0], m_dofU[1], m_dofU[2]);

	// Update the spatial nodal positions
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vec3d du(0, 0, 0);
		int nx = -node.m_ID[m_dofU[0]] - 2; if (nx >= 0) du.x = ui[nx];
		int ny = -node.m_ID[m_dofU[1]] - 2; if (ny >= 0) du.y = ui[ny];
		int nz = -node.m_ID[m_dofU[2]] - 2; if (nz >= 0) du.z = ui[nz];

		if (node.m_rid == -1)
		{
			vec3d rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]) + du;
			node.m_rt = rt;
		}
	}

	// update model state
	UpdateModel();
}
//-----------------------------------------------------------------------------
bool FELCESolver::InitStep(double time)
{
	FEModel& fem = *GetFEModel();
	
	// set time integration parameters
	FETimeInfo& tp = fem.GetTime();
	tp.alpha = m_alpha;
	tp.beta = m_beta;
	tp.gamma = m_gamma;
	tp.alphaf = m_alphaf;
	tp.alpham = m_alpham;

    // evaluate load curve values at current (or intermediate) time
	double t = tp.currentTime;
//	double dt = tp.timeIncrement;
//	double ta = (t > 0) ? t - (1-m_alpha)*dt : m_alpha*dt;
//	return FESolver::InitStep(ta);
    return FESolver::InitStep(t);
	

}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration. 
void FELCESolver::PrepStep()
{
	FEModel& fem = *GetFEModel();
	
    FETimeInfo& tp = fem.GetTime();
	double dt = tp.timeIncrement;
	tp.augmentation = 0;
    
	// zero total displacements
	zero(m_Ui);

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEMesh& mesh = fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		ni.m_rp = ni.m_rt;
		//ni.m_vp = ni.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		ni.m_ap = ni.m_at;
        ni.m_dp = ni.m_dt;
		ni.UpdateValues();

        // initial guess at start of new time step
        // solid
        ni.m_at = ni.m_ap*(1-0.5/m_beta) - ni.m_vp/(m_beta*dt);
        vec3d vs = ni.m_vp + (ni.m_at*m_gamma + ni.m_ap*(1-m_gamma))*dt;
        //ni.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vs);
        
        // solid shell
        //vec3d aqp = ni.get_vec3d_prev(m_dofSA[0], m_dofSA[1], m_dofSA[2]);
        //vec3d vqp = ni.get_vec3d_prev(m_dofSV[0], m_dofSV[1], m_dofSV[2]);
        //vec3d aqt = aqp*(1-0.5/m_beta) - vqp/(m_beta*dt);
        //ni.set_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2], aqt);
        //vec3d vqt = vqp + (aqt*m_gamma + aqp*(1-m_gamma))*dt;
        //ni.set_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2], vqt);
    }

    // apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
//	vector<double> dummy(m_neq, 0.0);
//	zero(m_Fn);
//	FEResidualVector Fn(*GetFEModel(), m_Fn, dummy);
//	NodalLoads(Fn, tp);

	// apply boundary conditions
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_ui;
	zero(ui);
	int nbc = fem.BoundaryConditions();
	for (int i=0; i<nbc; ++i)
	{
		FEBoundaryCondition& dc = *fem.BoundaryCondition(i);
		if (dc.IsActive()) dc.PrepStep(ui);
	}


	// initialize rigid bodies
	//m_rigidSolver.PrepStep(tp, ui);

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
		if (dom.IsActive()) dom.PreSolveUpdate(tp);
	}

	// update model state
	UpdateModel();
	
}

//-----------------------------------------------------------------------------
// Performs the quasi-newton iterations.
bool FELCESolver::Quasin()
{
	vector<double> u0(m_neq);
	vector<double> Rold(m_neq);

	// convergence norms
	double	normR1;		// residual norm
	double	normE1;		// energy norm
	double	normU;		// displacement norm
	double	normu;		// displacement increment norm
	double	normRi;		// initial residual norm
	double	normEi;		// initial energy norm
	double	normEm;		// max energy norm
	double	normUi;		// initial displacement norm

	FEModel& fem = *GetFEModel();

	// Get the current step
	FEAnalysis* pstep = fem.GetCurrentStep();

	// set the time information
	FETimeInfo& tp = fem.GetTime();

	// initialize arc length stuff
	if (m_arcLength > 0)
	{
		m_al_inc = 0.0;

		// store the total increment from the last time step
		// we need it later to decide in which direction to proceed.
		m_Uip = m_Ui;
	}

	// prepare for the first iteration
	PrepStep();
	
	// Initialize the QN-method
	if (QNInit() == false) return false;
	
	// loop until converged or when max nr of reformations reached
	bool bconv = false;		// convergence flag
	do
	{
		feLog(" %d\n", m_niter+1);

		// assume we'll converge. 
		bconv = true;
		
		// solve the equations
		SolveEquations(m_ui, m_R0);
		
		// apply arc-length method
		if (m_arcLength > 0) DoArcLength();

		// do the line search
		double s = DoLineSearch();

		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_R0*m_R0);
			normEi = fabs(m_ui*m_R0);
			normUi = fabs(m_ui*m_ui);
			normEm = normEi;
		}
		
		// calculate actual displacement increment
		// NOTE: We don't apply the line search directly to m_ui since we need the unscaled search direction for the QN update below
		int neq = (int)m_Ui.size();
		vector<double> ui(m_ui);
		for (int i = 0; i<neq; ++i) ui[i] *= s;

		// update total displacements
		UpdateIncrements(m_Ui, ui, false);

		// calculate norms
		normR1 = m_R1*m_R1;
		normu  = ui*ui;
		normU  = m_Ui*m_Ui;
		normE1 = fabs(ui*m_R1);

		// check for nans
		if (ISNAN(normR1)) throw NANInResidualDetected();
		if (ISNAN(normu)) throw NANInSolutionDetected();

		// check residual norm
		if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;	

		// check displacement norm
		if ((m_Dtol > 0) && (normu  > (m_Dtol*m_Dtol)*normU )) bconv = false;

		// check energy norm
		if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;

		// check linestep size
		if ((m_lineSearch->m_LStol > 0) && (s < m_lineSearch->m_LSmin)) bconv = false;

		// check energy divergence
		if (m_bdivreform)
		{
			if (normE1 > normEm) bconv = false;
		}

		// Written by rabin
		// Find largest residual and and largest correction in current iteration and largest correcton of whole solution
		// for both displacement and angle
		FEMesh& mesh = fem.GetMesh();
		int n;
		vector<double> UiDis(neq,0.0), UiAng(neq,0.0), uiDis(neq,0.0), uiAng(neq,0.0);
		for (int i = 0; i < mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);

			// displacement dofs
			// current position = initial + total at prev conv step + total increment so far + current increment
			if ((n = node.m_ID[m_dofU[0]]) >= 0) {

				UiDis[3 * i] = fabs(m_Ui[n]); uiDis[3*i] = fabs(m_ui[n]);

			}
			if ((n = node.m_ID[m_dofU[1]]) >= 0) {
				UiDis[3 * i + 1] = fabs(m_Ui[n]); uiDis[3 * i + 1] = fabs(m_ui[n]);
			}
			if ((n = node.m_ID[m_dofU[2]]) >= 0) {
				UiDis[3 * i + 2] = fabs(m_Ui[n]); uiDis[3 * i + 2] = fabs(m_ui[n]);
			}
			if ((n = node.m_ID[m_dofDA[0]]) >= 0) {
				UiAng[i] = fabs(m_Ui[n]); uiAng[i] = fabs(m_ui[n]);
			}
		}

		double MaxUiDis = *max_element(UiDis.begin(), UiDis.end());
		double MaxuiDis = *max_element(uiDis.begin(), uiDis.end());
		double MaxUiAng = *max_element(UiAng.begin(), UiAng.end());
		double MaxuiAng = *max_element(uiAng.begin(), uiAng.end());

		// print convergence summary
		feLog(" Nonlinear solution status: time= %lg\n", tp.currentTime);
		// feLog("\tstiffness updates             = %d\n", m_qnstrategy->m_nups);
		// feLog("\tright hand side evaluations   = %d\n", m_nrhs);
		// feLog("\tstiffness matrix reformations = %d\n", m_nref);
		if (m_lineSearch->m_LStol > 0) feLog("\tstep from line search         = %lf\n", s);
		feLog("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		// feLog("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
		// feLog("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
		feLog("\t   displacement     %15le %15le %15le \n", normUi, normu ,(m_Dtol*m_Dtol)*normU );
		// feLog("\t   displacement     %15le %15le %15le \n", normUi, normu ,(m_Dtol*m_Dtol)*normU );
		// feLog("\t   Max dis cor Ui   %15le %15le %15le \n", 0, MaxUiDis, 0);
		// feLog("\t   Max dis inc ui   %15le %15le %15le \n", 0, MaxuiDis, 0);
		// feLog("\t   Max Ang cor Ui   %15le %15le %15le \n", 0, MaxUiAng, 0);
		// feLog("\t   Max Ang inc ui   %15le %15le %15le \n", 0, MaxuiAng, 0);
		// see if we may have a small residual
		if ((bconv == false) && (normR1 < m_Rmin))
		{
			// check for almost zero-residual on the first iteration
			// this might be an indication that there is no force on the system
			feLogWarning("No force acting on the system.");
			bconv = true;
		}

		// see if we have exceeded the max residual
		if ((bconv == false) && (m_Rmax > 0) && (normR1 >= m_Rmax))
		{
			// doesn't look like we're getting anywhere, so let's retry the time step
			throw MaxResidualError();
		}

		// check if we have converged. 
		// If not, calculate the BFGS update vectors
		if (bconv == false)
		{
			// do additional checks that may trigger a stiffness reformation
			if (s < m_lineSearch->m_LSmin)
			{
				// check for zero linestep size
				feLogWarning("Zero linestep size. Stiffness matrix will now be reformed");
				QNForceReform(true);
			}
			else if ((normE1 > normEm) && m_bdivreform)
			{
				// check for diverging
				feLogWarning("Problem is diverging. Stiffness matrix will now be reformed");
				normEm = normE1;
				normEi = normE1;
				normRi = normR1;
				QNForceReform(true);
			}

			// Do the QN update (This may also do a stiffness reformation if necessary)
			bool bret = QNUpdate();

			// something went wrong with the update, so we'll need to break
			if (bret == false) break;
		}
		else if (m_baugment)
		{
			// do the augmentations
			bconv = DoAugmentations();
		}
	
		// increase iteration number
		m_niter++;

		// do minor iterations callbacks
		fem.DoCallback(CB_MINOR_ITERS);
	}
	while (bconv == false);

	// if converged we update the total displacements
	if (bconv)
	{
        UpdateIncrements(m_Ut, m_Ui, true);

		// TODO: To zero or not to zero. That is the question!
		//       The arc-length method requires that we do NOT zero
		//       here, so that m_Uip gets initialized properly.
		//       However, the "jfnk-tangent test" requires that we zero this
		//       Otherwise the displacement increment can get counted twice in
		//       Update2 (from m_Ui and m_Ut since we just added it)
		if (m_arcLength == 0)
			zero(m_Ui);
	}

	return bconv;
}

//-----------------------------------------------------------------------------
// Exception that is thrown when the arc-length method has failed
class ArcLengthFailed : public FEException
{
public:
	ArcLengthFailed() : FEException("arc-length update has failed") {}
};

//-----------------------------------------------------------------------------
bool quadratic_solve_LCE(double a, double b, double c, double x[2])
{
	double D2 = b*b - 4.0*a*c;
	if (D2 < 0) return false;
	double D = sqrt(D2);
	if (b >= 0)
	{
		x[0] = (-b - D) / (2.0*a);
		x[1] = 2.0*c / (-b - D);
	}
	else
	{
		x[0] = 2.0*c / (-b + D);
		x[1] = (-b + D) / (2.0*a);
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Apply arc-length
void FELCESolver::DoArcLength()
{
	// auxiliary displacement
	vector<double> uF(m_neq, 0.0);

	// the arc-length scale factor
	double psi = m_al_scale;
	assert(psi == 0.0);

	const FETimeInfo& tp = GetFEModel()->GetTime();

	m_al_gamma = 0.0;
	if (m_arcLength == ARC_LENGTH_METHOD::CRISFIELD)
	{
		// solve for auxiliary displacement
		SolveEquations(uF, m_Fext);

		// if this is the first time step, we pick a special gamma
		if (m_niter == 0)
		{
			if (GetFEModel()->GetCurrentStep()->m_ntimesteps == 0)
			{
				// The first time we get here, we simply pick the arc-length step size
				// from the solution 
				m_al_gamma = tp.timeIncrement;
				m_al_ds = m_al_gamma*sqrt(uF*uF + psi*(m_Fext*m_Fext));
			}
			else
			{
				m_al_gamma = m_al_ds / sqrt(uF*uF + psi*(m_Fext*m_Fext));
				double uFdx = uF*m_Uip;
				if (uFdx < 0.0) m_al_gamma = -m_al_gamma;
			}
		}
		else
		{
			// The general case requires solving a quadratic equation
			// setup quadratic equation
			double Fe_norm2 = m_Fext*m_Fext;
			double a = uF*uF + (psi*psi)*Fe_norm2;
			double b = 2.0*(uF*(m_Ui + m_ui)) + 2 * m_al_inc*(psi*psi)*Fe_norm2;
			double c = m_ui*(m_Ui*2.0 + m_ui) + m_Ui*m_Ui - m_al_ds*m_al_ds + (psi*psi)*(m_al_inc*m_al_inc)*Fe_norm2;

			// solve quadratic equation
			double g[2];
			if (quadratic_solve_LCE(a, b, c, g) == false)
			{
				m_al_ds *= 0.5;
				throw ArcLengthFailed();
			}

			// two possible solution vectors
			vector<double> u1 = m_ui + uF*g[0];
			vector<double> u2 = m_ui + uF*g[1];

			// calculate two s-vectors
			vector<double> sk(2 * m_neq, 0.0), s1(2 * m_neq, 0.0), s2(2 * m_neq, 0.0);
			for (int i = 0; i < m_neq; ++i)
			{
				sk[i] = m_Ui[i];
				sk[i + m_neq] = m_al_inc*psi*m_Fext[i];

				s1[i] = m_Ui[i] + u1[i];
				s1[i + m_neq] = (m_al_inc + g[0])*psi*m_Fext[i];

				s2[i] = m_Ui[i] + u2[i];
				s2[i + m_neq] = (m_al_inc + g[1])*psi*m_Fext[i];
			}

			// calculate the norms of these vectors
			// NOTE: Should be Ds!!
			double norm_sk = 0, norm_s1 = 0, norm_s2 = 0;
			for (int i = 0; i < 2 * m_neq; ++i)
			{
				norm_sk += sk[i] * sk[i];
				norm_s1 += s1[i] * s1[i];
				norm_s2 += s2[i] * s2[i];
			}
			norm_sk = sqrt(norm_sk);
			norm_s1 = sqrt(norm_s1);
			norm_s2 = sqrt(norm_s2);

			// see which one produces the closest angle to the current path
			double c1 = (sk*s1) / (norm_sk*norm_s1);
			double c2 = (sk*s2) / (norm_sk*norm_s2);

			if (c1 >= c2)
			{
				m_al_gamma = g[0];
			}
			else
			{
				m_al_gamma = g[1];
			}
		}

		m_al_inc += m_al_gamma;
		m_al_lam += m_al_gamma;
		m_ui = m_ui + uF*m_al_gamma;
	}

	// evaluate the arc-length equation. 
	double sk2 = (m_Ui + m_ui)*(m_Ui + m_ui) + (psi*psi)*m_al_inc*m_al_inc*(m_Fext*m_Fext);
	double sk = sqrt(sk2);
	double serr = fabs((sk - m_al_ds) / m_al_ds);
	feLog("\tarc-length increment : %lg (%lg)\n", m_al_inc, m_al_gamma);
	feLog("\tarc-length factor    : %lg\n", m_al_lam);
	feLog("\tarc-length constraint: %lg (%lg, err = %lg)\n", sk, m_al_ds, serr);
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.

bool FELCESolver::StiffnessMatrix()
{
	FEModel& fem = *GetFEModel();

	const FETimeInfo& tp = fem.GetTime();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// setup the linear system
	FESolidLinearSystem LS(this, &m_rigidSolver, *m_pK, m_Fd, m_ui, (m_msymm == REAL_SYMMETRIC), m_alpha, m_nreq);

	// calculate the stiffness matrix for each domain
	for (int i=0; i<mesh.Domains(); ++i) 
	{
		if (mesh.Domain(i).IsActive()) 
		{
			FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
			dom.StiffnessMatrix(LS);
		}
	}

	// calculate the body force stiffness matrix for each non-rigid domain
	for (int j = 0; j<fem.ModelLoads(); ++j)
	{
		FEModelLoad* pml = fem.ModelLoad(j);
		if (pml->IsActive()) pml->StiffnessMatrix(LS);
	}

	// calculate stiffness matrices for surface loads
	// for arclength method we need to apply the scale factor to all the 
	// external forces stiffness matrix. 
	if (m_arcLength > 0) LS.StiffnessAssemblyScaleFactor(m_al_lam);
/*	int nsl = fem.SurfaceLoads();
	for (int i = 0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = fem.SurfaceLoad(i);
		if (psl->IsActive())
		{
			psl->StiffnessMatrix(LS, tp);
		}
	}*/
	if (m_arcLength > 0) LS.StiffnessAssemblyScaleFactor(1.0);

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FELCESolver::Residual(vector<double>& R)
{
	// get the time information
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// zero nodal reaction forces
	zero(m_Fr);

	// setup the global vector
	zero(R);
	FEResidualVector RHS(fem, R, m_Fr);

	// calculate the internal (stress) forces
	InternalForces(RHS);
	
	// calculate nodal reaction forces
	for (int i = 0; i < m_neq; ++i) m_Fr[i] -= R[i];

	// extract the internal forces
	// (only when we really need it, below)
	if (m_logSolve && fem.GetCurrentStep()->m_ntimesteps > 0)
	{
		m_Fint = R;
	}

	if (m_arcLength > 0)
	{
		// Note the negative sign. This is because during residual assembly
		// a negative sign is applied to the internal force vector. 
		// The model loads assume the residual is Fe - Fi (i.e. -R)
		m_Fint = -R;
	}
	
	// calculate external forces
	ExternalForces(RHS);

	// For arc-length we need the external loads
	if (m_arcLength > 0)
	{
		// extract the external force
		m_Fext = R + m_Fint;

		// we need to apply the arc-length factor to the external loads
        for (int i=0; i<R.size();++i) R[i] = m_Fext[i]* m_al_lam - m_Fint[i];
	}

	// apply the residual transformation
	// NOTE: This is an implementation of Ankush Aggarwal method to accelerate the Newton convergence
	if (m_logSolve && fem.GetCurrentStep()->m_ntimesteps > 0)
	{
		double TOL = 1.e-8;
		bool logused = false;
		vector<double> RHSlog;
		RHSlog.resize(R.size());
		for (int i = 0; i<m_Fint.size(); ++i)
		{
			if (fabs(RHS[i] - m_Fint[i])>TOL && fabs(m_Fint[i])>TOL && (m_Fint[i] - RHS[i]) / m_Fint[i]>0)
			{
				RHSlog[i] = -m_Fint[i] * log((m_Fint[i] - RHS[i]) / m_Fint[i]);
				logused = true;
			}
			else
			{
				RHSlog[i] = RHS[i];
			}
		}
		for (int i = 0; i<m_Fint.size(); ++i) R[i] = RHSlog[i];
		if (logused)
			feLog("Log method used\n");
	}

	// increase RHS counter
	m_nrhs++;
	return true;
}

//-----------------------------------------------------------------------------
//! Internal forces
void FELCESolver::InternalForces(FEGlobalVector& R)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
		if (edom) edom->InternalForces(R);
	}
	
}

//-----------------------------------------------------------------------------
//! external forces
void FELCESolver::ExternalForces(FEGlobalVector& RHS)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	FEMesh& mesh = fem.GetMesh();

	// apply loads
	for (int j = 0; j<fem.ModelLoads(); ++j)
	{
		FEModelLoad* pml = fem.ModelLoad(j);
		if (pml->IsActive()) pml->LoadVector(RHS);
	}

	// calculate inertial forces for dynamic problems
	if (fem.GetCurrentStep()->m_nanalysis == FESolidAnalysis::DYNAMIC)
	{
		// allocate F
		vector<double> F;

		// calculate the inertial forces for all elastic domains
		for (int nd = 0; nd < mesh.Domains(); ++nd)
		{
			FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(nd));
			if (edom) edom->InertialForces(RHS, F);
		}

		// update rigid bodies
		m_rigidSolver.InertialForces(RHS, tp);
	}

	// calculate forces due to surface loads
/*	int nsl = fem.SurfaceLoads();
	for (int i = 0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->LoadVector(RHS, tp);
	}
*/



	// forces due to point constraints
	//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->Residual(this, R);

	// set the nodal reaction forces
	// TODO: Is this a good place to do this?
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		node.set_load(m_dofU[0], 0);
		node.set_load(m_dofU[1], 0);
		node.set_load(m_dofU[2], 0);

		int n;
		if ((n = node.m_ID[m_dofU[0]]) >= 0) node.set_load(m_dofU[0], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[0]] - 2) >= 0) node.set_load(m_dofU[0], -m_Fr[n]);

		if ((n = node.m_ID[m_dofU[1]]) >= 0) node.set_load(m_dofU[1], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[1]] - 2) >= 0) node.set_load(m_dofU[1], -m_Fr[n]);

		if ((n = node.m_ID[m_dofU[2]]) >= 0) node.set_load(m_dofU[2], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[2]] - 2) >= 0) node.set_load(m_dofU[2], -m_Fr[n]);

		// add nodal loads
		double s = (m_arcLength>0 ? m_al_lam : 1.0);
//		if ((n = node.m_ID[m_dofU[0]]) >= 0) node.set_load(m_dofU[0], -m_Fn[n]*s);
//		if ((n = node.m_ID[m_dofU[1]]) >= 0) node.set_load(m_dofU[1], -m_Fn[n]*s);
//		if ((n = node.m_ID[m_dofU[2]]) >= 0) node.set_load(m_dofU[2], -m_Fn[n]*s);
	}
}

