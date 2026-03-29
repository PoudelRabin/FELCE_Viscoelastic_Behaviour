#include "stdafx.h"
#include "FELCEDomain.h"
#include "FELCEMaterial444.h"
#include "FEElasticMaterial.h"
#include "FEBodyForce.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioMech.h"
#include <FECore/FELinearSystem.h>
#include "FELCEInitialAngle.h"



//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FELCEDomain::FELCEDomain(FEModel* pfem) : FESolidDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem), m_dofDA(pfem), m_dof(pfem)
{
	m_pMat = 0;
    m_alphaf = m_beta = 1;
    m_alpham = 2;
	m_update_dynamic = true; // default for backward compatibility
    
	m_secant_stress = false;
	m_secant_tangent = false;
    
	// TODO: Can this be done in Init, since  there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
        m_dofDA.AddVariable(FEBioMech::GetVariableName(FEBioMech::DIRECTOR_ANGLE));
	}
}
//-----------------------------------------------------------------------------
// get the total dof list
const FEDofList& FELCEDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FELCEDomain& FELCEDomain::operator = (FELCEDomain& d)
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
//! Set flag for update for dynamic quantities
void FELCEDomain::SetDynamicUpdateFlag(bool b)
{
	m_update_dynamic = b;
}

//-----------------------------------------------------------------------------
//! Assign material
void FELCEDomain::SetMaterial(FEMaterial* pmat)
{
    
	FEDomain::SetMaterial(pmat);
	if (pmat)
	{
		m_pMat = dynamic_cast<FELCEMaterial444*>(pmat);
		assert(m_pMat);
	}
	else m_pMat = 0;
}

//-----------------------------------------------------------------------------
void FELCEDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.set_active(m_dofU[0]);
				node.set_active(m_dofU[1]);
				node.set_active(m_dofU[2]);
                node.set_active(m_dofDA[0]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! serialization
void FELCEDomain::Serialize(DumpStream& ar)
{
	//erialize the base class, which instantiates the elements
	FESolidDomain::Serialize(ar);
	if (ar.IsShallow()) return;

	// serialize class variables
	ar & m_alphaf;
	ar & m_alpham;
	ar & m_beta;
	ar & m_update_dynamic;
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FELCEDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    m_alphaf = timeInfo.alphaf;
    m_alpham = timeInfo.alpham;
    m_beta = timeInfo.beta;
    bool ifError = false;
    int NE = Elements();
    //#pragma omp parallel for shared(NE, ifError)
	for (int i=0; i<NE; ++i)
	{
        //try
        {
            FESolidElement& el = m_Elem[i];
            if (el.isActive())
            {
                
                
                int n = el.GaussPoints();

                for (int j = 0; j < n; ++j)
                {

                    FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
                    pt.m_Wp = pt.m_Wt;

                    mp.Update(timeInfo);
                }
            }
        }
        //catch (FEException e)
        //{
        //    #pragma omp critical
        //    {
        //        ifError = true;
        //        //feLogError(e.what());
        //    }
        //}
    
	}

   //if (ifError) throw FEException("Iternal Variable Calculation Error");
}


//-----------------------------------------------------------------------------
void FELCEDomain::InitMaterialPoints()
{
    FEMesh& m = *GetMesh();
    int Tnodes = m.Nodes();

    vector< vector<double> > initangle(Tnodes, vector<double>(8));
    
    for (int j = 0; j < (int)m_Elem.size(); ++j)
    {
        // get the solid element
        FESolidElement& el = m_Elem[j];

        // get the number of nodes
        int neln = el.Nodes();
        // get initial values of fluid pressure and solute concentrations
        for (int i = 0; i < neln; ++i)
        {
            FENode& ni = m.Node(el.m_node[i]);
            for (int itnodes = 0; itnodes < Tnodes; ++itnodes)
            {
                initangle[itnodes][i] = ni.get(m_dofDA[0]);
            }
        }

        // get the number of integration points
        int nint = el.GaussPoints();

        //initialize shape function
        matrix N_t(1, Tnodes);
        

        // loop over the integration points
        for (int n = 0; n < nint; ++n)
        {
            
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEElasticMaterialPoint& pda = *(mp.ExtractData<FEElasticMaterialPoint>());
            FELCEMaterialPoint& Lpt = *(mp.ExtractData<FELCEMaterialPoint>());
                        

            //Get shape function           
            NshapeFunction(el, n, N_t);
            
            //director angle theta0 at integraton point
            Lpt.m_theta0 = 0.0;  
            
            for (int itnodes = 0; itnodes < Tnodes; ++itnodes)
            { 
                pda.m_DA0 = el.Evaluate(initangle[itnodes], n);
                Lpt.m_theta0 = Lpt.m_theta0 +  N_t(0, itnodes)* pda.m_DA0;               
            }
            Lpt.m_theta = Lpt.m_theta0;
        }
    } 
}

//-----------------------------------------------------------------------------
void FELCEDomain::InternalForces(FEGlobalVector& R)
{
	int NE = Elements();
	#pragma omp parallel for shared (NE)
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
        
		if (el.isActive()) {
			// element force vector
			vector<double> fe;
			vector<int> lm;

			// get the element force vector and initialize it to zero
			int ndof = 4 * el.Nodes();
			fe.assign(ndof, 0);
            
			// calculate internal force vector
			ElementInternalForce(el, fe);
            
			// get the element's LM vector
			UnpackLM(el, lm);
   
			// assemble element 'fe'-vector into global R vector
			R.Assemble(el.m_node, lm, fe);
            
		}
	}
}

//-----------------------------------------------------------------------------
void FELCEDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// repeat over all solid elements
	int NE = Elements();
    
	#pragma omp parallel for shared (NE)
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		if (el.isActive()) {

			// get the element's LM vector
			vector<int> lm;
			UnpackLM(el, lm);
   
			// element stiffness matrix
			FEElementMatrix ke(el, lm);

			// create the element's stiffness matrix
			int ndof = 4 * el.Nodes();
			ke.resize(ndof, ndof);
			ke.zero();
            
			// calculate element stiffness
			ElementStiffness(el, ke);
           
			LS.Assemble(ke);         
        }
        
	}
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FELCEDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*4); //Changed from 6 to 4 by Rabin
	for (int i=0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		vector<int>& id = node.m_ID;
//todo should i be accomodating for the new varaible director angle?
		// first the displacement dofs
		lm[4*i  ] = id[m_dofU[0]];
		lm[4*i+1] = id[m_dofU[1]];
		lm[4*i+2] = id[m_dofU[2]];

        //Added by Rabin 31.05.2023
        lm[4 * i + 3] = id[m_dofDA[0]];
	}
    
    // substitute interface dofs for solid-shell interfaces
	//FESolidElement& sel = static_cast<FESolidElement&>(el);
	//for (int i = 0; i<sel.m_bitfc.size(); ++i)
    //{
    //    if (sel.m_bitfc[i]) {
    //        FENode& node = m_pMesh->Node(el.m_node[i]);
    //        vector<int>& id = node.m_ID;
    //        
    //        // first the displacement dofs
    //        lm[3*i  ] = id[m_dofSU[0]];
    //        lm[3*i+1] = id[m_dofSU[1]];
    //        lm[3*i+2] = id[m_dofSU[2]];
    //    }
    //}
}

//-----------------------------------------------------------------------------
void FELCEDomain::Update(const FETimeInfo& tp)
{

	bool berr = false;
	int NE = Elements();
	#pragma omp parallel for shared(NE, berr)
	for (int i=0; i<NE; ++i)
	{
		try
		{
			FESolidElement& el = Element(i);
			if (el.isActive())
			{
                UpdateElementStress(i, tp);
            
                //Added by rabin
                UpdateLCEquantities(i, tp);

                //update the InternalVAriable
                // should this be done after update of LCE quantities
                UpdateInternalVariable_Av(i, tp);
                UpdateInternalVariable_Fv(i, tp);
			}
        } 
		//catch (NegativeJacobian e)
        catch (FEException e)
		{
			#pragma omp critical
			{
				// reset the logfile mode
				berr = true;
				//if (e.DoOutput()) feLogError(e.what());
			}
		}
	}
	if (berr) throw NegativeJacobianDetected();
}

//----------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
//! \todo Remove the remodeling solid stuff
void FELCEDomain::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double dt =tp.timeIncrement;

	// get the solid element
	FESolidElement& el = m_Elem[iel];

	// get the number of integration points
	int nint = el.GaussPoints();

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
    const int NELN = 8;
    vec3d r[NELN], v[NELN], a[NELN];
	GetCurrentNodalCoordinates(el, r, m_alphaf);

	// update dynamic quantities
	if (m_update_dynamic)
	{
		for (int j = 0; j<neln; ++j)
		{
			FENode& node = m_pMesh->Node(el.m_node[j]);
			//v[j] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2])*m_alphaf + node.m_vp*(1 - m_alphaf);
			a[j] = node.m_at*m_alpham + node.m_ap*(1 - m_alpham);
		}
	}

	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// material point coordinates
		mp.m_rt = el.Evaluate(r, n);

		// get the deformation gradient and determinant at intermediate time
        double Jt;
        mat3d Ft, Fp;
        Jt = defgrad(el, Ft, n);
        defgradp(el, Fp, n);
        
		if (m_alphaf == 1.0)
		{
			pt.m_F = Ft;
            pt.m_J = Jt;
		}
		else
		{
			pt.m_F = Ft*m_alphaf + Fp*(1-m_alphaf);
            pt.m_J = pt.m_F.det();
		}

        mat3d Fi = pt.m_F.inverse();
        pt.m_L = (Ft - Fp)*Fi / dt;
		if (m_update_dynamic)
		{
			pt.m_v = el.Evaluate(v, n);
			pt.m_a = el.Evaluate(a, n);
		}
        

        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, tp);
		// calculate the stress at this material point
//		pt.m_s = m_pMat->Stress(mp);
		pt.m_s = (m_secant_stress ? m_pMat->SecantStress(mp) : m_pMat->Stress(mp));
        
        // adjust stress for strain energy conservation
        if (m_alphaf == 0.5) 
		{
			FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(m_pMat);

			// evaluate strain energy at current time
			mat3d Ftmp = pt.m_F;
			double Jtmp = pt.m_J;
			pt.m_F = Ft;
			pt.m_J = Jt;
			pt.m_Wt = pme->StrainEnergyDensity(mp);
			pt.m_F = Ftmp;
			pt.m_J = Jtmp;

            mat3ds D = pt.RateOfDeformation();
            double D2 = D.dotdot(D);
            if (D2 > 0)
                pt.m_s += D*(((pt.m_Wt-pt.m_Wp)/(dt*pt.m_J) - pt.m_s.dotdot(D))/D2);
        }
    } 
}

//-----------------------------------------------------------------------------
void FELCEDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    // define some parameters that will be passed to lambda
    FEBodyForce* bodyForce = &BF;

    // TODO: a remaining issue here is that dofU does not consider the shell displacement
    // dofs for interface nodes (see UnpackLM). Is that an issue?

    // evaluate the residual contribution
    LoadVector(R, m_dofU, [=](FEMaterialPoint& mp, int node_a, std::vector<double>& fa) {

        // evaluate density
        double density = m_pMat->Density(mp);

        // get the force
        vec3d f = bodyForce->force(mp);

        // get element shape functions
        double* H = mp.m_shape;

        // get the initial Jacobian
        double J0 = mp.m_J0;

        // set integrand
        //fa[0] = -H[node_a] * density * f.x * J0;
        //fa[1] = -H[node_a] * density * f.y * J0;
        //fa[2] = -H[node_a] * density * f.z * J0;
        fa[0] =0.0; 
        fa[1] =0.0;
        fa[2] =0.0;
        });

    
}

//-----------------------------------------------------------------------------
void FELCEDomain::MassMatrix(FELinearSystem& LS, double scale)
{
    // TODO: a remaining issue here is that dofU does not consider the shell displacement
    // dofs for interface nodes (see UnpackLM). Is that an issue?

    // evaluate body force stiffness
    LoadStiffness(LS, m_dofU, m_dofU, [=](FEMaterialPoint& mp, int node_a, int node_b, matrix& Kab) {

        // density
        double density = m_pMat->Density(mp);

        // shape functions
        double* H = mp.m_shape;

        // Jacobian
        double J0 = mp.m_J0;

        // mass
        double kab = scale * density * H[node_a] * H[node_b] * J0;
        Kab.zero();
        //Kab[0][0] = kab;
        //Kab[1][1] = kab;
        //Kab[2][2] = kab;
        });
}

//-----------------------------------------------------------------------------
void FELCEDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
    // define some parameters that will be passed to lambda
    FESolidMaterial* mat = m_pMat;
    FEBodyForce* bodyForce = &bf;

    // TODO: a remaining issue here is that dofU does not consider the shell displacement
    // dofs for interface nodes (see UnpackLM). Is that an issue?

    // evaluate body force stiffness
    LoadStiffness(LS, m_dofU, m_dofU, [=](FEMaterialPoint& mp, int node_a, int node_b, matrix& Kab) {

        // loop over integration points
        double detJ = mp.m_J0 * m_alphaf;

        // density
        double dens_n = mat->Density(mp);
        
        // get the stiffness
        mat3d K = bodyForce->stiffness(mp);

        // shape functions
        double* H = mp.m_shape;

        // put it together
        //Kab = K * (-H[node_a] * H[node_b] * dens_n * detJ);
        Kab.zero();
        });
}

//-----------------------------------------------------------------------------
// Calculate inertial forces \todo Why is F no longer needed?
void FELCEDomain::InertialForces(FEGlobalVector& R, vector<double>& F)
{
    int NE = Elements();
#pragma omp parallel for shared(R, F)
    for (int i = 0; i < NE; ++i)
    {
        // get the element
        FESolidElement& el = m_Elem[i];

        if (el.isActive()) {
            // element force vector
            vector<double> fe;
            vector<int> lm;

            // get the element force vector and initialize it to zero
            int ndof = 4 * el.Nodes();
            fe.assign(ndof, 0);

            // calculate internal force vector
            //ElementInertialForce(el, fe);
            
            // get the element's LM vector
            UnpackLM(el, lm);

            // assemble element 'fe'-vector into global R vector
            R.Assemble(el.m_node, lm, fe);
        }
    }
}
//-----------------------------------------------------------------------------
//                Below mostly written by Rabin Poudel @ Cardiff University
//-----------------------------------------------------------------------------
// ! Update the current director angle value at integration point m_theta
void FELCEDomain::UpdateLCEquantities(int iel, const FETimeInfo& tp)
{
    
    // get the solid element
    FESolidElement& el = m_Elem[iel];

    // get the number of integration points
    int nint = el.GaussPoints();
    
    // number of nodes
    int neln = el.Nodes();
    
    //Get current DOF values
    matrix u_elem(3, neln);
    matrix t_elem(1, neln);
    CurrentDofValues(el, u_elem, t_elem);
    //Initialize shape function 
    matrix N_t(1, neln);

    // loop over the integration points
    for (int n = 0; n < nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FELCEMaterialPoint& Lpt = *(mp.ExtractData<FELCEMaterialPoint>());
        //Get shape function
        NshapeFunction(el, n, N_t);
        //director angle theta0 at integraton point
        matrix t_temp(1, 1);
        t_temp = N_t * t_elem.transpose();
        Lpt.m_theta = t_temp(0, 0);       
    }

}
//----------------------------------------------------------------------------------------------
// ! Calculate and update internal variable Av
void FELCEDomain::UpdateInternalVariable_Av(int iel, const FETimeInfo& tp)
{
    
    double dt = tp.timeIncrement;

    
 
    // get the solid element
    FESolidElement& el = m_Elem[iel];

    // get the number of integration points
    int nint = el.GaussPoints();

    // number of nodes
    int neln = el.Nodes();
    

    // loop over the integration points
    for (int n = 0; n < nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& Ept = *(mp.ExtractData<FEElasticMaterialPoint>());
        FELCEMaterialPoint& Lpt = *(mp.ExtractData<FELCEMaterialPoint>());

        //Get material Property
        //Get material point data as mp of the gaussian point n
        FELCEMaterial444& mat = dynamic_cast<FELCEMaterial444&>(*m_pMat);
        double mu_neq =  mat.m_mu_neq_nc(mp);
        double eta_u =  mat.m_eta_u(mp);
        double a = mat.m_a(mp);





        //Get all deformation gradient
        mat3d F, G0, G, Ginv, Avp, Av, Ae, Av_inv, A_DG;
        mat3d eye;
        eye.unit();
              
        F = Ept.m_F;
        G0 = Lpt.InitialSpontaniousDG(a);
        G = Lpt.SpontaniousDG(a);
        Ginv = G.inverse();
        Avp = Lpt.m_Avp;
        A_DG = (Ginv * F) * G0;
        double detF = F.det();


        // linearly varying a and eta_u based on strain
    // double lval[3];
    // vec3d rvec[3];
    // mat3ds Csym = (F.transpose() * F).sym();
    // Csym.eigen(lval, rvec);
    // double lam_x = pow(lval[0],0.5);
    // double final_total_strain = 3.5;
    // double a0 = 4.0;
    // double ax = 3.0;
    // a = a0 + (ax-a0)*(lam_x-1.0)/(final_total_strain-1);
    //// below not needed
    //double eta_t0 = 800;
   // double eta_tx = 100.0;
   // eta_t = eta_t0 + (eta_tx-eta_t0)*(lam_x-1.0)/(final_total_strain-1);
       
        // Calculate new Av for this time step
        
        //Initial Guess for Av
        Av = Avp;
        
        //Newton Raphson Iteration
        double restol = 1.0e-12;
        int maxiter = 500; 
        for (int iter = 0; iter < maxiter; ++iter)
        {
             Av_inv = Av.inverse();
             Ae = A_DG * Av_inv;
            
             mat3d dWvdAe = Ae * mu_neq * pow(detF, -2.0 / 3.0);
             mat3d tempmat1 = (Ae.transpose() * dWvdAe) * Av;
             
             //changing to voigt notation
             matrix tempvgt1(9, 1), Avp_vgt(9, 1), Av_vgt(9, 1);            
             tens33_9(tempmat1, tempvgt1);             
             tens33_9(Avp, Avp_vgt);
             tens33_9(Av, Av_vgt);

             //calculating residual and its norm
             matrix residual(9, 1);
             residual = Av_vgt - (Avp_vgt + (tempvgt1 * dt * (2.0 / eta_u)));
             matrix tempressq(1, 1);
             tempressq = residual.transpose() * residual;
             double resnormsquared = tempressq(0, 0);
             double resnorm = fabs(sqrt(resnormsquared));               
             
             //terminate NR iteration if norm is less than tol
             if (iter>0 && resnorm <= restol) {  
                 
                 break;
             }
             
             //calculating the tangents
             mat3d eye;
             eye.unit();
             tens4d TempA = dyad2(eye, eye.transpose());
             tempmat1  = (Ae.transpose() * dWvdAe) * Av;
             tens4d TempB = dyad3(-Av_inv.transpose(), tempmat1.transpose());
             tens4d tangent = TempA - (TempB * dt * (2.0 / eta_u));
             
             //fourth order tens tangent to 9by9 matrix
             matrix tangent_vgt(9, 9);
             double tangent_vgt_temp[9][9];
             tangent.extract(tangent_vgt_temp);
             matrix tangent_vgt_temp2(9, 9);
             for (int i = 0; i < 9; ++i) {
                 for (int j = 0; j < 9; ++j) {
//
                     tangent_vgt_temp2(i, j) = tangent_vgt_temp[i][j];
                 }
             }
             //tangent_vgt = tangent_vgt_temp2.transpose();
             tens3333_99(tangent, tangent_vgt);
             //tangent_vgt = tangent_vgt.transpose();

             //call solve for solving the equation
             vector<double> dAv_vec(9);
             vector<double> negresidual_vec(9);
             
             for (int i = 0; i < 9; ++i) {
                negresidual_vec[i] = -residual(i, 0);
             }
             tangent_vgt.solve(dAv_vec, negresidual_vec);
             //matrix dAv(9,1);
             //tangent_vgt = tangent_vgt.transpose();
             //dAv = tangent_vgt.inverse() * residual * -1.0;
             
        
             
             for (int i = 0; i < 9; ++i) {
                // dAv_vec[i] = dAv(i,0);
                 Av_vgt(i, 0) = Av_vgt(i, 0) + dAv_vec[i];
             }
             
             //change voigt to matrix for easy calculation of Ae

             int indx[9] = { 0,3,8,6,1,4,5,7,2 };
             int u = 0;
             for (int i = 0; i < 3; i++) {
                 for (int j = 0; j < 3; j++) {
                     int v = indx[u];
                     Av(i, j) = Av_vgt(v, 0);
                     u = u + 1; 
                 }
             }
             
             if (iter == (maxiter-1)) {
                 throw FEException("No Convergence of Internal Variable after Max iteration in Av");
             }  
        } 

        if (Av.det() <= 0.0) {
            throw FEException("Negetive Jacobian of Internal Variable Av detected");
        }
        
        Lpt.m_Av = Av;        

    }
    
}
// -------------------------------------------------------------------------------------------------- 
// //----------------------------------------------------------------------------------------------
// ! Calculate and update internal variable Fv
void FELCEDomain::UpdateInternalVariable_Fv(int iel, const FETimeInfo& tp)
{
    //return;
    double dt = tp.timeIncrement;

    // get the solid element
    FESolidElement& el = m_Elem[iel];

    // get the number of integration points
    int nint = el.GaussPoints();

    // number of nodes
    int neln = el.Nodes();


    // loop over the integration points
    for (int n = 0; n < nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& Ept = *(mp.ExtractData<FEElasticMaterialPoint>());
        FELCEMaterialPoint& Lpt = *(mp.ExtractData<FELCEMaterialPoint>());

        //Get material Property
        //Get material point data as mp of the gaussian point n
        FELCEMaterial444& mat = dynamic_cast<FELCEMaterial444&>(*m_pMat);
        double mu_neq_nh = mat.m_mu_neq_nh(mp);
        double eta_t = mat.m_eta_t(mp);
        double a = mat.m_a(mp);
        //Get all deformation gradient and previous viscous angle
        mat3d F, Fvp, Fv, Febar, Fv_inv, Fvdot; // Fv is at integration point
        mat3d eye;
        eye.unit();

        F = Ept.m_F;
		Fvp = Lpt.m_Fvp;
        double detF = F.det();
		mat3d Fbar = F * pow(detF, -1.0 / 3.0);







        // linearly varying a and eta_u based on strain
    // double lval[3];
    // vec3d rvec[3];
    // mat3ds Csym = (F.transpose() * F).sym();
    // Csym.eigen(lval, rvec);
    // double lam_x = pow(lval[0],0.5);
    // double final_total_strain = 3.5;
    // double a0 = 4.0;
    // double ax = 3.0;
    // a = a0 + (ax-a0)*(lam_x-1.0)/(final_total_strain-1);
    // double eta_t0 = 800;
    // double eta_tx = 100.0;
    // eta_t = eta_t0 + (eta_tx-eta_t0)*(lam_x-1.0)/(final_total_strain-1);

        // Calculate new Fv for this time step

        //Initial Guess for Av
        Fv = Fvp;
        

        //Newton Raphson Iteration
        double restol = 1.0e-12;
        int maxiter = 300;
        for (int iter = 0; iter < maxiter; ++iter)
        {

            Fv_inv = Fv.inverse();
            Febar = Fbar * Fv_inv;
            Fvdot = (Fv-Fvp)/dt;

            // Residual Calculation
            mat3d tempmat1 = Febar.transpose() * Febar;
            mat3d tempmat2 = Fvdot * Fv_inv;

            mat3d resmat = tempmat1*mu_neq_nh - tempmat2*0.50*eta_t;

            //changing to voigt notation
            matrix residual(9, 1), Fv_vgt(9, 1);
            tens33_9(resmat, residual);
            tens33_9(Fv, Fv_vgt);

            // residual norm
            matrix tempressq(1, 1);
            tempressq = residual.transpose() * residual;
            double resnormsquared = tempressq(0, 0);
            double resnorm = fabs(sqrt(resnormsquared));

            //terminate NR iteration if norm is less than tol
            if (iter > 0 && resnorm <= restol) {

                break;
            }

            //calculating the tangents
            tempmat1 = Febar.transpose() * Febar;
            tempmat2 = Fvdot * Fv_inv;
            tens4d TempA = dyad2(tempmat1, Fv_inv.transpose());
            tens4d TempB = dyad3(Fv_inv.transpose(), tempmat1.transpose());
            tens4d TempC = dyad2(tempmat2, Fv_inv.transpose());
            tens4d TempD = dyad2(eye, Fv_inv.transpose());

            tens4d tangent = (TempA + TempB) * (-mu_neq_nh) +
                             (TempC - TempD/dt) * 0.50 * eta_t;

            //fourth order tens tangent to 9by9 matrix
            matrix tangent_vgt(9, 9);
            tens3333_99(tangent, tangent_vgt);

            //call solve for solving the equation
            vector<double> dFv_vec(9);
            vector<double> negresidual_vec(9);

            for (int i = 0; i < 9; ++i) {
                negresidual_vec[i] = -residual(i, 0);
            }
            tangent_vgt.solve(dFv_vec, negresidual_vec);

            for (int i = 0; i < 9; ++i) {
                Fv_vgt(i, 0) = Fv_vgt(i, 0) + dFv_vec[i];
            }

            //change voigt to matrix for easy calculation of Febar

            int indx[9] = { 0,3,8,6,1,4,5,7,2 };
            int u = 0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    int v = indx[u];
                    Fv(i, j) = Fv_vgt(v, 0);
                    u = u + 1;
                }
            }

            if (iter == (maxiter - 1)) {
                throw FEException("No Convergence of Internal Variable after Max iteration in Fv");
            }


        }

        if (Fv.det() <= 0.0) {
            throw FEException("Negetive Jacobian of Internal Variable Fv detected");
        }

        Lpt.m_Fv = Fv;

    }

}
// -------------------------------------------------------------------------------------------------- 
//! calculates the internal equivalent nodal forces for solid elements
void FELCEDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
      
    // Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 4*neln;
    double*	gw = el.GaussWeights();
    
    //Initialize Relem
    matrix Relem(32,1);
    Relem.zero();

    //Get Initial Angles
    //todo get initial angles from the initial condition not specify here
    matrix t0_elem(1, 8);  
    //for (int q = 0; q < 8; ++q) { t0_elem(0, q) = 1.5707963267948966;}
    //for (int q = 0; q < 8; ++q) { t0_elem(0, q) = acos(-1.0)/2.0; }


    for (int j = 0; j < 8; ++j) 
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(j);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        t0_elem(0, j) = pt.m_DA0;
        // feLog("t0 %f\n", t0_elem(0,j)); 
    }
    
    
    //Get current dof values of displacement and Angle in el
    matrix u_elem(3, 8);
    matrix t_elem(1, 8);
    CurrentDofValues(el,u_elem,t_elem);

    
    //Get the Value of the Pressure of the mixed formualtion for the element el
    double pressmix = PressureMix(el);
    
    //-------------------Gaussian Loop for Material---------------------
    //Declatration and initialization      
    
    //Declare and Initialize unassembled residual matrices
    matrix ru(24, 1);
    matrix rt(8, 1);
    ru.zero(); rt.zero();
    
    //Start of the looop
    for (int ngauss=0; ngauss<8; ++ngauss)
    {
        //Declare and Initialize the stresses
        matrix PK1_vgt(9, 1);
        matrix Pi0(3, 1);
        matrix gt(1, 1);
        PK1_vgt.zero();  Pi0.zero(); gt.zero(); 
        
        //GEt material Point and data
        FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        
        //Get Bt matrix that interpolate Grad(theta)
        matrix Bt(3, 8);
        Btmatrix(el, ngauss, Bt);
        
        //Get Bu matrix that interpoalte Grad(u)
        matrix Bu(9, 24);
        Bumatrix(Bt, Bu);
               
        //Get deformation Gradient F
        mat3d F;
        //F = pt.m_F;
        F = DeformationGradientF(u_elem, Bt);

        //Get Shape function
        matrix N_t(1,8);
        NshapeFunction(el,ngauss,N_t);
        
         
        matrix Gradt(3, 1);
        Gradt = Bt * t_elem.transpose();
        
  
        //Get value of stresses from material class function
        // get the material
        FELCEMaterial444& mat = dynamic_cast<FELCEMaterial444&>(*m_pMat);
        
        mat.StressLCE(el, ngauss, F, Gradt, pressmix,
                      PK1_vgt, Pi0,gt); 
        
        //Get unassembled residual matrices
        residualmatrices(el, ngauss, N_t, Bu, Bt, PK1_vgt, Pi0, gt, ru, rt);  

    }
      
    //------------------------Elemetal Assembly------------------------- 
    //Position of displacement and angle dof
    int pos_u[24] = {0,1,2,4,5,6,8,9,10,12,13,14,16,17,18,20,21,22,24,
                     25,26,28,29,30};
    int pos_t[8] = {3,7,11,15,19,23,27,31};
       
    for (int i = 0; i < 24 ; ++i) { for (int j = 0; j < 1 ; ++j) {
            Relem(pos_u[i],j) = -ru(i,j);}}

    for (int i = 0; i < 8 ; ++i) { for (int j = 0; j < 1 ; ++j) {
            Relem(pos_t[i],j) = -rt(i,j);}}
    
    
    for (int i=0; i<32; ++i)
	{
       fe[i] = Relem(i,0);
    }

    //Check for the values output - Print dof values and the matrices and input in matlab to check with material443---
        //double U_check[32] = { 0.0 };
        //for (int i = 0; i < 8; ++i)
        //{
        //    FEModel& fem_check = *GetFEModel();
        //    FEMesh& mesh_check = fem_check.GetMesh();
        //    FENode Noderef_check = mesh_check.Node(el.m_node[i]);
        //    U_check[4 * i] = Noderef_check.get(m_dofU[0]);
        //    U_check[4 * i + 1] = Noderef_check.get(m_dofU[1]);
        //    U_check[4 * i + 2] = Noderef_check.get(m_dofU[2]);
        //    U_check[4 * i + 3] = Noderef_check.get(m_dofDA[0]);
        //}
        //feLog("U=\n");
        //for (int nr = 0; nr < 32; ++nr) {
        //        feLog("%.15le\n", U_check[nr]);
        //} feLog("\n***************************************************\n");
        //
        //
        //feLog("ru=\n");
        //for (int nr = 0; nr < 24; ++nr) {
        //    for (int nc = 0; nc < 1; ++nc) {
        //        feLog("%.15le\t", ru(nr, nc));
        //    } feLog("\n");
        //} feLog("\n****************************************************\n");
    
}

//-----------------------------------------------------------------------------
//! Calculates element stiffness element matrix
void FELCEDomain::ElementStiffness(FESolidElement &el, matrix& ke)
{
    
    // Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 4*neln;
    double*	gw = el.GaussWeights();
    
    //Initialize Kelem
    matrix Kelem(32,32);
    Kelem.zero();

    //Get Initial Angles
    //todo get initial angles from the initial condition not specify here
    matrix t0_elem(1, 8);
    //for (int q = 0; q < 8; ++q) { t0_elem(0, q) = 1.5707963267948966; }
    //for (int q = 0; q < 8; ++q) { t0_elem(0, q) = acos(-1.0)/2.0; }

    for (int j = 0; j < 8; ++j)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(j);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        t0_elem(0, j) = pt.m_DA0;
        //feLog("t0 %f\n", t0_elem(0, j));
    }

    //Get current dof values of displacement and Angle in el
    matrix u_elem(3, 8);
    matrix t_elem(1, 8);
    CurrentDofValues(el,u_elem,t_elem);

    
    //Get the Value of the Pressure of the mixed formualtion for the element el
    double pressmix = PressureMix(el);
    
    //-------------------Gaussian Loop for Material---------------------
    //Declatration and initialization      
    
    //Declare and Initialize unassembled Stiffness matrices
    matrix K11(24, 24), K_uu(24, 24);
    double K22=0.0, K23=0.0, K32=0.0;
    matrix K13(24, 1);
    matrix K31(1, 24);
    matrix K_ut(24, 8);
    matrix K_tu(8, 24);
    matrix K_tt(8, 8);
    K11.zero(); K_uu.zero(); K13.zero(); K31.zero(); K_ut.zero();
    K_tu.zero(); K_tt.zero();
    
    //Start of the looop
    for (int ngauss=0; ngauss<8; ++ngauss)
    {
        //Declare and Initialize the stress derivatives
        matrix  dgdF_vgt(9, 1), dJdF_vgt(9, 1), dPdt_vgt(9, 1);
        matrix  dgdt(1, 1);
        matrix  dPdF_vgt(9, 9);
        mat3d dPidGradt;
        dgdF_vgt.zero(); dJdF_vgt.zero();
        dgdt.zero(); dPdF_vgt.zero(); dPdt_vgt.zero();
        dPidGradt.zero();

        //Get Material Point and data
        FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        
        //Get Bt matrix that interpolate Grad(theta)
        matrix Bt(3, 8);
        Btmatrix(el, ngauss, Bt);
        

        //Get Bu matrix that interpoalte Grad(u)
        matrix Bu(9, 24);
        Bumatrix(Bt, Bu);
               
        //Get deformation Gradient F
        mat3d F;
        //F = pt.m_F;
        F = DeformationGradientF(u_elem, Bt);
        
        //Get Shape function
        matrix N_t(1,8);
        NshapeFunction(el,ngauss,N_t);
        
        
        matrix Gradt(3, 1);
        Gradt = Bt * t_elem.transpose();
        
  
        //Get value of stresses from material class function
        // get the material
        FELCEMaterial444& mat = dynamic_cast<FELCEMaterial444&>(*m_pMat);
        
        mat.TangentLCE(el, ngauss, F, Gradt, pressmix, 
                       dPdF_vgt, dPdt_vgt, dgdF_vgt, dPidGradt, dgdt,
                       dJdF_vgt); 

        double bulkmod =1000.0*( mat.m_mu_nh(mp) + mat.m_mu_nc(mp));
        
        //Calculate unassembled element matrices
        stiffnessmatrices(el, ngauss,N_t, Bu, Bt, bulkmod, dPdF_vgt, dPdt_vgt,
                        dgdF_vgt, dPidGradt, dgdt, dJdF_vgt, K11, K22,
                         K23, K32, K13, K31, K_ut, K_tu, K_tt);  

    }
    
    K_uu = K11 + ((K13*K31) * (K22/(K23*K32)));
    
    //------------------------Elemetal Assembly------------------------- 
    //Position of displacement and angle dof
    int pos_u[24] = {0,1,2,4,5,6,8,9,10,12,13,14,16,17,18,20,21,22,24,
                     25,26,28,29,30};
    int pos_t[8] = {3,7,11,15,19,23,27,31};
                   
    for (int i = 0; i < 24 ; ++i) { for (int j = 0; j < 24 ; ++j) {
            Kelem(pos_u[i],pos_u[j]) = K_uu(i,j);}}

    for (int i = 0; i < 24 ; ++i) { for (int j = 0; j < 8 ; ++j) {
            Kelem(pos_u[i],pos_t[j]) = K_ut(i,j);}}

    for (int i = 0; i < 8 ; ++i) { for (int j = 0; j < 24 ; ++j) {
            Kelem(pos_t[i],pos_u[j]) = K_tu(i,j);}}

    for (int i = 0; i < 8 ; ++i) { for (int j = 0; j < 8 ; ++j) {
            Kelem(pos_t[i],pos_t[j]) = K_tt(i,j);}}
    
    //------------------------------------------------------------------    
    ke = Kelem; 
    
}


//------------------Function definations--------------------------------
//----------------------------------------------------------------------
double FELCEDomain::PressureMix(FESolidElement& el)
{

    //Get material point data as mp of the gaussian point n
    FELCEMaterial444& mat = dynamic_cast<FELCEMaterial444&>(*m_pMat);
    double bulkmod = 1000.0 * (mat.m_mu_nh(0) + mat.m_mu_nc(0));; //Here assumed that bulkmod is uniform in all gauss point

    //Get current dof values of displacement and Angle in el
    matrix u_elem(3, 8);
    matrix t_elem(1, 8);
    CurrentDofValues(el, u_elem, t_elem);
    //---------------Gaussian loop for volume change--------------------
    double vol0 = 0.0;
    double vol = 0.0;
    for (int ngauss=0; ngauss<8; ++ngauss)
    {
        
        //Get determinant of Jacobi in reference (From FESolidDomain) 
        vec3d N_X[8];
        double detJac = ShapeGradient0(el, ngauss, N_X);
        double* gw = el.GaussWeights();

        //Get Bt matrix that interpolate Grad(theta)
        matrix Bt(3, 8);
        Btmatrix(el, ngauss, Bt);

        //get deformation Gradient F
        mat3d F;
        FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        //F = pt.m_F;
        F = DeformationGradientF(u_elem, Bt);
      

        //Get determinat of deformation Gradient F
        double detF = F.det();
        vol0 = vol0 + detJac*gw[ngauss];
        vol = vol + detF*detJac*gw[ngauss];
    }
    
    //double vol0 = Volume(el);
    //double vol = CurrentVolume(el);
    double jbar = vol/vol0;
    double pressmix = bulkmod*(jbar-1);
    
    return pressmix;
}

//----------------------------------------------------------------------
void FELCEDomain::CurrentDofValues(FESolidElement& el,matrix& u_elem,matrix& t_elem)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh1 = fem.GetMesh();
    
    for (int j = 0; j < 8; ++j) {
        FENode Noderef = mesh1.Node(el.m_node[j]);
        u_elem(0, j) = Noderef.get(m_dofU[0]);       
        u_elem(1, j) = Noderef.get(m_dofU[1]);
        u_elem(2, j) = Noderef.get(m_dofU[2]);       
        t_elem(0, j) = Noderef.get(m_dofDA[0]);        
    }
}

//----------------------------------------------------------------------
void FELCEDomain::NshapeFunction(FESolidElement& el,int ngauss,matrix& N_t)
{
    double xi = el.gr(ngauss);
    double eta = el.gs(ngauss);
    double zeta = el.gt(ngauss);
    double Nshape[8];
    
    el.shape_fnc(Nshape, xi, eta, zeta);
    for (int i = 0; i < 8; ++i)
    {
       N_t(0,i) = Nshape[i]; 
       
    }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
void FELCEDomain::Btmatrix(FESolidElement el, int ngauss, matrix &Bt)
{
    
    vec3d N_X[8];
    
    double detJac = ShapeGradient0(el, ngauss, N_X);
    Bt.zero();    
    for (int i = 0; i < 8; ++i)
		{
			Bt[0][i] = N_X[i].x;
			Bt[1][i] = N_X[i].y;
			Bt[2][i] = N_X[i].z;
		}
    
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
mat3d FELCEDomain::DeformationGradientF(matrix u_elem, matrix Bt)
{
    matrix eye(3,3);
    eye.zero();
    eye[0][0] = 1.0;
    eye[1][1] = 1.0;
    eye[2][2] = 1.0;

    matrix Ftemp;
    Ftemp = eye + (u_elem * Bt.transpose());

    mat3d F;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                F(i, j) = Ftemp(i, j);
            }
        }

    return F;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
void FELCEDomain::Bumatrix(matrix Bt, matrix &Bu)
{
    matrix dNu_dX(3, 3);
    Bu.zero();
    int pos1[3] = { 0 };
    int pos2[3] = { 0 };
    
    for (int i = 0; i < 8; ++i)
    {
        dNu_dX.zero();
        dNu_dX(0, 0) = Bt(0, i);
        dNu_dX(1, 1) = Bt(1, i);
        dNu_dX(2, 2) = Bt(2, i);
        
        pos1[0]=0; pos1[1] = 1; pos1[2] = 2;
        pos2[0] = 3 * i; pos2[1] = 3 * i + 1; pos2[2] = 3 * i + 2;
        for (int ii = 0; ii < 3 ; ++ii) {
            for (int j = 0; j < 3 ; ++j) {
            Bu(pos1[ii],pos2[j]) = dNu_dX(ii,j);
            }
        }
               
        dNu_dX.zero();
        dNu_dX(0, 0) = Bt(1, i);
        dNu_dX(1, 1) = Bt(2, i);
        dNu_dX(2, 2) = Bt(0, i);
            pos1[0] = 3; pos1[1] = 4; pos1[2] = 5;
        pos2[0] = 3 * i; pos2[1] = 3 * i + 1; pos2[2] = 3 * i + 2;
        for (int ii = 0; ii < 3 ; ++ii) {
            for (int j = 0; j < 3 ; ++j) {
            Bu(pos1[ii],pos2[j]) = dNu_dX(ii,j);
            }
        }

        dNu_dX.zero();
        dNu_dX(0, 1) = Bt(0, i);
        dNu_dX(1, 2) = Bt(1, i);
        dNu_dX(2, 0) = Bt(2, i);
        pos1[0] = 6; pos1[1] = 7; pos1[2] = 8;
        pos2[0] = 3 * i; pos2[1] = 3 * i + 1; pos2[2] = 3 * i + 2;
        for (int ii = 0; ii < 3 ; ++ii) { 
            for (int j = 0; j < 3 ; ++j) {
            Bu(pos1[ii],pos2[j]) = dNu_dX(ii,j);
            }
        }
        
    } 
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
void FELCEDomain::stiffnessmatrices(FESolidElement el, int ngauss, 
                     matrix N_t,matrix Bu, matrix Bt,double bulkmod,
                     matrix dPdF_vgt, matrix dPdt_vgt,matrix dgdF_vgt,
                     mat3d dPidGradt, matrix dgdt, matrix dJdF_vgt,
                     matrix& K11, double& K22, double& K23,double& K32,
                     matrix& K13, matrix& K31, matrix& K_ut, matrix& K_tu,
                     matrix& K_tt)
{
    
    vec3d N_X[8];
    double detJac = ShapeGradient0(el, ngauss, N_X);
    double*	gw = el.GaussWeights();
    
    double dvol = detJac*gw[ngauss];
    
    K11 = K11 + (Bu.transpose()*dPdF_vgt)*Bu *dvol;
    
    K13 = K13 + (Bu.transpose()*dJdF_vgt) *dvol;
    
    K22 = K22 + bulkmod*dvol;  
    K23 = K23 - dvol;
    
    K31 = K13.transpose();
    K32 = K23;
    
    K_ut = K_ut + ((Bu.transpose()*dPdt_vgt)*N_t) * dvol;
   
    K_tu = K_tu - ((N_t.transpose()*dgdF_vgt.transpose()))*Bu * dvol;
    
    K_tt = K_tt + ((Bt.transpose()*dPidGradt)*Bt) * dvol 
            - ((N_t.transpose()*dgdt)*N_t) * dvol;  
    
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
void FELCEDomain::residualmatrices(FESolidElement el, int ngauss, 
                     matrix N_t, matrix Bu, matrix Bt,matrix PK1_vgt, 
                     matrix Pi0, matrix gt, matrix& ru,matrix& rt)
                        
{
    vec3d N_X[8];
    double detJac = ShapeGradient0(el, ngauss,  N_X);
    double*	gw = el.GaussWeights();
       
    double dvol = detJac*gw[ngauss];
    
    //residual matrices
    ru = ru + (Bu.transpose()*PK1_vgt) * dvol;
    rt = rt + (Bt.transpose()*Pi0) * dvol - (N_t.transpose()*gt) * dvol;
    
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
mat3d FELCEDomain::LCEAveragePK1Stress(FESolidElement el)
{
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    const int ndof = 4 * neln;
    double* gw = el.GaussWeights();

    //Initialize Average PK1;
    mat3d AvgPK1;
    AvgPK1.zero();
    

    //Get Initial Angles
    //todo get initial angles from the initial condition not specify here
    matrix t0_elem(1, 8);
    for (int j = 0; j < 8; ++j)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(j);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        t0_elem(0, j) = pt.m_DA0;
        //feLog("t0 %f\n", t0_elem(0, j));
    }

    //Get current dof values of displacement and Angle in el
    matrix u_elem(3, 8);
    matrix t_elem(1, 8);
    CurrentDofValues(el, u_elem, t_elem);


    //Get the Value of the Pressure of the mixed formualtion for the element el
    double pressmix = PressureMix(el);

    //-------------------Gaussian Loop for Material---------------------
    //Declatration and initialization      

    //Start of the looop
    for (int ngauss = 0; ngauss < 8; ++ngauss)
    {
        //Declare and Initialize the stresses
        matrix PK1_vgt(9, 1);
        matrix Pi0(3, 1);
        matrix gt(1, 1);
        PK1_vgt.zero();  Pi0.zero(); gt.zero();

        //GEt material Point and data
        FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

        //Get Bt matrix that interpolate Grad(theta)
        matrix Bt(3, 8);
        Btmatrix(el, ngauss, Bt);

        //Get Bu matrix that interpoalte Grad(u)
        matrix Bu(9, 24);
        Bumatrix(Bt, Bu);

        //Get deformation Gradient F
        mat3d F;
        //F = pt.m_F;
        F = DeformationGradientF(u_elem, Bt);

        //Get Shape function
        matrix N_t(1, 8);
        NshapeFunction(el, ngauss, N_t);

        //Calculate theta, theta0 and Grad of theta at the gauss point
        matrix t_temp(1, 1);
        t_temp = N_t * t_elem.transpose();
        double theta = t_temp(0, 0);

        matrix t0_temp(1, 1);
        t0_temp = N_t * t0_elem.transpose();
        double theta0 = t0_temp(0, 0);

        matrix Gradt(3, 1);
        Gradt = Bt * t_elem.transpose();


        //Get value of stresses from material class function
        // get the material
        FELCEMaterial444& mat = dynamic_cast<FELCEMaterial444&>(*m_pMat);

        mat3d PK1;
        PK1.zero();
        PK1 = mat.StressLCE(el, ngauss, F, Gradt, pressmix,
            PK1_vgt, Pi0, gt);
        AvgPK1 = AvgPK1 + PK1;
    }
    AvgPK1 = AvgPK1 * (1.0 / 8.0);

    return AvgPK1;
}
//----------------------------Function definations----------------------
//----------------------------------------------------------------------
void FELCEDomain::tens33_9(mat3d A, matrix& B)
{
    //Notice the voigt order is different than I used in matlab tens33_9

    B(0, 0) = A(0, 0);
    B(1, 0) = A(1, 1);
    B(2, 0) = A(2, 2);
    B(3, 0) = A(0, 1);
    B(4, 0) = A(1, 2);
    B(5, 0) = A(2, 0);
    B(6, 0) = A(1, 0);
    B(7, 0) = A(2, 1);
    B(8, 0) = A(0, 2);
}
//----------------------------------------------------------------------
void FELCEDomain::tens3333_99(tens4d A, matrix& B)
{
    //Notice the voigt order is different than I used in matlab tens3333_99

    double ind1[9] = {0,1,2,0,1,2,1,2,0};
    double ind2[9] = {0,1,2,1,2,0,0,1,2};

    for (int m=0; m<9; m++)
    {
        int i = ind1[m];
        int j = ind2[m];
        for (int n=0; n<9; n++)
        {
            int k = ind1[n];
            int l = ind2[n];
            B(m,n) = A(i,j,k,l);
            
        }
    }
}
//=================================================================================================
BEGIN_FECORE_CLASS(FEStandardLCEDomain, FELCEDomain)
	ADD_PARAMETER(m_elemType, "elem_type", FE_PARAM_ATTRIBUTE, "$(solid_element)\0");
END_FECORE_CLASS();

FEStandardLCEDomain::FEStandardLCEDomain(FEModel* fem) : FELCEDomain(fem)
{

}
