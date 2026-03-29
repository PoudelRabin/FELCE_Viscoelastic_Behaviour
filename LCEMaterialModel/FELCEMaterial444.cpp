//********************NOTE BY RABIN POUDEL******************************
/*Note for the fourth order operation
dyad operation of the tens4d are diffrent to the one done in ...
matlab (ie Mikhail Itskov Book). The equivalent are

TensorProd_AoB(A,B) = dyad1(A,B) 
TensorProd_AxB(A,B) = dyad2(A,B.transpose())
TensorProd_AxBt(A,B) = dyad3(A, B.transpose())

AoB(i,j,k,l) = A(i,j)*B(k,l)----
AxB(i,j,k,l) = A(i,k)*B(l,j)    |--Defination in Mikhail Itskov Book 
AxBt(i,j,k,l) = A(i,l)*B(k,j)---

A(dyad1)B(i,j,k,l) = A(i,j)*B(k,l)----
A(dyad2)B(i,j,k,l) = A(i,k)*B(j,l)    |--Defination in tens4d.hpp
A(dyad3)B(i,j,k,l) = A(i,l)*B(j,k)----


The 9*9 matrix obtaine by tensor3333_99 from 4th order is not same
extract 9*9 used in tens4d. The order of indices for voigt notation
are not in same order (see tens4d.h and tens4d.hpp). Also 9*9 
matrix function is to be first transposed for our use(since tens4d.h
first fills the single column). we use the voigt order as in tens4d.h*/
//**********************************************************************




#include "stdafx.h"
#include "FELCEMaterial444.h"
#include<FECore/FESolidDomain.h>
#include "FECore/log.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FELCEMaterial444, FEElasticMaterial)
	ADD_PARAMETER(m_mu_nh, FE_RANGE_GREATER(0.0), "mu_nh");
	ADD_PARAMETER(m_kfrank, FE_RANGE_GREATER(0.0), "K");
    ADD_PARAMETER(m_a, FE_RANGE_GREATER(0.0), "a");
    ADD_PARAMETER(m_mu_nc, FE_RANGE_GREATER(0.0), "mu_nc");
    ADD_PARAMETER(m_mu_neq_nc, FE_RANGE_GREATER(0.0), "mu_neq_nc");
    ADD_PARAMETER(m_eta_u, FE_RANGE_GREATER(0.0), "eta_u");
	ADD_PARAMETER(m_mu_neq_nh, FE_RANGE_GREATER(0.0), "mu_neq_nh");
	ADD_PARAMETER(m_eta_t, FE_RANGE_GREATER(0.0), "eta_t");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//constructer
FELCEMaterial444::FELCEMaterial444(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_Base = 0;
    
}
//-----------------------------------------------------------------------------
//! create material point data for this material
FEMaterialPointData* FELCEMaterial444:: CreateMaterialPointData()
{
    //return new FELCEMaterialPoint(m_Base->CreateMaterialPointData());
    return new FELCEMaterialPoint(new FEElasticMaterialPoint);
    
}

//-----------------------------------------------------------------------------
mat3ds FELCEMaterial444::Stress(FEMaterialPoint& mp)
{
    mat3ds s;
    s.zero();
    return s;
}
//-----------------------------------------------------------------------------
tens4ds FELCEMaterial444::Tangent(FEMaterialPoint& mp)
{
    tens4ds T;
    T.zero();
    return T;
}

//-----------------------------------------------------------------------------
mat3d FELCEMaterial444::StressLCE(FESolidElement el, int ngauss, 
        mat3d F, matrix Gradt, double pressmix, matrix &PK1_vgt, 
        matrix &Pi0, matrix &gt)
            
{
    //Declare some temporary double, matrix and tensors 
    double temp1,temp2,temp3,temp4,temp5; 
    mat3d tempmat1, tempmat2, tempmat3, tempmat4, tempmat5;
    tens4d tempA, tempB, tempC, tempD, tempE;

    
    
    //-----------------------------------------------------
    //Get material point
    FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
    FEElasticMaterialPoint& Ept = *(mp.ExtractData<FEElasticMaterialPoint>());
    FELCEMaterialPoint& Lpt = *(mp.ExtractData<FELCEMaterialPoint>());
    
    //Get Current time
    double cTime = GetFEModel()->GetTime().currentTime;

    //GEt material parameter value
    //todo is writing (0) good? Change it for material point, here assumed parameter is same everywhere
    double mu_nh = m_mu_nh(mp);
    double mu_nc = m_mu_nc(mp);
    double mu_neq_nc = m_mu_neq_nc(mp);
    double eta_u = m_eta_u(mp);
	double mu_neq_nh = m_mu_neq_nh(mp);
	double eta_t = m_eta_t(mp);

    //mu = ((1.0e-10 - 2.0e-2)/(1.20-1.0))*(cTime+1.0-1.20) + 1.0e-10;
    //if (cTime>0.20) {mu = 1.0e-10;}
    double kfrank = m_kfrank(mp);
    double a = m_a(mp);

    //Get deformation gradient and its determinant
    double detF = F.det();
    mat3d Finv = F.inverse();
    mat3d Fbar = F* pow(detF, -1.0 / 3.0);

    mat3d Cbar = Fbar.transpose() * Fbar;
    
    // For neckign uncomment below block on both stress and Tangent
   // double lval[3];
   // vec3d rvec[3];
   // mat3ds Cbarsym = Cbar.sym();
   // Cbarsym.eigen(lval, rvec);
   // double lam_x = pow(lval[0],0.5);
   // //double lam_x = cTime+1.0;
   // double t1 = exp(16.0);
   // double t2 = pow(t1, (lam_x-1.45));
   // double Qorder = 0.35*((t2-1)/(t2+1) - (t1-1)/(t1+1) + 2.0);
   // a = (1.0+2.0*Qorder)/(1.0-Qorder);
    

    double a0 = a;
    double bulkmod = 1000.0*(mu_nh+mu_nc);

    //--------------------------


    // linearly varying a and eta_u based on strain
    // double lval[3];
    // vec3d rvec[3];
    // mat3ds Csym = (F.transpose() * F).sym();
    // Csym.eigen(lval, rvec);
    // double lam_x = pow(lval[0],0.5);
    // double final_total_strain = 3.5;
    // a0 = 4.0;
    // double ax = 3.0;
    // a = a0 + (ax-a0)*(lam_x-1.0)/(final_total_strain-1);
    // double eta_t0 = 800;
    // double eta_tx = 100.0;
    // eta_t = eta_t0 + (eta_tx-eta_t0)*(lam_x-1.0)/(final_total_strain-1);
    
    
    
    //-------------
    mat3d eye;
    eye.unit();
    
    //-------------
    //get theta0
    double theta0 = Lpt.m_theta0;

    //get initial spontanious deformation gradient
    mat3d G0 = Lpt.InitialSpontaniousDG(a);
    
    mat3d G02 = G0*G0;
    //-------------------

    //-------------
    //get theta0
    double theta = Lpt.m_theta;

    //get current spontanious deformation gradient
    mat3d G = Lpt.SpontaniousDG(a);
    mat3d Ginv = G.inverse();
    mat3d Ginv2 = Ginv * Ginv;


    mat3d dnTndt; //derivative of nTn wrt theta
    dnTndt.zero();
    dnTndt[0][0] = -sin(theta*2.0 );
    dnTndt[0][1] = cos(theta * 2.0);
    dnTndt[1][0] = cos(theta * 2.0);
    dnTndt[1][1] = sin(theta * 2.0);
  
    temp1 = pow(a, -2.0/3.0);
    temp2 = pow(a, 2.0/6.0);
    mat3d dGinv2dt = dnTndt * (temp1-temp2);
    //----------------------------------------------
    //Calculation of other deformation gradient
    mat3d A_DG = (Ginv * F) * G0;
    mat3d Av = Lpt.m_Av;
    mat3d Av_inv = Av.inverse();
    mat3d Ae = A_DG * Av_inv;
    mat3d Fv = Lpt.m_Fv;
    mat3d Fv_inv = Fv.inverse();
    mat3d Fe = F * Fv_inv;

    
    
    //Derivative of J=detF
    mat3d dJdF = Finv.transpose() * detF;
    
    //--------First order derivatives of free energies-----------------
    temp1 = (F*G02).dotdot(Ginv2*F); //doublecontraction(F*G02,Ginv2*F) 
    tempmat1 =  (Finv.transpose()*temp1)*(-1.0 / 3.0);
    tempmat2 = (Ginv2*F)*G02;
    mat3d delas2dF = (tempmat1 + tempmat2) * mu_nc * pow(detF,(-2.0/3.0));

    //Semi-Soft elasticity, Pure neohookean
    temp1 = (F ).dotdot( F); //doublecontraction(F,F) 
    tempmat1 = (Finv.transpose() * temp1) * (-1.0 / 3.0);
    tempmat2 = F;
    mat3d delas1dF = (tempmat1 + tempmat2) * mu_nh * pow(detF, (-2.0 / 3.0));

    mat3d delasdF = delas1dF + delas2dF;

    //derivative of neq2 enerry wrt  F
    temp1 = Ae.dotdot(Ae); //double contraction (Ae,Ae)
    tempmat1 = (Finv.transpose() * temp1) * (-1.0 / 3.0);
    tempmat2 = (Ginv * Ae) * (Av_inv.transpose() * G0);
    mat3d dneq2dF = (tempmat1 + tempmat2) * mu_neq_nc * pow(detF, (-2.0 / 3.0));

    //derivative of neq1 energy wrt to F
    temp1 = Fe.dotdot(Fe); //doublecontraction(Fe:Fe)
    tempmat1 = (Finv.transpose() * temp1) * (-1.0 / 3.0);
    tempmat2 = Fe * Fv_inv.transpose();
    mat3d dneq1dF = (tempmat1 + tempmat2) * mu_neq_nh * pow(detF, (-2.0 / 3.0));

    mat3d dneqdF = dneq1dF + dneq2dF;

 
    //------
    matrix dfrankdGradt(3, 1);
    dfrankdGradt = Gradt * kfrank;
    
    //-----
    matrix delasdt(1, 1);
    tempmat1 = (G02*Fbar.transpose()) ;
    tempmat2 = (dGinv2dt * Fbar);
    tempmat3 = tempmat1 * tempmat2;
    temp1 = tempmat3.trace();
    delasdt(0, 0) = 0.5 * mu_nc * temp1;

    //derivative of neq energy (only neq2) wrt theta
    matrix dneqdt(1, 1);
    tempmat1 = (Av_inv.transpose() * G0) * Fbar.transpose();
    tempmat2 = (dGinv2dt * Fbar) * (G0 * Av_inv);
    tempmat3 = tempmat1 * tempmat2;
    temp1 = tempmat3.trace();
    dneqdt(0,0) = 0.5 * mu_neq_nc * temp1;
     
    
    //--------Calculation of Stresses---------------
    mat3d PK1 = delasdF + dJdF * pressmix + dneqdF;
    Pi0 = dfrankdGradt;

    ////--------Later Addition to check------///
    //TODO Testing only as gt = skewsymmetric part of cauchy stress
    // mat3d Cauchy = PK1*F.transpose();
    // matrix dndt(3, 1); //director vector at current gauss point
    // dndt(0,0) = -sin(theta);
    // dndt(1,0) = cos(theta);
    // dndt(2,0) = 0.0;
    // matrix nTdndt_temp(3,3);
    // nTdndt_temp = n*dndt.transpose();
    // mat3d nTdndt;
    //  for (int i = 0; i < 3; ++i) {
    //     for (int j = 0; j < 3; ++j) {
    //         nTdndt(i, j) = nTdndt_temp(i, j);
    //     }
    // }
    // temp1 = (Cauchy.transpose() - Cauchy).dotdot(nTdndt);
    // delasdt(0,0) = -temp1;
    ///--------Later addition to be end here------////

    gt = (delasdt + dneqdt) * (-1.0);



    
    //Get PK1 in voigt notation
    tens33_9(PK1,PK1_vgt);
    
    return PK1;
}



//==============================================================================================
//==============================================================================================
void FELCEMaterial444::TangentLCE(FESolidElement el, int ngauss,
           mat3d F, matrix Gradt, double pressmix, matrix &dPdF_vgt, 
          matrix &dPdt_vgt, matrix &dgdF_vgt, mat3d &dPidGradt, matrix &dgdt,
          matrix &dJdF_vgt)          
{
    //Declare some temporary double, matrix and tensors 
    double temp1,temp2,temp3,temp4,temp5; 
    mat3d tempmat1, tempmat2, tempmat3, tempmat4, tempmat5;
    tens4d tempA, tempB, tempC, tempD, tempE;
    
    

    //-----------------------------------------------------
    //Get material point
    FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
    FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
    FELCEMaterialPoint& Lpt = *(mp.ExtractData<FELCEMaterialPoint>());

    //Get Current time
    double cTime = GetFEModel()->GetTime().currentTime;


    //GEt material parameter value
    //todo is writing (0) good? Change it for material point, here assumed parameter is same everywhere
    double mu_nh = m_mu_nh(mp);
    double mu_nc = m_mu_nc(mp);
    double mu_neq_nc = m_mu_neq_nc(mp);
    double eta_u = m_eta_u(mp);
    double mu_neq_nh = m_mu_neq_nh(mp);
    double eta_t = m_eta_t(mp)*1000;
    
    double kfrank = m_kfrank(mp);
    double a = m_a(mp);

    //Get deformation gradient and its determinant
    double detF = F.det();
    mat3d Finv = F.inverse();
    mat3d Fbar = F* pow(detF, -1.0 / 3.0);

    mat3d Cbar = Fbar.transpose() * Fbar;


    // For neckign uncomment below block on both stress and Tangent
   // double lval[3];
   // vec3d rvec[3];
   // mat3ds Cbarsym = Cbar.sym();
   // Cbarsym.eigen(lval, rvec);
   // double lam_x = pow(lval[0],0.5);
   // //double lam_x = cTime+1.0;
   // double t1 = exp(16.0);
   // double t2 = pow(t1, (lam_x-1.45));
   // double Qorder = 0.35*((t2-1)/(t2+1) - (t1-1)/(t1+1) + 2.0);
   // a = (1.0+2.0*Qorder)/(1.0-Qorder);


   


    

    double a0 = a;
    double bulkmod = 1000.0 * (mu_nh + mu_nc);
    // pertubating the current theta
    //theta = theta + theta0 - 90.0;
    //--------------------------
    
    // linearly varying a and eta_u based on strain
    // double lval[3];
    // vec3d rvec[3];
    // mat3ds Csym = (F.transpose() * F).sym();
    // Csym.eigen(lval, rvec);
    // double lam_x = pow(lval[0],0.5);
    // double final_total_strain = 3.5;
    // a0 = 4.0;
    // double ax = 3.0;
    // a = a0 + (ax-a0)*(lam_x-1.0)/(final_total_strain-1);
    // double eta_t0 = 800;
    // double eta_tx = 100.0;
    // eta_t = eta_t0 + (eta_tx-eta_t0)*(lam_x-1.0)/(final_total_strain-1);
//
    
    
    //-------------
    mat3d eye;
    eye.unit();

    //-------------
    //get theta0
    double theta0 = Lpt.m_theta0;

    //get initial spontanious deformation gradient
    mat3d G0 = Lpt.InitialSpontaniousDG(a);
    mat3d G02 = G0 * G0;
    //-------------------

    //-------------
    //get theta0
    double theta = Lpt.m_theta;

    //get current spontanious deformation gradient
    mat3d G = Lpt.SpontaniousDG(a);
    mat3d Ginv = G.inverse();
    mat3d Ginv2 = Ginv * Ginv;

    
      
    mat3d dnTndt; //derivative of nTn wrt theta
    dnTndt.zero();
    dnTndt[0][0] = -sin(theta*2.0 );
    dnTndt[0][1] = cos(theta * 2.0);
    dnTndt[1][0] = cos(theta * 2.0);
    dnTndt[1][1] = sin(theta * 2.0);
    
    mat3d ddnTnddt; //second derivative of nTn wrt theta   
    ddnTnddt.zero();
    ddnTnddt[0][0] = -2.0 * cos(theta * 2.0);
    ddnTnddt[0][1] = -2.0 * sin(theta * 2.0);
    ddnTnddt[1][0] = -2.0 * sin(theta * 2.0);
    ddnTnddt[1][1] = 2.0*cos(theta * 2.0);
   
    
    temp1 = pow(a, -2.0/3.0);
    temp2 = pow(a, 2.0/6.0);
    mat3d dGinv2dt = dnTndt * (temp1-temp2);
    mat3d ddGinv2ddt = ddnTnddt * (temp1-temp2);

    //--------------------------------------------
    //Calculation of other deformation gradient
    mat3d A_DG = (Ginv * F) * G0;
    mat3d Av = Lpt.m_Av;
    mat3d Av_inv = Av.inverse();
    mat3d Ae = A_DG * Av_inv;
    mat3d Fv = Lpt.m_Fv;
    mat3d Fv_inv = Fv.inverse();
    mat3d Fe = F * Fv_inv;
    
    //Derivative of J=detF
    mat3d dJdF = Finv.transpose() * detF;
    tempA = dyad1(Finv.transpose(),dJdF);
    tempB = dyad3(-Finv.transpose(),Finv);
    tens4d ddJddF = tempA + tempB*detF;
       
    //---------Second order derivate of free energies----------------
    //First, we need delasdF to calculate ddelasddF 
    temp1 = (F*G02).dotdot(Ginv2*F); //doublecontraction(F*G02,Ginv2*F) 
    tempmat1 =  (Finv.transpose()*temp1)*(-1.0 / 3.0);
    tempmat2 = (Ginv2*F)*G02;
    mat3d delas2dF = (tempmat1 + tempmat2) * mu_nc * pow(detF,(-2.0/3.0)); 

    //Semi-Soft elasticity, Pure neohookean
    temp1 = (F).dotdot(F); //doublecontraction(F,F) 
    tempmat1 = (Finv.transpose() * temp1) * (-1.0 / 3.0);
    tempmat2 = F;
    mat3d delas1dF = (tempmat1 + tempmat2) * mu_nh * pow(detF, (-2.0 / 3.0));

    mat3d delasdF = delas1dF + delas2dF;

    //derivative of neq2 enerfy wt to F
    temp1 = Ae.dotdot(Ae); //double contraction (Ae,Ae)
    tempmat1 = (Finv.transpose() * temp1) * (-1.0 / 3.0);
    tempmat2 = (Ginv * Ae) * (Av_inv.transpose() * G0);
    mat3d dneq2dF = (tempmat1 + tempmat2) * mu_neq_nc * pow(detF, (-2.0 / 3.0));

	//derivative of neq1 energy wrt to F
    temp1 = Fe.dotdot(Fe); //doublecontraction(Fe:Fe)
    tempmat1 = (Finv.transpose() * temp1) * (-1.0 / 3.0);
    tempmat2 = Fe * Fv_inv.transpose();
    mat3d dneq1dF = (tempmat1 + tempmat2) * mu_neq_nh * pow(detF, (-2.0 / 3.0));


    
    //neo classical 
    tempA = dyad2(Ginv2,G02.transpose());
    tempB= dyad1(Finv.transpose(),(Ginv2*(F*G02)))*2.0;
    temp1 = (F*G02).dotdot(Ginv2*F); //doublecontraction(F*G02,Ginv2*F);
    tempC = dyad3(-Finv.transpose(),Finv) * temp1;
    tempD = dyad1(delas2dF,Finv.transpose())*(2.0/3.0);
    tens4d ddelas2ddF = ((tempB + tempC)*(-1.0/3.0) + tempA)* 
                             mu_nc*pow(detF,(-2.0/3.0)) - tempD;

    //Semi-Soft elasticity, Pure neohookean
    tempA = dyad2(eye, eye.transpose());
    tempB = dyad1(Finv.transpose(), F) * 2.0;
    temp1 = (F ).dotdot( F); //doublecontraction(F,F);
    tempC = dyad3(-Finv.transpose(), Finv) * temp1;
    tempD = dyad1(delas1dF, Finv.transpose()) * (2.0 / 3.0);
    tens4d ddelas1ddF = ((tempB + tempC) * (-1.0 / 3.0) + tempA) *
        mu_nh * pow(detF, (-2.0 / 3.0)) - tempD;
   
    tens4d ddelasddF = ddelas1ddF + ddelas2ddF;

    // derivative of dneq2dF energy wrt to F
    tempmat1 = (G0 * Av_inv) * (Av_inv.transpose() * G0);
    tempA = dyad2(Ginv2, tempmat1.transpose());
    tempmat2 = (Ginv * Ae) * (Av_inv.transpose() * G0);
    tempB = dyad1(Finv.transpose(), tempmat2) * 2.0;
    temp1 = Ae.dotdot(Ae);
    tempC = dyad3(-Finv.transpose(), Finv) * temp1;
    tempD = dyad1(dneq2dF, Finv.transpose()) * (2.0 / 3.0);
    tens4d ddneq2ddF = ((tempB + tempC) * (-1.0 / 3.0) + tempA) *
        mu_neq_nc * pow(detF, (-2.0 / 3.0)) - tempD;

	// derivative of dneq1dF energy wrt to F
	tempmat1 = (Fv_inv * Fv_inv.transpose());
	tempA = dyad2(eye, tempmat1.transpose());
	tempmat2 = Fe * Fv_inv.transpose();
	tempB = dyad1(Finv.transpose(), tempmat2) * 2.0;
	temp1 = Fe.dotdot(Fe); //doublecontraction(Fe:Fe)
	tempC = dyad3(-Finv.transpose(), Finv) * temp1;
	tempD = dyad1(dneq1dF, Finv.transpose()) * (2.0 / 3.0);
	tens4d ddneq1ddF = ((tempB + tempC) * (-1.0 / 3.0) + tempA) * 
		mu_neq_nh * pow(detF, (-2.0 / 3.0)) - tempD;

	tens4d ddneqddF = ddneq1ddF + ddneq2ddF;
                             
    //----
    temp1 = (F*G02).dotdot(dGinv2dt*F); //doublecontraction(F*G02,dGinv2dt*F);
    tempmat1 = Finv.transpose()*(1.0/3.0)*temp1;
    mat3d ddelasdFdt = ((dGinv2dt*(F*G02))- tempmat1)*mu_nc*pow(detF,(-2.0/3.0));

    //derivative of dneqdF (only neq2) wrt to theta
    tempmat1 = (Av_inv.transpose() * G0);
    tempmat2 = (dGinv2dt * F) * (G0 * Av_inv);
    tempmat3 = (tempmat1 * F.transpose()) * tempmat2;
    temp1 = eye.dotdot(tempmat3);
    tempmat4 = (tempmat2 * tempmat1);
    tempmat5 = Finv.transpose() * (1.0 / 3.0) * temp1;
    temp2 = mu_neq_nc * pow(detF, (-2.0 / 3.0));
    mat3d ddneqdFdt = (tempmat4 - tempmat5) * temp2;

 
    //----
    mat3d ddelasdtdF = ddelasdFdt;

    // derivative of dneqdt wrt F
    mat3d ddneqdtdF = ddneqdFdt;
    
    //-----
    mat3d ddfrankddGradt = eye*kfrank;
    
    //-----
    matrix ddelasdtdt(1, 1);
    tempmat1 = (G02*Fbar.transpose()) * (ddGinv2ddt*Fbar);
    temp1 = tempmat1.trace();
    ddelasdtdt(0, 0) = 0.5 * mu_nc * temp1;

    //derivative of dneqdt (only neq2) wrt theta
    matrix ddneqdtdt(1, 1);
    tempmat1 = (Av_inv.transpose() * G0) * Fbar.transpose();
    tempmat2 = (ddGinv2ddt * Fbar) * (G0 * Av_inv);
    tempmat3 = tempmat1 * tempmat2;
    temp1 = tempmat3.trace();
    ddneqdtdt(0, 0) = 0.5 * mu_neq_nc * temp1;

    
    //--------Calculation of Stress Tangents---------------
    tens4d dPdF = ddJddF*pressmix + ddelasddF + ddneqddF;
    mat3d dPdt = ddelasdFdt + ddneqdFdt;
    mat3d dgdF = -ddelasdtdF - ddneqdtdF;
    dPidGradt = ddfrankddGradt;

    ////--------Later Addition to check------///
    //TODO Testing only as gt = skewsymmetric part of cauchy stress
    // mat3d PK1 = delasdF + dJdF*pressmix;
    // mat3d Cauchy = PK1*F.transpose();
    // mat3d dCauchydt = dPdt * F.transpose();
    // matrix dndt(3, 1); //director vector at current gauss point
    // dndt(0,0) = -sin(theta);
    // dndt(1,0) = cos(theta);
    // dndt(2,0) = 0.0;
    // matrix nTdndt_temp(3,3);
    // nTdndt_temp = n*dndt.transpose();
    // mat3d nTdndt;
    //  for (int i = 0; i < 3; ++i) {
    //     for (int j = 0; j < 3; ++j) {
    //         nTdndt(i, j) = nTdndt_temp(i, j);
    //     }
    // }
    // matrix dndtTdndt_temp(3,3);
    // dndtTdndt_temp = dndt*dndt.transpose();
    // mat3d dndtTdndt;
    //  for (int i = 0; i < 3; ++i) {
    //     for (int j = 0; j < 3; ++j) {
    //         dndtTdndt(i, j) = dndtTdndt_temp(i, j);
    //     }
    // }
    // temp1 = (dCauchydt.transpose() - dCauchydt).dotdot(nTdndt);
    // temp2 = (Cauchy.transpose() - Cauchy).dotdot(dndtTdndt);
    // ddelasdtdt(0,0) = -temp1 - temp2;
    ///--------Later addition to be end here------////

    dgdt = (ddelasdtdt+ddneqdtdt)*(-1.0);
    
    //------Transforming Tensor to voigt notation
    //matrix(9,9) dPdF_vgt = tens3333_99(dPdF);
    //double dPdF_vgt_temp[9][9];
    //dPdF.extract(dPdF_vgt_temp);
    //matrix dPdF_vgt_temp2(9, 9);
    //for (int i = 0; i < 9; ++i) {
    //    for (int j = 0; j < 9; ++j) {
//
    //        dPdF_vgt_temp2(i, j) = dPdF_vgt_temp[i][j];
    //    }
    //}
    //dPdF_vgt = dPdF_vgt_temp2;//.transpose(); //I dont think we need to transpose

    tens3333_99(dPdF, dPdF_vgt);
    //dPdF_vgt = dPdF_vgt.transpose();
    tens33_9(dPdt,dPdt_vgt);
    tens33_9(dgdF,dgdF_vgt);
    tens33_9(dJdF,dJdF_vgt);
    
}

//----------------------------Function definations----------------------
//----------------------------------------------------------------------
void FELCEMaterial444::tens33_9(mat3d A, matrix& B)
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
void FELCEMaterial444::tens3333_99(tens4d A, matrix& B)
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
   
