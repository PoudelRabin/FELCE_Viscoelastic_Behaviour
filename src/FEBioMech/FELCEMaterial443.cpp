#include "stdafx.h"
#include "FELCEMaterial443.h"
#include<FECore/FESolidDomain.h>
#include "FECore/log.h"
#include <FECore/FEModel.h>




//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FELCEMaterial443, FEElasticMaterial)
	ADD_PARAMETER(m_mu, FE_RANGE_GREATER(0.0), "mu");
	ADD_PARAMETER(m_kfrank, FE_RANGE_GREATER(0.0), "K");
    ADD_PARAMETER(m_a, FE_RANGE_GREATER(0.0), "a");
    ADD_PARAMETER(m_bulk, FE_RANGE_GREATER(0.0), "B");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FELCEMaterial443::FELCEMaterial443(FEModel* pfem) : FEElasticMaterial(pfem) {}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
mat3ds FELCEMaterial443::Stress(FEMaterialPoint& mp)
{
    mat3ds s;
    s.zero();
    return s;
}
//-----------------------------------------------------------------------------
tens4ds FELCEMaterial443::Tangent(FEMaterialPoint& mp)
{
    tens4ds T;
    T.zero();
    return T;
}

//-----------------------------------------------------------------------------
void FELCEMaterial443::MaterialRoutine(FESolidElement el, int ngauss, matrix Bt, mat3d F, matrix t_e,
                                  double pressmix, matrix& PK1_vgt, matrix& Pi0,
                                  matrix& gt, matrix& dPdF_vgt, matrix& dPdt_vgt,
                                  matrix& dgdF_vgt, mat3d& dPidGradt, matrix& dgdt,
                                  matrix& dJdF_vgt)
                     
{
    
    //Get material point
    FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
    FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
    //GEt material parameter value
    //todo is writing (0) good? Change it for material point, here assumed parameter is same everywhere
    double mu = m_mu(mp);
    double kfrank = m_kfrank(mp);
    double a = m_a(mp);
    double bulkmod = m_bulk(mp);
    
    
    //Get deformation gradient and its determinant
    double detF = F.det();
    mat3d Finv = F.inverse();
    mat3d B = F*F.transpose();
    mat3d Bbar = B* pow(detF, -2.0 / 3.0);
    
    //Get Shape function in matrix form
    //Get integration point from ngauss (See FESolidElement.h)
    double xi = el.gr(ngauss);
    double eta = el.gs(ngauss);
    double zeta = el.gt(ngauss);
    double Nshape[8];
    
    el.shape_fnc(Nshape, xi, eta, zeta);
    matrix N_t(1, 8);
    for (int i = 0; i < 8; ++i)
    {
       N_t(0,i) = Nshape[i]; 
       
    }
    //Value of angle theta at current gaussian point
    //todo better to pass theta and Gradt values from domain
    matrix theta(1, 1);
    theta = N_t * t_e.transpose(); 
    matrix Gradt(3, 1);
    Gradt = Bt * t_e.transpose(); //Grad of theta at Gauss Point
    
    
    matrix n(3, 1); //director vector at current gauss point
    n(0,0) = cos(theta(0,0));
    n(1,0) = sin(theta(0,0));
    n(2,0) = 0.0;
    
    matrix nTn_temp(3, 3);
    mat3d nTn;
    nTn_temp = n*n.transpose(); //dyadic product of n and n 
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            nTn(i, j) = nTn_temp(i, j);
        }
    }
    
    matrix dndt(3, 1); //derivatives of n wrt theta   
    dndt[0][0] = -sin(theta(0, 0));
    dndt[1][0] = cos(theta(0, 0));
    dndt[2][0] = 0.0;
    
    matrix ddnddt(3, 1);
    ddnddt= n*(-1.0); //second derivate of n wrt theta
    
    mat3d dnTndt; //derivative of nTn wrt theta
    dnTndt.zero();
    dnTndt[0][0] = -sin(theta(0,0)*2.0 );
    dnTndt[0][1] = cos(theta(0, 0) * 2.0);
    dnTndt[1][0] = cos(theta(0, 0) * 2.0);
    dnTndt[1][1] = sin(theta(0, 0) * 2.0);

    mat3d ddnTnddt; //second derivative of nTn wrt theta   
    ddnTnddt.zero();
    ddnTnddt[0][0] = -2.0 * cos(theta(0, 0) * 2.0);
    ddnTnddt[0][1] = -2.0 * sin(theta(0, 0) * 2.0);
    ddnTnddt[1][0] = -2.0 * sin(theta(0, 0) * 2.0);
    ddnTnddt[1][1] = 2.0*cos(theta(0, 0) * 2.0);

    mat3d eye;
    eye.zero();
    eye[0][0] = 1.0;
    eye[1][1] = 1.0;
    eye[2][2] = 1.0;
    
    //Declare some temporary double, matrix and tensors
    double temp1,temp2,temp3,temp4,temp5;
    mat3d tempmat1, tempmat2, tempmat3, tempmat4, tempmat5;
    tens4d tempA, tempB, tempC, tempD, tempE;
    
    temp1 = pow(a, 1.0/3.0);
    temp2 = pow(a, -1.0/6.0);
    mat3d G = nTn * (temp1-temp2)  + eye*temp2; //Spontanious deformation tensor
    
    //temp1 = pow(m_a, -1.0/3.0);
    //temp2 = pow(m_a, 1.0/6.0);
    //matrix (3,3) Ginv = nTn * (temp1-temp2)  + eye*temp2; //inverse of G
    mat3d Ginv = G.inverse();
    
    //temp1 = pow(m_a, -2.0/3.0);
    //temp2 = pow(m_a, 2.0/6.0);
    //matrix (3,3) Ginv2 = nTn * (temp1-temp2)  + eye*temp2; //Ginv*Ginv
    
    mat3d Ginv2 = Ginv*Ginv;
    
    temp1 = pow(a, -2.0/3.0);
    temp2 = pow(a, 2.0/6.0);
    mat3d dGinv2dt = dnTndt * (temp1-temp2);
    mat3d ddGinv2ddt = ddnTnddt * (temp1-temp2);
    
    /*Note for the fourth order operation
    dyad operation of the tens4d are diffrent to the one done in ...
    matlab (ie Mikhail Itskov Book). The equivalent are

    TensorProd_AoB(A,B) = dyad1(A,B)
    TensorProd_AxB(A,B) = dyad2(A,B.transpose())
    TensorProd_AxBt(A,B) = dyad3(A, B.transpose())

    The 9*9 matrix obtaine by tensor3333_99 from 4th order is not same
    extract 9*9 used in tens4d. The order of indeces for voigt notation
    are not in same order (see tens4d.h and tens4d.hpp). Also 9*9 
    matrix function is to be first transposed for our use(since tens4d.h
    first fills the single column). we use the voight order as in tens4d.h*/
    
    //Derivative of J=detF
    mat3d dJdF = Finv.transpose() * detF;
    tempA = dyad1(Finv.transpose(),dJdF);
    tempB = dyad3(-Finv.transpose(),Finv);
    tens4d ddJddF = tempA + tempB*detF;
    
    //--------First order derivatives of free energies-----------------
    temp1 = Ginv2.dotdot(B); //doublecontraction(Ginv2,B); 
    tempmat1 =  (Finv.transpose()*temp1)*(-1.0 / 3.0);
    mat3d delasdF = (tempmat1+ Ginv2*F) * mu * pow(detF,(-2.0/3.0));
    
    //------
    matrix dfrankdGradt(3, 1);
    dfrankdGradt = Gradt * kfrank;
    
    //-----
    matrix delasdt(1, 1);
    delasdt(0, 0) = 0.5 * mu * ((Bbar * dGinv2dt).trace());
    
    //---------Second order derivate of free energies----------------
    tempA = dyad2(Ginv2,eye.transpose());
    tempB= dyad1(Finv.transpose(),(Ginv2*F))*2.0;
    temp1 = Ginv2.dotdot(B); //doublecontraction(Ginv,B);
    tempC = dyad3(-Finv.transpose(),Finv) * temp1;
    tempD = dyad1(delasdF,Finv.transpose())*(2.0/3.0);
    tens4d ddelasddF = ((tempB+tempC)*(-1.0/3.0) + tempA)*mu*pow(detF,(-2.0/3.0)) - tempD;
    
    //----
    temp1 = dGinv2dt.dotdot(B); //doublecontraction(dGinv2dt,B);
    tempmat1 = Finv.transpose()*(1.0/3.0)*temp1;
    mat3d ddelasdFdt = ((dGinv2dt*F) - tempmat1)*mu*pow(detF,(-2.0/3.0));
    
    //----
    mat3d ddelasdtdF = ddelasdFdt;
    
    //-----
    mat3d ddfrankddGradt = eye*kfrank;
    
    //-----
    matrix ddelasdtdt(1, 1);
    ddelasdtdt(0, 0) = 0.5 * mu * ((Bbar * ddGinv2ddt).trace());
    
    //--------Calculation of Stresses and its derivatives---------------
    mat3d PK1 = delasdF + dJdF*pressmix;
    Pi0 = dfrankdGradt;
    gt = delasdt*(-1.0);
    
    tens4d dPdF = ddJddF*pressmix + ddelasddF;
    mat3d dPdt = ddelasdFdt;
    mat3d dgdF = -ddelasdtdF;
    dPidGradt = ddfrankddGradt;
    dgdt = ddelasdtdt*(-1.0);
    
    //------Transforming Tensor to voigt notation
    //matrix(9,9) *dPdF_vgt = tens3333_99(dPdF);
    double dPdF_vgt_temp[9][9];
    dPdF.extract(dPdF_vgt_temp);
    matrix dPdF_vgt_temp2(9, 9);
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {

            dPdF_vgt_temp2(i, j) = dPdF_vgt_temp[i][j];
        }
    }
    dPdF_vgt = dPdF_vgt_temp2.transpose();
    
    // todo find if using &* is valid and legal
    tens33_9(PK1,PK1_vgt);
    tens33_9(dPdt,dPdt_vgt);
    tens33_9(dgdF,dgdF_vgt);
    tens33_9(dJdF,dJdF_vgt);
    
    //feLog("theta=\n");
    //for (int nr = 0; nr < 1; ++nr) {
    //    for (int nc = 0; nc < 1; ++nc) {
    //        feLog("%.15le\t", theta(nr, nc));
    //    } feLog("\n");
    //} feLog("\n*****************************************************************\n");
}

//----------------------------Function definations----------------------
//----------------------------------------------------------------------
void FELCEMaterial443::tens33_9(mat3d A, matrix &B)
{
   //Notice the voigt order is different than I used in matlab tens33_9

    B(0,0) = A(0,0);
    B(1,0) = A(1,1);
    B(2,0) = A(2,2);
    B(3,0) = A(0,1);
    B(4,0) = A(1,2);
    B(5,0) = A(2,0);
    B(6,0) = A(1,0);
    B(7,0) = A(2,1);
    B(8,0) = A(0,2); 
   
}
