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
#include "FELCEMaterialPoint.h"



//-----------------------------------------------------------------------------
FELCEMaterialPoint::FELCEMaterialPoint(FEMaterialPointData* mp) : FEMaterialPointData(mp)
{
    
	//Added by Rabin on 15.08.2023
    m_theta0 = 0.0;
	m_theta = m_theta0;
	m_Avp.unit();
	m_Av = m_Avp;
    m_Fvp.unit();
	m_Fv = m_Fvp;
    
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FELCEMaterialPoint::Copy()
{
	FELCEMaterialPoint* pt = new FELCEMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FELCEMaterialPoint::Init()
{
    
    m_theta0 = 0.0;
    m_theta = m_theta0;
	m_Avp.unit();
	m_Av = m_Avp;
    m_Fvp.unit();
    m_Fv = m_Fvp;
    
	// don't forget to initialize the base class
	FEMaterialPointData::Init();
}
//-----------------------------------------------------------------------------
//! Update material point data.
void FELCEMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    m_Avp = m_Av;
    m_Fvp = m_Fv;
    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}
//-----------------------------------------------------------------------------
void FELCEMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	
	ar & m_theta0 & m_theta & m_Av & m_Avp & m_Fv & m_Fvp;
}

//-----------------------------------------------------------------------------
//! calculates the Initial Spontanious Deformation Gradient

mat3d FELCEMaterialPoint::InitialSpontaniousDG(double a0) const
{
	//get the initial director angle at the nodes
	double theta0 = m_theta0;
    // Get the material properties a0
	
    
    mat3d eye;
    eye.unit();

    //-------------
    matrix n0(3, 1); //intial director vector at current gauss point
    n0(0, 0) = cos(theta0);
    n0(1, 0) = sin(theta0);
    n0(2, 0) = 0.0;

    matrix n0Tn0_temp(3, 3);
    mat3d n0Tn0;
    n0Tn0_temp = n0 * n0.transpose(); //dyadic product of n0 and n0 

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            n0Tn0(i, j) = n0Tn0_temp(i, j);
        }
    }

    double temp1 = pow(a0, 1.0 / 3.0);
    double temp2 = pow(a0, -1.0 / 6.0);
    mat3d G0 = n0Tn0 * (temp1 - temp2) + eye * temp2; //Initial Spontanious deformation tensor
    
    return G0;
}

//-----------------------------------------------------------------------------
//! calculates the current Spontanious Deformation Gradient

mat3d FELCEMaterialPoint::SpontaniousDG(double a) const
{
    //get the initial director angle at the nodes
    double theta = m_theta;
    // Get the material properties a0
    

    mat3d eye;
    eye.unit();

    //-------------
    matrix n(3, 1); //intial director vector at current gauss point
    n(0, 0) = cos(theta);
    n(1, 0) = sin(theta);
    n(2, 0) = 0.0;

    matrix nTn_temp(3, 3);
    mat3d nTn;
    nTn_temp = n * n.transpose(); //dyadic product of n0 and n0 

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            nTn(i, j) = nTn_temp(i, j);
        }
    }

    double temp1 = pow(a, 1.0 / 3.0);
    double temp2 = pow(a, -1.0 / 6.0);
    mat3d G = nTn * (temp1 - temp2) + eye * temp2; //current Spontanious deformation tensor

    return G;
}
//-----------------------------------------------------------------------------



