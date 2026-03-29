/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FELCEInitialAngle.h"
#include <ctime>


//=============================================================================
BEGIN_FECORE_CLASS(FELCEInitialAngle, FEInitialCondition)
	ADD_PARAMETER(m_data, "value");  
END_FECORE_CLASS();

FELCEInitialAngle::FELCEInitialAngle(FEModel* fem) : FEInitialDOF(fem)
{
}

bool FELCEInitialAngle::Init()
{
	if (SetDOF("vart") == false) return false; //vart defination in FESolidModule.cpp
	double av_initangle = 45.0; // initial angle in degrees desired in average sense 
	double perturb = 0.5; //perturbation angle 
	double rand_t = -1.0 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2.0)));
	double rand1 = (acos(-1.0)/180.0 * av_initangle) + (perturb * acos(-1.0)/180.0 * rand_t);
	FEInitialDOF::SetValue(rand1);
	return FEInitialDOF::Init();
}
