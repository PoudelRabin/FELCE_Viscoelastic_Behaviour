//This file contains the value of the quantities at the material point for the LCE material

#pragma once
#include <FECore/FEMaterialPoint.h>
#include <FECore/tens4d.h>
#include <FECore/FEMaterialPointProperty.h>
#include "febiomech_api.h"



//-----------------------------------------------------------------------------
//! This class defines material point data for elastic materials.
class FEBIOMECH_API FELCEMaterialPoint : public FEMaterialPointData
{
public:
    //! constructor
    FELCEMaterialPoint(FEMaterialPointData* mp = nullptr);
    
    //! Initialize material point data
    void Init() ;

    //! Update material point data.
    void Update(const FETimeInfo& timeInfo);

    //! create a shallow copy
    FEMaterialPointData* Copy() ;

    //! serialize material point data
    void Serialize(DumpStream& ar) ;

public:
    //below functions take material property stretch parameter as an argument
    mat3d InitialSpontaniousDG(double a0) const;
    mat3d SpontaniousDG(double a) const;


public:
    //inital director angle
    double m_theta0;
    //current director angle
    double m_theta;
    //current Internal variable
    mat3d m_Av;  //! nonequilibrium deformation gradient A = AeAv
	mat3d m_Fv; // nonequilibrium deformation gradient F = FeFv
    //previous internal variable
    mat3d m_Avp;
	mat3d m_Fvp;
    
};
