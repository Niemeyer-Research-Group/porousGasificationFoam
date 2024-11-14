/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    medium.H

Description
    template for porous medium shape
\*---------------------------------------------------------------------------*/


// volVectorField containig coordinates of cell centers
const volVectorField& Ccells = porosityF.mesh().C();

double h = 0.01;
//double h = 0.0025;
double dp = h/20;
scalar ep = 1e-6;

// loop over cells containing logical conditions for porosity shape and value
forAll(Ccells, celli)
{
  porosityF[celli] = 0.875;
// Example: const value in whole space.
   //porosityF[celli]=0.46;
   //porosityF[0]=0.46;

// Example: a half plane x < 0.
    //if ((Ccells[celli].x() < 0.1) && (Ccells[celli].x() > -0.1) && (Ccells[celli].y() < 0.1)
    //    && (Ccells[celli].y() > -0.1) && (Ccells[celli].z() < 0.0) && (Ccells[celli].z() > -0.1))
    //{
    //    porosityF[celli] = 0.6;
    //}
    //else
    //{
      //  porosityF[celli] = 1.0;
    //}
  
  if ((Ccells[celli].y() > (h-3.0*dp)-ep) && (Ccells[celli].y() < h))
  {
	porosityF[celli] = (1.0-0.875)/(3.0*dp)*(Ccells[celli].y()-h)+1.0;
  }
  //else if (Ccells[celli].y() < (h-3.0*dp))
  //{
	//porosityF[celli] = 0.875;
  //}
  //else
  //{
	//porosityF[celli] = 1.0;
  //}
  
  if (Ccells[celli].y() > h){
	  porosityF[celli] = 1.0;
  }
  
  if (Ccells[celli].y() < 0){
	  porosityF[celli] = .875;
  }
  
}



// Base vectors of the local coordinate system in which porosity tensor is diagonal
// v1 - x, v2 - y, v3 = v1xv2


// Porosity tensor in local coordinates (chosen for porosity tensor to be diagonal).
tensor D_loc(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
//tensor anisotropyK_loc(0, 0, 0, 0, 0, 0, 0, 0, 0);
scalar vv = 1/(Foam::pow((1.0-Foam::pow((1.0-.875),(1.0/3))),3.0)*(1.0+Foam::pow((1.0-.875),(1.0/3)))/(11.4*(1.0-.875))*Foam::pow(dp,2.0));
forAll(Ccells,celli){
    //D_loc.xx() = 16e8;
    //D_loc.yy() = 16e8;
    //D_loc.zz() = 16e8;
	if ((porosityF[celli] < 1.0) && (porosityF[celli] > 0.0)){
	if (Ccells[celli].y() < h)
	{
		double a = 1/(Foam::pow((1.0-Foam::pow((1.0-porosityF[celli]),(1.0/3))),3.0)*(1.0+Foam::pow((1.0-porosityF[celli]),(1.0/3)))/(11.4*(1.0-porosityF[celli]))*Foam::pow(dp,2.0))*1.5;
		D_loc.xx() = a;
		D_loc.yy() = a;
		D_loc.zz() = a;
		Df[celli] = (D_loc);
		//Df[celli] = tensor (a, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, a);
	}
	
    //anisotropyK_loc.xx() = 1;
    //anisotropyK_loc.yy() = 1;
    //anisotropyK_loc.zz() = 1;

    // Transformation of porosity tensor to normal base in cartesian coordinates.
    //Df[celli] = (T_loc & D_loc & T_loc.T());
	
    //anisotropyK[celli] = (T_loc & anisotropyK_loc & T_loc.T());
	//anisotropyK[celli] = anisotropyK_loc;
}
}

//const fvPatchScalarField& tbp = porosityF.boundaryField()[porosityF.mesh().boundaryMesh().findPatchID("top")];
//forAll(tbp, celli){
//	tbp[celli] = 1.0;
//}

porosityF.correctBoundaryConditions();
Df.correctBoundaryConditions();
