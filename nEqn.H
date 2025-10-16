//solve ntPA_F 
 fvScalarMatrix ntPA_FEqn
    (
         fvm::ddt(ntPA_F)
       - tPA_Fsource
     );
// ntPA_FEqn.relax();
 fvOptions.constrain(ntPA_FEqn);
 ntPA_FEqn.solve();
 fvOptions.correct(ntPA_F);

//solve nPLG_F 
 fvScalarMatrix nPLG_FEqn
    (
         fvm::ddt(nPLG_F)
       - PLG_Fsource
     );
 //nPLG_FEqn.relax();
 fvOptions.constrain(nPLG_FEqn);
 nPLG_FEqn.solve();
 fvOptions.correct(nPLG_F);

//solve nPLS_F 
 fvScalarMatrix nPLS_FEqn
    (
         fvm::ddt(nPLS_F)
       - PLS_Fsource
     );
 //nPLS_FEqn.relax();
 fvOptions.constrain(nPLS_FEqn);
 nPLS_FEqn.solve();
 fvOptions.correct(nPLS_F);

//solve nPLS_Fw 
 fvScalarMatrix nPLS_FwEqn
    (
         fvm::ddt(nPLS_Fw)
       - PLS_Fwsource
     );
 //nPLS_FwEqn.relax();
 fvOptions.constrain(nPLS_FwEqn);
 nPLS_FwEqn.solve();
 fvOptions.correct(nPLS_Fw);

//solve ntot 
 fvScalarMatrix ntotEqn
    (
         fvm::ddt(ntot)
       - ntotsource
     );
// ntotEqn.relax();
 fvOptions.constrain(ntotEqn);
 ntotEqn.solve();
 fvOptions.correct(ntot);

nfree=ntot-ntPA_F-nPLG_F-nPLS_F;
EL=1-ntot/ntot0;
//epss=epss0*(1-EL);
#include "epssEqn.H"
epsf=1-epss;
epsf=min(epsf,0.99999);


// Step 5: Update kclot
kclot = 4.0 * pow(Rf0, 2.0) / (70 * pow(1.0 - epsf, 1.5) * (1.0 + 52 * pow(1.0 - epsf, 1.5)));

forAll(EL, celli)
{
    if (EL[celli] >= 0.8)
    {
        Solid[celli] = 0; // Set isplasma to 0 for Solid == 1
    }
    // Optional: handle other cases
    // else {isplasma[celli] = 0;}
}

forAll(EL, celli)
{
    if (EL[celli] >= 0.8)
    {
        isplasma[celli] = 1; // Set isplasma to 0 for Solid == 1
    }
    // Optional: handle other cases
    // else {isplasma[celli] = 0;}
}

forAll(mesh.boundary(), patchi) 
{
	forAll(Solid.boundaryField()[patchi],facei)
	{
		if(EL.boundaryField()[patchi][facei]<0.8) 
		{Solid.boundaryFieldRef()[patchi][facei]=1; }//ifplasma.boundaryFieldRef()[patchi][facei]=1;}
		else 
		{Solid.boundaryFieldRef()[patchi][facei]=0; }//ifplasma.boundaryFieldRef()[patchi][facei]=0;}
	}
	forAll(isplasma.boundaryField()[patchi],facei)
	{
		if(EL.boundaryField()[patchi][facei]<0.8) 
		{isplasma.boundaryFieldRef()[patchi][facei]=0; }//ifplasma.boundaryFieldRef()[patchi][facei]=1;}
		else 
		{isplasma.boundaryFieldRef()[patchi][facei]=1; }//ifplasma.boundaryFieldRef()[patchi][facei]=0;}
	}
}

Solidf = fvc::interpolate(Solid);
isplasmaf = fvc::interpolate(isplasma);
//Eliminating intermediate values of Solid Indicator Function
forAll(Solidf,facei)
{
	if(Solidf[facei]<1) // if <1 ==0.  
	{	
	Solidf[facei]=0; 
	}
}

forAll(isplasmaf,facei)
{
	if(isplasmaf[facei]>0) // if <1 ==0.  
	{	
	isplasmaf[facei]=1; 
	}
}


