//Update Solid Speed from Displacement Field
Us = fvc::ddt(D);
Us.correctBoundaryConditions();
phis = fvc::interpolate(Us) & mesh.Sf();  
Dsup = D*epss;
krxn = -epss0*kdeg*gamma*nPLS_F/ntot0;
forAll (krxn, celli)
{   
    if (Solid[celli] == 0)
    {
        krxn[celli] = 0;
    }
}

Info<< "\nCalculating epss" << endl;
fvScalarMatrix epssEqn
(
	fvm::ddt(epss)
      //+ fvm::div(phis,epss)
      - krxn //units of 1/(seconds * control Vol)
);

epssEqn.solve();
epss.correctBoundaryConditions();
/*epss.max(0.);
epss.min(0.9999);
epsf = 1- epss;
epsf.max(0.);
epsf.min(0.9999);*/

//#include "updateSolidState.H"
