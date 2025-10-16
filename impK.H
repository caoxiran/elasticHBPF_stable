   volScalarField impK
   (
       IOobject
       (
           "impK",
           runTime.timeName(),
           mesh,
           IOobject::READ_IF_PRESENT,
           IOobject::NO_WRITE
       ),
       mesh,
       dimensionedScalar("impK", dimPressure/dimDensity, 0.0),
       "zeroGradient"
   );
    if (min(nu_).value() == 0.5 || mesh.foundObject<volScalarField>("p"))
    {
        impK = 2.0*mu_/rho_s*epss;
    }
    else
    {
        impK = (2.0*mu_ + (4.0/3.0)*K_)/rho_s*epss;    
    }
    
    // Add numerical stability check
    impK = max(impK, dimensionedScalar("impK_min", dimPressure/dimDensity, SMALL));
    surfaceScalarField impKf("impKf", fvc::interpolate(impK));
