Info<< "\nCalculating Solid Displacement Field" << endl;

//Reading nCorrectors and stress type (thermal and/or plane or neither)
#include "readSolidDisplacementFoamControls.H"

//Initializing number of corrections and residuals
int iCorr = 0; 
scalar initialResidual = 0;

//Calculating Source Terms
Terzaghi = -(epss*rhof/rhos)*fvc::grad(p);  
volTensorField gradD(fvc::grad(D));
//Calculating Displacement
        do
        {
            if (thermalStress)
            {
                volScalarField& T = Tptr();
                solve
                (
                    fvm::ddt(T) == fvm::laplacian(DT, T)
                );
            }

    { //This is always done
        // Create viscoNeoHookeanElastic model directly
                #include "viscoNeoHookeanElastic.H"
                #include "impK.H"                
                fvVectorMatrix DEqn  //[0 2 -2 0 0 0 0]
                (
		    //transient term, D is displacement vector
                    fvm::d2dt2(epss,D) 
                 ==
		    //implicit part of the stress tensor (for numerical stability)
                    fvm::laplacian(impKf, D, "laplacian(DD,D)")
                    
                    //subtract implicit part to avoid double counting
                    - fvc::laplacian(impKf, D, "laplacian(DD,D)")
                    
                    //add full stress tensor (including hyperelastic/viscoelastic)
                    + epss*fvc::div(sigmaD, "div(sigmaD)")

		    //Terzaghi Term
		    + Terzaghi

                    //Transport Momentum Source 
                    + fvc::Sp(epsf*nuf*kinv*rhof/rhos, U)
                    - fvc::Sp(epss*epsf*nuf*kinv*rhof/rhos, Us)                  
                );

                if (mechanicalProperties.found("stabilization"))
                {
                    #include "stabilization.H"
                    DEqn += stabilizationTerm*epss;
                }
                if (thermalStress)
                {
                    const volScalarField& T = Tptr();
                    DEqn += fvc::grad(threeKalpha*T);
                }

                initialResidual = DEqn.solve().max().initialResidual();



                //Calculating Explicit Stress Tensor using ViscoNeoHookeanElastic
                // Calculate stress using the material model
                #include "correct.H"
            }

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

	//here we reconstruct the final stress tensor from everything to print
        #include "calculateStress.H"

// ************************************************************************* //
