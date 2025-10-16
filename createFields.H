#include "createFluidMechFields.H"
#include "createPhi.H"
#include "createSolidMechFields.H"
#include "createTransportRxnFields.H"
#include "createFvOptions.H"
//Initializing Terzaghi Stress 
volVectorField Terzaghi= -(sigmacoe*rhof/rho_s)*epss*fvc::grad(p);

//Pressure reference cells
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

//Reading Max Courant Number for Porous Media
scalar maxCoUs = runTime.controlDict().lookupOrDefault<scalar>("maxCoUs", 1.0);
maxCoUs=maxCoUs;

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
#include "createMRF.H"


