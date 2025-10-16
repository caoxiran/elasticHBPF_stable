//Creating Fields for Solid Mechanics 
#include "createControlsSolid.H"
#include "readMechanicalProperties.H"
#include "readThermalProperties.H"

Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field superficial displacement Dsup\n" << endl;
volVectorField Dsup
(
    IOobject
    (
        "Dsup",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    epss*D
);

Info<< "Reading Solid Velocity Field\n" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::ddt(D)
);

#include "createSurfaceFluxes.H"

autoPtr<volScalarField> Tptr(nullptr);

if (thermalStress)
{
    Info<< "Reading field T\n" << endl;
    Tptr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );
}

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimPressure/dimDensity, symmTensor::zero)
);

/*Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

mesh.setFluxRequired(D.name());*/

// Reading Solid Mechanics activators
scalar activateSolidMechanics (transportProperties.lookupOrDefault<scalar>("activateSolidMechanics",0.));
activateSolidMechanics = activateSolidMechanics; //prevents unused variable warning

scalar pCoeff (transportProperties.lookupOrDefault<scalar>("pCoeff",1.));
pCoeff=pCoeff;

// Material model will be created in solidMechanics.H
Info<< "Material model parameters will be read by viscoNeoHookeanElastic\n" << endl;  
