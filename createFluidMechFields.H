#include "readTransportProperties.H"

// Creating Fields //

Info<< "Reading field fluid pressure field 'p'\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fluid velocity field 'U'\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading solid field 'Solid'\n" << endl;
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid", dimensionSet(0,0,0,0,0,0,0), 0) 
);

Info<< "Creating solid fraction field 'epss' from epss0*Solid\n" << endl;
volScalarField epss
(
    IOobject
    (
        "epss",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



//volVectorField SolidBorder(fvc::grad(Solid));
Info<< "Creating plasma field 'isplasma' from 1-Solid\n" << endl;
volScalarField isplasma
(
    IOobject
    (
        "isplasma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("isplasma", dimensionSet(0,0,0,0,0,0,0), 1) 
);
//volScalarField nfree("nfree",ntot-ntPA_F-nPLS_F-nPLG_F);
volScalarField EL
(
    IOobject
    (
        "EL",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("EL", dimensionSet(0,0,0,0,0,0,0), 1.0) 
);

volScalarField epsf
(
    IOobject
    (
        "epsf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (1-epss)
);
 
