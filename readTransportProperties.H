Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Reading Model Parameters //
Info<< "Reading the fluid's kinematic viscocity 'nuf'\n" << endl;
dimensionedScalar nuf
(
    "nuf",
    transportProperties.lookup("nuf")
);

Info<< "Reading density of fluid 'rhof'\n" << endl;
dimensionedScalar rhof
(
    "rhof",
    transportProperties.lookup("rhof")
);
dimensionedScalar Rf0
(
    transportProperties.lookup("Rf0")
);
Info<< "Reading density of solid 'rhos'\n" << endl;
dimensionedScalar rhos
(
    "rhos",
    transportProperties.lookup("rhos")
);
// rhos is now defined in readMechanicalProperties.H as rho_s (volScalarField)

Info<< "Reading inverse of intrinsic permeability 'k0inv'\n" << endl;
dimensionedScalar k0inv
(
    "k0inv",
    transportProperties.lookup("k0inv")
);

Info<< "Reading Ferrari Drag term\n" << endl;
dimensionedScalar Fkinv
(
    "Fkinv",
    transportProperties.lookup("Fkinv")
);

Info<< "Reading Terzaghi Coeff. 'sigma'\n" << endl;
dimensionedScalar sigmacoe
(
    "sigmacoe",
    transportProperties.lookup("sigmacoe")
);

dimensionedScalar epss0 //1-eps0  0.05
(
    transportProperties.lookup("epss0")
);
dimensionedScalar epsf0("epsf0",1-epss0);
