Info<< "Reading mechanical properties\n" << endl;

IOdictionary mechanicalProperties
(
    IOobject
    (
        "mechanicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//- Solid's density
Info<< "Reading the solid's density rho_s" << endl;
dimensionedScalar rho_s_value(mechanicalProperties.lookup("rho_s")); //dimMass/dimVolume
volScalarField rho_s
(
    IOobject
    (
        "rho_s",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rho_s_value
);

//- List of modulii of elasticity for each Maxwell model
//Info<< "Reading the Elastic modulus rhoE" << endl;
//List<dimensionedScalar> rhoE_(mechanicalProperties.lookup("rhoE")); //dimMass/dimLength/sqr(dimTime)

//- List of normalized modulii of elasticity for each Maxwell model
//Info<< "Normalising E : E/rho\n" << endl;
//List<dimensionedScalar> E_("E",rhoE_/rho_s.value());

//- List of modulii of elasticity for each Maxwell model
Info<< " Reading modulii of elasticity: E\n" << endl;
List<scalar> E_(mechanicalProperties.lookup("E"));

//- Poission ratio
Info<< "Reading the Poission ratio nu" << endl;
dimensionedScalar nu_value(mechanicalProperties.lookup("nu")); //dimless
volScalarField nu_
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    nu_value
);

//- infinite Elastic modulii
Info<< "Reading the Elastic modulii EInfinity" << endl;
dimensionedScalar EInf_value(mechanicalProperties.lookup("EInfinity")); //dimPressure
volScalarField EInf_
(
    IOobject
    (
        "EInf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    EInf_value
);


//- Relaxed (long term) second Lame parameter (shear modulus)
volScalarField mu_
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    EInf_/(2.0*(1.0 + nu_))
);

//- Relaxed (long term) first Lame parameter
volScalarField lambda_
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    nu_*EInf_/((1.0 + nu_)*(1.0 - 2.0*nu_))
);

//
volScalarField threeK_
(
    IOobject
    (
        "threeK",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    EInf_/(1.0 - 2.0*nu_)
);

//- List of relaxation times for each Maxwell model
Info<< "Reading the relaxation times for each Maxwell model tau" << endl;
List<scalar> tau_(mechanicalProperties.lookup("relaxationTimes")); //dimTime

//- Relaxed (long term) nondimensional relative modulus
Info<< "Reading the nondimensional relative modulus gammaInf" << endl;
scalar gammaInf_(0.0); //dimless

//- List of nondimensional relative modulii
Info<< "Reading the nondimensional relative modulii gamma" << endl;
List<scalar> gamma_(E_.size(), 0.0); //dimless

//- Shear modulii in Ogden model
Info<< "Reading the hyperelastic model parameter mu0, mu1, mu2" << endl;
dimensionedScalar mu0_E(mechanicalProperties.lookup("mu0")); //dimPressure
dimensionedScalar mu1_E(mechanicalProperties.lookup("mu1")); //dimPressure
dimensionedScalar mu2_E(mechanicalProperties.lookup("mu2")); //dimPressure

volScalarField mu0_(mu0_E/rho_s); //dimPressure/dimMass*DimLength^3
volScalarField mu1_(mu1_E/rho_s); //dimPressure/dimMass*DimLength^3
volScalarField mu2_(mu2_E/rho_s); //dimPressure/dimMass*DimLength^3
//- Dimensionless material constants
Info<< "Reading the hyperelastic model parameter alpha" << endl;
List<scalar> alpha_(mechanicalProperties.lookup("alpha")); //dimless

//- An empirical coefficient
Info<< "Reading the hyperelastic model parameter beta" << endl;
dimensionedScalar beta_(mechanicalProperties.lookup("beta")); //dimless

//- Bulk modulus
/*Info<< "Reading the hyperelastic model parameter K" << endl;
dimensionedScalar K_value(mechanicalProperties.lookup("K")); //dimPressure
volScalarField K_
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    K_value
);*/

Switch planeStress(mechanicalProperties.lookup("planeStress"));


if (planeStress)
{
    Info<< "Plane Stress\n" << endl;

    lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - nu_));
    threeK_ = EInf_/(1.0 - nu_);
}
else
{
    lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    Info<< "Plane Strain\n" << endl;
}

volScalarField K_
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    lambda_ + (2.0/3.0)*mu_
);





