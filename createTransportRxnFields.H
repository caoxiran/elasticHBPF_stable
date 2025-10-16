// Scalar Transport Fields //
#include "readTransportRxnProperties.H"

/*Info<< "Reading Corrected Ion Diffusivity 'Dion'\n" << endl;
volScalarField Dion
(
    IOobject
    (
        "Dion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Df*pow((1-epss),(n-1))
);

Info<< "Reading Ionic Concentration 'C'\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading Superficial Ionic Concentration 'C'\n" << endl;
volScalarField Cepsf
(
    IOobject
    (
        "Cepsf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    C*epsf
);*/

Info<< "Reading Corrected Ion Diffusivity 'Dion'\n" << endl;

volScalarField DTtPAion
(
    IOobject
    (
        "DTtPAion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DTtPA*pow((1-epss),(n-1))
);

volScalarField DTPLSion
(
    IOobject
    (
        "DTPLSion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DTPLS*pow((1-epss),(n-1))
);

volScalarField DTPLGion
(
    IOobject
    (
        "DTPLGion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DTPLG*pow((1-epss),(n-1))
);

volScalarField DTFBGion
(
    IOobject
    (
        "DTFBGion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DTFBG*pow((1-epss),(n-1))
);

volScalarField DTAPion
(
    IOobject
    (
        "DTAPion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DTAP*pow((1-epss),(n-1))
);

volScalarField DTMGion
(
    IOobject
    (
        "DTMGion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DTMG*pow((1-epss),(n-1))
);

volScalarField DTPAIion
(
    IOobject
    (
        "DTPAIion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DTPAI*pow((1-epss),(n-1))
);

volScalarField DTPLS_APion
(
    IOobject
    (
        "DTPLS_APion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DTPLS_AP*pow((1-epss),(n-1))
);


Info<< "Reading field CtPA\n" << endl;
volScalarField CtPA
(
    IOobject
    (
        "CtPA",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field CPLS\n" << endl;
volScalarField CPLS
(
    IOobject
    (
        "CPLS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field CPLG\n" << endl;
volScalarField CPLG
(
    IOobject
    (
        "CPLG",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field CFBG\n" << endl;
volScalarField CFBG
(
    IOobject
    (
        "CFBG",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field CAP\n" << endl;
volScalarField CAP
(
    IOobject
    (
        "CAP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field CMG\n" << endl;
volScalarField CMG
(
    IOobject
    (
        "CMG",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field CPAI\n" << endl;
volScalarField CPAI
(
    IOobject
    (
        "CPAI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field CPLS_AP\n" << endl;
volScalarField CPLS_AP
(
    IOobject
    (
        "CPLS_AP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field ntPA_F\n" << endl;
volScalarField ntPA_F
(
    IOobject
    (
        "ntPA_F",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field nPLG_F\n" << endl;
volScalarField nPLG_F
(
    IOobject
    (
        "nPLG_F",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field nPLS_F\n" << endl;
volScalarField nPLS_F
(
    IOobject
    (
        "nPLS_F",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field nPLS_Fw\n" << endl;
volScalarField nPLS_Fw
(
    IOobject
    (
        "nPLS_Fw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Creating field ntot from ntot0*Solid\n" << endl;
volScalarField ntot
(
    IOobject
    (
        "ntot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


//Info << "ntot dimensions: " << ntot.dimensions() << endl;
//Info << "ntPA_F dimensions: " << ntPA_F.dimensions() << endl;
//Info << "nPLG_F dimensions: " << nPLG_F.dimensions() << endl;
//Info << "nPLS_F dimensions: " << nPLS_F.dimensions() << endl;
volScalarField nfree
(
    IOobject
    (
        "nfree",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);





volScalarField kclot("kclot",4*pow(Rf0,2.)/(70*pow(1-epsf,1.5)*(1.+52*pow(1.-epsf,1.5))));
volScalarField Drag("Drag", nuf/kclot);
//Labeling the places that have porous media
forAll(epsf,celli) 
{
	if(epsf[celli]>0.9999) 
	{	
	Solid[celli]=0; 
	isplasma[celli]=1;
	}
	else
	{
	Solid[celli]=1; 
	isplasma[celli]=0;
	}
}

forAll(mesh.boundary(), patchi) 
{
	forAll(Solid.boundaryField()[patchi],facei)
	{
		if(epsf.boundaryField()[patchi][facei]>0.9999) 
		{Solid.boundaryFieldRef()[patchi][facei]=0;}
		else
		{Solid.boundaryFieldRef()[patchi][facei]=1;}
	}
}

forAll(mesh.boundary(), patchi) 
{
	forAll(isplasma.boundaryField()[patchi],facei)
	{
		if(epsf.boundaryField()[patchi][facei]>0.9999) 
		{isplasma.boundaryFieldRef()[patchi][facei]=1;}
		else
		{isplasma.boundaryFieldRef()[patchi][facei]=0;}
	}
}

surfaceScalarField Solidf(fvc::interpolate(Solid));
surfaceScalarField isplasmaf(fvc::interpolate(isplasma));
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
Info<< "Reading inverse permeability field 'kinv'\n" << endl;
volScalarField kinv
(
    IOobject
    (
        "kinv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE 
    ),
    //(1/kclot)
    (k0inv*(epss*epss)/(max((1-epss)*(1-epss)*(1-epss),SMALL)))
    //(k0inv*(pow(epss,1.5)*(1+56*pow(epss,3.0)))/(max(pow(epss0,1.5)*(1+56*pow(epss0,3.0)),dimensionedScalar("SMALL",dimless, SMALL))))
);
/*Info<< "Reading equilibrium porosity\n" << endl;
dimensionedScalar epssEq
(
    "epssEq",
    transportProperties.lookup("epssEq")
);*/

// Declare krxn as volScalarField
volScalarField krxn
(
    IOobject
    (
        "krxn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (-epss0*kdeg*gamma*nPLS_F/ntot0)
);


//volScalarField EL("EL", dimensionedScalar("EL", dimless, 1.0) - ntot / ntot0); // Extent of lysis
//volScalarField phif("phif", phif0.value() * (1 - EL));
