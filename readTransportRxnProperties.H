// Rxn Fields //
Info<< "Reading diffusivity 'Df'\n" << endl;
dimensionedScalar Df
(
    "Df",
    transportProperties.lookup("Df")
);

Info<< "Reading diffusivity tourtuosity coefficient n\n" << endl;
dimensionedScalar n
(
    "n",
    transportProperties.lookup("n")
);

Info<< "Define if we want to solve Convection Eqn.\n" << endl;
scalar SolveConvectionEqn(transportProperties.lookupOrDefault<scalar>("SolveConvectionEqn",1));
SolveConvectionEqn =  SolveConvectionEqn;

Info<< "Reaction Rate Constant\n" << endl;

dimensionedScalar DTtPA //5e-11 m^2 s^-1
(
    transportProperties.lookup("DTtPA")
);

dimensionedScalar DTPLG //5e-11 m^2 s^-1
(
    transportProperties.lookup("DTPLG")
);

dimensionedScalar DTPLS //5e-11 m^2 s^-1
(
    transportProperties.lookup("DTPLS")
);

dimensionedScalar DTFBG //5e-11 m^2 s^-1
(
    transportProperties.lookup("DTFBG")
);

dimensionedScalar DTAP //5e-11 m^2 s^-1
(
    transportProperties.lookup("DTAP")
);

dimensionedScalar DTMG //5e-11 m^2 s^-1
(
    transportProperties.lookup("DTMG")
);

dimensionedScalar DTPAI //5e-11 m^2 s^-1
(
    transportProperties.lookup("DTPAI")
);

dimensionedScalar DTPLS_AP //5e-11 m^2 s^-1
(
    transportProperties.lookup("DTPLS_AP")
);

dimensionedScalar ntot0 //1.74 muM
(
    transportProperties.lookup("ntot0")
);


//adsorption and desorption of fibrinolytic proteins

dimensionedScalar katPA // 0.01 1/(muM s)
(
    transportProperties.lookup("katPA")
);

dimensionedScalar kdtPA // 0.0058 1/s
(
    transportProperties.lookup("kdtPA")
);

dimensionedScalar kaPLG // 0.1 1/(muM s)
(
    transportProperties.lookup("kaPLG")
);

dimensionedScalar kdPLG // 3.8 1/s
(
    transportProperties.lookup("kdPLG")
);

dimensionedScalar kaPLS // 0.1 1/(muM s)
(
    transportProperties.lookup("kaPLS")
);

dimensionedScalar kdPLS // 0.05 1/s
(
    transportProperties.lookup("kdPLS")
);

//lysis related parameters

dimensionedScalar k_MM // 0.3 1/s
(
    transportProperties.lookup("k_MM")
);

dimensionedScalar gamma // 0.1
(
    transportProperties.lookup("gamma")
);

dimensionedScalar KM // 0.16 muM
(
    transportProperties.lookup("KM")
);

dimensionedScalar kdeg // 2.178 1/s
(
    transportProperties.lookup("kdeg")
);

//Kinetic parameters
/*********Reaction 1. PLG-->PLS by tPA, Michaelis Menten****************/
dimensionedScalar k1f // 10 1/(muM s)
(
    transportProperties.lookup("k1f")
);

dimensionedScalar k1r // 280 1/s
(
    transportProperties.lookup("k1r")
);

dimensionedScalar k2 // 0.3 1/s
(
    transportProperties.lookup("k2")
);
//------------------------------------------------------------
dimensionedScalar k_PLG_cat("k_PLG_cat",k2);
dimensionedScalar KM_PLG("KM_PLG",(k1r+k2)/k1f);
//-------------------------------------------------------

/******Reaction 2. PLS inactivation by AP, second order Michaelis-Menten****************/
dimensionedScalar k_AP_r // 0.0021 1/s
(
    transportProperties.lookup("k_AP_r")
);

dimensionedScalar k_AP_f // 10 1/(muM s)
(
    transportProperties.lookup("k_AP_f")
);

dimensionedScalar k_AP_cat // 0.004 1/s
(
    transportProperties.lookup("k_AP_cat")
);

/*********Reaction 3. FBG-->FDP by PLS, Michaelis-Menten************/
dimensionedScalar k5f // 10 1/(muM s)
(
    transportProperties.lookup("k5f")
);

dimensionedScalar k5r // 300 1/s
(
    transportProperties.lookup("k5r")
);

dimensionedScalar k6 // 25 or 250 1/s
(
    transportProperties.lookup("k6")
);
//-------------------------------------------------------------
dimensionedScalar k_FBG_cat("k_FBG_cat",k6);
dimensionedScalar KM_FBG("KM_FBG",(k5r+k6)/k5f);
//------------------------------------------------------------

/******Reaction 4. PLS inactivation by MG, second order Michaelis-Menten****************/
dimensionedScalar kMG // 0.35 1/(muM s)
(
    transportProperties.lookup("kMG")
);

/******Reaction 5. tPA inhibiation by PAI, second order Michaelis-Menten****************/
dimensionedScalar kPAI // 37 1/(muM s)
(
    transportProperties.lookup("kPAI")
);


dimensionedScalar ELcrit //0.995
(
    transportProperties.lookup("ELcrit")
);


