// Convection-Diffusion Equation //
#include "SourceTerm.H"

Info<< "\nCalculating Scalar Transport " << endl;

/*fvScalarMatrix ConvectionEqn //(C=Concentration=rho*wf/mw)
(
         fvm::ddt(epsf,C) 
        + fvm::div(phi,C) 
        - fvm::laplacian((Dion*epsf),C)
        //+ fvm::div(phic,C) //solid movement correction
);

ConvectionEqn.solve();
Cepsf = C*epsf;
C.max(0.); */
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
//solve CtPA 
fvScalarMatrix CtPAEqn
    (
         fvm::ddt(CtPA)
       + fvm::div(phi, CtPA)
       - fvm::laplacian(DTtPA, CtPA)
       - R_tPA_wtClot
       - R_tPA_plasma*isplasma
     );
 CtPAEqn.relax();
 fvOptions.constrain(CtPAEqn);
 CtPAEqn.solve();
 fvOptions.correct(CtPA);
 
forAll(CtPA, celli)
{
    if (CtPA[celli] < 0)
    {
        CtPA[celli] = 0;
    }
}
//solve CPLG 
 fvScalarMatrix CPLGEqn
    (
         fvm::ddt(CPLG)
       + fvm::div(phi, CPLG)
       - fvm::laplacian(DTPLG, CPLG)
       - R_PLG_wtClot
       - R_PLG_plasma*isplasma
     );
 CPLGEqn.relax();
 fvOptions.constrain(CPLGEqn);
 CPLGEqn.solve();
 fvOptions.correct(CPLG);

forAll(CPLG, celli)
{
    if (CPLG[celli] < 0)
    {
        CPLG[celli] = 0;
    }
}
//solve CPLS 
 fvScalarMatrix CPLSEqn
    (
         fvm::ddt(CPLS)
       + fvm::div(phi, CPLS)
       - fvm::laplacian(DTPLS, CPLS)
       - R_PLS_wtClot
       - R_PLS_plasma*isplasma
     );
 CPLSEqn.relax();
 fvOptions.constrain(CPLSEqn);
 CPLSEqn.solve();
 fvOptions.correct(CPLS);

forAll(CPLS, celli)
{
    if (CPLS[celli] < 0)
    {
        CPLS[celli] = 0;
    }
}

//solve CFBG 
 fvScalarMatrix CFBGEqn
    (
         fvm::ddt(CFBG)
       + fvm::div(phi, CFBG)
       - fvm::laplacian(DTFBG, CFBG)
       - R_FBG_plasma*isplasma
     );
 CFBGEqn.relax();
 fvOptions.constrain(CFBGEqn);
 CFBGEqn.solve();
 fvOptions.correct(CFBG);
forAll(CFBG, celli)
{
    if (CFBG[celli] < 0)
    {
        CFBG[celli] = 0;
    }
}
//solve CAP 
 fvScalarMatrix CAPEqn
    (
         fvm::ddt(CAP)
       + fvm::div(phi, CAP)
       - fvm::laplacian(DTAP, CAP)
       - R_AP_plasma*isplasma
     );
 CAPEqn.relax();
 fvOptions.constrain(CAPEqn);
 CAPEqn.solve();
 fvOptions.correct(CAP);
forAll(CAP, celli)
{
    if (CAP[celli] < 0)
    {
        CAP[celli] = 0;
    }
}
//solve CMG 
 fvScalarMatrix CMGEqn
    (
         fvm::ddt(CMG)
       + fvm::div(phi, CMG)
       - fvm::laplacian(DTMG, CMG)
       - R_MG_plasma*isplasma
     );
 CMGEqn.relax();
 fvOptions.constrain(CMGEqn);
 CMGEqn.solve();
 fvOptions.correct(CMG);
forAll(CMG, celli)
{
    if (CMG[celli] < 0)
    {
        CMG[celli] = 0;
    }
}
//solve CPAI 
 fvScalarMatrix CPAIEqn
    (
         fvm::ddt(CPAI)
       + fvm::div(phi, CPAI)
       - fvm::laplacian(DTPAI, CPAI)
       - R_PAI_plasma*isplasma
     );
 CPAIEqn.relax();
 fvOptions.constrain(CPAIEqn);
 CPAIEqn.solve();
 fvOptions.correct(CPAI);
forAll(CPAI, celli)
{
    if (CPAI[celli] < 0)
    {
        CPAI[celli] = 0;
    }
}
//solve CPLS_AP 
 fvScalarMatrix CPLS_APEqn
    (
         fvm::ddt(CPLS_AP)
       + fvm::div(phi, CPLS_AP)
       - fvm::laplacian(DTPLS_AP, CPLS_AP)
       - R_PLS_AP_plasma*isplasma
     );
 CPLS_APEqn.relax();
 fvOptions.constrain(CPLS_APEqn);
 CPLS_APEqn.solve();
 fvOptions.correct(CPLS_AP);
forAll(CPLS_AP, celli)
{
    if (CPLS_AP[celli] < 0)
    {
        CPLS_AP[celli] = 0;
    }
}



#include "nEqn.H"

