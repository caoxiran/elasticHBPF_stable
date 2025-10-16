#include <fstream>
// Calculating the Drag Coefficient
// Drag is initialized as nu / kclot
Drag = nuf / kclot;
kinv = k0inv*(epss*epss)/(max((1-epss)*(1-epss)*(1-epss),SMALL));

// Loop through all cells in isplasma
forAll(isplasma, celli)
{
    if (Solid[celli] == 1)
    {
        isplasma[celli] = 0; // Set isplasma to 0 for Solid == 1
    }
    // Optional: handle other cases
    // else {isplasma[celli] = 0;}
}

// Loop through Drag and set Drag to 0 for cells where Solid == 0
forAll(Drag, celli)
{
    if (Solid[celli] == 0)
    {
        Drag[celli] = 0;
    }
}

// Loop through other variables and reset values where Solid == 0
forAll(ntot, celli)
{
    if (Solid[celli] == 0)
    {
        ntot[celli] = 0;
    }
}

forAll(nPLS_Fw, celli)
{
    if (Solid[celli] == 0)
    {
        nPLS_Fw[celli] = 0;
    }
}

forAll(nPLS_F, celli)
{
    if (Solid[celli] == 0)
    {
        nPLS_F[celli] = 0;
    }
}

forAll(ntPA_F, celli)
{
    if (Solid[celli] == 0)
    {
        ntPA_F[celli] = 0;
    }
}

forAll(nPLG_F, celli)
{
    if (Solid[celli] == 0)
    {
        nPLG_F[celli] = 0;
    }
}

forAll(EL, celli)
{
    if (Solid[celli] == 0)
    {
        EL[celli] = 1;
    }
}

forAll(epss, celli)
{
    if (Solid[celli] == 0)
    {
        epss[celli] = 1e-6;
    }
}
