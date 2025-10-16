//Reactions in the clot
//fibrinolytic reactions within the clot
volScalarField r1("r1", katPA*CtPA*nfree-kdtPA*ntPA_F);
volScalarField r2("r2", kaPLG*CPLG*nfree-kdPLG*nPLG_F);
volScalarField r3("r3", kaPLS*CPLS*nfree-kdPLS*nPLS_F);
volScalarField r4("r4", k_MM*ntPA_F*nPLG_F/(KM*(1-epsf)+nPLG_F));
volScalarField r5("r5", kdeg*gamma*nPLS_F);
volScalarField r6("r6", kdPLS*nPLS_Fw);
forAll(r4, celli)
{
    if (Solid[celli] == 0)
    {
        r4[celli] = 0;
    }
}
/**************Bounde phase species reaction rates***************/
volScalarField tPA_Fsource("tPA_Fsource", r1);
volScalarField PLG_Fsource("PLG_Fsource", r2-r4);
volScalarField PLS_Fsource("PLS_Fsource", r3+r4-r5);
volScalarField PLS_Fwsource("PLS_Fwsource", r5-r6);
volScalarField ntotsource("ntotsource", -r5);
/***************Free phase species reaction rates**************/
//tPA, PLG, and PLS react with the fibre
//AP, FBG, MG, PAI and intermediate complex do not react with the fibre
volScalarField R_tPA_wtClot("R_tPA_wtClot", -r1/epsf);
volScalarField R_PLG_wtClot("R_PLG_wtClot", -r2/epsf);
volScalarField R_PLS_wtClot("R_PLS_wtClot", (-r3 -  r6)/epsf);

//Reactions in the plasma
//fibrinolytic reactions in PLASMA
// 1: tPA + PLG <--> tPA*PLG --> PLS + tPA
volScalarField R1 ("R1", k_PLG_cat*CtPA*CPLG/(KM_PLG+CPLG));
// 2: PLS + FBG <--> PLS*FBG --> PLS + FDP
volScalarField R2 ("R2", k_FBG_cat*CPLS*CFBG/(KM_FBG+CFBG));
// 3: PLS + AP <--> PLS*AP
volScalarField R3 ("R3", k_AP_f*CPLS*CAP-k_AP_r*CPLS_AP);
// 4: PLS*AP --> inactive
volScalarField R4 ("R4", k_AP_cat*CPLS_AP);
// 5: PLS + MG --> inactive
volScalarField R5 ("R5", kMG*CPLS*CMG);
// 6: tPA + PAI --> inactive
volScalarField R6 ("R6", kPAI*CtPA*CPAI);

//Only free phase species are involved in the reactions, no bound phase
//tPA, PLG and PLS reactions occur in both the clot and plasma.
volScalarField R_tPA_plasma ("R_tPA_plasma", -R6);
volScalarField R_PLG_plasma ("R_PLG_plasma", -R1);
volScalarField R_PLS_plasma ("R_PLS_plasma", R1 - R3 - R5);
volScalarField R_FBG_plasma ("R_FBG_plasma", -R2);
volScalarField R_AP_plasma ("R_AP_plasma", -R3);
volScalarField R_MG_plasma ("R_MG_plasma", -R5);
volScalarField R_PAI_plasma ("R_PAI_plasma", -R6);
volScalarField R_PLS_AP_plasma ("R_PLS_AP_plasma", R3 - R4);

/*****************************************************/
/*scalar r1Sum = 0.0;
scalar r2Sum = 0.0;
scalar r3Sum = 0.0;
scalar r4Sum = 0.0;
scalar r5Sum = 0.0;
scalar r6Sum = 0.0;
scalar R_tPA_wtClotSum = 0.0;
scalar R_PLG_wtClotSum = 0.0;
scalar R_PLS_wtClotSum = 0.0;
label countSolidCells = 0; // Counts the number of cells where Solid == 1*/
forAll(isplasma, celli)
{
    if (Solid[celli] == 1)
    {
        isplasma[celli] = 0; // Set isplasma to 0 for Solid == 1
        //kclotSum += kclot[celli]; // Accumulate kclot values
       /* r1Sum += r1[celli];
        r2Sum += r2[celli];
        r3Sum += r3[celli];
        r4Sum += r4[celli];
        r5Sum += r5[celli];
        r6Sum += r6[celli];
        R_tPA_wtClotSum += R_tPA_wtClot[celli];
        R_PLG_wtClotSum += R_PLG_wtClot[celli];
        R_PLS_wtClotSum += R_PLS_wtClot[celli];        countSolidCells++; // Increment the counter for Solid == 1 cells*/
    }
    // Optional: handle other cases
    // else {isplasma[celli] = 0;}
}
/*scalar r1Average = (countSolidCells > 0) ? r1Sum / countSolidCells : 0.0;
scalar r2Average = (countSolidCells > 0) ? r2Sum / countSolidCells : 0.0;
scalar r3Average = (countSolidCells > 0) ? r3Sum / countSolidCells : 0.0;
scalar r4Average = (countSolidCells > 0) ? r4Sum / countSolidCells : 0.0;
scalar r5Average = (countSolidCells > 0) ? r5Sum / countSolidCells : 0.0;
scalar r6Average = (countSolidCells > 0) ? r6Sum / countSolidCells : 0.0;
scalar R_tPA_wtClotAverage = (countSolidCells > 0) ? R_tPA_wtClotSum / countSolidCells : 0.0;
scalar R_PLG_wtClotAverage = (countSolidCells > 0) ? R_PLG_wtClotSum / countSolidCells : 0.0;
scalar R_PLS_wtClotAverage = (countSolidCells > 0) ? R_PLS_wtClotSum / countSolidCells : 0.0;
Info << "Average r3 in Solid == 1 region: " << r3Average << endl;
Info << "Average r5 in Solid == 1 region: " << r5Average << endl;
Info << "Average r6 in Solid == 1 region: " << r6Average << endl;
Info << "Average R_PLS_wtClot in Solid == 1 region: " << R_PLS_wtClotAverage << endl;

std::ofstream outFile1("r3Average.txt", std::ios::app);
if (outFile1.is_open())
{
    // 
    outFile1 << "Time: " << runTime.timeName() 
            << ", Average r3 in Solid == 1 region: " << r3Average << std::endl;
    outFile1.close(); 
}
else
{
    FatalErrorInFunction
        << "Unable to open file for appending data!"
        << exit(FatalError);
}

std::ofstream outFile2("r5Average.txt", std::ios::app);
if (outFile2.is_open())
{
    // 
    outFile2 << "Time: " << runTime.timeName() 
            << ", Average r5 in Solid == 1 region: " << r5Average << std::endl;
    outFile2.close(); 
}
else
{
    FatalErrorInFunction
        << "Unable to open file for appending data!"
        << exit(FatalError);
}

std::ofstream outFile3("r6Average.txt", std::ios::app);
if (outFile3.is_open())
{
    // 
    outFile3 << "Time: " << runTime.timeName() 
            << ", Average r6 in Solid == 1 region: " << r6Average << std::endl;
    outFile3.close(); 
}
else
{
    FatalErrorInFunction
        << "Unable to open file for appending data!"
        << exit(FatalError);
}

std::ofstream outFile4("R_PLS_wtClotAverage.txt", std::ios::app);
if (outFile4.is_open())
{
    // 
    outFile4 << "Time: " << runTime.timeName() 
            << ", Average R_PLS_wtClot in Solid == 1 region: " << R_PLS_wtClotAverage << std::endl;
    outFile4.close(); 
}
else
{
    FatalErrorInFunction
        << "Unable to open file for appending data!"
        << exit(FatalError);
}*/
// Open the CSV file for appending
/*std::ofstream outFile1("averages.csv", std::ios::app);

// Check if the file is open
if (outFile1.is_open())
{
    // If it's the first time writing, write headers
    if (runTime.timeIndex() == 0)
    {
        outFile1 << "Time, Average_r3, Average_r5, Average_r6, Average_R_PLS_wtClot" << std::endl;
    }

    // Write the data for the current time step
    outFile1 << runTime.timeName() << ", "
         << r1Average << ", "
         << r2Average << ", "
         << r3Average << ", "
         << r4Average << ", "
         << r5Average << ", "
         << r6Average << ", "
         << R_tPA_wtClotAverage << ", "
         << R_PLG_wtClotAverage << ", "
         << R_PLS_wtClotAverage << std::endl;

    outFile1.close(); // Close the file after writing
}
else
{
    FatalErrorInFunction
        << "Unable to open file for appending data!"
        << exit(FatalError);
}
*/

