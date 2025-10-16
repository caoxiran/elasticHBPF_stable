#include <cmath>

// Get displacement field
const volVectorField& D = mesh.lookupObject<volVectorField>("D"); //m

// Calculate deformation gradient F = I + grad(D)
volTensorField F = I + fvc::grad(D); //dimless

// Add numerical stability checks for F
volScalarField detF(det(F));
volScalarField detF_safe(max(detF, SMALL));

// Check for extreme deformations
scalar maxStretchValue = max(mag(F)).value();
if (maxStretchValue > 10.0)
{
    WarningIn("correct.H")
        << "Large deformation detected: max stretch = " 
        << maxStretchValue << endl;
}
    
// Calculate the Jacobian of the deformation gradient
volScalarField J(detF_safe); //dimless - use safe determinant
    
// Calculate the volume preserving right Cauchy Green tensor
volTensorField C(F.T() & F); //dimless

// Force symmetry to eliminate numerical errors that cause complex eigenvalues
// This handles tiny asymmetries like (1, 1e-6, 0, 1e-6, 1, 0, 0, 0, 1)
volTensorField C_symmetric = 0.5*(C + C.T());
C = C_symmetric;
    
// Define Fbar := J^(-1/3)*F
volTensorField Fbar(pow(J, -1.0/3.0)*F); //dimless
    
// Define Cbar := J^(-2/3)*C
volTensorField Cbar(pow(J, -2.0/3.0)*C);//dimless
    
// Take references to the internal fields for efficiency
const tensorField& CI = C.primitiveField();//dimless
const scalarField& JI = J.primitiveField();//dimless
symmTensorField& transformNeededI = transformNeeded_.primitiveFieldRef();//dimPressure/dimDensity
    
// Partial derive of WBar(CBar), internalField operation
vector eigenValues = vector::zero;//dimless
vector eigenValuesBar = vector::zero;//dimless
    
forAll(transformNeededI, cellI)
{
    // Calculate principal stretches: lambda^2, which are eigenvalues of C
    // Use identity matrix for very small deformations to avoid complex eigenvalues
    tensor C_cell = CI[cellI];
    
    // Check if deformation is very small (close to identity)
    scalar maxOffDiagonal = max(max(mag(C_cell.xy()), mag(C_cell.xz())), mag(C_cell.yz()));
    scalar maxDiagonalDeviation = max(max(mag(C_cell.xx() - 1.0), mag(C_cell.yy() - 1.0)), mag(C_cell.zz() - 1.0));
    
    if (maxOffDiagonal < 1e-8 && maxDiagonalDeviation < 1e-8)
    {
        // Use identity matrix for very small deformations
        eigenValues = vector(1.0, 1.0, 1.0);
    }
    else
    {
        // Force exact symmetry
        tensor C_sym = C_cell;
        C_sym.xy() = C_sym.yx() = 0.5*(C_cell.xy() + C_cell.yx());
        C_sym.xz() = C_sym.zx() = 0.5*(C_cell.xz() + C_cell.zx());
        C_sym.yz() = C_sym.zy() = 0.5*(C_cell.yz() + C_cell.zy());
        
        eigenValues = Foam::eigenValues(C_sym);//dimless
    }
    
    // Check for negative eigenvalues (should not happen for C = F^T * F)
    if (eigenValues.x() < 0.0) eigenValues.x() = SMALL;
    if (eigenValues.y() < 0.0) eigenValues.y() = SMALL;
    if (eigenValues.z() < 0.0) eigenValues.z() = SMALL;
        
    // Calculate modified principal stretches
    eigenValuesBar.x() = std::pow(double(JI[cellI]), -1.0/3)*std::sqrt(double(eigenValues.x()));//dimless
    eigenValuesBar.y() = std::pow(double(JI[cellI]), -1.0/3)*std::sqrt(double(eigenValues.y()));//dimless
    eigenValuesBar.z() = std::pow(double(JI[cellI]), -1.0/3)*std::sqrt(double(eigenValues.z()));//dimless
        
    transformNeededI[cellI].xx() =
        mu0_[cellI]*std::pow(eigenValuesBar.x(), double(alpha_[0] - 1))
      + mu1_[cellI]*std::pow(eigenValuesBar.x(), double(alpha_[1] - 1))
      + mu2_[cellI]*std::pow(eigenValuesBar.x(), double(alpha_[2] - 1));
        
    transformNeededI[cellI].yy() =
        mu0_[cellI]*std::pow(eigenValuesBar.y(), double(alpha_[0] - 1))
      + mu1_[cellI]*std::pow(eigenValuesBar.y(), double(alpha_[1] - 1))
      + mu2_[cellI]*std::pow(eigenValuesBar.y(), double(alpha_[2] - 1));
        
    transformNeededI[cellI].zz() =
        mu0_[cellI]*std::pow(eigenValuesBar.z(), double(alpha_[0] - 1))
      + mu1_[cellI]*std::pow(eigenValuesBar.z(), double(alpha_[1] - 1))
      + mu2_[cellI]*std::pow(eigenValuesBar.z(), double(alpha_[2] - 1));
}
    
    // Partial derive of WBar(CBar), boundaryField operation
    forAll(C.boundaryField(), patchI)
    {
        // Take references to the patch fields for efficiency
        const tensorField& CP = C.boundaryField()[patchI];
        const scalarField& JP = J.boundaryField()[patchI];
        symmTensorField& transformNeededP = transformNeeded_.boundaryFieldRef()[patchI]; //dimPressure/dimDensity
        
        forAll(CP, faceI)
        {
            // Calculate principal stretches: lambda^2, which are eigenvalues of C
            // Use identity matrix for very small deformations to avoid complex eigenvalues
            tensor C_face = CP[faceI];
            
            // Check if deformation is very small (close to identity)
            scalar maxOffDiagonal = max(max(mag(C_face.xy()), mag(C_face.xz())), mag(C_face.yz()));
            scalar maxDiagonalDeviation = max(max(mag(C_face.xx() - 1.0), mag(C_face.yy() - 1.0)), mag(C_face.zz() - 1.0));
            
            if (maxOffDiagonal < 1e-8 && maxDiagonalDeviation < 1e-8)
            {
                // Use identity matrix for very small deformations
                eigenValues = vector(1.0, 1.0, 1.0);
            }
            else
            {
                // Force exact symmetry
                tensor C_sym = C_face;
                C_sym.xy() = C_sym.yx() = 0.5*(C_face.xy() + C_face.yx());
                C_sym.xz() = C_sym.zx() = 0.5*(C_face.xz() + C_face.zx());
                C_sym.yz() = C_sym.zy() = 0.5*(C_face.yz() + C_face.zy());
                
                eigenValues = Foam::eigenValues(C_sym);
            }
            
            // Check for negative eigenvalues (should not happen for C = F^T * F)
            if (eigenValues.x() < 0.0) eigenValues.x() = SMALL;
            if (eigenValues.y() < 0.0) eigenValues.y() = SMALL;
            if (eigenValues.z() < 0.0) eigenValues.z() = SMALL;
            
            // Calculate modified principal stretches
            eigenValuesBar.x() = std::pow(double(JP[faceI]), -1.0/3)*std::sqrt(double(eigenValues.x()));
            eigenValuesBar.y() = std::pow(double(JP[faceI]), -1.0/3)*std::sqrt(double(eigenValues.y()));
            eigenValuesBar.z() = std::pow(double(JP[faceI]), -1.0/3)*std::sqrt(double(eigenValues.z()));
            
            transformNeededP[faceI].xx() = 
                mu0_[faceI]*std::pow(eigenValuesBar.x(), double(alpha_[0] - 1))
              + mu1_[faceI]*std::pow(eigenValuesBar.x(), double(alpha_[1] - 1))
              + mu2_[faceI]*std::pow(eigenValuesBar.x(), double(alpha_[2] - 1));
            
            transformNeededP[faceI].yy() =
                mu0_[faceI]*std::pow(eigenValuesBar.y(), double(alpha_[0] - 1))
              + mu1_[faceI]*std::pow(eigenValuesBar.y(), double(alpha_[1] - 1))
              + mu2_[faceI]*std::pow(eigenValuesBar.y(), double(alpha_[2] - 1));
            
            transformNeededP[faceI].zz() =
                mu0_[faceI]*std::pow(eigenValuesBar.z(), double(alpha_[0] - 1))
              + mu1_[faceI]*std::pow(eigenValuesBar.z(), double(alpha_[1] - 1))
              + mu2_[faceI]*std::pow(eigenValuesBar.z(), double(alpha_[2] - 1));
        }
    }
    
    // 2 * FBar * (Partial Cbar of Partial WBar) * FBar.T()
    // Manual transform: F & S & F.T()
    transformFbar_ = 2*symm(Fbar & transformNeeded_ & Fbar.T());//dimPressure/dimDensity
    
    // dev operation
    Info << "sigmaD dimensions: " << sigmaD.dimensions() << endl;
    Info << "dev(transformFbar_) dimensions: " << dev(transformFbar_)().dimensions() << endl;
    sigmaD = dev(transformFbar_);//dimPressure/dimDensity
    
    // Calculate initial stress
    volTensorField invFbar(inv(Fbar));
    // Manual transform: invFbar & sigmaD & invFbar.T()
    s_ = symm(invFbar & sigmaD & invFbar.T());//dimPressure/dimDensity
    
    // Update internal stress variables, representing stress relaxations for
    // each Maxwell model
    scalar deltaT = mesh.time().deltaTValue();
    
    forAll(H_, MaxwellModelI)
    {
        scalar tau_safe = max(tau_[MaxwellModelI], SMALL);
        H_[MaxwellModelI] =
            Foam::exp(-deltaT/tau_safe)*h_[MaxwellModelI].oldTime()
          - Foam::exp(-deltaT/(2.0*tau_safe))*s_.oldTime(); //dimPressure/dimDensity
    }
    
    forAll(h_, MaxwellModelI)
    {
        scalar tau_safe = max(tau_[MaxwellModelI], SMALL);
        h_[MaxwellModelI] =
            H_[MaxwellModelI]
          + Foam::exp(-deltaT/(2.0*tau_safe))*s_; //dimPressure/dimDensity
    }
    
    // Calculate the current total stress, where the volumetric term is
    // elastic and the deviatoric term is viscoelastic
    scalar gRelax = gammaInf_;
    
    forAll(gamma_, MaxwellModelI)
    {
        scalar tau_safe = max(tau_[MaxwellModelI], SMALL);
        gRelax += gamma_[MaxwellModelI]*Foam::exp(-deltaT/(2*tau_safe));
    }
    
    // Calculate hydrostatic pressure, defined in (G. A. HOLZAPFEL,1996)
    // Add numerical stability checks
    volScalarField J_safe(max(J, SMALL));
    dimensionedScalar beta_safe("beta_safe", beta_.dimensions(), max(beta_.value(), SMALL));
    
    Info << "K_ dimensions: " << K_.dimensions() << endl;
    Info << "rho_s dimensions: " << rho_s.dimensions() << endl;
    volScalarField pressure
    (
        K_*(1.0/(beta_safe*J_safe))*(1 - pow(J_safe, -beta_safe))/rho_s //dimPressure/dimDensity
    );
    Info << "pressure dimensions: " << pressure.dimensions() << endl;
    
    Info << "J*pressure*I dimensions: " << (J*pressure*I)().dimensions() << endl;
    Info << "gRelax*sigmaD dimensions: " << (gRelax*sigmaD)().dimensions() << endl;
    sigmaD += J*pressure*I + gRelax*sigmaD; //dimPressure/dimDensity
    
    forAll(H_, MaxwellModelI)
    {
        // Manual transform: Fbar & H & Fbar.T()
        transformH_[MaxwellModelI] = symm(Fbar & H_[MaxwellModelI] & Fbar.T());//dimPressure/dimDensity
    }
    
    forAll(transformH_, MaxwellModelI)
    {
        Info << "gamma_[MaxwellModelI]*dev(transformH_[MaxwellModelI]) dimensions: " << (gamma_[MaxwellModelI]*dev(transformH_[MaxwellModelI]))().dimensions() << endl;
        sigmaD += gamma_[MaxwellModelI]*dev(transformH_[MaxwellModelI]);//dimPressure/dimDensity
    }


//void Foam::viscoNeoHookeanElastic::correct(surfaceSymmTensorField& sigmaD)
//{
    // Simplified surface implementation
 //   sigmaD = dimensionedSymmTensor("zero", dimPressure, symmTensor::zero);
//}
