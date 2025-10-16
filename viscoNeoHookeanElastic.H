/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    viscoNeoHookeanElastic

Description
    Neo-visco-elastic Hookean mechanical law based on Simo & Hughes (2000):
    
        tau_ = J*p*I + g(deltaT)*tau_+hBar

    The material is represented as a generalised Maxwell model (Prony series)
    i.e. a spring in parallel with multiple Maxwell models.

Author
    Based on solids4foam implementation, adapted for OpenFOAM-7

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

//- List of internal variables for each Maxwell model
PtrList<volSymmTensorField> h_;

//- List of internal variables for each Maxwell model (surface)
PtrList<surfaceSymmTensorField> hf_;

//- List of intermediate variables for each Maxwell model
PtrList<volSymmTensorField> H_;

//- List of intermediate variables for each Maxwell model (surface)
PtrList<surfaceSymmTensorField> Hf_;

//- List of transformed variables for each Maxwell model
PtrList<volSymmTensorField> transformH_;

//- List of transformed variables for each Maxwell model (surface)
PtrList<surfaceSymmTensorField> transformHf_;

//- Current deviatoric stress
volSymmTensorField s_
(
    IOobject
    (
        "s",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimPressure/dimDensity, symmTensor::zero)
);

//- Current deviatoric stress (surface)
surfaceSymmTensorField sf_
(
    IOobject
    (
        "sf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero",  dimPressure/dimDensity, symmTensor::zero)
);

//- Transform needed field
volSymmTensorField transformNeeded_
(
    IOobject
    (
        "transformNeeded",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("I", dimPressure/dimDensity, I)
);

//- Transform needed field (surface)
surfaceSymmTensorField transformNeededf_
(
    IOobject
    (
        "transformNeededf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("I",  dimPressure/dimDensity, I)
);

//- Transform Fbar field
volSymmTensorField transformFbar_
(
    IOobject
    (
        "transformFbar",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("I", dimPressure/dimDensity, I)
);

//- Transform Fbar field (surface)
surfaceSymmTensorField transformFbarf_
(
    IOobject
    (
        "transformFbarf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("I", dimPressure/dimDensity, I)
);

// Initialize viscoNeoHookeanElastic model
{
// Check for physical Poisson's ratio
if (min(nu_).value() < -1.0 || max(nu_).value() > 0.5)
{
    FatalErrorIn
    (
        "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
        "(\n"
        "    const word& name,\n"
        "    const fvMesh& mesh,\n"
        "    dictionary& dict\n"
        ")"
    )   << "Unphysical Poisson's ratio: nu should be >= -1.0 and <= 0.5"
        << abort(FatalError);
}

// Check for incompressibility
if (min(nu_).value() == 0.5)
{
    Info<< "Material is incompressible: make sure to use a hybrid"
        << " approach solid model" << endl;

    // Set lambda and k to GREAT
    lambda_ = dimensionedScalar("lambda", dimPressure, GREAT);
    K_ = dimensionedScalar("K", dimPressure, GREAT);

     //Note: dict variable not available in this scope
     if (mechanicalProperties.found("K"))
     {
         FatalErrorIn
         (
                 "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
                 "(\n"
                 "    const word& name,\n"
                 "    const fvMesh& mesh,\n"
                 "    dictionary& dict\n"
                 ")"
          )   << "We cannot specify k and nu independently as well as "
              << "Young's/shear modulii!" << abort(FatalError);

             // PC: what is going on here?
          K_ = dimensionedScalar(mechanicalProperties.lookup("K"));
      }
}
    else
    {
        if (planeStress)
        {
            lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - nu_));
        }
        else
        {
            lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - 2.0*nu_));
        }

        K_ = lambda_ + (2.0/3.0)*mu_;

        // Note: dict variable not available in this scope
        // if (dict.found("k"))
        // {
        //     WarningIn
        //     (
        //         "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
        //         "(\n"
        //         "    const word& name,\n"
        //         "    const fvMesh& mesh,\n"
        //         "    dictionary& dict\n"
        //         ")"
        //     )   << "k is directly looked up from the dictionary; "
        //         << "it is not calculated from mu and lambda" << endl;

        //     // PC: what is going on here?
        //     K_ = dimensionedScalar(dict.lookup("k"));
        // }
    }

    // Check E_ and tau_ are the same length
    if (E_.size() != tau_.size())
    {
        FatalErrorIn
        (
            "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    dictionary& dict\n"
            ")"
        )   << "The E and tau lists should have the same length!"
            << abort(FatalError);
    }

    // Calculate relative modulii
    // E_ is dimensionless, EInf_ has pressure dimensions
    // Calculate total modulus as scalar first
    scalar E0_scalar = max(EInf_).value();
    forAll(E_, i)
    {
        E0_scalar += E_[i];
    }

    gammaInf_ = max(EInf_).value() / E0_scalar;

    forAll(gamma_, i)
    {
        gamma_[i] = E_[i] / E0_scalar;
    }

    // Check all the relaxation times are positive
    if (min(tau_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    dictionary& dict\n"
            ")"
        )   << "All relaxation times should be positive!"
            << abort(FatalError);
    }

    // Check all the E values are positive
    if (min(E_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    dictionary& dict\n"
            ")"
        )   << "All values of stiffness E should be positive!"
            << abort(FatalError);
    }

    // Print out the relative module

    Info<< "Relative modulii" << nl
        << "    gammaInfinity: " << gammaInf_ << nl;

    forAll(gamma_, i)
    {
        Info<< "    gamma[" << i << "] : " << gamma_[i] << nl;
    }

    Info<< endl;

    transformH_.setSize(gamma_.size());
    transformHf_.setSize(gamma_.size());

    forAll(transformH_, MaxwellModelI)
    {
        transformH_.set
        (
            MaxwellModelI,
            new volSymmTensorField
            (
                IOobject
                (
                    "transformH" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero",  dimPressure/dimDensity, symmTensor::zero)
            )
        );

        transformHf_.set
        (
            MaxwellModelI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "transformHf" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero",  dimPressure/dimDensity, symmTensor::zero)
            )
        );
    }


    // Create the internal stress variables for each Maxwell model

    h_.setSize(gamma_.size());
    hf_.setSize(gamma_.size());

    forAll(h_, MaxwellModelI)
    {
        h_.set
        (
            MaxwellModelI,
            new volSymmTensorField
            (
                IOobject
                (
                    "h" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero",  dimPressure/dimDensity, symmTensor::zero)
            )
        );

        hf_.set
        (
            MaxwellModelI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "hf" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero",  dimPressure/dimDensity, symmTensor::zero)
            )
        );

        // We need to store the old time field
        h_[MaxwellModelI].storeOldTime();
        hf_[MaxwellModelI].storeOldTime();
    }


    H_.setSize(gamma_.size());
    Hf_.setSize(gamma_.size());

    forAll(H_, MaxwellModelI)
    {
        H_.set
        (
            MaxwellModelI,
            new volSymmTensorField
            (
                IOobject
                (
                    "H" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero",  dimPressure/dimDensity, symmTensor::zero)
            )
        );

        Hf_.set
        (
            MaxwellModelI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "Hf" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero",  dimPressure/dimDensity, symmTensor::zero)
            )
        );
    }

    // Store the old time s field
    s_.storeOldTime();
    sf_.storeOldTime();

} // End of viscoNeoHookeanElastic initialization


