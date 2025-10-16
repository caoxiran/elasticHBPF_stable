/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Momentum stabilisation for solid mechanics equations

\*---------------------------------------------------------------------------*/

#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Read stabilization parameters


const dictionary& stabilizationDict(mechanicalProperties.subDict("stabilization"));
// Lookup method
word method = word(stabilizationDict.lookup("type"));
volTensorField gradD(fvc::grad(D));
// Merge scale factor into impK field
volScalarField gammaMod("gammaStabilisation", impK);
scalarField& gammaModI = gammaMod;

// Optional: cell zone scale factors
if (method != "none")
{
    if (stabilizationDict.lookupOrDefault<Switch>("useCellZones", false))
    {
        forAll(mesh.cellZones(), czI)
        {
            cellZone& cz = mesh.cellZones()[czI];

            // Lookup the scale factor for this cell zone
            scalar czScaleFac =
                readScalar(stabilizationDict.lookup(cz.name() + "ScaleFactor"));

            forAll(cz, cI)
            {
                label cellI = cz[cI];
                gammaModI[cellI] *= czScaleFac;
            }
        }

        // Correct processor boundaries
        gammaMod.correctBoundaryConditions();
    }
    else
    {
        // Read scale factor
        scalar scaleFactor = readScalar(stabilizationDict.lookup("scaleFactor"));

        gammaMod *= scaleFactor;
    }
}

// Interpolate gamma to faces
surfaceScalarField gammaf(fvc::interpolate(gammaMod, "interpolate(impK)"));

// Specify different scale factor at material interfaces
scalar interfaceScaleFactor
(
    stabilizationDict.lookupOrDefault<scalar>("interfaceScaleFactor", 0.01)
);

const labelList& own = mesh.owner();
const labelList& nei = mesh.neighbour();
const scalarField& gammaI = impK;

#ifdef OPENFOAM_NOT_EXTEND
scalarField& gammafI = gammaf.ref();
#else
scalarField& gammafI = gammaf.ref();
#endif

forAll(gammafI, faceI)
{
    scalar gOwn = gammaI[own[faceI]];
    scalar gNei = gammaI[nei[faceI]];

    if (mag(gOwn - gNei) > SMALL)
    {
        // Arithmetric average
        gammafI[faceI] = interfaceScaleFactor*0.5*(gOwn + gNei);
    }
}

// Correct processor patches
forAll(gammaf.boundaryField(), patchI)
{
    if (D.boundaryField()[patchI].coupled())
    {
        const scalarField pif
        (
            impK.boundaryField()[patchI].patchInternalField()
        );
        const scalarField pnf
        (
            impK.boundaryField()[patchI].patchNeighbourField()
        );

#ifdef OPENFOAM_NOT_EXTEND
        scalarField& gammafP = gammaf.boundaryFieldRef()[patchI];
#else
        scalarField& gammafP = gammaf.boundaryFieldRef()[patchI];
#endif

        if (D.boundaryField()[patchI].type() == "processor")
        {
            forAll(gammafP, faceI)
            {
                scalar gOwn = pif[faceI];
                scalar gNei = pnf[faceI];

                if (mag(gOwn - gNei) > SMALL)
                {
                    gammafP[faceI] = interfaceScaleFactor*0.5*(gOwn + gNei);
                }
            }
        }
        else
        {
            gammafP = interfaceScaleFactor*0.5*(pif + pnf);
        }
    }
    else
    {
        // Set stabilisation to zero on non-coupled boundaries
#ifdef OPENFOAM_NOT_EXTEND
        gammaf.boundaryFieldRef()[patchI] = 0.0;
#else
        gammaf.boundaryFieldRef()[patchI] = 0.0;
#endif
    }
}

// Create stabilization term
volVectorField stabilizationTerm
(
    IOobject
    (
        "stabilizationTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "zero",
        impK.dimensions()*gradD.dimensions()/dimLength,
        vector::zero
    )
);

// Calculate stabilisation term
if (method == "RhieChow")
{
    stabilizationTerm =
    (
       fvc::laplacian(gammaf, D, "laplacian(DD,D)")
     - fvc::div(gammaf*(mesh.Sf() & fvc::interpolate(gradD)))
    );
}
else if (method == "JamesonSchmidtTurkel")
{
    stabilizationTerm = -fvc::laplacian
    (
        mesh.magSf(),
        fvc::laplacian(gammaf, D, "JSTinner"),
        "JSTouter"
    );
}
else if (method == "alpha")
{
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const vectorField& nI = n.internalField();
    
    // Calculate delta vectors manually (equivalent to deltaVectors function)
    tmp<surfaceVectorField> tdelta
    (
        new surfaceVectorField
        (
            IOobject
            (
               "deltaCoeffs",
               mesh.pointsInstance(),
               mesh
            ),
            mesh,
            dimLength
        )
    );
    surfaceVectorField& delta = tdelta.ref();
    
    const volVectorField& C = mesh.C();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    
    // Internal faces
    forAll(owner, faceI)
    {
        delta[faceI] = C[neighbour[faceI]] - C[owner[faceI]];
    }
    
    // Boundary faces
    forAll(delta.boundaryField(), patchI)
    {
        const labelList faceCells(mesh.boundary()[patchI].faceCells());
        
        if (delta.boundaryField()[patchI].coupled())
        {
            delta.boundaryFieldRef()[patchI] =
                C.boundaryField()[patchI].patchNeighbourField()
              - C.boundaryField()[patchI].patchInternalField();
        }
        else
        {
            delta.boundaryFieldRef()[patchI] =
                C.boundaryField()[patchI]
              - C.boundaryField()[patchI].patchInternalField();
        }
    }
    
    const vectorField& deltaI = delta.internalField();
    
    const scalarField& magSfI = mesh.magSf().internalField();
    const vectorField& DI = D.internalField();
    const tensorField& gradDI = gradD.internalField();

    vectorField& resultI = stabilizationTerm;

    forAll(gammafI, faceI)
    {
        label ownCellID = own[faceI];
        label neiCellID = nei[faceI];

        const vector& ownD = DI[ownCellID];
        const vector& neiD = DI[neiCellID];

        const tensor& ownGradD = gradDI[ownCellID];
        const tensor& neiGradD = gradDI[neiCellID];

        const vector& n = nI[faceI];
        const vector& d = deltaI[faceI];

        vector extrapOwnD = ownD + 0.5*(d & ownGradD);
        vector extrapNeiD = neiD - 0.5*(d & neiGradD);

        scalar curMagSf = magSfI[faceI];
        scalar curGamma = gammafI[faceI];

        vector faceDamping =
            curMagSf*curGamma*(extrapNeiD - extrapOwnD)/mag(n & d);

        resultI[ownCellID] += faceDamping;
        resultI[neiCellID] -= faceDamping;
    }

    if (Pstream::parRun())
    {
        notImplemented
        (
            "Parallel boundaries for alpha scheme have to be implemented"
        );
    }

    // Divide by the volume
    resultI /= mesh.V();
}
else if (method == "Laplacian")
{
    stabilizationTerm = fvc::laplacian(gammaf, D);
}
else if (method != "none")
{
    FatalErrorIn("stabilization.H")
        << "Unknown method = " << method << nl
        << "Methods are: none, RhieChow, JamesonSchmidtTurkel, alpha and "
        << "Laplacian" <<  abort(FatalError);
}

// ************************************************************************* //
