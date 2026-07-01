/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 DHI
    Copyright (C) 2017 OpenCFD Ltd.
    Copyright (C) 2018 Johan Roenby
    Copyright (C) 2024 Dr. Meng Zhang
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

Application
    compressibleInterIsoLptFoam

Description
    Solver derived from interFoam for two compressible, isothermal immiscible
    fluids using the isoAdvector phase-fraction based interface capturing
    approach, with optional mesh motion and mesh topology changes including
    adaptive re-meshing.


    The solver integrates a compressible isoAdvector scheme with Lagrangian particle 
    tracking (LPT) to establish a VoF–LPT framework tailored for spray atomization. 
    The algorithm continuously interrogates the VoF field; micro-scale liquid structures 
    are instantaneously converted into Lagrangian droplets once prescribed geometric criteria
    (size and sphericity) are satisfied. Consequently, diminutive droplets are efficiently
    represented as discrete particles, eliminating the prohibitive computational cost 
    associated with fully resolved VoF meshes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"  
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"
#include "isoAdvection.H"
#include "twoPhaseMixtureThermo.H"


#include "dropletCloud.H"
#include "phaseCoupling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUf.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
        
    while (runTime.run())
    {    

        #include "readDyMControls.H"

        // Store divU and divUp from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));

        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"


        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
                
       // Store the particle positions
        cloud.storeGlobalPositions();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                if (isA<dynamicRefineFvMesh>(mesh))
                {
                    advector.surf().reconstruct();
                }
                

                mesh.update();

                if (mesh.changing())
               {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (isA<dynamicRefineFvMesh>(mesh))
                    {
                        advector.surf().mapAlphaField();
                        alpha2 = 1.0 - alpha1;
                        alpha2.correctBoundaryConditions();
                        rho == alpha1*rho1 + alpha2*rho2;
                        rho.correctBoundaryConditions();
                        rho.oldTime() = rho;
                        alpha2.oldTime() = alpha2;
                    }

                    MRF.update();

                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

                }

                if ((mesh.changing() && correctPhi))
                {
                    // Calculate absolute flux from the mapped surface velocity
                    phi = mesh.Sf() & Uf;

                    #include "correctPhi.H"

                    // Make the fluxes relative to the mesh motion
                    fvc::makeRelative(phi, U);

                    mixture.correct();
                }
                


                if (mesh.changing() && checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
            
           

            #include "alphaControls.H"
            #include "compressibleAlphaEqnSubCycle.H"

            turbulence.correctPhasePhi();
            

         
            #include "UEqn.H"
            volScalarField divUp("divUp", fvc::div(fvc::absolute(phi, U), p));
            #include "TEqn.H"   
                     
            // --- Pressure corrector loop
         
            while (pimple.correct())
            {
             
                #include "pEqn.H"
            }
           
            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }

             // Inject droplets
            coupling.update();

            // Move droplets
            cloud.move();


        runTime.write();

        runTime.printExecutionTime(Info);
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
