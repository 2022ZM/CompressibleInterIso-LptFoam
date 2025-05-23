/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "dropletCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
     
namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<droplet>, 0);
}


bool Foam::droplet::move
(
    dropletCloud& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;
    
    
    while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
    {
        if (debug)
        {
            Info<< "Time = " << mesh().time().timeName()
                << " trackTime = " << trackTime
                << " steptFraction() = " << stepFraction() << endl;
        }
    
     
        const scalar sfrac = stepFraction();

        const scalar f = 1 - stepFraction();
        trackToAndHitFace(f*trackTime*U_, f, cloud, td);

        const scalar dt = (stepFraction() - sfrac)*trackTime;

        const tetIndices tetIs = this->currentTetIndices();
        scalar rhoc = td.rhoInterp().interpolate(this->coordinates(), tetIs);
        vector Uc = td.UInterp().interpolate(this->coordinates(), tetIs);
        scalar nuc = td.muInterp().interpolate(this->coordinates(), tetIs) / rhoc;
    
       
        scalar rhop = cloud.rhop();  
        scalar magUr = mag(Uc - U_);
        label cellI = cell();
        scalar ReFunc = 1.0;
        scalar Re = magUr*d_/nuc;


        if (Re > 0.01)
        {
            ReFunc += 0.15*pow(Re, 0.687);
        }

        scalar Dc = (24.0*nuc/d_)*ReFunc*(3.0/4.0)*(rhoc/(d_*rhop));
        
        // calculate momentum
        scalar m = rhop*(4.0/3.0)*constant::mathematical::pi*pow(d_/2.0, 3.0);
        vector i1 = U_*m;

        U_ = (U_ + dt*(Dc*Uc + (1.0 - rhoc/rhop)*td.g()))/(1.0 + dt*Dc);


        
        // update momentum source for coupled simulation
        cloud.source()[cellI] += U_*m - i1;

        // Check if particle passed given faceZones
        if (face() > -1)
        {
            const faceZoneMesh& fzm = mesh().faceZones();
            
            const labelList fzIDs = cloud.faceZoneIDs();
            
            forAll(fzIDs, i)
            {
                const faceZone& fz = fzm[fzIDs[i]];

                forAll(fz, j)
                {                
                    if (fz[j] == face())
                    {
                        cloud.dataDiameter()[Pstream::myProcNo()][i].append(d_);
                        cloud.dataPosition()[Pstream::myProcNo()][i].append(position());
                        cloud.dataNParticle()[Pstream::myProcNo()][i].append(nParticle_);
                    }
                }
            }
        }
    }
    return td.keepParticle;
}


bool Foam::droplet::hitPatch(dropletCloud&, trackingData&)
{
    return false;
}


void Foam::droplet::hitProcessorPatch
(
    dropletCloud&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::droplet::hitWallPatch(dropletCloud& cloud, trackingData&)
{
    const vector nw = normal();

    scalar Un = U_ & nw;
    vector Ut = U_ - Un*nw;

    if (Un > 0)
    {
        U_ -= (1.0 + 0.97)*Un*nw;
    }

    U_ -= 0.09*Ut;
}



void Foam::droplet::transformProperties(const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}


void Foam::droplet::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


// ************************************************************************* //
