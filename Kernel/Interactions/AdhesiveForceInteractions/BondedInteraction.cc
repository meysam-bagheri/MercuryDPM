//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "BondedInteraction.h"
#include "Species/AdhesiveForceSpecies/BondedSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
BondedInteraction::BondedInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    bonded_ = false;
    bondForce_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"BondedInteraction::BondedInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
BondedInteraction::BondedInteraction(const BondedInteraction& p)
        : BaseInteraction(p)
{
    ///\todo tw check if the parameters are valid when inserting the species into the handler
    bonded_ = p.bonded_;
    bondForce_ = p.bondForce_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"BondedInteraction::BondedInteraction(const BondedInteraction &p finished"<<std::endl;
#endif
}

BondedInteraction::BondedInteraction()
{
#ifdef MERCURYDPM_USE_MPI
    logger(FATAL,"ChargedBondedInteractions are currently not implemented in parallel MercuryDPM");
#endif
}

/*!
 *
 */
BondedInteraction::~BondedInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"BondedInteraction::~BondedInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void BondedInteraction::write(std::ostream& os) const
{
    os << " bonded " << bonded_;
    os << " bondForce " << bondForce_;
}

/*!
 * \param[in,out] is
 */
void BondedInteraction::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> bonded_;
    is >> dummy >> bondForce_;
}

/*!
 * \details Uses the most basic adhesion contact model.
 */
void BondedInteraction::computeAdhesionForce()
{
    if (bonded_ && getOverlap() >= 0)
    {
        const BondedSpecies* species = getSpecies();
        addForce(getNormal() * (-getBondForce()
                                - species->getBondDissipation() * getNormalRelativeVelocity()));
    }
}

/*!
 * \return const pointer of BondedSpecies*
 */
const BondedSpecies* BondedInteraction::getSpecies() const
{
    return static_cast<const BondedSpecies*>(getBaseSpecies()->getAdhesiveForce());
;
}

/*!
 * \return std::string
 */
std::string BondedInteraction::getBaseName() const
{
    return "Bonded";
}

/*!
 * \details Elastic (=potential) energy is defined as the energy gained by separating two interactables.
 * As it costs energy to separate adhesive interactables, the elastic energy is negative.
 * \return the elastic energy stored in the adhesive spring. 
 */
Mdouble BondedInteraction::getElasticEnergy() const
{
    if (!bonded_)
        return 0.0;
    else
        return -getBondForce() * getOverlap();
}

bool BondedInteraction::getBonded() const
{
    return bonded_;
}

void BondedInteraction::setBonded(bool bonded)
{
    bonded_ = bonded;
}


void BondedInteraction::bond()
{
    bondForce_=getSpecies()->getBondForceMax();
    bonded_ = true;
}

void BondedInteraction::bondInPlace()
{
    bondForce_= getForce().getLength();
    bonded_ = true;
}

void BondedInteraction::unbond()
{
    bonded_ = false;
}
