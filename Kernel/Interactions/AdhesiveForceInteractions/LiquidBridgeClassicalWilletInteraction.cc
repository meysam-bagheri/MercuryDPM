//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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


#include "LiquidBridgeClassicalWilletInteraction.h"
#include "Species/AdhesiveForceSpecies/LiquidBridgeClassicalWilletSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
LiquidBridgeClassicalWilletInteraction::LiquidBridgeClassicalWilletInteraction(BaseInteractable* P, BaseInteractable* I,
                                                             unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    wasInContact_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeClassicalWilletInteraction::LiquidBridgeClassicalWilletInteraction() finished"<<std::endl;
#endif
}

//used for mpi
LiquidBridgeClassicalWilletInteraction::LiquidBridgeClassicalWilletInteraction()
        : BaseInteraction()
{
    wasInContact_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeClassicalWilletInteraction::LiquidBridgeClassicalWilletInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
LiquidBridgeClassicalWilletInteraction::LiquidBridgeClassicalWilletInteraction(const LiquidBridgeClassicalWilletInteraction& p)
        : BaseInteraction(p)
{
    wasInContact_ = p.wasInContact_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeClassicalWilletInteraction::LiquidBridgeClassicalWilletInteraction(const LiquidBridgeClassicalWilletInteraction &p finished"<<std::endl;
#endif
}

/*!
 * 
 */
LiquidBridgeClassicalWilletInteraction::~LiquidBridgeClassicalWilletInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LiquidBridgeClassicalWilletInteraction::~LiquidBridgeClassicalWilletInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void LiquidBridgeClassicalWilletInteraction::write(std::ostream& os) const
{
    os << " wasInContact " << wasInContact_;
}

/*!
 * \param[in,out] is
 */
void LiquidBridgeClassicalWilletInteraction::read(std::istream& is)
{
    // we read in wasInContact_ like this because an early version did not initialize it.
    std::string dummy;
    is >> dummy >> wasInContact_;
}

/*!
 * 
 */
// full model from Apendix Eq. A1, Willet et al. 2000, Langmuir
void LiquidBridgeClassicalWilletInteraction::computeAdhesionForce()
{
    const LiquidBridgeClassicalWilletSpecies* species = getSpecies();
    if (getOverlap() >= 0)
    {
        wasInContact_ = true;
        Mdouble effectiveRadius = 2.0 * getEffectiveRadius();
        Mdouble V_star = species->getLiquidBridgeVolume() / std::pow(effectiveRadius,3);

        Mdouble f1 = (-0.44507 + 0.050832*species->getContactAngle() - 1.1466*std::pow(species->getContactAngle(),2)) + \
        (-0.1119 - 0.000411*species->getContactAngle() - 0.1490*std::pow(species->getContactAngle(),2))*std::log(V_star) + \
        (-0.012101 - 0.0036456*species->getContactAngle() - 0.01255*std::pow(species->getContactAngle(),2))*std::pow(std::log(V_star),2) + \
        (-0.0005 - 0.0003505*species->getContactAngle() - 0.00029076*std::pow(species->getContactAngle(),2))*std::pow(std::log(V_star),3);

        Mdouble F_star = std::exp(f1);

        Mdouble fdotn = -2.0 * constants::pi * effectiveRadius * species->getSurfaceTension() * F_star;
        addForce(getNormal() * fdotn);
    }
    else if (wasInContact_)
    {
        Mdouble effectiveRadius = 2.0 * getEffectiveRadius();
        Mdouble s_c = -getOverlap() * std::sqrt(effectiveRadius / species->getLiquidBridgeVolume())/2;
        Mdouble V_star = species->getLiquidBridgeVolume() / std::pow(effectiveRadius,3);

        Mdouble f1 = (-0.44507 + 0.050832*species->getContactAngle() - 1.1466*std::pow(species->getContactAngle(),2)) + \
        (-0.1119 - 0.000411*species->getContactAngle() - 0.1490*std::pow(species->getContactAngle(),2))*std::log(V_star) + \
        (-0.012101 - 0.0036456*species->getContactAngle() - 0.01255*std::pow(species->getContactAngle(),2))*std::pow(std::log(V_star),2) + \
        (-0.0005 - 0.0003505*species->getContactAngle() - 0.00029076*std::pow(species->getContactAngle(),2))*std::pow(std::log(V_star),3);

        Mdouble f2 = (1.9222 - 0.57473*species->getContactAngle() - 1.2918*std::pow(species->getContactAngle(),2)) + \
        (-0.0668 - 0.1201*species->getContactAngle() - 0.22574*std::pow(species->getContactAngle(),2))*std::log(V_star) + \
        (-0.0013375 - 0.0068988*species->getContactAngle() - 0.01137*std::pow(species->getContactAngle(),2))*std::pow(std::log(V_star),2);

        Mdouble f3 = (1.268 - 0.01396*species->getContactAngle() - 0.23566*std::pow(species->getContactAngle(),2)) + \
        (0.198 +0.092*species->getContactAngle() - 0.06418*std::pow(species->getContactAngle(),2))*log(V_star) + \
        (0.02232 + 0.02238*species->getContactAngle() - 0.009853*std::pow(species->getContactAngle(),2))*std::pow(std::log(V_star),2) + \
        (0.0008585 + 0.001318*species->getContactAngle() - 0.00053*std::pow(species->getContactAngle(),2))*std::pow(std::log(V_star),3);

        Mdouble f4 = (-0.010703 + 0.073776*species->getContactAngle() - 0.34742*std::pow(species->getContactAngle(),2)) + \
        (0.03345 + 0.04543*species->getContactAngle() - 0.09056*std::pow(species->getContactAngle(),2))*std::log(V_star) + \
        (0.0018574 + 0.004456*species->getContactAngle() - 0.006257*std::pow(species->getContactAngle(),2))*std::pow(std::log(V_star),2);

        Mdouble F_star = std::exp(f1 - f2*std::exp(f3*std::log(s_c)+f4*std::pow(std::log(s_c),2)));

        Mdouble fdotn = -2.0 * constants::pi * effectiveRadius * species->getSurfaceTension() * F_star;
        addForce(getNormal() * fdotn);
    }
}

/*!
 * \return Mdouble
 */
Mdouble LiquidBridgeClassicalWilletInteraction::getElasticEnergy() const
{
    ///\todo TW
    return 0.0;
}

/*!
 * \return const LiquidBridgeClassicalWilletSpecies*
 */
const LiquidBridgeClassicalWilletSpecies* LiquidBridgeClassicalWilletInteraction::getSpecies() const
{
    return static_cast<const LiquidBridgeClassicalWilletSpecies*>(getBaseSpecies()->getAdhesiveForce());
;
}

/*!
 * \return std::string
 */
std::string LiquidBridgeClassicalWilletInteraction::getBaseName() const
{
    return "LiquidBridgeClassicalWillet";
}

bool LiquidBridgeClassicalWilletInteraction::getWasInContact() const
{
    return wasInContact_;
}

void LiquidBridgeClassicalWilletInteraction::setWasInContact(bool wasInContact)
{
    wasInContact_ = wasInContact;
}
