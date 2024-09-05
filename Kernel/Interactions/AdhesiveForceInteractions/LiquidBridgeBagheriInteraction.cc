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


#include "LiquidBridgeBagheriInteraction.h"
#include "Species/AdhesiveForceSpecies/LiquidBridgeBagheriSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
LiquidBridgeBagheriInteraction::LiquidBridgeBagheriInteraction(BaseInteractable* P, BaseInteractable* I,
                                                             unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    wasInContact_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeBagheriInteraction::LiquidBridgeBagheriInteraction() finished"<<std::endl;
#endif
}

//used for mpi
LiquidBridgeBagheriInteraction::LiquidBridgeBagheriInteraction()
        : BaseInteraction()
{
    wasInContact_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeBagheriInteraction::LiquidBridgeBagheriInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
LiquidBridgeBagheriInteraction::LiquidBridgeBagheriInteraction(const LiquidBridgeBagheriInteraction& p)
        : BaseInteraction(p)
{
    wasInContact_ = p.wasInContact_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeBagheriInteraction::LiquidBridgeBagheriInteraction(const LiquidBridgeBagheriInteraction &p finished"<<std::endl;
#endif
}

/*!
 * 
 */
LiquidBridgeBagheriInteraction::~LiquidBridgeBagheriInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LiquidBridgeBagheriInteraction::~LiquidBridgeBagheriInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void LiquidBridgeBagheriInteraction::write(std::ostream& os) const
{
    os << " wasInContact " << wasInContact_;
}

/*!
 * \param[in,out] is
 */
void LiquidBridgeBagheriInteraction::read(std::istream& is)
{
    // we read in wasInContact_ like this because an early version did not initialize it.
    std::string dummy;
    is >> dummy >> wasInContact_;
}

/*!
 * 
 */
void LiquidBridgeBagheriInteraction::computeAdhesionForce()
{
    const LiquidBridgeBagheriSpecies* species = getSpecies();
    if (getOverlap() >= 0)
    {
        wasInContact_ = true;
        Mdouble effectiveRadius = 2 * getEffectiveRadius();
        Mdouble V = species->getLiquidBridgeVolume()/std::pow(effectiveRadius,3);
        
        // full model according to Eq. 24-26 from Bagheri et al. 2024, CPM
        Mdouble f0 = 2*constants::pi*(1 - 0.3823*std::pow(V,0.2586));
        // std::cout << species->getMeanRadius() << std::endl;
        Mdouble a = 0.4158*std::pow(V,0.2835)+0.6474;
        Mdouble b = -0.2087*std::pow(V,0.3113)+2.267;
        Mdouble ft = (1 - a*std::pow(std::sin(species->getContactAngle()),b))*f0;

        Mdouble force_full_fit = ft;
        
        Mdouble fdotn = -force_full_fit * effectiveRadius * species->getSurfaceTension();
        
        addForce(getNormal() * fdotn);
    }
    else if (wasInContact_)
    {
        Mdouble effectiveRadius = 2 * getEffectiveRadius();
        Mdouble V = species->getLiquidBridgeVolume()/std::pow(effectiveRadius,3);
        Mdouble S_rupture = (1.0 + species->getContactAngle()*0.5)*(std::pow(V,1.0/3.0)+0.1*std::pow(V,2.0/3.0));        
        Mdouble S = -getOverlap()/(effectiveRadius*S_rupture);
        Mdouble log_V = std::log(V);

        // full model according to Eq. 27-30 from Bagheri et al. 2024, CPM
        Mdouble f0 = 2*constants::pi*(1 - 0.3823*std::pow(V,0.2586));

        Mdouble a = 0.4158*std::pow(V,0.2835)+0.6474;
        Mdouble b = -0.2087*std::pow(V,0.3113)+2.267;
        Mdouble ft = (1 - a*std::pow(std::sin(species->getContactAngle()),b))*f0;
        
        Mdouble a_s = -0.3319*std::pow(V,0.4974) + 0.6717*std::pow(V,0.1995);
        Mdouble b_s = 13.84*std::pow(V,-0.3909) - 12.11*std::pow(V,-0.3945);

        Mdouble a_t_a = -0.007815*std::pow(log_V,2)-0.2105 *log_V -1.426;
        Mdouble b_t_a = -1.78*std::pow(V,0.8351) + 0.6669*std::pow(V,-0.01391);
        Mdouble a_t = a_t_a*std::pow(species->getContactAngle(),3)+b_t_a*species->getContactAngle()+1.0;

        Mdouble force_full_fit = ft*(1+a_s*S)/(1+a_t*b_s*a_s*S+a_t*b_s*std::pow(S,2));

        Mdouble fdotn = -force_full_fit * effectiveRadius * species->getSurfaceTension();

        addForce(getNormal() * fdotn);
    }
}

/*!
 * \return Mdouble
 */
Mdouble LiquidBridgeBagheriInteraction::getElasticEnergy() const
{
    ///\todo TW
    return 0.0;
}

/*!
 * \return const LiquidBridgeBagheriSpecies*
 */
const LiquidBridgeBagheriSpecies* LiquidBridgeBagheriInteraction::getSpecies() const
{
    return static_cast<const LiquidBridgeBagheriSpecies*>(getBaseSpecies()->getAdhesiveForce());
;
}

/*!
 * \return std::string
 */
std::string LiquidBridgeBagheriInteraction::getBaseName() const
{
    return "LiquidBridgeBagheri";
}

bool LiquidBridgeBagheriInteraction::getWasInContact() const
{
    return wasInContact_;
}

void LiquidBridgeBagheriInteraction::setWasInContact(bool wasInContact)
{
    wasInContact_ = wasInContact;
}
