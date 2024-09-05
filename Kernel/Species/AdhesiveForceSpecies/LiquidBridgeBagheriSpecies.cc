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

#include "LiquidBridgeBagheriSpecies.h"
#include "Logger.h"
#include "Species/BaseSpecies.h"

LiquidBridgeBagheriSpecies::LiquidBridgeBagheriSpecies()
{
    liquidBridgeVolume_ = 0;
    cbrtLiquidBridgeVolume_ = 0;
    surfaceTension_ = 0;
    contactAngle_ = 0;
    // setInteractionDistance();
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeBagheriSpecies::LiquidBridgeBagheriSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
LiquidBridgeBagheriSpecies::LiquidBridgeBagheriSpecies(const LiquidBridgeBagheriSpecies& s)
{
    liquidBridgeVolume_ = s.liquidBridgeVolume_;
    cbrtLiquidBridgeVolume_ = s.cbrtLiquidBridgeVolume_;
    surfaceTension_ = s.surfaceTension_;
    contactAngle_ = s.contactAngle_;

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeBagheriSpecies::LiquidBridgeBagheriSpecies(const LiquidBridgeBagheriSpecies &p) finished"<<std::endl;
#endif
}

LiquidBridgeBagheriSpecies::~LiquidBridgeBagheriSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LiquidBridgeBagheriSpecies::~LiquidBridgeBagheriSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void LiquidBridgeBagheriSpecies::write(std::ostream& os) const
{
    os << " liquidBridgeVolume " << liquidBridgeVolume_;
    os << " surfaceTension " << surfaceTension_;
    os << " contactAngle " << contactAngle_;
}

/*!
 * \param[in] input stream (typically the restart file)
 */
void LiquidBridgeBagheriSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> liquidBridgeVolume_;
    is >> dummy >> surfaceTension_;
    is >> dummy >> contactAngle_;
    setInteractionDistance();
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string LiquidBridgeBagheriSpecies::getBaseName() const
{
    return "LiquidBridgeBagheri";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void LiquidBridgeBagheriSpecies::mix(LiquidBridgeBagheriSpecies* const S, LiquidBridgeBagheriSpecies* const T)
{
    setLiquidBridgeVolume(BaseSpecies::average(S->getLiquidBridgeVolume(), T->getLiquidBridgeVolume()));
    surfaceTension_ = BaseSpecies::average(S->getSurfaceTension(), T->getSurfaceTension());
    contactAngle_ = BaseSpecies::average(S->getContactAngle(), T->getContactAngle());
}

///\return the maximum separation distance below which adhesive forces can occur (needed for contact detection)
void LiquidBridgeBagheriSpecies::setInteractionDistance()
{
    getBaseSpecies()->setInteractionDistance((1.0 + 0.5 * contactAngle_) * cbrtLiquidBridgeVolume_);
}

/*!
 * \param[in] liquidBridgeVolume the volume of the liquid bridge.
 */
void LiquidBridgeBagheriSpecies::setLiquidBridgeVolume(Mdouble liquidBridgeVolume)
{
    if (liquidBridgeVolume >= 0)
    {
        liquidBridgeVolume_ = liquidBridgeVolume;
        cbrtLiquidBridgeVolume_ = cbrt(liquidBridgeVolume);
        setInteractionDistance();
    }
    else
    {
        logger(ERROR, "Error in setLiquidBridgeVolume(%)", liquidBridgeVolume);
    }
}

/*!
 * \return the volume of the liquid bridge.
 */
Mdouble LiquidBridgeBagheriSpecies::getLiquidBridgeVolume() const
{
    return liquidBridgeVolume_;
}

/*!
 * \param[in] surfaceTension the surface tension of the liquid.
 */
void LiquidBridgeBagheriSpecies::setSurfaceTension(Mdouble surfaceTension)
{
    if (surfaceTension >= 0)
        surfaceTension_ = surfaceTension;
    else
    {
        logger(ERROR, "Error in setSurfaceTension");
    }
}

/*!
 * \return the surface tension of the liquid.
 */
Mdouble LiquidBridgeBagheriSpecies::getSurfaceTension() const
{
    return surfaceTension_;
}

/*!
 * \param[in] contactAngle the contact angle between particle and liquid bridge surface.
 */
void LiquidBridgeBagheriSpecies::setContactAngle(Mdouble contactAngle)
{
    logger.assert_debug(contactAngle >= 0,"Error in setContactAngle");
    contactAngle_ = contactAngle;
    setInteractionDistance();
}

/*!
 * \return the contact angle between particle and liquid bridge surface.
 */
Mdouble LiquidBridgeBagheriSpecies::getContactAngle() const
{
    return contactAngle_;
}