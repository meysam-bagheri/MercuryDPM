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

#ifndef LiquidBridgeBagheriSpecies_H
#define LiquidBridgeBagheriSpecies_H

#include "Species/AdhesiveForceSpecies/BaseAdhesiveForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/AdhesiveForceInteractions/LiquidBridgeBagheriInteraction.h"

/*!
 * \brief LiquidBridgeBagheriSpecies contains the parameters used to describe a short-range force caused by liquid bridges.
 * \details See LiquidBridgeBagheriInteraction::computeForce for a description of the force law.
 */
class LiquidBridgeBagheriSpecies : public BaseAdhesiveForce
{
public:
    ///\brief The correct Interaction type for this AdhesiveForceSpecies
    typedef LiquidBridgeBagheriInteraction InteractionType;
    
    ///\brief The default constructor.
    LiquidBridgeBagheriSpecies();
    
    ///\brief The default copy constructor.
    LiquidBridgeBagheriSpecies(const LiquidBridgeBagheriSpecies& s);
    
    ///\brief The default destructor.
    ~LiquidBridgeBagheriSpecies();
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;
    
    ///\brief creates default values for mixed species
    void mix(LiquidBridgeBagheriSpecies* S, LiquidBridgeBagheriSpecies* T);
    
    void setInteractionDistance();
    
//setters and getters
    /*!
     * \brief used to set the Volume of the liquid bridge.
     */
    void setLiquidBridgeVolume(Mdouble liquidBridgeVolume);
    
    /*!
     * \brief used to access the Volume of the liquid bridge.
     */
    Mdouble getLiquidBridgeVolume() const;
    
    /*!
     * \brief used to set the surface tension of the liquid.
     */
    void setSurfaceTension(Mdouble surfaceTension);
    
    /*!
     * \brief used to access the surface tension of the liquid.
     */
    Mdouble getSurfaceTension() const;
    
    /*!
     * \brief used to set the contact angle between particle and liquid bridge surface.
     */
    void setContactAngle(Mdouble contactAngle);
    
    /*!
     * \brief used to access the contact angle between particle and liquid bridge surface.
     */
    Mdouble getContactAngle() const;

    /*!
     * \brief used to access the mean radius of the particle.
     */
    void setMeanRadius(Mdouble meanRadius);
    
    /*!
     * \brief used to access the mean radius of the particle.
     */
    Mdouble getMeanRadius() const;
private:
    /*! 
     * \brief Volume of the liquid bridge. 
     * \details Currently, this is a Species property (i.e. all particles have 
     * the same liquid bridge volume), but this should probably become a 
     * Particle or Interaction property.
     */
    Mdouble liquidBridgeVolume_;
    
    /*! 
     * \brief Cube root of the liquid bridge. 
     * \details This is a slave variable of liquidBridgeVolume_ (i.e. it is 
     * always set together) but is useful, because the cube root has to be 
     * calculated to obtain the liquid bridge adhesive force.
     */
    Mdouble cbrtLiquidBridgeVolume_;
    
    /*! 
     * \brief contact angle between particle and liquid bridge surface
     */
    Mdouble contactAngle_;
    
    /*! 
     * \brief surface tension of the liquid
     */
    Mdouble surfaceTension_;

     /*! 
     * \brief mean radius of the particles
     */
    Mdouble meanRadius_;   
    
};

#endif
