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

#ifndef LiquidMigrationLSSpecies_H
#define LiquidMigrationLSSpecies_H

#include "Species/AdhesiveForceSpecies/BaseAdhesiveForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/AdhesiveForceInteractions/LiquidMigrationLSInteraction.h"

/*!
 * \brief LiquidMigrationLSSpecies contains the parameters used to describe a short-range force caused by liquid bridges.
 * \details See LiquidMigrationLSInteraction::computeForce for a description of the force law.
 */
class LiquidMigrationLSSpecies : public BaseAdhesiveForce
{
public:
    ///\brief The correct Interaction type for this AdhesiveForceSpecies
    typedef LiquidMigrationLSInteraction InteractionType;
    
    ///\brief The default constructor.
    LiquidMigrationLSSpecies();
    
    ///\brief The default copy constructor.
    LiquidMigrationLSSpecies(const LiquidMigrationLSSpecies& s);
    
    ///\brief The default destructor.
    ~LiquidMigrationLSSpecies();
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;
    
    ///\brief creates default values for mixed species
    void mix(LiquidMigrationLSSpecies* S, LiquidMigrationLSSpecies* T);

//adhesion-specific functions
    
    ///\brief returns the largest separation distance at which adhesive short-range forces can occur.
    void setInteractionDistance();

    /*!
     * \brief Returns a pointer to the adhesive force interaction LiquidMigrationLSInteraction.
     */
    //const LiquidMigrationLSInteraction*;

//setters and getters
    /*!
     * \brief used to set the Volume of the liquid bridge.
     */
    void setLiquidBridgeVolumeMax(Mdouble liquidBridgeVolume);
    
    /*!
     * \brief used to access the Volume of the liquid bridge.
     */
    Mdouble getLiquidBridgeVolumeMax() const;

    /*!
      * \brief used to set the Volume of the liquid bridge.
      */
    void setLiquidBridgeVolumeMin(Mdouble);

    /*!
     * \brief used to access the Volume of the liquid bridge.
     */
    Mdouble getLiquidBridgeVolumeMin() const;

    /*!
     * \brief used to set the distribution coefficient of the liquid.
     */
    void setDistributionCoefficient(Mdouble distributionCoefficient);
    
    /*!
     * \brief used to access the surface tension of the liquid.
     */
    Mdouble getDistributionCoefficient() const;
    
    /*!
     * \brief used to set the surface tension of the liquid.
     */
    void setSurfaceTension(Mdouble surfaceTension);
    
    /*!
     * \brief used to access the surface tension of the liquid.
     */
    Mdouble getSurfaceTension() const;
    
    /*!
     * \brief used to set the contact angle between particle and liquid bridge surface (in radians)
     */
    void setContactAngle(Mdouble contactAngle);
    
    /*!
     * \brief used to access the contact angle between particle and liquid bridge surface.
     */
    Mdouble getContactAngle() const;

    /*!
     * \brief used to set the Viscosity of the liquid.
     */
    void setViscosity(Mdouble viscosity);

    /*!
     * \brief used to access the Viscosity of the liquid.
     */
    Mdouble getViscosity() const;

private:
    /*! 
     * \brief Max volume a liquid bridge gets during formation; any excess remains in the liquid film of the particle.
     */
    Mdouble liquidBridgeVolumeMax_;

    /*!
     * \brief Minimum volume needed to form a liquid bridge.
     */
    Mdouble liquidBridgeVolumeMin_;

    /*! 
     * \brief distribution coefficient of the liquid
     * \details This is a variable that determines the fraction of the liquid that 
     * is to be distributed to the neighboring contacts of the particles after a 
     * liquid bridge rupture.
     */
    Mdouble distributionCoefficient_;
    
    /*! 
     * \brief contact angle between particle and liquid bridge surface (in radians)
     */
    Mdouble contactAngle_;
    
    /*! 
     * \brief surface tension of the liquid
     */
    Mdouble surfaceTension_;

    /*!
     * \brief viscosity of the liquid
     */
    Mdouble viscosity_;
    
};

#endif
