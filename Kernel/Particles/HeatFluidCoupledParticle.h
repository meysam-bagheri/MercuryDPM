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

#ifndef HeatFluidCoupledParticle_H
#define HeatFluidCoupledParticle_H

#include "ThermalParticle.h"
#include "LiquidFilmParticle.h"

/*!
 * \class HeatFluidCoupled
 * \brief Class of particles that store both temperature and liquid volume,
 *        which is adapted for the CFD-DEM studies.
 * \details In some CFD-DEM studies, a drying process needs to be simulated.
 *          Therefore, we selected the drying/liquid evaporation model of Azmir et al. (2018),
 *          which considers heat and mass transfer between the liquid film on the particles and the surrounding air.
 *          Their model considers convection, conduction, radiation of heat, and evaporation of liquid.
 *          Conduction is already implemented in ThermalParticles, from which this class derives;
 *          convection and radiation is still missing.
 *          This class implements evaporation model.
 */
template<class Particle>
class HeatFluidCoupled : public Thermal<LiquidFilm<Particle>>
{
public:
    /*!
     * \brief HeatFluidCoupled constructor creates a HeatFluidCoupled at (0,0,0) with radius,
     *        mass and inertia equal to 1.
     */
    HeatFluidCoupled() = default;

    /*!
     * \brief HeatFluidCoupled copy constructor, which accepts as input a reference to a HeatFluidCoupled.
     *        It creates a copy of this HeatFluidCoupled and all it's information.
     *        Usually it is better to use the copy() function for polymorphism.
     * \details Constructor that copies most of the properties of the given particle.
     *          Please note that not everything is copied, for example the position
     *          in the HGrid is not determined yet by the end of this constructor.
     *          It also does not copy the interactions and the pointer to the handler
     *          that handles this particle. Use with care.
     * \param[in,out] p  Reference to the HeatFluidCoupled this one should become a copy of.
     */
    HeatFluidCoupled(const HeatFluidCoupled& p) = default;

    /*!
     * \brief HeatFluidCoupled destructor, needs to be implemented and checked if it removes tangential spring information.
     * \details Destructor that asks the ParticleHandler to check if this was the smallest or largest particle and adjust itself accordingly.
     */
    ~HeatFluidCoupled() override = default;

    /*!
     * \brief HeatFluidCoupled copy method. Use copy constructor of this HeatFluidCoupled to create a copy on the heap,
     *        useful for polymorphism.
     *\return pointer to the particle's copy.
     */
    HeatFluidCoupled* copy() const override
    {
        return new HeatFluidCoupled<Particle>(*this);
    }

    /*!
     * \brief Returns the name of the object; in this case "HeatFluidCoupledParticle".
     * \return The object name.
     */
    std::string getName() const override
    {
        return "HeatFluidCoupled" + Particle::getName();
    }
    /// Tells the vtkWriter how many fields should be written for this particle type.
    unsigned getNumberOfFieldsVTK() const override { return 5; }

    /// Tells the vtkWriter the type of each field written for this particle type.
    std::string getTypeVTK(unsigned) const override { return "Float32"; }

    /// Tells the vtkWriter the name of each field written for this particle type.
    std::string getNameVTK(unsigned i) const override {
        if (i==0)
            return "fullLiquidVolume";
        else if (i==1)
            return "liquidFilmVolume";
        else if (i==2)
            return "liquidBridgeVolume";
        else if (i==3)
            return "totalEvaporatedLiquidVolume";
        else /* i=4 */
            return "temperature";
    }

    /// Tells the vtkWriter the value of each field written for this particle type.
    std::vector<Mdouble> getFieldVTK(unsigned i) const override {
        if (i==0)
            return { this->getFullLiquidVolume() };
        else if (i==1)
            return { this->liquidVolume_ };
        else if (i==2)
            return { this->getLiquidBridgeVolume() };
        else if (i==3)
            return { this->totalEvaporatedLiquidVolume_ };
        else /* i=4 */
            return { this->temperature_ };
    }

    /// The actionAfterTimeStep is defined in the species, as we cannot extract the species properties of a HeatFluidCoupled*Species
    void actionsAfterTimeStep() override {
        this->getSpecies()->actionsAfterTimeStep(this);
    }
};

/// Template specialisation of HeatFluidCoupled<Particle> for spherical particles.
typedef HeatFluidCoupled<SphericalParticle> HeatFluidCoupledParticle;
#endif
