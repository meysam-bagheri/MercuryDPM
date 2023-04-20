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


#ifndef SUBCRITICALMASERBOUNDARY_H
#define SUBCRITICALMASERBOUNDARY_H

#include <map>

#include "Particles/BaseParticle.h"
#include "Boundaries/BaseBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Interactions/BaseInteraction.h"
#include "Math/Vector.h"

class ParticleSpecies;

/*!
 * \brief Variation on the PeriodicBoundary which also has an outflow part
 * \details Creates a boundary which divides the domain into two parts: a
 * periodic domain and an outflow domain.
 * Any particle flowing through the right of the periodic domain is
 * moved to both the left side of the periodic domain (as with a PeriodicBoundary), 
 * and also copied the left side of the outflow domain. Furthermore, the
 * particles near the right side of the periodic domain also exert forces on the
 * left side of the outflow domain, and vice versa. 
 * The left side of the periodic domain does not influence the right side of the
 * periodic domain.
 * When a particle in the outflow domain tries to enter the periodic domain, it
 * gets removed.
 * If a particle tries to leave the periodic domain to the left, it gets
 * removed.
 *
 * The difference between SubcriticalMaserBoundary (formerly known as
 * MaserBoundaryOldStyle) and ConstantMassFlowMaserBoundary (formerly known
 * simply as MaserBoundary) is that in ConstantMassFlowMaserBoundary, the left
 * hand side of the periodic box does not have any influence on the rest of the
 * flow, but the right side of the periodic box and the left side of the outflow
 * domain interact.  The ConstantMassFlowMaserBoundary is most useful for fast
 * (supercritical) flows, and for flows for which the flux across the boundary
 * needs to be controlled. The SubcriticalMaserBoundary is more useful for slow
 * flows, as the ConstantMassFlowMaserBoundary might generate "pulse-waves" in
 * those cases. 
 *
 * For a compact overview of the behaviour of SubcriticalMaserBoundary, please
 * look at the output of SubcriticalMaserBoundarySelfTest.
 *
 * \todo Which Maser is it used in Denissen2019?
 * To cite the Maser:
 *   I. F. C. Denissen, T. Weinhart, A. Te Voortwis, S. Luding, J. M. N. T. Gray
 *   and A. R. Thornton, Bulbous head formation in bidisperse shallow granular
 *   flow over an inclined plane. Journal of Fluid Mechanics, 866:263--297,
 *   mar 2019.
 */
class SubcriticalMaserBoundary : public BaseBoundary
{
public:
    /*!
     * \brief MaserBoundary constructor
     */
    SubcriticalMaserBoundary();
    
    /*!
     * \brief Maserboundary constructor that takes a periodic boundary, and converts it to a maser boundary
     */
    explicit SubcriticalMaserBoundary(const PeriodicBoundary& periodicBoundary);
    
    /*!
     * \brief Creates a copy of this maser on the heap.
     */
    SubcriticalMaserBoundary* copy() const override;
    
    /*!
     * \brief Sets all boundary properties at once and adds particles of the handler to the maser.
     * This deactivates the Maser.
     */
    void set(Vec3D normal, Vec3D planewiseShift, Mdouble distanceLeft, Mdouble distanceRight);

    /*!
     * \brief Sets the Maser's normal and positions, and sets the planewise
     * shift to zero. This deactivates the Maser. 
     */
    void set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight);
    
    /*!
     * \brief Sets a planewise direction to the shift. Doesn't change the normal
     * or the positions.
     */
    void setPlanewiseShift(Vec3D planewiseShift);

    /*!
     * \brief Sets the shift of the Maser. Usually don't use this directly, use
     * set() or setPlanewiseShift() instead.
     */
    void setShift(Vec3D shift);

    /*!
     * \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief writes boundary properties to ostream
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;
    
    /*!
     * \brief Creates periodic particles when the particle is a maser particle
     * and is sufficiently close to one of the boundary walls.
     */
    void createPeriodicParticle(BaseParticle* p, ParticleHandler& pH) override;
    
    void createPeriodicParticles(ParticleHandler& pH) override;
    
    /*!
     * \brief Shifts the particle to its 'periodic' position if it is a maser particle
     * and has crossed either of the walls. Creates a 'normal' particle at its current
     * position if it is a maser particle which crossed the RIGHT boundary wall.
     */
    bool checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH) const;
    
    /*!
     * \brief Evaluates what the particles have to do after they have changed position
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    /*!
     * \brief Converts a 'normal' particle into a maser particle.
     */
    void addParticleToMaser(BaseParticle* p);
    
    /*!
     * \brief Convert a maser particle into a 'normal' particle
     */
    void removeParticleFromMaser(BaseParticle* p);
    
    /*!
     * \brief Returns true if the particle is a Maser particle, and false otherwise.
     */
    bool isMaserParticle(BaseParticle* p) const;
    
    /*!
     * \brief Returns true if the particle is a Normal particle, and false otherwise.
     */
    bool isNormalParticle(BaseParticle* p) const;
    
    /*!
     * \brief Does everything that needs to be done for this boundary between setupInitialConditions and the time loop,
     * in this case, it activates the maser.
     */
    void actionsBeforeTimeLoop() override;
    
    /*!
     * \brief Opens the gap, and transforms particles to maser particles. Also calls turnOnCopying().
     */
    void activateMaser();
    
    /*!
     * \brief Stops copying particles (and act merely as a chute)
     */
    void deactivateMaser();
    
    Mdouble getDistanceLeft() const;
    
    Mdouble getDistanceRight() const;

protected:
    
    /*!
     * \brief Shifts the particle to its 'periodic' position
     */
    void shiftPosition(BaseParticle* p) const;
    
    /*!
     * \brief Creates a copy of the input particle, that gets removed again in DPMBase::removeDuplicatePeriodicParticles
     */
    BaseParticle* createGhostCopy(BaseParticle* p) const;
    
    /*!
     * \brief Returns whether the given particle is closer to the right boundary of the periodic part.
     * \param[in] p Particle for which we would like to know whether it is closest to the right boundary
     * \return      True if p is closer to the right boundary, false otherwise
     */
    bool isClosestToRightBoundary(const BaseParticle* const p) const
    {
        const Mdouble distance = Vec3D::dot(p->getPosition(), normal_);
        return (distanceRight_ - distance < distance - distanceLeft_);
    }
    
    /*!
     * \brief Returns the distance of the wall to the particle
     * \param[in] p     Pointer to the particle of which we want to know the distance to the wall to
     * \return          Distance of the particle to the boundary: positive means the particle is inside the periodic
     *                  part of the boundary, negative means it's outside.
     */
    Mdouble getDistance(BaseParticle* p) const
    {
        const Mdouble distance = Vec3D::dot(p->getPosition(), normal_);
        return std::min(distance - distanceLeft_, distanceRight_ - distance);
    }
    
    /*!
     * \brief Normal unit vector of both maser walls. Points in the flowing direction.
     */
    Vec3D normal_;
    
    /*!
     * \brief position of left boundary wall, s.t. normal*x=position_left
     */
    Mdouble distanceLeft_;
    
    /*!
     * \brief position of right boundary wall, s.t. normal*x=position_right
     */
    Mdouble distanceRight_;
    
    /*!
     * \brief Direction in which particles are to be shifted when they cross the boundary.
     * \details I.e., the vector pointing from a point the left boundary wall to the equivalent point
     * on the right one.
     * \details By default this is equal to normal_ * (distanceRight_ - distanceLeft_)
     * but you can also set a planewise direction to the shift.
     */
    Vec3D shift_;
    
    /*!
     * \brief List of 'normal' particles' species, and their maser counterparts
     */
    std::map<const ParticleSpecies*, const ParticleSpecies*> speciesConversionNormalToMaser_;
    
    /*!
     * \brief List of 'maser' particles' species, and their normal counterparts
     */
    std::map<const ParticleSpecies*, const ParticleSpecies*> speciesConversionMaserToNormal_;
    
    /*!
     * \brief Flag whether or not the gap is created and particles transformed already.
     */
    bool maserIsActivated_;
    
};

#endif
