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

#ifndef CircularPeriodicBoundary_H
#define CircularPeriodicBoundary_H

#include "BaseBoundary.h"

class ParticleHandler;

class BaseParticle;

/*!
 * \class CircularPeriodicBoundary
 * \brief used to create a circular periodic boundary
 */
class CircularPeriodicBoundary : public BaseBoundary
{
public:
    /*!
     * \brief default constructor
     */
    CircularPeriodicBoundary();
    
    /*!
     * \brief Constructor
     */
    explicit CircularPeriodicBoundary(double innerRadius);
    
    /*!
     * \brief destructor
     */
    ~CircularPeriodicBoundary() override;
    
    /*!
     * \brief
     */
    CircularPeriodicBoundary* copy() const override;
    
    /*!
     * \brief
     * \param[in]
     */
    void rotateParticle(BaseParticle* P, double angle);
    
    /*!
     * \brief 
     * \param[in]
     */
    void createPeriodicParticle(BaseParticle* p, ParticleHandler& pH) override;
    
    void createPeriodicParticles(ParticleHandler& pH) override;
    
    /*!
     * \brief
     */
    bool checkBoundaryAfterParticleMoved(BaseParticle* P, ParticleHandler& pH);
    
    /// \todo MX: When implementing this function, I realised there might not be a unit test for this boundary
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    /*!
     * \brief reads the CircularPeriodicBoundary
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief
     */
    void oldRead(std::istream& is);
    
    /*!
     * \brief outputs the CircularPeriodicBoundary
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     * \return string type
     */
    std::string getName() const override;

private:
    
    /*!
     * \brief
     */
    double innerRadius_;
    
    //A particle is between to Radii in plane (i) and has two (straight) walls (plane 0 defines centre)
    //If it is close to its inner Radius it should be copied once with a positive rotation of 2*pi/2^i
    //If it is close to one of its straight walls is should be rotative with +/- 2*pi/2^i
    
    //Particles can only apply forces to real particles
    
    //If a particle crosses a straight wall it should simply be shifted
    //If a particle crosses its inner Radius it should be copied
    //If a particle crosses its outer Radius it may need to be deleted
    
};

#endif
