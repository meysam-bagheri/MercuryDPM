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


#ifndef PARTICLE_HANDLER_H
#define PARTICLE_HANDLER_H

#include "BaseHandler.h"
#include "Particles/BaseParticle.h"
#include "Math/PSD.h"

class SpeciesHandler;

class BaseSpecies;

enum AttributeType
{
    MIN = 0, MAX = 1
};

/*!
 * \class ParticleHandler
 * \brief Container to store pointers to all particles.
 * \details The ParticleHandler is a container that stores pointers to all particles belonging to a simulation.
 * It is implemented as a vector of pointers of type BaseParticle*.
 */
class ParticleHandler : public BaseHandler<BaseParticle>
{
public:
    /*!
     * \brief Default constructor, it creates an empty ParticleHandler.
     */
    ParticleHandler();
    
    /*!
     * \brief Constructor that copies all BaseParticle it contains and then sets the smallest and largest particle.
     */
    ParticleHandler(const ParticleHandler& PH);
    
    /*!
     * \brief Assignment operator. 
     */
    ParticleHandler& operator=(const ParticleHandler& rhs);
    
    /*!
     * \brief Destructor, it destructs the ParticleHandler and all BaseParticle it contains.
     */
    ~ParticleHandler() override;
    
    /*!
     * \brief Adds a BaseParticle to the ParticleHandler.
     */
    void addExistingObject(BaseParticle* P) override;
    
    /*!
     * \brief Adds a BaseParticle to the ParticleHandler. 
     */
    void addObject(BaseParticle* P) override;
    
    /*!
     * \brief Adds a BaseParticle located at processor fromProcessor to toProcessor
     */
    void addObject(int fromProcessor, BaseParticle* P);
    
    /*!
     * \brief Adds a ghost particle located at fromProcessor to toProcessor
     */
    void addGhostObject(int fromProcessor, int toProcessor, BaseParticle* p);
    
    /*!
     * \brief Adds a BaseParticle to the ParticleHandler. 
     */
    void addGhostObject(BaseParticle* P) override;
    
    /*!
     * \brief Removes a BaseParticle from the ParticleHandler. 
     */
    void removeObject(unsigned int index) override;
    
    /*!
     * \brief Removes a BaseParticle from the ParticleHandler without a global check, this is only to be done for mpi routines
     * 
     */
    void removeGhostObject(unsigned int index);
    
    /*!
     * \brief Removes the last BaseParticle from the ParticleHandler.
     */
    void removeLastObject();
    
    /*!
     * \brief Computes the smallest particle (by interaction radius) and sets it in smallestParticle_
     */
    void computeSmallestParticle();
    
    /*!
     * \brief Computes the largest particle (by interaction radius) and sets it in largestParticle_
     */
    void computeLargestParticle();
    
    /*!
     * \brief Gets a pointer to the smallest BaseParticle (by interactionRadius) in this ParticleHandler of the local domain.
     */
    BaseParticle* getSmallestParticleLocal() const;
    
    /*!
     * \brief Gets a pointer to the smallest BaseParticle (by interactionRadius) in this ParticleHandler of the local domain.
     * When mercury is running in parallel this function throws an error since a pointer to another domain is useless
     */
    BaseParticle* getSmallestParticle() const;
    
    /*!
     * \brief Gets a pointer to the largest BaseParticle (by interactionRadius) in the ParticleHandler of the local domain
     */
    BaseParticle* getLargestParticleLocal() const;
    
    /*!
     * \brief Returns the pointer of the largest particle in the particle handler.
     * When mercury is running in parallel this function throws an error since a pointer to another domain is useless
     */
    BaseParticle* getLargestParticle() const;
    
    /*!
     * \brief Returns the smallest interaction radius of the current domain
     */
    Mdouble getSmallestInteractionRadiusLocal() const;
    
    /*!
     * \brief Returns the smallest interaction radius
     */
    Mdouble getSmallestInteractionRadius() const;
    
    /*!
     * \brief Returns the largest interaction radius of the current domain
     */
    Mdouble getLargestInteractionRadiusLocal() const;
    
    /*!
     * \brief Returns the largest interaction radius
     */
    Mdouble getLargestInteractionRadius() const;
    
    /*!
     * \brief Gets a pointer to the fastest BaseParticle in this ParticleHandler.
     */
    BaseParticle* getFastestParticleLocal() const;
    
    BaseParticle* getFastestParticle() const;
    
    Mdouble getKineticEnergy() const;
    
    Mdouble getRotationalEnergy() const;
    
    Mdouble getMass() const;
    
    Vec3D getMassTimesPosition() const;
    
    Vec3D getCentreOfMass() const;
    
    Vec3D getMomentum() const;
    
    Vec3D getAngularMomentum() const;
    
    Mdouble getVolume() const;
    
    Mdouble getMeanRadius() const;
    
    /*!
     * \brief Gets a pointer to the particle with the lowest coordinates in direction i in this ParticleHandler.
     */
    BaseParticle* getLowestPositionComponentParticleLocal(int i) const;
    
    /*!
     * \brief Gets a pointer to the particle with the lowest coordinates in direction i in this ParticleHandler.
     */
    BaseParticle* getLowestPositionComponentParticle(int i) const;
    
    /*!
     * \brief Gets a pointer to the particle with the highest coordinates in direction i in this ParticleHandler.
     */
    BaseParticle* getHighestPositionComponentParticleLocal(int i) const;
    
    /*!
     * \brief Gets a pointer to the particle with the highest coordinates in direction i in this ParticleHandler. 
     */
    BaseParticle* getHighestPositionComponentParticle(int i) const;
    
    /*!
     * \brief Gets a pointer to the particle with the lowest velocity in direction i in this ParticleHandler.
     */
    BaseParticle* getLowestVelocityComponentParticleLocal(int i) const;
    
    /*!
     * \brief Gets a pointer to the particle with the lowest velocity in direction i in this ParticleHandler. 
     */
    BaseParticle* getLowestVelocityComponentParticle(int i) const;
    
    /*!
     * \brief Gets a pointer to the particle with the highest velocity in direction i in this ParticleHandler.
     */
    BaseParticle* getHighestVelocityComponentParticleLocal(int i) const;
    
    /*!
     * \brief Gets a pointer to the particle with the highest velocity in direction i in this ParticleHandler. 
     */
    BaseParticle* getHighestVelocityComponentParticle(int i) const;
    
    /*!
     * \brief Computes an attribute type (min/max/..) of a particle attribute (position/velocity) in a local domain
     */
    //Mdouble getParticleAttributeLocal(std::function<Mdouble (BaseParticle*)> attribute, AttributeType type) const;
    
    /*!
     * \brief Computes an attribute type (min/max/..) of a particle attribute (position/velocity) in a local domain
     * \details Many functions the particleHandler return a pointer to a BaseParticle, i.e. the fastest particle, however
     * in parallel this is not useful as pointers can't be send across processes. This function gives the flexibility to
     * find a global particle extremum such as the fastest particle Velocity or lowest particle position.
     * \param[in] attribute A function that obtains a scalar from a particle. i.e. position or radius
     * \param[in] type An AttributeType tells this function what to do with the attribute, take the maximum or the minimum
     * \return Returns a scalar value representing the AttributeType of the Attribute: In easier terms returns the max/min/... of a particle attribute such as radius
      */
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, T>::type
    getParticleAttributeLocal(std::function<T(BaseParticle*)> attribute, AttributeType type) const
    {
        T attributeParticle;
        T attributeFinal;
        
        if ((*this).getSize() == 0)
        {
            logger(WARN, "ParticleHandler is empty: returning 0.0");
            attributeFinal = 0;
            return 0.0;
        }
        else
        {
            //Initialise with the first particle found
            attributeParticle = attribute(objects_[0]);
            attributeFinal = attributeParticle;
            
            //Find the final attribute
            for (BaseParticle* particle : (*this))
            {
                //Obtain the attribute
                if (!(particle->isMPIParticle() || particle->isPeriodicGhostParticle() || particle->isFixed()))
                {
                    attributeParticle = attribute(particle);
                }
                
                //Decide what to do with the magnitude of the attribute
                switch (type)
                {
                    case AttributeType::MIN :
                        if (attributeParticle < attributeFinal)
                        {
                            attributeFinal = attributeParticle;
                        }
                        break;
                    case AttributeType::MAX :
                        if (attributeParticle > attributeFinal)
                        {
                            attributeFinal = attributeParticle;
                        }
                        break;
                    default :
                        logger(ERROR, "Attribute type is not recognised");
                        break;
                }
            }
            return attributeFinal;
        }
    }
    
    /*!
     * \brief Computes an attribute type (min/max/..) of a particle attribute(position/velocity) in the global domain
     * \details Many functions the particleHandler return a pointer to a BaseParticle, i.e. the fastest particle, however
     * in parallel this is not useful as pointers can't be send across processes. This function gives the flexibility to
     * find a global particle extremum such as the fastest particle Velocity or lowest particle position.
     * \param[in] attribute A function that obtains a scalar from a particle. i.e. position or radius
     * \param[in] type An AttributeType tells this function what to do with the attribute, take the maximum or the minimum
     * \return Returns a scalar value representing the AttributeType of the Attribute: In easier terms returns the max/min/... of a particle attribute such as radius
     */
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, T>::type
    getParticleAttribute(std::function<T(BaseParticle*)> attribute, AttributeType type) const
    {
#ifdef MERCURYDPM_USE_MPI
        T particleAttributeLocal = getParticleAttributeLocal(attribute, type);
        T particleAttributeGlobal;

        //Communicate local to global using the approriate rule
        MPIContainer& communicator = MPIContainer::Instance();
        switch(type)
        {
            case AttributeType::MIN :
                communicator.allReduce(particleAttributeLocal,particleAttributeGlobal,MPI_MIN);
                break;
            case AttributeType::MAX :
                communicator.allReduce(particleAttributeLocal,particleAttributeGlobal,MPI_MAX);
                break;
            default :
                logger(ERROR,"Attribute type is not recognised");
                break;
        }
        return particleAttributeGlobal;
#else
        return getParticleAttributeLocal(attribute, type);
#endif
    }
    
    /*!
     * \brief Function returns the highest position in the x-direction
     */
    Mdouble getHighestPositionX() const;
    
    /*!
     * \brief Empties the whole ParticleHandler by removing all BaseParticle.
     */
    void clear() override;
    
    /*!

     * \brief Gets the number of particles that are fixed.
     */
    unsigned int getNumberOfFixedParticles() const;
    
    /*!
     * \brief Gets the number of particles that are not fixed.
     */
    unsigned int getNumberOfUnfixedParticles() const;
    
    /*!
     * \brief Reads BaseParticle into the ParticleHandler from restart data.
     */
    static BaseParticle* createObject(const std::string& type);
    
    /*!
     * \brief Create a new particle, based on the information provided in a restart file.
     */
    BaseParticle* readAndCreateObject(std::istream& is);

//    /*!
//     * \brief Create a new particle, based on the information from old-style restart data.
//     */
//    void readAndCreateOldObject(std::istream &is, const std::string& type);
    
    /*!
     * \brief Create a new particle in the WallHandler, based on the information provided in a restart file.
     */
    void readAndAddObject(std::istream& is) override;
    
    void write(std::ostream& os) const;
    
    /*!
     * \brief Checks if the extrema of this ParticleHandler needs updating. 
     */
    void checkExtrema(BaseParticle* P);
    
    /*!
     * \brief Checks if the extrema of this ParticleHandler needs updating when a particle is deleted.
     */
    void checkExtremaOnDelete(BaseParticle* P);
    
    /*!
     * \brief Computes the mass for all BaseParticle of the given species in this ParticleHandler.
     */
    void computeAllMasses(unsigned int indSpecies);
    
    /*!
     * \brief Computes the mass for all BaseParticle in this ParticleHandler.
     */
    void computeAllMasses();
    
    /*!
     * \brief Increment of the number of fixed particles.
     * \todo MX: For Jonny, is this still required, keeping the parallel code in mind?
     */
    void addedFixedParticle();
    
    /*!
     * \brief Decrement of the number of fixed particles.
     * \todo MX: For Jonny, is this still required, keeping the parallel code in mind?
     */
    void removedFixedParticle();
    
    /*!
     *  \brief Returns the name of the handler, namely the string "ParticleHandler".
     */
    std::string getName() const override;
    
    /*!
     * \brief Returns the number of real objects (on all processors)
     */
    unsigned int getNumberOfRealObjects() const;
    
    /*!
     * \brief Returns the number of objects in the container. In parallel code this practice is forbidden to avoid confusion
     * with real and fake particles. If the size of the container is wished, call size()
     */
    unsigned int getNumberOfObjects() const override;
    
    /*!
     * \brief Computes the number of Fixed particles on a local domain
     */
    unsigned int getNumberOfFixedObjectsLocal() const;
    
    /*!
     * \brief Computes the number of fixed particles in the whole simulation
     */
    unsigned int getNumberOfFixedObjects() const;
    
    /*!
     * \brief Returns the number of real objects on a local domain. MPI particles and periodic particles are neglected
     */
    unsigned int getNumberOfRealObjectsLocal() const;

    void actionsBeforeTimeStep();

    void actionsAfterTimeStep();

    double getLiquidFilmVolume() const;

    /*!
     * \brief Writes the PSD of the particles currently in the handler to a CSV file, with type PROBABILITYDENSITY_NUMBER_DISTRIBUTION.
     */
    void saveNumberPSDtoCSV(std::string csvFileName, std::vector<double> diameterBins = {}) const;

    /*!
     * \brief Gets the PSD from the particles currently in the handler.
     */
    PSD getPSD(std::vector<Mdouble> bins = {}, Mdouble scaleFactor = 1.0) const;

private:
    
    Mdouble getKineticEnergyLocal() const;
    
    Mdouble getRotationalEnergyLocal() const;
    
    Mdouble getMassLocal() const;
    
    Mdouble getVolumeLocal() const;
    
    Vec3D getMassTimesPositionLocal() const;
    
    Mdouble getSumRadiusLocal() const;
    
    /*!
     * \brief A pointer to the largest BaseParticle (by interactionRadius) in this ParticleHandler
     * \todo TW: note that checkExtrema gets called if a particle gets created 
     * and its Species and Radius gets set, even if it's not yet included in the 
     * particleHandler! This is necessary to check a not included particle for 
     * overlaps before inserting it into the handler. Not sure if this is a 
     * sensible structure; to be discussed.
     * Note: these statistic now include mpi and ghost particles as well
     */
    BaseParticle* largestParticle_;
    
    /*!
     * \brief A pointer to the smallest BaseParticle (by interactionRadius) in this ParticleHandler
     */
    BaseParticle* smallestParticle_;
    
    /*!
 * \brief Number of fixed particles
 */
    unsigned int NFixedParticles_;
};

#endif

