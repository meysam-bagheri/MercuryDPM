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

#include "InsertionBoundary.h"
#include "DPMBase.h"
#include "Particles/BaseParticle.h"
#include<iostream>

#ifdef MERCURYDPM_USE_MPI
#include "MpiDataClass.h"
#endif

/*!
 * \details Default constructor, sets all data members to 0 or default.
 */
InsertionBoundary::InsertionBoundary() : BaseBoundary()
{
    numberOfParticlesInserted_ = 0;
    massInserted_ = 0;
    volumeInserted_ = 0;
    maxFailed_ = 0;
    isActivated_ = true;
    volumeFlowRate_ = constants::inf;
    initialVolume_ = 0;
    samplingInterval_ = 0;
    checkParticleForInteraction_ = true;
    velMin_ = Vec3D(0.0, 0.0, 0.0);
    velMax_ = Vec3D(0.0, 0.0, 0.0);
    particleSizeDistributionVector_.resize(1);
    isManuallyInserting_ = false;
    chosenSpecies_ = 0;
}

/*!
 * \details Copy constructor
 */
InsertionBoundary::InsertionBoundary(const InsertionBoundary& other)
        : BaseBoundary(other)
{
    numberOfParticlesInserted_ = other.numberOfParticlesInserted_;
    massInserted_ = other.massInserted_;
    volumeInserted_ = other.volumeInserted_;
    maxFailed_ = other.maxFailed_;
    isActivated_ = other.isActivated_;
    volumeFlowRate_ = other.volumeFlowRate_;
    initialVolume_ = other.initialVolume_;
    samplingInterval_ = other.samplingInterval_;
    variableCumulativeVolumeFlowRate_ = other.variableCumulativeVolumeFlowRate_;
    checkParticleForInteraction_ = other.checkParticleForInteraction_;
    particleSizeDistributionVector_ = other.particleSizeDistributionVector_;
    probability_ = other.probability_;
    velMin_ = other.velMin_;
    velMax_ = other.velMax_;
    isManuallyInserting_ = other.isManuallyInserting_;
    chosenSpecies_ = other.chosenSpecies_;

    for (int i = 0; i < other.particleToCopy_.size(); i++)
    {
        particleToCopy_.resize(other.particleToCopy_.size());
        particleToCopy_[i] = other.particleToCopy_[i]->copy();
    }
}

/*!
 * \details Destructor that deletes the BaseParticle that is copied and inserted
 *          at every insertion.
 */
InsertionBoundary::~InsertionBoundary()
{
    for (auto& particleToCopy: particleToCopy_)
    {
        delete particleToCopy;
    }

}

/*!
 * \details The default behaviour will be to return particleToCopy_, but this
 * can be overridden by the children (to get size or species dispersity).
 */
BaseParticle* InsertionBoundary::generateParticle(RNG& random)
{
    double check = random.getRandomNumber(0, 1);
    for (int i = 0; i < probability_.size(); i++)
    {
        if (check < probability_[i])
        {
            chosenSpecies_ = i;
            break;
        }
        else
        {
            check -= probability_[i];
        }
    }
    BaseParticle* P = particleToCopy_[chosenSpecies_]->copy();
    if (particleSizeDistributionVector_[chosenSpecies_].getParticleSizeDistribution().empty())
    {
        particleSizeDistributionVector_[chosenSpecies_].setDistributionUniform(P->getRadius(), P->getRadius(), 50);
        logger(INFO, "Assembling a discrete uniform particle size distribution (50 bins) with the radius set for the "
                     "InsertionBoundary.");
    }
    // manual insertion routine to insert PSDs as accurate as possible into a given volume
    if (isManuallyInserting_)
    {
        logger.assert_debug(initialVolume_ > 0.0, "Use setInitialVolume to define the particle insertion volume");
        Mdouble radius;
        // getVolumeFlowRate() * time + initialVolume_ - volumeInserted_ lead to more inaccurate results, therfore
        // -volumeInserted was removed.
        radius = particleSizeDistributionVector_[chosenSpecies_].insertManuallyByVolume(getVolumeFlowRate() *
                                                                                       getHandler()->getDPMBase()->getTime() +
                                                                                       initialVolume_);
        P->setRadius(radius);
        return P;
    }
    Mdouble radius;
    radius = particleSizeDistributionVector_[chosenSpecies_].drawSample();
    P->setRadius(radius);
    return P;
}

/*!
 * \details Checks the inserted total volume by the flowrate and initialVolume and returns if a particle
 * is allowed to be inserted.
 * \return TRUE if the inserted volume is lower than the total volume and FALSE if the inserted volume is higher than
 * the total volume.
 */
bool InsertionBoundary::insertParticle(Mdouble time)
{
    // check if the flow rate limit has been reached
    if (variableCumulativeVolumeFlowRate_.empty())
    {
        return volumeInserted_ < getVolumeFlowRate() * time + initialVolume_;
    }
    else
    {
        const Mdouble iMax = (Mdouble) variableCumulativeVolumeFlowRate_.size() - 2;
        const Mdouble i = std::min(time / samplingInterval_, iMax);
        if (i == iMax)
        {
            static unsigned count = 0;
            if (count == 0)
            {
                logger(WARN, "Reached end of volume flowrate function");
            }
            ++count;
        }
        const size_t id = i;
        const Mdouble allowedVolume = variableCumulativeVolumeFlowRate_[id] +
                (variableCumulativeVolumeFlowRate_[id + 1] -
                        variableCumulativeVolumeFlowRate_[id]) * (i - id);
        return volumeInserted_ < allowedVolume;
    }
}

/*!
 * \details Is used to fill the insides of the boundary with particles until 
 * it is filled up. 
 * \param[in,out] md    the problem's DPMBase object
 * \todo rename to something like "insertUntilMaxFailed"?
 */
void InsertionBoundary::checkBoundaryBeforeTimeStep(DPMBase* md)
{
    logger(VERBOSE, "In InsertionBoundary::checkBoundaryBeforeTimeStep\n");

    if (!isActivated_)
    {
        return;
    }

    /* Each timestep, the InsertionBoundary attempts to fill up a region with
     * particles.
     *
     * It first calls generateParticle() to get a randomised particle, subject
     * to a specified distribution over sizes and species. (The basic class
     * supports size dispersity only but PolydisperseInsertionBoundary will
     * support species dispersity.) 
     * Then it repeatedly calls placeParticle(), which gives the particle a
     * random location (and possibly velocity) in a specified,
     * geometry-dependent bound. Each time, it checks whether the new particle
     * would have an interaction with another particle or a wall. 
     *
     * If it manages to do that within maxFailed_ tries, then:
     *   * the new particle is inserted, 
     *   * the failure counter is reset, and
     * the processes is repeated with a new generateParticle().
     *
     * Otherwise, the processes terminates for this timestep.
     * */

    // Keep count of how many successive times we have failed to place a new
    // particle. 
    unsigned int failed = 0;
    while (failed <= maxFailed_ && insertParticle(md->getNextTime())) // 'generating' loop
    {
        /* Generate random *intrinsic* properties for the new particle. */
        logger(VERBOSE, "about to call generateParticle\n");

        //generate a particle the first time
        if (!p0) p0 = generateParticle(md->random);
        // Important for particle generation with a particle size distribution as it generates a particle with zero
        // radius. If a particle is not allowed to be inserted by the PSD criteria it will generate a particle with
        // zero diameter. This if statement prevents inserting particles with zero radius, which would else be a problem.
        /// \todo create a definition for zero-size particles. Right now the zero-size particle is only used as a stop criterion for the manual PSD insertion
        if (p0->getRadius() == 0)
        {
            logger(VERBOSE, "The PSD for the specified volume is fully set");
            // free up memory space
            delete p0;
            // out of the 'placing' loop
            failed = maxFailed_ + 1;
            continue;
        }
        logger(VERBOSE, "generated a particle with intrinsics %", p0);

        while (true) // 'placing' loop
        {
            /* Generate extrinsic properties (position and velocity) for this
             * new particle. */
            placeParticle(p0, md->random);
            logger(VERBOSE, "attempting to place particle with radius % at %, vel %", p0->getMaxInteractionRadius(), p0->getPosition(), p0->getVelocity());

#ifdef MERCURYDPM_USE_MPI
            /* Communicate the new particle's properties by setHandler (note
                 * that this doesn't actually add the particle to the handler). */
                if (NUMBER_OF_PROCESSORS > 1)
                {
                    MPIParticle particle;
                    // 	//Every domain generates a particle (to get the species right etc)

                    //Send particle data from root to other processors to sync the particle properties
                    if (PROCESSOR_ID == 0)
                    {
                        particle.copyDataFromParticleToMPIParticle(p0);
                    }

                    MPIContainer::Instance().broadcast(&particle,MercuryMPIType::PARTICLE);

                    //Process the received data
                    if (PROCESSOR_ID != 0)
                    {
                        copyDataFromMPIParticleToParticle(&particle, p0, &(md->particleHandler));
                    }
                }
#endif
            p0->setHandler(&md->particleHandler);
            /* Check whether the particle has any interactions. */
            if (!checkParticleForInteraction_ || md->checkParticleForInteraction(*p0))
            {
                //Note: in parallel only one of the domains will actually add the particle
                md->particleHandler.copyAndAddObject(p0);
                failed = 0;
                
                ++numberOfParticlesInserted_;
                const double volume = p0->getVolume();
                volumeInserted_ += volume;
                massInserted_ += p0->getSpecies()->getDensity() * volume;
                logger(VERBOSE, "successfully placed a particle %, with position: % after % fails.", p0,
                       p0->getPosition(), failed);
                /* JMFT: The generateParticle() routine allocates memory, so we should
                 * free it here. (Don't worry, the particle will have been copied to the
                 * particleHandler by this point iff we want it.) */
                delete p0;
                // generate a new particle to be inserted only if the last particle could be successfully places
                p0 = generateParticle(md->random);
                break; // out of the 'placing' loop
            }
            else
            {
                failed++;
                logger(VERBOSE, "failed to place a particle; have failed % times", failed);
            }

            if (failed > maxFailed_)
            {
                logger(VERBOSE, "failed too many times; giving up");
                if (isManuallyInserting_)
                {
                    particleSizeDistributionVector_[chosenSpecies_].decrementNParticlesPerClass();
                    particleSizeDistributionVector_[chosenSpecies_].decrementVolumePerClass(p0->getVolume());
                }
                break; // out of the 'placing' loop (and will leave the 'generating' loop too
            }
        }
        logger(VERBOSE, "failed % times, so breaking out of InsertionBoundary loop for this timestep.", failed);
    }
    // logger(INFO, "volumeInserted_ = %", volumeInserted_);
}

/*!
 * \detail calls the function checkBoundaryBeforeTimeStep() to fill a domain with particles; also reports how many
 * particles where inserted at the end of the routine.
 */
void InsertionBoundary::insertParticles(DPMBase* md)
{
    checkBoundaryBeforeTimeStep(md);
    logger(INFO, "Inserted % particles", getNumberOfParticlesInserted());
}

/*!
 * \details Returns the number of particles inserted in the boundary
 * \return the number of particles inserted
 */
unsigned int InsertionBoundary::getNumberOfParticlesInserted() const
{
    return numberOfParticlesInserted_;
}

/*!
 * \details Returns the mass of particles inserted in the boundary
 * \return the mass of particles inserted
 */
Mdouble InsertionBoundary::getMassOfParticlesInserted() const
{
    return massInserted_;
}

/*!
 * \details Returns the volume of particles inserted in the boundary
 * \return the volume of particles inserted
 */
Mdouble InsertionBoundary::getVolumeOfParticlesInserted() const
{
    return volumeInserted_;
}

/*! 
 * \details set all particle property counter variables to zero. reset() does not activate or deactivate the
 * InsertionBoundary.
 */
void InsertionBoundary::reset()
{
    numberOfParticlesInserted_ = 0;
    massInserted_ = 0;
    volumeInserted_ = 0;
}

/*!
 * \details Turns on the InsertionBoundary by setting the boolean to TRUE.
 */
void InsertionBoundary::activate()
{
    isActivated_ = true;
}

/*!
 * \details Turns off the InsertionBoundary by setting the boolean to FALSE.
 */
void InsertionBoundary::deactivate()
{
    isActivated_ = false;
}

/*!
 * \details checks the activation status of the InsertionBoundary by checking the respective boolean variable.
 * \return TRUE for an activated InsertionBoundary and FALSE for a deactivated InsertionBoundary
 */
bool InsertionBoundary::isActivated()
{
    return isActivated_;
}

/*!
 * \details Return maxFailed_ (see InsertionBoundary::set).
 * \return the maximum number of particle insertion trials
 */
unsigned int InsertionBoundary::getMaxFailed() const
{
    return maxFailed_;
}

/*!
 * \details Sets the vector of pointers to particles which will be inserted by the insertion boundary. This is mainly
 * used to insert particles with different intrinsic properties (such as PSD, mechanical properties, etc.)
 * \param[in] particleToCopy    vector of pointers to the particles which are to be inserted
 */
void InsertionBoundary::setParticleToCopy(std::vector<BaseParticle*> particleToCopy)
{
    if (particleToCopy.empty())
    {
        logger(ERROR, "Setting particleToCopy to be empty?");
    }
    if (!particleToCopy_.empty())
    {
        for (auto& ParticleToCopy: particleToCopy_)
        {
            delete ParticleToCopy;
        }
    }
    for (int i = 0; i < particleToCopy.size(); i++)
    {
        particleToCopy_.resize(particleToCopy.size());
        particleToCopy_[i] = particleToCopy[i]->copy();
    }
}

void InsertionBoundary::shiftBoundary(Vec3D shift)
{

}

void InsertionBoundary::rotateBoundary(Vec3D angle)
{

}

/*!
 * \details Sets the vector of pointers to particles which will be inserted by the insertion boundary.
 * \param[in] particleToCopy    pointer to the particle to be inserted
 */
void InsertionBoundary::setParticleToCopy(BaseParticle* particleToCopy)
{
    if (particleToCopy == nullptr)
    {
        logger(ERROR, "Setting particleToCopy to be a null pointer?");
    }
    else
    {
        if (!particleToCopy_.empty())
        {
            for (auto& particleToCopy: particleToCopy_)
            {
                delete particleToCopy;
            }
        }
        particleToCopy_.resize(1);
        particleToCopy_[0] = particleToCopy->copy();
    }
}

/*!
 * \details returns pointer to the particle copies which are to be inserted
 */
std::vector<BaseParticle*> InsertionBoundary::getParticleToCopy()
{
    if (particleToCopy_.empty())
    {
        logger(ERROR, "particleToCopy not set");
    }
    return particleToCopy_;
}

/*!
 * \details reads the boundary's id_ and maxFailed_ from the given istream
 * \param[in,out] is    stream the data members are read from
 */
void InsertionBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy, type;
    is >> dummy >> maxFailed_ >> dummy;
    if (dummy == "volumeFlowRate")
    {
        is >> volumeFlowRate_ >> dummy;
    }
    is >> massInserted_;
    is >> dummy >> volumeInserted_;
    is >> dummy >> numberOfParticlesInserted_;
    is >> dummy >> isActivated_;
    size_t psdVectorSize;
    is >> dummy >> psdVectorSize;
    particleSizeDistributionVector_.clear();
    particleSizeDistributionVector_.resize(psdVectorSize);
    for (auto& particleSizeDistributionVector: particleSizeDistributionVector_)
    {
        size_t psdSize;
        DistributionElements radiusAndProbability{};
        std::vector<DistributionElements> particleSizeDistribution{};
        is >> dummy >> psdSize;
        particleSizeDistribution.clear();
        particleSizeDistribution.reserve(psdSize);
        for (size_t i = 0; i < psdSize; i++)
        {
            is >> radiusAndProbability.internalVariable;
            is >> radiusAndProbability.probability;
            particleSizeDistribution.push_back(radiusAndProbability);
        }
        particleSizeDistributionVector.setPSDFromVector(particleSizeDistribution, PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION);
    }
    if (psdVectorSize > 1)
    {
        is >> dummy;
        Mdouble psdRatio;
        probability_.clear();
        probability_.reserve(psdVectorSize);
        for (size_t i = 0; i < psdVectorSize; i++)
        {
            is >> psdRatio;
            probability_.push_back(psdRatio);
        }
    }

    ///\todo make theses reads non-optional
    helpers::readOptionalVariable(is, "checkParticleForInteraction", checkParticleForInteraction_);
    helpers::readOptionalVariable(is, "initialVolume", initialVolume_);
    if (helpers::readOptionalVariable(is, "samplingInterval", samplingInterval_))
    {
        size_t n;
        Mdouble flowRate;
        is >> dummy >> n;
        //variableCumulativeVolumeFlowRate_.clear();
        variableCumulativeVolumeFlowRate_.reserve(n);
        for (size_t i = 0; i < n; ++i)
        {
            is >> flowRate;
            variableCumulativeVolumeFlowRate_.push_back(flowRate);
        }
    }
    is >> dummy;
    if (dummy != "noParticleToCopy")
    {
        size_t particleToCopySize;
        for (auto& particleToCopy: particleToCopy_)
        {
            delete particleToCopy;
        }
        is >> particleToCopySize;
        particleToCopy_.resize(particleToCopySize);
        for (auto& particleToCopy: particleToCopy_)
        {
            particleToCopy = getHandler()->getDPMBase()->particleHandler.readAndCreateObject(is);
            // The .restart file records the index of the particle's species, but
            // doesn't record the pointer, i.e. the memory address of the species within
            // the speciesHandler. The latter needs to be reset now.
            particleToCopy->setSpecies(getHandler()->getDPMBase()->speciesHandler.getObject(
                    particleToCopy->getIndSpecies()));
        }
    }
}

/*!
 * \details adds the boundary's id_ and maxFailed_ to the given ostream
 * \param[in,out] is    stream the data members are to be added to
 */
void InsertionBoundary::write(std::ostream& os) const
{
    logger(VERBOSE, "In InsertionBoundary::write\n");
    BaseBoundary::write(os);
    os << " maxFailed " << maxFailed_;
    if (std::isfinite(volumeFlowRate_))
        os << " volumeFlowRate " << volumeFlowRate_;
    os << " massInserted " << massInserted_;
    os << " volumeInserted " << volumeInserted_;
    os << " numberOfParticlesInserted " << numberOfParticlesInserted_;
    os << " isActivated " << isActivated_;
    os << " psdCount " << particleSizeDistributionVector_.size();
    for (auto& particleSizeDistribution: particleSizeDistributionVector_)
    {
        os << " psd " << particleSizeDistribution.getParticleSizeDistribution().size();
        for (auto p: particleSizeDistribution.getParticleSizeDistribution())
        {
            os << " " << p.internalVariable
               << " " << p.probability;
        }
    }
    if (!probability_.empty())
    {
        os << " psdRatio";
        for (auto& psdRatio: probability_)
        {
            os << " " << psdRatio;
        }
    }
    if (!checkParticleForInteraction_)
    {
        os << " checkParticleForInteraction " << checkParticleForInteraction_;
    }
    os << " initialVolume " << initialVolume_;
    os << " samplingInterval " << samplingInterval_;
    os << " variableCumulativeVolumeFlowRate " << variableCumulativeVolumeFlowRate_.size();
    for (const auto flowRate: variableCumulativeVolumeFlowRate_)
    {
        os << ' ' << flowRate;
    }
    if (!particleToCopy_.empty())
    {
        os << " particleToCopy " << particleToCopy_.size();
        for (auto& particleToCopy: particleToCopy_)
        {
            os << " " << *particleToCopy;
        }
    }
    else
    {
        os << " noParticleToCopy";
    }
}

/*!
 * \details Gets the volumetric flow rate of the insertion routine.
 * \return A double which corresponds to the flow rate of the insertion routine
 */
Mdouble InsertionBoundary::getVolumeFlowRate() const
{
    return volumeFlowRate_;
}

/*!
 * \details Sets the volumetric flow rate of the insertion routine.
 */
void InsertionBoundary::setVolumeFlowRate(Mdouble volumeFlowRate)
{
    volumeFlowRate_ = volumeFlowRate;
}

/*!
 * \details Gets the volume to be inserted by the insertion routine.
 * \return A double which corresponds to volume to be inserted by the insertion routines
 */
Mdouble InsertionBoundary::getInitialVolume() const
{
    return initialVolume_;
}

/*!
 * \details Sets the volume to be inserted by the insertion routine. (Total volume to insert = getVolumeFlowRate() *
 * time + initialVolume_)
 */
void InsertionBoundary::setInitialVolume(Mdouble initialVolume)
{
    initialVolume_ = initialVolume;
    if (!std::isfinite(volumeFlowRate_))
    {
        volumeFlowRate_ = 0;
    }
}

/*!
 * \details Sets a variable volume flow rate taken at fixed sampling intervals; the values are cumulative; thus,
 * we need to ensure the volume inserted before time t=n*samplingInterval is less than
 * variableCumulativeVolumeFlowRate[n].
 * \see variableCumulativeVolumeFlowRate_
 */
void InsertionBoundary::setVariableVolumeFlowRate(const std::vector<Mdouble>& variableCumulativeVolumeFlowRate,
                                                  Mdouble samplingInterval)
{
    logger.assert_debug(samplingInterval > 0, "sampling interval needs to be positive");
    const Mdouble endTime = variableCumulativeVolumeFlowRate.size() * samplingInterval;
    logger(INFO, "variable flowrate is defined up to %", endTime);
    logger.assert_always(getHandler()->getDPMBase()->getTimeMax() < endTime,
                         "variable flowrate is defined up to %, but tMax is set to %", endTime,
                         getHandler()->getDPMBase()->getTimeMax());
    variableCumulativeVolumeFlowRate_ = variableCumulativeVolumeFlowRate;
    samplingInterval_ = samplingInterval;
}

/*!
 * \details Sets the range of particle radii that may be generated to a PSD defined by the user.
 */
/// \todo TP: Consider std::move instead of a set function. This would result in a speedup + no one has to set the PSD for insertionBoundaries again as it is set when a PSD is inserted.
void InsertionBoundary::setPSD(const PSD psd)
{
    particleSizeDistributionVector_.resize(1);
    particleSizeDistributionVector_[0] = psd;
}

/*!
 * \details Sets the ranges of particle radii that may be generated by different PSDs defined by the user.
 */
/// \todo TP: Consider std::move instead of a set function. This would result in a speedup + no one has to set the PSD for insertionBoundaries again as it is set when a PSD is inserted.
void InsertionBoundary::setPSD(std::vector<PSD> psd, std::vector<Mdouble> psdRatio)
{
    particleSizeDistributionVector_.resize(psd.size());
    particleSizeDistributionVector_ = psd;
    logger.assert_always(std::accumulate(psdRatio.begin(), psdRatio.end(), 0.0) == 1.0, "Please make sure that the sum of "
                                                                                 "psdRatios adds up to unity");
    probability_.resize(psdRatio.size());
    probability_ = psdRatio;

}

/*!
 * \details gets the user defined particle size distributions
 * \return a vector of PSD class objects containing containing the differnt user defined particle size distributions
 */
std::vector<PSD> InsertionBoundary::getPSD()
{
    return particleSizeDistributionVector_;
}

/*!
 * \details sets the isManuallyInserting_ to TRUE, resulting in a top-down class-by-class insertion routine to insert PSDs
 * as accurate as possible.
 */
void InsertionBoundary::setManualInsertion(bool isManuallyInserting)
{
    isManuallyInserting_ = isManuallyInserting;
}

/*!
 * \details Gets the variable that checks if a particle has an interaction with a wall or another particle.
 * \return if TRUE the particle has an interaction and if FALSE the particle has no interaction.
 */
bool InsertionBoundary::getCheckParticleForInteraction() const
{
    return checkParticleForInteraction_;
}

/*!
 * \details Sets the distribution type from the Distribution class.
 * \param[in] checkParticleForInteraction      boolean which determines if a particle has an interaction.
 */
void InsertionBoundary::setCheckParticleForInteraction(bool checkParticleForInteraction)
{
    checkParticleForInteraction_ = checkParticleForInteraction;
}

// /*!
// * \details write distribution type variables to file.
// */
//std::ostream& operator<<(std::ostream& os, InsertionBoundary::Distribution type)
//{
//    os << static_cast<unsigned>(type);
//    return os;
//}
//
// /*!
// * \details write Distribution type variables from file.
// */
//std::istream& operator>>(std::istream& is, InsertionBoundary::Distribution& type)
//{
//    unsigned uType;
//    is >> uType;
//    type = static_cast<InsertionBoundary::Distribution>(uType);
//    return is;
//}

