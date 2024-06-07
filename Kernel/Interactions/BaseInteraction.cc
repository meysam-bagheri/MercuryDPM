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


#include "BaseInteraction.h"
#include "InteractionHandler.h"
#include "DPMBase.h"

/*!
 * \details This is the constructor which creates a new interactions between two
 * BaseInteractable objects. The timeStamp is time the interactions is created
 * and is used to check if the interations is current or has ended.
 * It adds
 * \param[in] P         BaseInteractable pointer which is the first object involved in the interaction normally a particle.
 * \param[in] I         BaseInteractable pointer which is the second object involved in the interaction often a wall or particle.
 * \param[in] timeStamp Mdouble which is the time the interaction starts.
 */
BaseInteraction::BaseInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseObject()
{
    P_ = P;
    I_ = I;
    normal_.setZero();
    contactPoint_.setNaN();
    overlap_ = 0;
    timeStamp_ = timeStamp;
    /// \bug Why is the species set to zero here and not the correct mixed type.
    species_ = nullptr;
    force_.setZero();
    torque_.setZero();
    if (P->getSpecies()->getHandler()->getDPMBase()->getInteractionFile().getFileType() == FileType::ONE_FILE)
    {
        writeInteraction(P->getSpecies()->getHandler()->getDPMBase()->getInteractionFile().getFstream(), true);
    }
    lagrangeMultiplier_ = 0;

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"BaseInteraction::BaseInteraction() finished"<<std::endl;
#endif
}


BaseInteraction::BaseInteraction()
        : BaseObject()
{
    P_ = nullptr;
    I_ = nullptr;
    normal_.setZero();
    overlap_ = 0;
    timeStamp_ = 0;
    species_ = nullptr;
    force_.setZero();
    torque_.setZero();
    lagrangeMultiplier_ = 0;
    contactPoint_.setNaN();
}

/*!
 * \details This an copy constructor for a BaseInteraction.
 * \param[in] p BaseInteraction
 */
BaseInteraction::BaseInteraction(const BaseInteraction& p)
        : BaseObject(p)
{
    P_ = p.P_;
    I_ = p.I_;
    normal_ = p.normal_;
    overlap_ = p.overlap_;
    force_ = p.force_;
    torque_ = p.torque_;
    species_ = p.species_;
    timeStamp_ = p.timeStamp_;
    lagrangeMultiplier_ = p.lagrangeMultiplier_;
    contactPoint_ = p.contactPoint_;
    
    ///\todo why are not all of the member variables copied?
//    InteractionHandler* handler_;
//    BaseInteractable* P_;
//    BaseInteractable* I_;
//    unsigned int identificationP_;
//    int identificationI_;
//    bool isWallInteraction_;
//    Vec3D contactPoint_;
//    Vec3D relativeVelocity_;
//    Mdouble normalRelativeVelocity_;
//    Mdouble absoluteNormalForce_;
//    Mdouble distance_;
//    Vec3D force_;
//    Vec3D torque_;
//    unsigned timeStamp_;
//    Vec3D normal_;
//    Mdouble overlap_;
//    const BaseSpecies* species_;
//    Mdouble lagrangeMultiplier_;
//    unsigned multiContactIdentifier_;

}

/*!
 * \details Destructor for BaseInteraction. Also removes the interaction from the
 *          list of interactions for both objects involved in the interaction.
 */
BaseInteraction::~BaseInteraction()
{
#if MERCURYDPM_USE_MPI
    if (P_ == nullptr)
    {
        logger(DEBUG,"Destroying a fictitious interaction used in MPI transmissions");
    }
    else
    {
#endif
    logger.assert_debug(P_ != nullptr, "Trying to destroy an interaction with P_ = nullptr");
    logger.assert_debug(I_ != nullptr, "Trying to destroy an interaction with I_ = nullptr");
    File& interactionFile = getHandler()->getDPMBase()->getInteractionFile();
    if (interactionFile.getFileType() == FileType::ONE_FILE)
    {
        writeInteraction(interactionFile.getFstream(), false);
    }
    
    P_->removeInteraction(this);
    I_->removeInteraction(this);
#if MERCURYDPM_USE_MPI
    }
#endif
}


/*!
 * \details BaseInteaction write function. Writes out all the information
 *          required to recreate this interaction.
 *          To write this interaction to the screen call write(std::cout).
 *          See also BaseInteraction::read
 * \param[in] os     std::ostream to which the information is written. Note, is
 *                   any ostram is can be file or screen.
 */
void BaseInteraction::write(std::ostream& os) const
{
    os << getName();
    if (dynamic_cast<BaseParticle*>(I_) != nullptr)
    {
        os << " particleIds " << P_->getId() << " " << I_->getId();
        ///\todo should we output id's here? os << " id " << getId() << " particleIds " << P_->getId() << " " << I_->getId();
    }
    else
    {
        os << " particleWallIds " << P_->getId() << " " << I_->getId();
    }
    os << " timeStamp " << timeStamp_;
    os << " contactPoint " << contactPoint_;
    os << " normal " << normal_;
    os << " distance " << distance_;
    os << " overlap " << overlap_;
    os << " force " << force_;
    os << " torque " << torque_;
}

/*!
 * \details BaseInteaction read functions. Reads in all the information about an
 *          interaction.
 *          Note, this can be from any istream but would normally be a file
 *          See also BaseInteaction::write
 * \param[in] is     std::istream to which the information is read from.
 */
void BaseInteraction::read(std::istream& is)
{
    //the rest gets read by the interaction handler
    std::string dummy;
    helpers::readOptionalVariable(is, "contactPoint", contactPoint_);
    helpers::readOptionalVariable(is, "normal", normal_);
    helpers::readOptionalVariable(is, "distance", distance_);
    helpers::readOptionalVariable(is, "overlap", overlap_);
    is >> dummy >> force_;
    is >> dummy >> torque_;
}

/*!
 * \details Functions which returns the name of the Interaction here is called
 *          BaseInteraction; but, it should be later overridden by the actual
 *          interaction classes.
 */
std::string BaseInteraction::getName() const
{
    return "BaseInteraction";
}

/*!
 * \details sets the normal of the interaction, in direction from I to P.
 * Must be a unit normal vector. This is not checked by the class.
 * \param[in] normal Vec3D which is the normal of the interaction.
 */
void BaseInteraction::setNormal(Vec3D normal)
{
    normal_ = normal;
}

/*!
 * \details set the distance of the interaction.
 * \param[in] distance  Mdouble which is the distance to set.
 */
void BaseInteraction::setDistance(Mdouble distance)
{
    distance_ = distance;
}

/*!
 * \details set the overlap between the two interactable object involved in the
 *          interactions.
 * \param[in] overlap   Mdouble which is the overlap to set.
 */
void BaseInteraction::setOverlap(Mdouble overlap)
{
    overlap_ = overlap;
}

/*
 * \details set the contact point between the two interactable objects involved
 * \param[in] contactPoint  Vec3D vector which will become the new contact point.
 */
void BaseInteraction::setContactPoint(Vec3D contactPoint)
{
    contactPoint_ = contactPoint;
}

/*!
 * \details Updated the time stamp of the interaction. The time stamp being old
 *          is how ended interactions are detected.
 * \param[in] timeStamp The new timeStamp for the interactions should be the
 *                      current time step.
 */
void BaseInteraction::setTimeStamp(unsigned timeStamp)
{
    timeStamp_ = timeStamp;
}

/*!
 * \details set the handler which this interactions belongs to.
 * \param[in] handler   InteractionHandler* pointer to the interaction handler,
 *                      this interaction will belong.
 */
void BaseInteraction::setHandler(InteractionHandler* handler)
{
    handler_ = handler;
}


InteractionHandler* BaseInteraction::getHandler() const
{
    return handler_;
}

/*!
    * \details Converts the linear overlap in the correct volume overlap i.e. this is the spherical cap approximation
    * Please see additional docs directory for the calculation of this and the mathematica file with the calculation in
    * @return The volume of the overlap between the two objects in the interaction
    */
Mdouble BaseInteraction::getOverlapVolume() const
{
    //logger(WARN,"This function is only tested for spheres and (may) give the wrong answer for other objects");
    //First get overlap and radii of each object
    Mdouble delta = getOverlap();
    Mdouble curvatureParticle = getP()->getCurvature(getP()->getPosition());
    Mdouble curvatureInteractable = getI()->getCurvature(getI()->getPosition());
    
    if (curvatureInteractable == 0)
    {
        logger(WARN, "This function is not implemented for walls yet and will return 0");
        return 0;
    }
    else
    {
        Mdouble r1 = 1.0 / curvatureInteractable;
        Mdouble r2 = 1.0 / curvatureParticle;
        
        if (delta > 2.0 * r1 || delta > 2.0 * r2)
            logger(WARN,
                   "Warning in BaseInteraction::getOverlapVolume: One of the particles is fully inside the other. The calculations are undefined for this case and will return nonsensical values.");
        
        //Eq (2.1) - see doc
//        logger(WARN,"Radius of first object % and second object %",r1,r2);
        Mdouble t1 = (r2 - delta / 2.0) / (r1 + r2 - delta) * delta;
        //Eq (2.2) - see doc
        Mdouble t2 = (r1 - delta / 2.0) / (r1 + r2 - delta) * delta;
        
        //Eq (2.3) - see doc
        Mdouble V = constants::pi * (
                std::pow(t1, 2.0) * r1
                - std::pow(t1, 3.0) / 3.0
                + std::pow(t2, 2.0) * r2
                - std::pow(t2, 3.0) / 3.0
        );
        
        return V;
    }
    
    return 0;
}

/*!
 * \details Removes the interaction from its InteractionHandler. Does no other
 *          cleaning up as it does not remove it from the particles.
 */
void BaseInteraction::removeFromHandler()
{
    handler_->removeObject(getIndex());
}

/*!
 * \details This functions copies the interactions of a real original particle.
 *          It first works out which of P and I is not the original particle.
 *          Then it create a new interactions between the new ghost copy and
 *          which every object is not the original particle from the P and I of
 *          the interaction.
 *          Note, at the end the ghost will be I in the new interaction and
 *          original item being interacted with will be P.
 * \todo    Can this be simpler if we replace the particle with the ghost.
 * \param[in] original  BaseInteractable* to the original particles who periodic
 *                      image is being created from.
 * \param[in] ghost     BaseInteractble* to the new ghost (periodic partcles)
 *                      that has just been created.
 */
void BaseInteraction::copySwitchPointer(const BaseInteractable* original, BaseInteractable* ghost) const
{
    //Copy the interaction of ghost
    BaseInteraction* C = copy();
    //Finds out which of P and I is that the particle from whom the ghost is begin created.
    //The object being interacted with is set to P
    if (C->getP() == original)
    {
        //Reverse some force history
        C->reverseHistory();
        //Set the P to the original particle
        C->P_ = C->getI();
    }
    //The new ghost particle is set to I in the interaction.
    C->I_ = ghost;
    
    //Add the the interaction to both original and the ghost
    handler_->addObject(C);
}

/*!
 * \details Returns a Mdouble to the current overlap for the objects involved in
 *          the interaction.
 * \return Mdouble which is the value of the overlap.
 */
Mdouble BaseInteraction::getContactRadius() const
{
    return getOverlap() < 0.0 ? 0.0 : sqrt(2.0 * getEffectiveRadius() * getOverlap());
}

/*!
 * \details Various variables in the force law need to be integrated. This is
 *          the place where this code goes.
 *          Note, it is empty at this point; it can be overriden in subclasses.
 *          For usage, see e.g. MindlinInteraction.cc.
 * \param[in] timeStep  The time-step dt.
 */
void BaseInteraction::integrate(Mdouble timeStep UNUSED)
{

}

/*!
 * \details Sets the species for the interactions.
 *          Note, this can be either a normal Species or a MixedSpecies;
 *          depending on if this interaction is between interactables of the
 *          same or different types.
 * \param[in] BaseSpecies* pointer to the actually species of the interaction.
 */
void BaseInteraction::setSpecies(const BaseSpecies* const species)
{
    species_ = species;
}

/*!
 * \details Changes the first object involved in the interaction; normally a
 *          particle.
 *          Note, set is slightly misleading as it removed the interaction from
 *          old particle and adds it to the new particle.
 * \param[in] P     BaseInteractable* An interactable object that is involved in
 *                  the interaction.
 */
void BaseInteraction::setP(BaseInteractable* P)
{
    P_->removeInteraction(this);
    P_ = P;
    P_->addInteraction(this);
    identificationP_ = P_->getId();
}

/*!
 * \details Set the first object involved in the interaction. This function differs from
 *          BaseInteraction::setP(BaseInteractable* P) because no interactable is removed:
 *          this is needed in DPMBase::importParticlesAs as in it there's no iteractable to remove
 *          and using BaseInteraction::setP(BaseInteractable* P) would led to an output printing
 *          overloaded with warnings.
 * \param[in] P     BaseInteractable* The particle involved in the interaction.
 */
void BaseInteraction::importP(BaseInteractable* P)
{
    P_->removeInteraction(this);
    P_ = P;
    P_->addInteraction(this);
    identificationP_ = P_->getId();
}

/*!
 * \details Changes the second object involved in the interaction; often a
 *          general interactable not always a particle.
 *          Note, set is slightly misleading as it removed the interaction from
 *          old interactable and adds it to the new interactable.
 * \param[in] I     BaseInteractable* An interactable object that is involved in
 *                  the interaction.
 */
void BaseInteraction::setI(BaseInteractable* I)
{
    I_->removeInteraction(this);
    I_ = I;
    I_->addInteraction(this);
    identificationI_ = I_->getId();
}

/*!
 * \details Changes the second object involved in the interaction. This function differs from
 *          BaseInteraction::setI(BaseInteractable* I) because no interactable is removed:
 *          this is needed in DPMBase::importParticlesAs as in it there's no iteractable to remove
 *          and using BaseInteraction::setI(BaseInteractable* I) would led to an output printing
 *          overloaded with warnings.
 * \param[in] I     BaseInteractable* The particle involved in the interaction.
 */
void BaseInteraction::importI(BaseInteractable* I)
{
    I_->removeInteraction(this);
    I_ = I;
    I_->addInteraction(this);
    identificationI_ = I_->getId();
}

Vec3D BaseInteraction::getIP() const
{
    //needs to be done by normal and distance or it fails in periodic systems
    //return getP()->getPosition() - getI()->getPosition();
    return getNormal() * getDistance();
}

Vec3D BaseInteraction::getIC() const
{
    ///\todo getDistance is not set in postprocessing
    return getNormal() * getDistance() + getContactPoint() - getP()->getPosition();
}

Vec3D BaseInteraction::getCP() const
{
    //works even in periodic systems becuase contact point is defined on the P-side
    return getP()->getPosition() - getContactPoint();
}

/*!
 * \details Writes the FStat information that is required for the coarse-
 *          graining package MercuryCG if you want stress and force information.
 *          Note, it takes a general ostream but is normally a file i.e.
 *          ofstream
 * \param[in] os    This is the ostream that the FStat information will be
 *                  written to. Normally, a file but could be a gerneral
 *                  ostream.
 */
void BaseInteraction::writeToFStat(std::ostream& os, Mdouble time) const
{
    ///\todo MX The documentation mentions that the first variable is the time - this is incorrect, is is the timeStamp the interaction started
    auto* IParticle = dynamic_cast<BaseParticle*>(I_);
    auto* PParticle = dynamic_cast<BaseParticle*>(P_);
    
    // do not write fstat output if the force is an internal bond
    if (isBonded()) return;
    
    Vec3D tangentialForce = getTangentialForce();
    Mdouble tangentialOverlap = getTangentialOverlap();
    
    Mdouble scalarNormalForce = Vec3D::dot(force_, getNormal());
    Mdouble scalarTangentialForce = tangentialForce.getLength();
    Vec3D tangential;
    if (scalarTangentialForce != 0.0)
        tangential = tangentialForce / scalarTangentialForce;
    else
        tangential = Vec3D(0.0, 0.0, 0.0);
    
    if (PParticle != nullptr && !PParticle->isFixed())
    {
        os << time << " " << P_->getIndex()
           << " " << static_cast<int>((IParticle == nullptr ? (-I_->getIndex() - 1) : I_->getIndex()))
           << " " << getContactPoint()
           << " " << getOverlap()
           << " " << tangentialOverlap
           << " " << scalarNormalForce
           << " " << scalarTangentialForce
           << " " << (IParticle == nullptr ? -normal_ : normal_)
           << " " << (IParticle == nullptr ? -tangential : tangential) << std::endl;
        ///\todo the flip in normal/tangential direction for walls should not be done; this is an old bug
    }
    if (IParticle != nullptr && !IParticle->isFixed() && IParticle->getPeriodicFromParticle() == nullptr)
    {
        os << time << " " << I_->getIndex()
           << " " << P_->getIndex()
           << " " << getContactPoint()
           << " " << getOverlap()
           << " " << tangentialOverlap
           << " " << scalarNormalForce
           << " " << scalarTangentialForce
           << " " << -normal_
           << " " << -tangential << std::endl;
    }
}

/*
 * \details Returns the distance between the two interactable objects involved
 *          in the interaction.
 * \return  Mdouble which is the distance between the two interacting objects.
 */
Mdouble BaseInteraction::getDistance() const
{
    return distance_;
}

/*!
 * \details Returns tangential overlap.
 *          Note, at this level there cannot be a tangential overlap
 *          hence by default it returns 0. This function will be overridden by
 *          interactions that have tangential components.
 * \return  Positive Mdouble that is the tangential overlap.
 */
Mdouble BaseInteraction::getTangentialOverlap() const
{
    return 0;
}

/*!
 * \details Returns the vector that is the tangential force
 *          Note, at this level there cannot be a tangential force therefore by
 *          default it return the zero vector. This function will be overridden
 *          by interactions that have tangential components.
 * \return  Vec3D that contains the current tangential force of the interaction.
 */
const Vec3D BaseInteraction::getTangentialForce() const
{
    return Vec3D(0.0, 0.0, 0.0);
}

/*!
 * \details Returns the relative velocity between the two interactable objects
 *          involved in the interaction.
 * \return  A reference to Vec3D that contains the relative velocity.
 */
const Vec3D& BaseInteraction::getRelativeVelocity() const
{
    return relativeVelocity_;
}

/*!
 * \details Returns the norm (length) of the relative normal velocity.
 *          Note, this can be negative or positive it is not a speed.
 *          \todo Ant: Check this comment.
 * \return  Mdouble that contains the norm (length) of the relative velocity.
 */
Mdouble BaseInteraction::getNormalRelativeVelocity() const
{
    return normalRelativeVelocity_;
}

/*!
 * \details Returns the absolute normal force. This is the magnitude of the normal
 *          force.
 *          \todo Ant: Check this comment.
 * \return  Mdouble that contains the absolute norm (length) of the normal
 *          force.
 */
Mdouble BaseInteraction::getAbsoluteNormalForce() const
{
    return absoluteNormalForce_;
}

/*!
 * \details add an increment to total force in the interaction. This is used by
 *          tangential and non-contact forces (e.g. adhesive forces) as this are
 *          'added' after the normal forces have been computed.
 */
void BaseInteraction::addForce(Vec3D force)
{
    force_ += force;
}

/*!
 * \details add an increment to total torque in the interaction. This is used by
 *          tangential and non-contact forces (e.g. adhesive forces) as this are
 *          'added' after the normal forces have been computed.
 */
void BaseInteraction::addTorque(Vec3D torque)
{
    torque_ += torque;
}

/*!
 * \details set the absolute values of the force. This is used by the normal
 *          forces as these are always called first and then the tangential and
 *          non-contact (e.g. adhesive forces) forces are added.
 *          See also BaseInteraction::addForce.
 */
void BaseInteraction::setForce(Vec3D force)
{
    force_ = force;
}

/*!
 * \details set the absolute values of the torque. This is used by the normal
 *          forces as these are always called first and then the tangential and
 *          non-contact (e.g. adhesive forces) forces/torques are added.
 *          See also BaseInteraction::addTorque.
 */
void BaseInteraction::setTorque(Vec3D torque)
{
    torque_ = torque;
}

/*!
 * \details set the relative velocity between the two particles involved in the
 *          interaction.
 * \param[in] relativeVelocity  This is Vec3D that contains the relative
 *                              velocity between the two interactable objects.
 */
void BaseInteraction::setRelativeVelocity(Vec3D relativeVelocity)
{
    relativeVelocity_ = relativeVelocity;
}

/*!
 * \details set the norm (length) of the normal relative velocity.
 * \param[in] normalRelativeVelocity    Mdouble containing the normal (length)
 *                                      of the normal velocity between the
 *                                      interactable objects.
 */
void BaseInteraction::setNormalRelativeVelocity(Mdouble normalRelativeVelocity)
{
    normalRelativeVelocity_ = normalRelativeVelocity;
}

/*!
 * \details set absolute normal force.
 * \param[in] absoluteNormalForce   Mdouble contain the value of the absolute
 *                                  normal force.
 */
void BaseInteraction::setAbsoluteNormalForce(Mdouble absoluteNormalForce)
{
    absoluteNormalForce_ = absoluteNormalForce;
}

/*!
 * \details Returns a BaseSpecies pointer to the current species.
 *          Note, this will be either a Species or a MixedSpecies done of which
 *          are derived from BaseSpecies.
 * \return  A BaseSpecies pointer to the species associated with this
 *          interaction.
 */
const BaseSpecies* BaseInteraction::getBaseSpecies() const
{
    return species_;
}

/*!
 * \details The children of this class will implement this function; however,
 *          it is blank.
 *          This function will do the actually force calculation for this
 *          interaction.
 *          Note, it is not virtual as it is not called from within this class.
 */
void BaseInteraction::computeForce()
{}

/*!
 * \details The children of this class will implement this function; however, it
 *          is blank.
 *          This function will contain the calculation for th elastic energy.
 *          Note, it is not virtual as it is not called from within this class.
 */
Mdouble BaseInteraction::getElasticEnergy() const
{
    return 0.0;
}

/*!
 * \details The children of this class will implement this function; however, it
 *          is blank.
 *          This function will contain code that changes some history
 *          information if a copy changes the P <-> I order; as can happen in
 *          creating a periodic particle.
 *          See also PeriodicBoundary::createPeriodicParticles and
 *          BaseInteraction::copySwitchPointer
 *          Note, it is not virtual as it is not called from within this class.
 */
void BaseInteraction::reverseHistory()
{
}

/*!
 * \details Method to write the interaction data to a stream, usually an interaction file.
 * It saves the IDs of the particles (in case there are two particles) or particle and wall.
 * Furthermore, the time when the interaction starts and ends is written to the stream.
 * \param[in/out] os The outputstream to which the interaction information should be written, usually a file
 * \param[in] created Whether or not this is the beginning of the interaction.
 */
void BaseInteraction::writeInteraction(std::ostream& os, bool created) const
{
    if (created)
    {
        os << "started ";
    }
    else
    {
        os << "ended   ";
    }
    
    if (dynamic_cast<BaseParticle*>(I_) != nullptr)
    {
        os << " particleIds     " << P_->getId() << " " << I_->getId() << " timeStamp ";
    }
    else
    {
        os << " particleWallIds " << P_->getId() << " " << I_->getId() << " timeStamp ";
    }
    
    if (created)
    {
        os << timeStamp_;
    }
    else
    {
        os << P_->getSpecies()->getHandler()->getDPMBase()->getTime();
    }
    
    os << std::endl;
    
}


unsigned int BaseInteraction::getMultiContactIdentifier() const
{
    return multiContactIdentifier_;
}

void BaseInteraction::setMultiContactIdentifier(unsigned int multiContactIdentifier)
{
    multiContactIdentifier_ = multiContactIdentifier;
}


void BaseInteraction::rotateHistory(Matrix3D& rotationMatrix)
{
    contactPoint_ = rotationMatrix * contactPoint_;
    relativeVelocity_ = rotationMatrix * relativeVelocity_;
    force_ = rotationMatrix * force_;
    torque_ = rotationMatrix * torque_;
    normal_ = rotationMatrix * normal_;
    ///\todo some of these might be unneccesary
}

/*!
 * \details Computes the effective radius of the two particles in the interaction. This is used by many of
 *          the later interaction models.
 *          This functions assumes P is the particle and I is either a particle
 *          or a wall.
 *          Effective Radius = \f$R_I*R_P/(R_I+R_P)\f$
 *          See also BaseInteraction::getEffectiveCorrectedRadius()
 * \return  A Mdouble which is the effective radius of the particles.
 */
Mdouble BaseInteraction::getEffectiveRadius() const
{
    Mdouble invEffectiveRadius = getP()->getCurvature(contactPoint_) + getI()->getCurvature(contactPoint_);
    logger.assert_debug(invEffectiveRadius > 0,
                        "getEffectiveRadius(): interaction % at % has infinite effective radius", getId(),
                        getContactPoint());
    return 1.0 / invEffectiveRadius;
}

/*!
 * \details Computes the effective mass of the particles in the interaction. This is used by many of
 *          the later interaction models.
 *          This functions assumes P is the particle and I is either a particle
 *          or a wall.
 *          Effective Radius = \f$R_I*R_P/(R_I+R_P)\f$
 *          See also BaseInteraction::getEffectiveCorrectedRadius()
 * \return  A Mdouble which is the effective radius of the particles.
 */
Mdouble BaseInteraction::getEffectiveMass() const
{
    Mdouble invEffectiveMass = getP()->getInvMass() + getI()->getInvMass();
    logger.assert_debug(invEffectiveMass > 0,
                        "getEffectiveMass(): interaction % at % has infinite effective mass", getId(),
                        getContactPoint());
    return 1.0 / invEffectiveMass;
}

void BaseInteraction::actionsAfterTimeStep()
{
}

/*!
 * \todo Thomas please document.
 */
void BaseInteraction::gatherContactStatistics()
{
    auto* IParticle = dynamic_cast<BaseParticle*>(I_);
    auto* PParticle = static_cast<BaseParticle*>(P_);
    
    Vec3D tangentialForce = getTangentialForce();
    Mdouble tangentialOverlap = getTangentialOverlap();
    
    Mdouble scalarNormalForce = Vec3D::dot(force_, getNormal());
    Mdouble scalarTangentialForce = tangentialForce.getLength();
    Vec3D tangential;
    if (scalarTangentialForce != 0.0)
        tangential = tangentialForce / scalarTangentialForce;
    else
        tangential = Vec3D(0.0, 0.0, 0.0);
    
    ///\todo TW centre is used just for backward compatibility; replace centre by contact law; we also have to correct it in StatisticsVector::gatherContactStatistics.
    ///There also seems to be an issue with the normal being defined differently for walls
    Vec3D centre;
    if (IParticle != nullptr)
        centre = getP()->getPosition() - normal_ * (PParticle->getRadius() + IParticle->getRadius() - overlap_) / 2.0;
    else
        centre = getP()->getPosition() - normal_ * (PParticle->getRadius() - overlap_);
    
    if (!PParticle->isFixed())
    {
        getHandler()->getDPMBase()->gatherContactStatistics(
                P_->getIndex(),
                static_cast<int>((IParticle == nullptr ? (-I_->getIndex() - 1) : I_->getIndex())),
                centre,
                getOverlap(),
                tangentialOverlap,
                scalarNormalForce,
                scalarTangentialForce,
                (IParticle == nullptr ? -normal_ : normal_),
                (IParticle == nullptr ? -tangential : tangential));
    }
    if (IParticle != nullptr && !IParticle->isFixed() && IParticle->getPeriodicFromParticle() == nullptr)
    {
        getHandler()->getDPMBase()->gatherContactStatistics(
                I_->getIndex(),
                static_cast<int>(P_->getIndex()),
                centre,
                getOverlap(),
                tangentialOverlap,
                scalarNormalForce,
                scalarTangentialForce,
                -normal_,
                -tangential);
        
    }
}

unsigned BaseInteraction::getNumberOfFieldsVTK() const
{
    return 0;
}

std::string BaseInteraction::getTypeVTK(unsigned i) const
{
    return "";
}

std::string BaseInteraction::getNameVTK(unsigned i) const
{
    return "";
}

std::vector<Mdouble> BaseInteraction::getFieldVTK(unsigned i) const
{
    return std::vector<Mdouble>();
}


void BaseInteraction::createMPIType()
{
}

void BaseInteraction::getMPIInteraction(void* historyDataArray, unsigned int index) const
{
}

void BaseInteraction::setMPIInteraction(void* historyDataArray, unsigned int index, const bool resetPointers)
{
}

void* BaseInteraction::createMPIInteractionDataArray(unsigned int numberOfInteractions) const
{
    logger(ERROR, "BaseInteraction::createMPIInteractionDataArray should never be called");
    return nullptr;
}

void BaseInteraction::deleteMPIInteractionDataArray(void* dataArray)
{
    logger(WARN, "Why on earth is this function called?");
}

void
BaseInteraction::getInteractionDetails(void* interactionDataArray, unsigned int index, unsigned int& identificationP,
                                       unsigned int& identificationI, bool& isWallInteraction, unsigned& timeStamp)
{
    logger(ERROR, "Something went wrong, this function should not be called");
}


void BaseInteraction::setIdentificationP(unsigned int identification)
{
    identificationP_ = identification;
}

void BaseInteraction::setIdentificationI(int identification)
{
    identificationI_ = identification;
}

unsigned int BaseInteraction::getIdentificationP()
{
    return identificationP_;
}

int BaseInteraction::getIdentificationI()
{
    return identificationI_;
}

void BaseInteraction::setWallInteraction(bool flag)
{
    isWallInteraction_ = flag;
}

bool BaseInteraction::isWallInteraction()
{
    return isWallInteraction_;
}

///\todo TW should P, I be of type unsigned?
void BaseInteraction::setBasicMPIInteractionValues(int P, int I, unsigned timeStamp, Vec3D force, Vec3D torque,
                                                   bool isWallInteraction, const bool resetPointers)
{
    this->setIdentificationP(P);
    this->setIdentificationI(I);
    this->setTimeStamp(timeStamp);
    this->setForce(force);
    this->setTorque(torque);
    this->setWallInteraction(isWallInteraction);
    if (resetPointers)
    {
        this->I_ = nullptr;
        this->P_ = nullptr;
    }
}

//note centre is unused
void BaseInteraction::setFStatData(std::fstream& fstat, BaseParticle* P, BaseParticle* I)
{
    Mdouble overlap, tangentialOverlap, scalarNormalForce, scalarTangentialForce;
    Vec3D centre, normal, tangential;
    fstat >> centre >> overlap >> tangentialOverlap >> scalarNormalForce >> scalarTangentialForce >> normal
          >> tangential;
    const Vec3D force = scalarNormalForce * normal + scalarTangentialForce * tangential;
    setForce(force);
    setNormal(normal);
    setOverlap(overlap);
    const Mdouble radius = P->getRadius();
    const Vec3D branch = (radius - 0.5 * getOverlap()) * getNormal();
    setContactPoint(P->getPosition() - branch);
    setDistance(radius + I->getRadius() - getOverlap());
}

void BaseInteraction::setFStatData(std::fstream& fstat, BaseParticle* P, BaseWall* I)
{
    Mdouble overlap, tangentialOverlap, scalarNormalForce, scalarTangentialForce;
    Vec3D centre, normal, tangential;
    fstat >> centre >> overlap >> tangentialOverlap >> scalarNormalForce >> scalarTangentialForce >> normal
          >> tangential;
    const Vec3D force = scalarNormalForce * normal + scalarTangentialForce * tangential;
    //note walls are defined different than particles (with an extra minus)
    setForce(-force);
    setNormal(-normal);
    setOverlap(overlap);
    const Mdouble radius = P->getRadius();
    const Vec3D branch = (radius - 0.5 * getOverlap()) * getNormal();
    setContactPoint(P->getPosition() - branch);
    setDistance(radius - getOverlap());
}
