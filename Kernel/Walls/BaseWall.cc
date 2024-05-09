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

#include <Particles/SuperQuadricParticle.h>
#include "BaseWall.h"
#include "DPMBase.h"

/*!
 * \details Default constructor for the \ref BaseWall class.
 * Simply creates an empty \ref BaseWall. Note that it also, by default, sets the
 * handler to a null pointer - i.e. does not automatically assign the current object
 * a given \ref WallHandler.
 */
BaseWall::BaseWall()
{
    handler_ = nullptr;
    logger(DEBUG, "BaseWall::BaseWall() finished");
}

/*!
 * \details An existing wall (i.e. a \ref BaseWall type object), w, is passed as an argument.
 * A copy of this object is then created. The \ref BaseWall class is relatively low-level, and this
 * copy constructor simply acts to provide a pointer to the \ref WallHandler belonging to w, i.e.
 * assigning the new wall to the same handler as the original. (A derived class' copy constructor
 * calls this, but does all of the other work.)
 * \param[in] w - The existing \ref BaseWall object to be copied.
 */
BaseWall::BaseWall(const BaseWall& w)
        : BaseInteractable(w)
{
    //sets the current handler to that belonging to the existing wall, w
    handler_ = w.handler_;
    renderedWalls_.reserve(w.renderedWalls_.capacity());
    vtkVisibility_ = w.vtkVisibility_;
    for (auto* const r : w.renderedWalls_)
        renderedWalls_.push_back(r->copy());
    logger(DEBUG, "BaseWall::BaseWall(const BaseWall &p) finished");
}

/*!
 * \details
 * Note that this is a <B>virtual destructor</B>, ensuring that derived classes can be deleted easily and safely.
 */
BaseWall::~BaseWall()
{
    logger(DEBUG, "BaseWall::~BaseWall() finished");
    removeRenderedWalls();
}

/*!
 * \details
 * The BaseWall takes no more information than for a \ref BaseInteractable.
 * (A derived class' read method does most of the work.)
 * \param[in] is - The input stream from which the BaseWall is read.
 */
void BaseWall::read(std::istream& is)
{
    BaseInteractable::read(is);
    unsigned size;
    std::string type;
    if (helpers::isNext(is,"vtkVisibility")) {
        is >> vtkVisibility_;
    }
    if (helpers::isNext(is,"renderedWalls")) {
        is >> size;
        for (unsigned i = 0; i < size; i++)
        {
            is >> type;
            BaseWall* wall = getHandler()->createObject(type);
            wall->setHandler(getHandler());
            wall->read(is);
            renderedWalls_.push_back(wall);
            renderedWalls_.back()->setId(renderedWalls_.size());
        }
    }
}

void BaseWall::write(std::ostream& os) const
{
    BaseInteractable::write(os);
    if (vtkVisibility_==false)
    {
        os << " vtkVisibility " << vtkVisibility_;
    }
    if (renderedWalls_.size()>0)
    {
        os << " renderedWalls " << renderedWalls_.size();
        for (auto w : renderedWalls_)
        {
            os << " ";
            w->write(os);
        }
    }
}

/*!
 * \details Setting the WallHandler also sets the DPMBase and therefore the SpeciesHandler for the
 * species. This wall's species pointer is updated to the new SpeciesHandler.
 * \param[in] handler - A pointer to the WallHandler that we want to handle this wall.
 *
 *
 */
void BaseWall::setHandler(WallHandler* handler)
{
    handler_ = handler;
    setSpecies(getHandler()->getDPMBase()->speciesHandler.getObject(getIndSpecies()));
}

/*!
 * \return A pointer to the WallHandler that manages this BaseWall.
 */
WallHandler* BaseWall::getHandler() const
{
    return handler_;
}

/*!
 * \param[in] indSpecies The index of the species of this BaseWall in the SpeciesHandler.
 */
void BaseWall::setIndSpecies(unsigned int indSpecies)
{
    if (handler_ != nullptr)
    {
        setSpecies(getHandler()->getDPMBase()->speciesHandler.getObject(getIndSpecies()));
    }
    else
    {
        BaseInteractable::setIndSpecies(indSpecies);
        logger(ERROR, "setIndSpecies called on a particle with no particle handler.\n"
                      "Therefore I can't request the given species from the species handler.\n"
                      " PartID = %", getId());
    }
}

/*!
 * \details Firstly, calls the "BaseInteractable" version of the setSpecies function (\ref  BaseInteractable::setSpecies())
 * which actually sets the \ref species_ of the particle (i.e. provides a pointer to the \ref ParticleSpecies which stores the relevant
 * material properties which we wish to assign to out wall) as well as assigning the relevant species index (\ref indSpecies_)
 * as well as performing relevant error checks.
 *
 * Secondly, sets a pointer to the relevant handler, which is needed to retrieve species information.
 *
 * \param[in] species - The pointer to the species whose properties we want to assign to the current wall - i.e. the species we want to "give" the wall.
 *
 * \todo TW: this function should also check if the particle is the correct particle for the species type.
 */
void BaseWall::setSpecies(const ParticleSpecies* species)
{
    BaseInteractable::setSpecies(species);
    //Checks if the current wall currently belongs to a defined wallHandler...
    if (getHandler() == nullptr)
    {
        //and if not:
        //Creates a pointer to the handler to which the chosen species belongs
        SpeciesHandler* sH = species->getHandler();
        //Creates a pointer back to the DPMBase class
        DPMBase* dB = sH->getDPMBase();
        //If this is not a null pointer, assigns the current wall to the relevant wallHandler.
        if (dB != nullptr)
            setHandler(&dB->wallHandler);
    }
}

bool BaseWall::isFixed() const
{
    return true;
}

/*!
 * \param[in] forceGoal - the desired force on the wall
 * \param[in] gainFactor - the rate at which the velocity of the wall should be adjusted
 * \param[in] baseVelocity - the velocity that the wall should travel at when the forceGoal is reached
 */
void BaseWall::setForceControl(Vec3D forceGoal, Vec3D gainFactor, Vec3D baseVelocity) {
    setPrescribedVelocity([this, forceGoal, gainFactor, baseVelocity] (double time){
        auto dForce = getForce()-forceGoal;
        return baseVelocity + gainFactor.multiplyElementwise(dForce);
    });
}

// returns the point intersecting a wall (x-p).n=0 and a line x=p0+t(p1-p0)
bool BaseWall::getLinePlaneIntersect(Vec3D& intersect, const Vec3D& p0, const Vec3D& p1, const Vec3D& n, const Vec3D& p)
{
    // t = (p-p0).n / (p1-p0).n
    //first compute the denominator
    Mdouble denominator = Vec3D::dot(p1 - p0, n);
    if (fabs(denominator) >= 1e-10)
    {
        Mdouble t = Vec3D::dot(p - p0, n) / denominator;
        if (t < 1 + 1e-12 && t > -1e-12)
            intersect = p0 + t * (p1 - p0);
        return true;
    }
    return false;
}

//checks if point is in wall (if close to the wall, the point is assumed out of the wall)
bool BaseWall::isInsideWallVTK(const Vec3D& point, const Vec3D& normal, const Vec3D& position) const
{
    return Vec3D::dot(position - point, normal) < -1e-12;
}

/*!
 * intersectVTK = point[i] + t*(point[i-1]-point[i])
 * and (intersectVTK - position) *normal_=0.
 * => (position-point[i]-t*(point[i-1]-point[i]))*normal_=0
 * t=(position-point[i]))*normal_ / (point[i-1]-point[i]))*normal_
 */
void BaseWall::projectOntoWallVTK(Vec3D& point0, const Vec3D& point1, const Vec3D& normal, const Vec3D& position) const
{
    Vec3D dPoint = point1 - point0;
    point0 += Vec3D::dot(position - point0, normal) / Vec3D::dot(dPoint, normal) * dPoint;
}

/*!
 * Checks if a set of VTK points is inside a half-space defined by position and normal; all points outside the half-space are projected onto the half-space boundary.
 * Thus, if the old set of points represented a wall object, the new set of points is represents the intersection of the wall with the half-space.
 */
void BaseWall::intersectVTK(std::vector<Vec3D>& points, const Vec3D normal, const Vec3D position) const
{
    // find first point in Wall
    std::vector<Vec3D>::iterator firstIn = points.begin();
    for (; firstIn != points.end(); firstIn++)
    {
        //stop if points[first] is in domain
        if (isInsideWallVTK(*firstIn, normal, position))
        {
            break;
        }
    }
    
    //if all points are out of the wall
    if (firstIn == points.end())
    {
        logger(DEBUG, "BaseWall::intersectVTK: all points out of wall");
        return;
    }
    
    // find first point out of the wall after firstIn
    std::vector<Vec3D>::iterator firstOut = firstIn + 1;
    for (; firstOut != points.end(); firstOut++)
    {
        if (!isInsideWallVTK(*firstOut, normal, position))
        {
            break;
        }
    }
    
    //if all points are in the wall
    if (firstOut == points.end() && firstIn == points.begin())
    {
        logger(DEBUG, "BaseWall::intersectVTK: points completely in wall; removing points");
        points.clear();
        return;
    }
    
    //if the sequence starts with a point out of the wall
    //Several cases have to be distinguished, multiple points in the wall the last point in or out of the wall: ooiiioo, ooioo, ooiii, or ooi
    //In addition, we add the case iiioo and ioo
    if (firstIn != points.begin() || !isInsideWallVTK(points.back(), normal, position))
    {
        // remove unnessesary points in the wall 
        // ooiiioo -> ooiioo, ooiii -> ooii, iiioo -> iioo
        if (firstOut - firstIn > 2)
        {
            logger(DEBUG, "BaseWall::intersectVTK: remove unnessesary points in the wall");
            //necessary reset of firstOut, as erase invalidates iterators after erase point
            points.erase(firstIn + 1, firstOut - 1); //note: erase does not delete the last point
            firstOut = firstIn + 2;
        }
        
        // if there is only one point in the wall, make it two
        // ooioo -> ooiioo, ooi -> ooii
        if (firstOut == firstIn + 1)
        {
            logger(DEBUG, "BaseWall::intersectVTK: there is only one point in the wall, make it two");
            //necessary reset of firstIn, firstOut, as insert invalidates iterators
            unsigned in = firstIn - points.begin();
            points.insert(firstIn + 1, *firstIn);
            firstIn = points.begin() + in;//necessary, unless capacity is set right
            firstOut = firstIn + 2;
        }
        
        // three cases remain: ooiioo, ooii, iioo
        
        //move both points onto the surface of the wall
        if (firstIn != points.begin())
        {
            logger(DEBUG, "BaseWall::intersectVTK: move first point onto the surface of the wall");
            projectOntoWallVTK(*firstIn, *(firstIn - 1), normal, position);
        }
        else
        {
            logger(DEBUG,
                   "BaseWall::intersectVTK: move first point (at the beginning of the list) onto the surface of the wall");
            projectOntoWallVTK(points.front(), points.back(), normal, position);
        }
        
        if (firstOut != points.end())
        {
            logger(DEBUG, "BaseWall::intersectVTK: move second point onto the surface of the wall");
            projectOntoWallVTK(*(firstOut - 1), *firstOut, normal, position);
        }
        else
        {
            logger(DEBUG,
                   "BaseWall::intersectVTK: move second point (at the end of the list) onto the surface of the wall");
            projectOntoWallVTK(points.back(), points.front(), normal, position);
        }
        //if sequence starts and ends with a point in the wall: iiiooiii
    }
    else
    {
        logger(DEBUG, "BaseWall::intersectVTK: sequence starts and ends with a point in the wall");
        
        // find first point in wall after firstOut
        for (firstIn = firstOut + 1; firstIn != points.end(); firstIn++)
        {
            if (isInsideWallVTK(*firstIn, normal, position))
            {
                break;
            }
        }
        
        // remove unnessesary points in the wall 
        // iiiooiii -> iooi
        points.erase(firstIn + 1, points.end());
        points.erase(points.begin(),
                     firstOut - 1); //note: erase does not delete the last point //note iterators are invalid now
        
        //move both points onto the surface of the wall: iooi
        projectOntoWallVTK(points.front(), *(points.begin() + 1), normal, position);
        projectOntoWallVTK(points.back(), *(points.end() - 2), normal, position);
    }
}

/*!
 * \param[in] p Pointer to the BaseParticle which we want to check the interaction for.
 * \param[in] timeStamp The time at which we want to look at the interaction.
 * \param[in] interactionHandler A pointer to the InteractionHandler in which the interaction can be found.
 * \return A pointer to the BaseInteraction that happened between this BaseWall
 * and the BaseParticle at the timeStamp.
 */
BaseInteraction*
BaseWall::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    Mdouble overlap;

    if (getDistanceNormalOverlap(*p, distance, normal, overlap))
    {
        #ifdef MERCURYDPM_TRIANGLE_WALL_CORRECTION
        // look for an existing interaction, or create a new one
        BaseInteraction *c = nullptr;
        // This if-statement deals with groups of walls. If a particle has multiple contacts with a group of walls, and if the contact areas of these contacts overlap, then we keep only the biggest of the overlapping contacts.
        if (getGroupId() > 0) {
            // if there is a contact with a group of walls, and if p had at least one previously detected contact (in the last timestep or the current)
            for (const auto i : p->getInteractions()) {
                const BaseInteractable* q = i->getI();
                if (q->getGroupId() == getGroupId() and i->getTimeStamp() == timeStamp) {
                    if (not(isFaceContact(normal) and q->isFaceContact(i->getNormal())))
                    {
                        ///\todo face-face contacts should be treated differently
                        if (overlap > i->getOverlap()) {
                            p->addForce(-i->getForce());
                            this->addForce(i->getForce());
                            if (getHandler()->getDPMBase()->getRotation()) {
                                p->addTorque(-i->getTorque() +
                                             Vec3D::cross(p->getPosition() - i->getContactPoint(), i->getForce()));
                                this->addTorque(i->getTorque() -
                                                Vec3D::cross(this->getPosition() - i->getContactPoint(),
                                                             i->getForce()));
                            }
                            i->setI(this);
                            c = i;
                        } else {
                            return nullptr;
                        }
                    }
                }
            }
        }
        if (c == nullptr) {
            // look for an existing interaction, or create a new one
            c = interactionHandler->getInteraction(p, this, timeStamp);
        }
        #else
        BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
        #endif

        c->setNormal(-normal);
        c->setDistance(distance);
        c->setOverlap(overlap);
        if (p->isSphericalParticle())
        {
            ///\todo{DK: What is the contact point for interactions with walls}
            c->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
        }
        else
        {
            Vec3D normalBodyFixed = normal;
            p->getOrientation().rotateBack(normalBodyFixed);
            auto furthestPoint = getFurthestPointSuperQuadric(normalBodyFixed, p->getAxes(),
                                                              p->getExponentEps1(), p->getExponentEps2());
            Vec3D overlapBody = overlap * normalBodyFixed;
            Vec3D contactPoint = furthestPoint - overlapBody / 2;
            p->getOrientation().rotate(contactPoint);
            contactPoint += p->getPosition();
            c->setContactPoint(contactPoint);
        }
        logger(DEBUG, "Particle contact with wall at %", c->getContactPoint());
        return c;
    }
    return nullptr;
}

void BaseWall::writeVTK(VTKContainer& vtk) const
{
    logger(WARN, "Walls % of type % have no vtk writer defined", getIndex(), getName());
}

void BaseWall::addToVTK(const std::vector<Vec3D>& points, VTKContainer& vtk)
{
    if (!points.empty())
    {
        //all all values in myPoints to points
        vtk.points.insert(vtk.points.end(), points.begin(), points.end());
        
        // create one cell object containing all indices of the added points (created a triangle strip connecting these points)
        std::vector<double> cell;
        cell.reserve(vtk.points.size() + 1);
        cell.push_back(vtk.points.size() - 1);
        for (unsigned i = vtk.points.size() - points.size(); i < vtk.points.size(); i++)
        {
            cell.push_back((double) i);
        }
        
        //add this triangle strip to the vtk file
        vtk.triangleStrips.push_back(cell);
    }
}


///\todo make it work with screw, coil and other weird walls
BaseInteraction* BaseWall::getInteractionWithSuperQuad(SuperQuadricParticle* p, unsigned timeStamp,
                                                                    InteractionHandler* interactionHandler)
{
    logger(ERROR, "Generic wall-superquad interactions not implemented yet.");
    return nullptr;
}

/// \details This functions returns a axis for a wall using it Quaternion descriptions. At the moment Quaternion are not implemented for a wall; so this is currently a workaround for the non-implementation Quaternion for the walls. In the future this functions will be replaced.
/// \return A Vec3D which is the axis of the wall
const Vec3D BaseWall::getAxis() const
{
    Quaternion Q = getOrientation();
    Vec3D axis;
    axis.X = Q.q1;
    axis.Y = Q.q2;
    axis.Z = Q.q3;
    return axis;
}

bool BaseWall::getVTKVisibility() const
{
    return vtkVisibility_;
}

void BaseWall::setVTKVisibility(const bool vtkVisibility)
{
    vtkVisibility_ = vtkVisibility;
}

bool BaseWall::getDistanceNormalOverlap(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return,
                                        Mdouble& overlap) const
{
    if (P.isSphericalParticle())
    {
        bool isInContact = getDistanceAndNormal(P, distance, normal_return);
        overlap = P.getRadius() - distance;
        return isInContact;
    }
    else
    {
        auto superQuadric = dynamic_cast<const SuperQuadricParticle*>(&P);
        return getDistanceNormalOverlapSuperquadric(*superQuadric, distance, normal_return, overlap);
    }
}

bool
BaseWall::getDistanceNormalOverlapSuperquadric(const SuperQuadricParticle& p, Mdouble& distance, Vec3D& normal_return,
                                               Mdouble& overlap) const
{
    logger(ERROR, "Generic wall-superquadric interactions not implemented yet.");
    return false;
}

Vec3D BaseWall::getFurthestPointSuperQuadric(const Vec3D& normalBodyFixed, const Vec3D& axes, Mdouble eps1, Mdouble eps2) const
{
    logger(ERROR, "Generic wall-superquadric interactions not implemented yet.");
    return {};
}

//todo how do i write a copy and add function?
void BaseWall::addRenderedWall(BaseWall* w)
{
    renderedWalls_.push_back(w);
}

void BaseWall::removeRenderedWalls() {
    for (auto* const r : renderedWalls_)
        delete r;
    renderedWalls_.clear();
}

BaseWall* BaseWall::getRenderedWall(size_t i) const
{
    return renderedWalls_[i];
}

std::vector<BaseWall*> BaseWall::getRenderedWalls() const
{
    return renderedWalls_;
}

void BaseWall::renderWall(VTKContainer& vtk)
{
    if (getVTKVisibility())
    {
        if (renderedWalls_.empty())
        {
            writeVTK(vtk);
        }
        else
        {
            const Mdouble time = getHandler()->getDPMBase()->getTime();
            const Mdouble timeStep = getHandler()->getDPMBase()->getTimeStep();
            const int saveCount = getHandler()->getDPMBase()->dataFile.getSaveCount();
            for (const auto& r: renderedWalls_)
            {
                // First write...
                r->writeVTK(vtk);

                // Then update the position and orientation for time steps that were skipped.
                // We have to use a loop, as doing it at once like so:
                // r->integrateBeforeForceComputation(time - timeStep * (saveCount - 1), timeStep * saveCount);
                // r->integrateAfterForceComputation(time - timeStep * (saveCount - 1), timeStep * saveCount);
                // gives small deviations which add up quickly to give very different results.
                for (int i = 0; i < saveCount; i++)
                {
                    r->integrateBeforeForceComputation(time - timeStep * (saveCount - i - 1), timeStep);
                    r->integrateAfterForceComputation(time - timeStep * (saveCount - i - 1), timeStep);
                }
            }
        }
    }
}

void BaseWall::setVelocityControl(Vec3D forceGoal, Vec3D gainFactor, Vec3D baseVelocity) {
    setPrescribedVelocity([this, forceGoal, gainFactor, baseVelocity] (double time){
        auto dForce = getForce()-forceGoal;
        return baseVelocity + gainFactor.multiplyElementwise(dForce);
    });
}

void BaseWall::addParticlesAtWall(unsigned numElements)
{
    auto& speciesHandler = getHandler()->getDPMBase()->speciesHandler;
    logger.assert_always(speciesHandler.getSize()>0,"addParticlesAtWall: You need to define at least one species");

    Vec3D max = getHandler()->getDPMBase()->getMax();
    Vec3D min = getHandler()->getDPMBase()->getMin();
    double h = Vec3D::min(max-min)/numElements;
    double r = 0.5*h;

    auto& particleHandler = getHandler()->getDPMBase()->particleHandler;
    double numParticles0 = particleHandler.getSize();
    SphericalParticle p;
    p.setSpecies(speciesHandler.getObject(0));
    Vec3D pos;
    for (pos.X = min.X; pos.X <= max.X; pos.X += h)
    for (pos.Y = min.Y; pos.Y <= max.Y; pos.Y += h)
    for (pos.Z = min.Z; pos.Z <= max.Z; pos.Z += h)
    {
        Vec3D normal;
        Mdouble distance;
        p.setRadius(2.0*r);
        p.setPosition(pos);
        //if touching the wall
        if (getDistanceAndNormal(p, distance, normal) && distance>=0)
        {
            p.setRadius(r);
            p.setPosition(pos+(distance-r)*normal);
            particleHandler.copyAndAddObject(p);
        }
    }
    logger(INFO,"Inserted % particles that touch wall %", particleHandler.getNumberOfObjects()-numParticles0, getIndex());
}
