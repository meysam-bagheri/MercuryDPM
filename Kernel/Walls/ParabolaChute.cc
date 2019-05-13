#include "ParabolaChute.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"

ParabolaChute::ParabolaChute()
{
    l_ = 1.0;
    ws_ = 1.0;
    logger(DEBUG, "ParabolaChute() constructor finished");
}

ParabolaChute::ParabolaChute(const ParabolaChute& other) : BaseWall(other)
{
    l_ = other.l_;
    ws_ = other.ws_;
    logger(DEBUG, "ParabolaChute copy constructor finished");
}


ParabolaChute::ParabolaChute(Mdouble length, Mdouble widthscale)
{
    l_ = length;
    ws_ = widthscale;
    logger(DEBUG, "ParabolaChute(params) constructor finished");
}

ParabolaChute::~ParabolaChute()
{
    logger(DEBUG, "ParabolaChute destructor finished");
}

void ParabolaChute::set(Mdouble length, Mdouble widthscale)
{
    l_ = length;
    ws_ = widthscale;
}

/*
ParabolaChute* ParabolaChute::copy() const {
    return new ParabolaChute(*this);
}
*/

bool ParabolaChute::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    /* Define some shortcuts */
    Mdouble y0 = p.getPosition().Y;
    Mdouble z0 = p.getPosition().Z;
    Mdouble a = p.getRadius();
    
    /* First, check whether the particle is definitely out of contact with the
     * chute. This will be so if the particle's interaction radius is small and
     * the particle's z-position is high. */
    // TODO do this properly
    if (false)
    {
        return false;
    }
    
    /* If not, then use Newton's method to minimise (half of the squared) distance */
    Mdouble R; //TODO
    Mdouble alpha; //TODO
    Mdouble dz = p.getPosition().Z; //TODO
    Mdouble q; // current guess
    Mdouble dd; // derivative of half of the squared distance at current guess
    Mdouble ddd; // second derivative at current guess
    Mdouble q0 = dz / l_; // minimum of the parabolic part
    
    // Iterate according to Newton
    do
    {
        dd = -2 * pow(q, 3) / pow(ws_, 2) + (1 - 2 * z0 / ws_) * q - y0;
        ddd = 6 * pow(q / ws_, 2) + (1 - 2 * z0 / ws_);
        q -= dd / ddd;
    } while (fabs(dd / ddd) > 1e-14);
    
    // Check whether this location is actually on the parabolic chute, otherwise
    // a point collision with the end of the chute calculated
    
    Mdouble distanceSquared = pow(q - y0, 2) + pow(pow(q, 2) / ws_ - z0, 2);
    // TODO
    if (distanceSquared >= mathsFunc::square(p.getWallInteractionRadius(this)))
    {
        return false;
    }
    else
    {
        Vec3D ContactPoint;
        distance = sqrt(distanceSquared);
        // ContactPoint.X = // TODO
        normal_return = ContactPoint - p.getPosition();
        normal_return /= normal_return.getLength();
        return true;
    }
}

/* TODO - getInteractionWith, read, write, getName */
