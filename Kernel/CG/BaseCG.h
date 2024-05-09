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

#ifndef BaseCG_H
#define BaseCG_H

#include <cstdlib>
#include <functional>
#include "GeneralDefine.h"
#include "BaseObject.h"
#include "File.h"
#include "Math/Vector.h"

class CGHandler;

class BaseInteractable;

/*!
 * \brief Base class of all CG objects, needed to store the various CG objects 
 * in the CGHandler.
 * \details
 * The CGHandler contains a set of BaseCG* pointers, which can contain different 
 * CG objects.
 * 
 * This class contains properties that any CG class needs:
 * - a pointer to a handler_ object with set and get functions
 * - the variables and functionality inherited from BaseObject
 * - initialize, evaluate, and finalize functions called by the CGHandler
 * - parameters that specify the options to initialize the mesh of CGPoints 
 * (min_, max_, nX_, nY_, nZ_) and properties of the CGFunction (width_)
 * - parameters that specify when evaluate will be called (timeMin_, timeMax_)
 * - the statFile into which statistical output is written.
 */
class BaseCG : public BaseObject
{
public:
    
    /*!
     * \brief Simple constructor, sets default values
     */
    BaseCG();
    
    /*!
    * \brief Default copy constructor, copies all values
    */
    BaseCG(const BaseCG& p) = default;
    
    /*!
     * \brief Default destructor, does nothing
     */
    ~BaseCG() override = default;
    
    /*!
     * \brief Currently, no read functions are implemented for the CGHandler, 
     * but the function is required for any derivative of BaseObject.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes class content into an output stream, usually a stat file
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief This class seems to have no use (?), but is required for any 
     * derivative of BaseObject.
     */
    void clear();
    
    /*!
     * \brief Copy operator. Required for BaseHandler::copyAndAddObject
     */
    virtual BaseCG* copy() const = 0;
    
    /*!
     * \brief Called at the beginning of the DPM simulation to initialise the cg 
     * evaluation and to open the statFile.
     */
    virtual void initialise() = 0;
    
    /*!
     * \brief Called after a given number of time steps (statFile::saveCount_)
     * to evaluate the CG fields.
     */
    virtual void evaluate() = 0;
    
    /*!
     * \brief Called at the end of the DPM simulation to finish the cg 
     * evaluation and to close the statFile.
     */
    virtual void finish() = 0;
    
    /*!
     * \brief Sets handler_, the pointer to the CGHandler.
     */
    void setHandler(CGHandler* handler);
    
    /*!
     * \brief Returns handler_, a pointer to the CGHandler.
     */
    CGHandler* getHandler() const;
    
    /*!
     * \brief Sets width_, the coarse-graining width.
     * \todo should be standard deviation, but is currently cutoff.
     */
    virtual void setWidth(Mdouble width)=0;
    
    /*!
     * \brief Returns width_, the coarse-graining width.
     */
    virtual Mdouble getWidth() const =0;
    
    /*!
     * \brief Sets nZ_, the number of spatial mesh points in the z-direction.
     */
    void setNZ(std::size_t nZ);
    
    /*!
     * \brief Returns nZ_, the number of spatial mesh points in the z-direction.
     */
    std::size_t getNZ() const;
    
    /*!
     * \brief Sets nY_, the number of spatial mesh points in the y-direction.
     */
    void setNY(std::size_t nY);
    
    /*!
     * \brief Returns nY_, the number of spatial mesh points in the y-direction.
     */
    std::size_t getNY() const;
    
    /*!
     * \brief Sets nX_, the number of spatial mesh points in the x-direction.
     */
    void setNX(std::size_t nX);
    
    /*!
     * \brief Returns nX_, the number of spatial mesh points in the x-direction.
     */
    std::size_t getNX() const;
    
    /*!
     * \brief Sets nX_, nY_, nZ_, the number of spatial mesh points in each cartesian direction.
     */
    void setN(std::size_t n);
    
    /*!
     * \brief Sets nX_, nY_, nZ_, the number of spatial mesh points in each cartesian direction.
     */
    void setN(std::array<std::size_t, 3> n);
    
    /*!
     * \brief Sets nX_, nY_, nZ_, the number of spatial mesh points in each cartesian direction.
     * However, instead of explicitly defining n, the mesh size h=(max-min)/n is defined
     */
    void setH(Mdouble h);
    
    /*!
     * \brief Sets nX_ the number of spatial mesh points in the X-direction.
     * Instead of explicitly defining nX, the mesh size hX=(max.X-min.X)/nX is defined
     */
    void setHX(Mdouble h);
    
    /*!
     * \brief Sets nX_ the number of spatial mesh points in the X-direction.
     * Instead of explicitly defining nX, the mesh size hX=(max.X-min.X)/nX is defined
     */
    void setHY(Mdouble h);
    
    /*!
     * \brief Sets nX_ the number of spatial mesh points in the X-direction.
     * Instead of explicitly defining nX, the mesh size hX=(max.X-min.X)/nX is defined
     */
    void setHZ(Mdouble h);
    
    /*!
     * \brief Sets timeMin_, the lower limit of the temporal domain.
     */
    void setTimeMin(Mdouble timeMin);
    
    /*!
     * \brief Sets timeMax_, the upper limit of the temporal domain.
     */
    void setTimeMax(Mdouble timeMax);
    
    /*!
     * \brief Returns timeMin_, the lower limit of the temporal domain.
     */
    Mdouble getTimeMin() const;
    
    /*!
     * \brief Returns timeMax_, the upper limit of the temporal domain.
     */
    Mdouble getTimeMax() const;
    
    /*!
     * \brief Sets max_, the lower limit of the spatial domain.
     */
    void setMin(Vec3D min);
    
    /*!
     * \brief Sets min_.X, max_.X, the limits of the spatial domain in X.
     */
    void setX(Mdouble min, Mdouble max);
    
    /*!
     * \brief Sets min_.Y, max_.Y, the limits of the spatial domain in Y.
     */
    void setY(Mdouble min, Mdouble max);
    
    /*!
     * \brief Sets min_.Z, max_.Z, the limits of the spatial domain in Z.
     */
    void setZ(Mdouble min, Mdouble max);
    
    void setXGrid(Mdouble min, Mdouble max, Mdouble h);
    
    void setYGrid(Mdouble min, Mdouble max, Mdouble h);
    
    void setZGrid(Mdouble min, Mdouble max, Mdouble h);
    
    void setGrid(Vec3D min, Vec3D max, Mdouble h);
    
    /*!
     * \brief Sets max_, the upper limit of the spatial domain.
     */
    void setMax(Vec3D max);
    
    /*!
     * \brief Returns min_, the lower limit of the spatial domain.
     */
    Vec3D getMin() const;
    
    /*!
     * \brief Returns max_, the upper limit of the spatial domain.
     */
    Vec3D getMax() const;
    
    /*!
     * sets selectedParticle_ such that only particles of a certain species are selected
     * @param speciesIndex
     */
    void selectSpecies(unsigned speciesIndex);
    
    void setSelectedParticle(const std::function<const bool(const BaseInteractable*)>& selectedParticle);
    
    void setEps(Mdouble eps);
    
    Mdouble getEps() const;

    void setAverageBeyondDomain(const bool val) {averageBeyondDomain_=val;}

    bool getAverageBeyondDomain() const {return averageBeyondDomain_;}
    
    void setVerbose(const bool verbose) {verbose_=verbose;}
    
    bool getVerbose() const {return verbose_;}
    
    virtual void setWidthTime(Mdouble widthTime);
    
    virtual Mdouble getWidthTime() const;
    
    virtual void setTimeStep(Mdouble timeStep);
    
    virtual Mdouble getTimeStep() const;
    
    /*
     * Sets width such that the CG function has a fixed standard deviation
     * See CGStandardDeviationUnitTest.
     */
    virtual void setStandardDeviation(Mdouble std) = 0;
    
    /*
     * Sets width such that the CG function has the same standard deviation as a spherical cg function.
     * See CGStandardDeviationUnitTest.
     */
    virtual void setRadius(Mdouble radius) = 0;

protected:
    
    /*!
     * the pointer to the CGHandler, used to get data from the CGHandler.
     */
    CGHandler* handler_; ///
    
    /*!
     * nX_, nY_, and nZ_ define the size of the mesh of CGPoints, i.e. the 
     * spatial positions at which the cg variables are evaluated. 
     * nX_, nY_, and nZ_ are the number of points in x-, y-, and z-direction.
     * the is used to set the width of the
     * These values are set to 1 by default, unless modified by the user.
     */
    std::size_t nX_;
    
    /*!
     * see nX_
     */
    std::size_t nY_;
    
    /*!
     * see nZ_
     */
    std::size_t nZ_;
    
    /*!
     * Finite difference step size  used to computed derivatives
     * of CG functions
     */
    Mdouble eps_;
    
    /*!
     * timeMin_ and timeMax_ define the temporal dimensions of the coarse-graining
     * volume; these parameters are set to \f$\pm\infty\f$ by default, unless 
     * modified by the user; otherwise, timeMin_ is set to DPMBase::time_ in 
     * CG::initialize.
     */
    Mdouble timeMin_;
    
    /*!
     * see timeMin_
     */
    Mdouble timeMax_;
    
    /*!
     * min_ and max_ define the spatial dimensions of the coarse-graining volume;
     * if these parameters are not defined by the user, they are set to the 
     * system dimensions (DPMBase::min_ and DPMBase::max_) in CG::initialize.
     */
    Vec3D min_;
    
    /*!
     * see min_.
     */
    Vec3D max_;
    
    /*!
     * A function returning true for each particle that should be included in the statistics (all by default).
     */
    std::function<bool(const BaseInteractable*)> selectedParticle_;

    /**
     * \brief Determines whether particles outside the domain are considered when computing the averaged fields
     * \details
     * If set to true, then the average of field f over a coordinate direction x is computed as
     *  \f$\frac{1}{x_{max}-x_{min}}\int_{-\infty}^{\infty} f dx\f$
     * If set to false, then the average of field f over a coordinate direction x is computed as
     *  \f$\frac{1}{x_{max}-x_{min}}\int_{x_{min}}^{x_{max}} f dx\f$
     *
     * \todo should the default be false?
     * \todo currently, the above description is not implemented; it simply ignores particles outside the domain.
     */
    bool averageBeyondDomain_ = true;
    
    bool verbose_ = false;

public:
    /*!
     * \brief File class to handle the output into a .stat file
     */
    File statFile;
};

#endif
