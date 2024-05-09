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

#include "BaseCG.h"
#include "DPMBase.h"

BaseCG::BaseCG()
{
    handler_ = nullptr;
    statFile.getFstream().precision(8);
    statFile.getFstream().setf(std::ios::left);
    statFile.setLastSavedTimeStep(NEVER);
    nX_ = 1;
    nY_ = 1;
    nZ_ = 1;
    timeMin_ = -constants::inf;
    timeMax_ = constants::inf;
    min_ = Vec3D(-constants::inf, -constants::inf, -constants::inf);
    max_ = Vec3D(constants::inf, constants::inf, constants::inf);
    selectedParticle_ = [](const BaseInteractable* p) { return true; };
    logger(DEBUG, "BaseCG::BaseCG() finished");
}

void BaseCG::clear()
{
    logger(WARN, "BaseCG::clear(), this function shouldn't be called");
}


void BaseCG::read(std::istream& is UNUSED)
{
}

/*!
 * \param[out] os output stream to which data is written
 */
void BaseCG::write(std::ostream& os) const
{
    //BaseObject::write(os);
    os << getName();
    os << " min " << min_;
    os << " max " << max_;
    if (!std::isinf(timeMin_)) os << " timeMin " << timeMin_;
    if (!std::isinf(timeMax_)) os << " timeMax " << timeMax_;
    os << " n " << nX_ << " " << nY_ << " " << nZ_;
    os << " width " << getWidth();
    if (verbose_) os << " verbose 1";
    //statFile
}

void BaseCG::setHandler(CGHandler* handler)
{
    handler_ = handler;
}

/*!
 * \return pointer to the CGHandler
 */
CGHandler* BaseCG::getHandler() const
{
#ifdef DEBUG_OUTPUT
    if (handler_ == nullptr)
    {
        std::cerr << "error: handler_==0 " << std::endl;
    }
#endif
    return handler_;
}

void BaseCG::setEps(Mdouble eps)
{
    eps_ = eps;
}

Mdouble BaseCG::getEps() const
{
    return eps_;
}

void BaseCG::setNZ(std::size_t nZ)
{
    nZ_ = nZ;
}

std::size_t BaseCG::getNZ() const
{
    return nZ_;
}

void BaseCG::setNY(std::size_t nY)
{
    nY_ = nY;
}

std::size_t BaseCG::getNY() const
{
    return nY_;
}

void BaseCG::setNX(std::size_t nX)
{
    nX_ = nX;
}

std::size_t BaseCG::getNX() const
{
    return nX_;
}

void BaseCG::setN(std::size_t n)
{
    nX_ = n;
    nY_ = n;
    nZ_ = n;
}

void BaseCG::setN(std::array<std::size_t, 3> n)
{
    nX_ = n[0];
    nY_ = n[1];
    nZ_ = n[2];
}

void BaseCG::setTimeMin(Mdouble timeMin)
{
    timeMin_ = timeMin;
}

void BaseCG::setTimeMax(Mdouble timeMax)
{
    timeMax_ = timeMax;
}

Vec3D BaseCG::getMin() const
{
    return min_;
}

Vec3D BaseCG::getMax() const
{
    return max_;
}

void BaseCG::setMin(Vec3D min)
{
    min_ = min;
}

void BaseCG::setMax(Vec3D max)
{
    max_ = max;
}

Mdouble BaseCG::getTimeMin() const
{
    return timeMin_;
}

Mdouble BaseCG::getTimeMax() const
{
    return timeMax_;
}

void BaseCG::setX(Mdouble min, Mdouble max)
{
    min_.X = min;
    max_.X = max;
}

void BaseCG::setY(Mdouble min, Mdouble max)
{
    min_.Y = min;
    max_.Y = max;
}

void BaseCG::setZ(Mdouble min, Mdouble max)
{
    min_.Z = min;
    max_.Z = max;
}

void BaseCG::setXGrid(Mdouble min, Mdouble max, Mdouble h)
{
    setX(min,max);
    setHX(h);
}

void BaseCG::setYGrid(Mdouble min, Mdouble max, Mdouble h)
{
    setY(min,max);
    setHY(h);
}

void BaseCG::setZGrid(Mdouble min, Mdouble max, Mdouble h)
{
    setZ(min,max);
    setHZ(h);
}

void BaseCG::setGrid(Vec3D min, Vec3D max, Mdouble h)
{
    setMin(min);
    setMax(max);
    setH(h);
}

void BaseCG::selectSpecies(unsigned speciesIndex)
{
    selectedParticle_ = [speciesIndex](const BaseInteractable* p)
    {
        return p->getIndSpecies() == speciesIndex;
    };
}

void BaseCG::setSelectedParticle(const std::function<const bool(const BaseInteractable*)>& selectedParticle)
{
    selectedParticle_ = selectedParticle;
}

void BaseCG::setH(Mdouble h)
{
    setHX(h);
    setHY(h);
    setHZ(h);
    logger(INFO, "min % max % h % nz %", min_, max_, h, nZ_);
}

void BaseCG::setHX(Mdouble h)
{
    logger.assert_debug(h > 0, "setHX(%): h has to be positive");
    logger.assert_debug(max_.X!=constants::inf && min_.X!=-constants::inf,
                  "setHX(%) can only be used after setting min and max values", h);
    setNX(static_cast<size_t>(std::ceil((max_.X - min_.X) / h)));
    logger.assert_always(getNX() > 0, "setHX(%) generated nX=% for %<x<%", h, getNX(), min_.X, max_.X);
}

void BaseCG::setHY(Mdouble h)
{
    logger.assert_debug(h > 0, "setHY(%): h has to be positive");
    logger.assert_debug(max_.Y!=constants::inf && min_.Y!=-constants::inf,
                  "setHY(%) can only be used after setting min and max values", h);
    setNY(static_cast<size_t>(std::ceil((max_.Y - min_.Y) / h)));
    logger.assert_always(getNY() > 0, "setHY(%) generated nY=% for %<y<%", h, getNY(), min_.Y, max_.Y);
}

void BaseCG::setHZ(Mdouble h)
{
    logger.assert_debug(h > 0, "setHZ(%): h has to be positive");
    logger.assert_debug(max_.Z!=constants::inf && min_.Z!=-constants::inf,
                  "setHZ(%) can only be used after setting min and max values", h);
    setNZ(static_cast<size_t>(std::ceil((max_.Z - min_.Z) / h)));
    logger.assert_always(getNZ() > 0, "setHZ(%) generated nZ=% for %<z<%", h, getNZ(), min_.Z, max_.Z);
}

void BaseCG::setWidthTime(Mdouble widthTime)
{
    logger(ERROR, "setWidthTime() should not be used, it is a dummy function for templating");
}

void BaseCG::setTimeStep(Mdouble timeStep)
{
    logger(ERROR, "setTimeStep() should not be used, it is a dummy function for templating");
}

Mdouble BaseCG::getWidthTime() const
{
    logger(ERROR, "getWidthTime() should not be used, it is a dummy function for templating");
    return 0;
}

Mdouble BaseCG::getTimeStep() const
{
    logger(ERROR, "getTimeStep() should not be used, it is a dummy function for templating");
    return 0;
}


