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
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearViscoelasticSpecies.h>
#include "CG/CGHandler.h"
#include "DPMBase.h"
#include "Walls/InfiniteWall.h"

/*!
 * \param[in] ch The CGHandler that has to be copied.
 */
CGHandler::CGHandler(const CGHandler& ch)
        : BaseHandler(ch)
{
    setDPMBase(ch.getDPMBase());
    //copyContentsFromOtherHandler(ch);
}

/*!
 * \param[in] rhs The WallHandler on the right hand side of the assignment.
 * \details This is not a copy assignment operator! It only copies the pointer to the
 *          DPMBase and the BaseWall in objects_, it sets the other data members
 *          to 0 or nullptr.
 */
CGHandler& CGHandler::operator=(const CGHandler& rhs)
{
    if (this != &rhs)
    {
        clear();
        setDPMBase(rhs.getDPMBase());
        copyContentsFromOtherHandler(rhs);
    }
    return *this;
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "CGHandler::operator =(const CGHandler&) finished" << std::endl;
#endif
}

/*!
 * \param[in] cg A pointer to the CG object that has to be added to the handler.
 */
void CGHandler::addObject(BaseCG* cg)
{
    //Puts the cg object in the list
    BaseHandler<BaseCG>::addObject(cg);
    //set the CGHandler pointer
    cg->setHandler(this);
}

/*!
 * \return a string containing the name.
 */
std::string CGHandler::getName() const
{
    return "CGHandler";
}

void CGHandler::readAndAddObject(std::istream& is UNUSED)
{

}

void CGHandler::write(std::ostream& os UNUSED) const
{

}

void CGHandler::initialise()
{
    //initialise all CG objects in the handler
    for (BaseCG* it : *this)
        it->initialise();
};

void CGHandler::evaluate()
{
    //evaluate all CG objects in the handler
    for (BaseCG* it : *this)
    {
        //if we are below timeMax and the next time step should be written
        if (it->getTimeMin() <= getDPMBase()->getTime() && it->getTimeMax() > getDPMBase()->getTime()
            && it->statFile.saveCurrentTimeStep(getDPMBase()->getNumberOfTimeSteps()))
        {
            //logger(INFO,"evaluate %, nt=%",it->statFile.getName(), getDPMBase()->getNumberOfTimeSteps());
            it->statFile.setLastSavedTimeStep(getDPMBase()->getNumberOfTimeSteps());
            it->evaluate();
        }
    }

    // Remember information from this read
    handlePreviousInformation();
};

void CGHandler::finish()
{
    //enforce that data gets written
    for (BaseCG* it : *this)
    {
        it->statFile.setLastSavedTimeStep(NEVER);
    }
    //evaluate all CG objects in the handler
    //evaluate();
    //finish all CG objects in the handler
    for (BaseCG* it : *this)
        it->finish();
}

void CGHandler::restart(std::string name)
{
    DPMBase* dpm = getDPMBase();
    
    // determines if the variable name contains the name of a restart file, or a problem name
    // (in which case the restart file is assumed to be name.restart)
    std::string::size_type endName = name.find(".restart");
    if (endName == std::string::npos)
    {
        dpm->setName(name);
    }
    else
    {
        dpm->setName(name.substr(0, endName));
        dpm->restartFile.setName(name);
    }
    
    // read restart file
    cgLogger(INFO, "Reading %", dpm->restartFile.getFullName());

    if (!dpm->readRestartFile(DPMBase::ReadOptions::ReadNoInteractions))
    {
        if (dpm->speciesHandler.getNumberOfObjects() == 0)
        {
            dpm->speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        }
        std::ifstream data(dpm->dataFile.getName());
        if (data.is_open())
        {
            Mdouble N = 0, t = 0;
            Vec3D min, max;
            data >> N >> t >> min >> max;
            data.close();
            dpm->setDomain(min, max);
            logger(INFO, "Reading domain size from %: min %, max % ", dpm->dataFile.getName(), min, max);
        }
        else
        {
            logger(ERROR, "Data file could not be opened");
        }
        //adding 10 walls by default
//        while (dpm->wallHandler.getNumberOfObjects()<10) {
//            dpm->wallHandler.copyAndAddObject(InfiniteWall(dpm->speciesHandler.getLastObject()));
//        }
        //what to do if restart file could not be loaded
        cgLogger(WARN, "No density is specified, thus density is assumed to be 1. Provide % file if you want to specify the density.", dpm->restartFile.getName());
    }
    else
    {
        cgLogger(INFO, "Successfully restarted from %, t=%, Np=%, Nc=%", dpm->restartFile.getFullName(), dpm->getTime(),
                 dpm->particleHandler.getSize(), dpm->interactionHandler.getNumberOfObjects());
        dpm->restartFile.decreaseCounter();
    }
}

void CGHandler::restartAndEvaluateRestartFiles(const std::string& name)
{
    restart(name);
    evaluateRestartFiles();
}

void CGHandler::restartAndEvaluateDataFiles(const std::string& name, bool evaluateFStatFiles)
{
    restart(name);
    evaluateDataFiles(evaluateFStatFiles);
}

void CGHandler::computeContactPoints()
{
    //recompute contact point (necessary for old restart files)
    for (BaseInteraction* const c : getDPMBase()->interactionHandler)
    {
        if (c->getContactPoint().isNaN())
        {
            const Vec3D& p = c->getP()->getPosition();
            const Vec3D& i = c->getI()->getPosition();
            const BaseParticle* const PParticle = dynamic_cast<BaseParticle*>(c->getP());
            const BaseParticle* const IParticle = dynamic_cast<BaseParticle*>(c->getI());
            if (IParticle != nullptr)
            {
                const Vec3D branchVector = p - i;
                const Mdouble distance = branchVector.getLength();
                c->setNormal(branchVector / distance);
                c->setOverlap(PParticle->getRadius() + IParticle->getRadius() - distance);
                c->setDistance(distance);
                c->setContactPoint(p - (PParticle->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
            }
            else
            {
                const InfiniteWall* const IWall = dynamic_cast<InfiniteWall*>(c->getI());
                if (IWall != nullptr)
                {
                    const Mdouble dist = Vec3D::dot(i - p, IWall->getNormal());
                    const Mdouble r = dynamic_cast<BaseParticle*>(c->getP())->getRadius();
                    const Vec3D branch = (0.5 * (dist + r)) * IWall->getNormal();
                    c->setNormal(-IWall->getNormal());
                    c->setContactPoint(p + branch);
                }
                else
                {
                    const Vec3D IP = p - i;
                    const Mdouble distance = Vec3D::getLength(IP);
                    c->setNormal(IP / distance);
                    c->setContactPoint(i);
                }
            }
            //logger(INFO,"c%",c->getContactPoint());
            static bool firstTime = true;
            if (firstTime)
            {
                cgLogger(WARN,
                         "recomputing contact point, as contact point information is not available\n"
                         "Contact point is placed assuming spherical particles\n"
                         "Complex walls, non-spherical particles, periodic walls can create errors");
                firstTime = false;
            }
        }
    }
}

bool CGHandler::evaluateRestartFiles()
{

#ifdef MERCURYDPM_USE_MPI
    //Make sure that the number of processors is equal to the number of processors used for the run
    MPIContainer& communicator = MPIContainer::Instance();
    std::vector<unsigned> numberOfDomains = this->getDPMBase()->getNumberOfDomains();
    int numberOfRequiredProcessors = numberOfDomains[0]*numberOfDomains[1]*numberOfDomains[2];
    if(!(numberOfRequiredProcessors == communicator.getNumberOfProcessors()))
    {
        if (communicator.getProcessorID() == 0)
        {
            logger(ERROR,"Please re-run the program with % cores",numberOfRequiredProcessors);
        }
        else
        {
            std::exit(-1);
        }
    }
#endif
    // reset counters so reading begins with the first data/fstat file
    DPMBase* const dpm = getDPMBase();
    
    //define and check time limits
    Mdouble timeMin = getTimeMin();
    Mdouble timeMax = getTimeMax();
    if (getDPMBase()->getTime() > timeMax)
    {
        logger(ERROR, "initial restart file (t=%) is beyond the maximum cg time (tMax=%)", dpm->getTime(), timeMax);
    }
    else
        while (getDPMBase()->getTime() < timeMin)
        {
            // Remember information from this read
            handlePreviousInformation();

            cgLogger(INFO, "Skipped %, t = %, because time is below tMin = %", dpm->restartFile.getFullName(),
                     dpm->getTime(), timeMin);
            //the particle and wall handler is cleared here, because  BaseSpecies doesn't delete particles belonging to it
            dpm->particleHandler.clear();
            dpm->wallHandler.clear();
            dpm->readRestartFile();
        }
    
    //call initialise to set up mesh, stat files, etc
    initialise();
    
    // evaluate restart files, starting with the one already read (since interactions where not read)
    while (dpm->readRestartFile() && getDPMBase()->getTime() < timeMax)
    {
        cgLogger(INFO, "Read %, t=%, Np=%, Nc=%", dpm->restartFile.getFullName(), dpm->getTime(),
                 dpm->particleHandler.getSize(), dpm->interactionHandler.getNumberOfObjects());
        
        //recompute contact point (necessary for old restart files)
        computeContactPoints();
        
        evaluate();
        
        //the particle and wall handler is cleared here, because  BaseSpecies doesn't delete particles belonging to it
        dpm->particleHandler.clear();
        dpm->wallHandler.clear();

        //continue if the next restart file can be read and the max time has not been reached
        if (dpm->restartFile.getFileType()==FileType::ONE_FILE) {
            logger(INFO,"Note: the evaluation is based on a single restart file. Please switch output to FileType::MULTIPLE_FILES if you want to coarse grain over multiple time steps");
            break;
        }
    }
    cgLogger(INFO, "Finished reading from %", dpm->restartFile.getFullName());
    
    finish();
    dpm->dataFile.close();
    return true;
}

bool CGHandler::evaluateDataFiles(bool evaluateFStatFiles)
{

#ifdef MERCURYDPM_USE_MPI
    //Make sure that the number of processors is equal to the number of processors used for the run
    MPIContainer& communicator = MPIContainer::Instance();
    std::vector<unsigned> numberOfDomains = this->getDPMBase()->getNumberOfDomains();
    int numberOfRequiredProcessors = numberOfDomains[0]*numberOfDomains[1]*numberOfDomains[2];
    if(!(numberOfRequiredProcessors == communicator.getNumberOfProcessors()))
    {
        if (communicator.getProcessorID() == 0)
        {
            logger(ERROR,"Please re-run the program with % cores",numberOfRequiredProcessors);
        }
        else
        {
            std::exit(-1);
        }
    }
#endif
    // reset counters so reading begins with the first data/fstat file
    DPMBase* const dpm = getDPMBase();
    //dpm->interactionHandler.clear();
    dpm->dataFile.setCounter(initialFileCounter);
    dpm->fStatFile.setCounter(initialFileCounter);
    
    initialise();
    
    // return false if no data file can be read
    ///\todo use ignore if time is out of bounds
    if (!dpm->readNextDataFile()) return false;
    
    //define and check time limits
    const Mdouble timeMin = getTimeMin();
    const Mdouble timeMax = getTimeMax();
    if (dpm->getTime() > timeMax)
    {
        logger(ERROR, "Time stamp of initial data file (%) is beyond the maximum cg time (tMax=%)", dpm->getTime(),
               timeMax);
    }
    else
        while (dpm->getTime() < timeMin)
        {
            // Remember information from this read
            handlePreviousInformation();

            if (evaluateFStatFiles) dpm->readNextFStatFile();
            cgLogger(INFO, "Skipped %, t = %, because time is below tMin = %", dpm->dataFile.getFullName(),
                     dpm->getTime(), timeMin);
            dpm->readNextDataFile();
        }
    
    // read data file
    do
    {
        if (evaluateFStatFiles) dpm->readNextFStatFile();
        
        cgLogger(INFO, "Read %, t=%, Np=%, Nc=%", dpm->dataFile.getFullName(), dpm->getTime(),
                 dpm->particleHandler.getSize(), dpm->interactionHandler.getNumberOfObjects());
        
        evaluate();
        
    } while (dpm->readNextDataFile() && getDPMBase()->getTime() <= timeMax);
    cgLogger(INFO, "Finished reading from %", dpm->dataFile.getFullName());
    
    finish();
    dpm->dataFile.close();
    return true;
}

Mdouble CGHandler::getTimeMin()
{
    Mdouble time = constants::inf;
    for (BaseCG* it : *this)
    {
        time = std::min(time, it->getTimeMin());
    }
    return time;
}

Mdouble CGHandler::getTimeMax()
{
    Mdouble time = -constants::inf;
    for (BaseCG* it : *this)
    {
        time = std::max(time, it->getTimeMax());
    }
    return time;
}

void CGHandler::setPreviousEvaluationTimeToCurrentTime() 
{
    previousEvaluationTime_ = getDPMBase()->getTime();
}

void CGHandler::setPreviousPositionToCurrentPosition() 
{
    for (BaseParticle* p : getDPMBase()->particleHandler)
    {
        p->setPreviousPosition(p->getPosition());
    }
}

void CGHandler::handlePreviousInformation()
{
    setPreviousEvaluationTimeToCurrentTime();
    setPreviousPositionToCurrentPosition();
}