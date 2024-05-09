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
#include "DPMBase.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <string>
#include <cstdio>
///\todo strcmp relies on this, should be changed to more modern version
#include <cstring>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include "Interactions/Interaction.h"
#include "Species/FrictionForceSpecies/SlidingFrictionSpecies.h"
#include "CMakeDefinitions.h"
#include "DPMBaseXBalls.icc" //This is part of this class and just separates out the stuff to do with xballs.
#include "Logger.h"
#include "Particles/SphericalParticle.h"
#include "Particles/NonSphericalParticle.h"
#include "Walls/BaseWall.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "VTKWriter/SuperQuadricParticleVtkWriter.h"
#include "VTKWriter/SphericalParticleVtkWriter.h"

//This is MPI related
#include "MpiContainer.h"
#include "MpiDataClass.h"
#include "Domain.h"
// this header for signal handler
#include <csignal>
// this header for memset function
#include <cstring>

#ifdef MERCURYDPM_USE_MPI
#include <mpi.h>
#endif

//This is OMP related
#ifdef MERCURYDPM_USE_OMP
#include <omp.h>
#endif

/*!
 * \details Warns the user of a fatal error and exits the program with a non-zero return
 * value to let the compiler know an error has occurred.
 * \param[in] module
 * \param[in] message
 * \todo Why is this here, and not in the logger?
 */
/**
* \deprecated
*/
MERCURYDPM_DEPRECATED
[[noreturn]] void logWriteAndDie(const std::string& module, std::string message)
{
    std::cerr << "A fatal   error has occured"
              << "\n Module  :" << module
              << "\n Message :" << message << std::endl;

    std::exit(-1);
}

/*!
* \details Overloads the insertion operator (<<) for use with Mercury objects. Utilises
* the write() function (see link for further information).
*
* \param[in] os - The output stream to which we want to 'insert' information relating to
* 					Mercury objects
* \param[in] md - An object (passed by reference) of the DPMBase class.
*/
std::ostream& operator<<(std::ostream& os, DPMBase& md)
{
md.write(os);
return os;
}

/*!
* \details A copy constructor which takes a "DPMBase"-type object and creates
* a "copy" - i.e. creates a new instance of a class possessing the same properties as the original. \n
* The argument "other" is the "original", i.e. the instance to be copied from. \n \n
* The first set of actions performed, which follow the general form: \n
*  <tt>[variable] = other.[variable]</tt>) \n
* simply copy the relevant variables (i.e. information such as particle details, system details, simulation details...)
* from the original ("other").
* \n \n
* The various <B>handlers</B> belonging to the original instance, however, are not directly copied,
* as this may cause problems (i.e. the handlers would still "point" to the original object,
* not the copy).
* \n
* Rather, new handlers are created (e.g. <tt>boundaryHandler.setDPMBase(this);</tt>), and the <B>contents</B>
* of the handlers is then passed over (e.g. <tt>boundaryHandler = other.boundaryHandler;</tt>).
* For each handler class, the assignment operator = has been overrided to copy the contents, not
* just the pointer.
* \param[in] other
*/
DPMBase::DPMBase(const DPMBase& other) : wallVTKWriter_(other.wallVTKWriter_),
                                         interactionVTKWriter_(other.interactionVTKWriter_),
                                         boundaryVTKWriter_(other.boundaryVTKWriter_),
                                         wallDetailsVTKWriter_(other.wallDetailsVTKWriter_)
{
    setName(other.getName());
    runNumber_ = other.runNumber_;
    systemDimensions_ = other.systemDimensions_;
    particleDimensions_ = other.particleDimensions_;
    gravity_ = other.gravity_;
    backgroundDrag_=other.backgroundDrag_;
/*    xMin_ = other.xMin_;
    xMax_ = other.xMax_;
    yMin_ = other.yMin_;
    yMax_ = other.yMax_;
    zMin_ = other.zMin_;
    zMax_ = other.zMax_;*/
    min_ = other.min_;
    max_ = other.max_;
    numberOfDomains_ = other.numberOfDomains_;
    time_ = other.time_;
    timeStep_ = other.timeStep_;
    numberOfTimeSteps_ = other.numberOfTimeSteps_;
    timeMax_ = other.timeMax_;
    restartVersion_ = other.restartVersion_; //to read new and old restart data
    restarted_ = other.restarted_; //to see if it was restarted or not
    append_ = other.append_;
    rotation_ = other.rotation_;
    xBallsColourMode_ = other.xBallsColourMode_; // sets the xballs argument cmode (see xballs.txt)
    xBallsVectorScale_ = other.xBallsVectorScale_; // sets the xballs argument vscale (see xballs.txt)
    xBallsScale_ = other.xBallsScale_; // sets the xballs argument scale (see xballs.txt)
    xBallsAdditionalArguments_ = other.xBallsAdditionalArguments_; // std::string where additional xballs argument can be specified (see xballs.txt)
    writeParticlesVTK_ = other.writeParticlesVTK_;
    readSpeciesFromDataFile_ = other.readSpeciesFromDataFile_;

//effectively saying "if there exists a CONTACT_LIST_HGRID, copy it, if not, ignore.
#ifdef CONTACT_LIST_HGRID
possibleContactList=other.possibleContactList;
#endif
    random = other.random;

    boundaryHandler.setDPMBase(this);
    particleHandler.setDPMBase(this);
    interactionHandler.setDPMBase(this);
    speciesHandler.setDPMBase(this);
    wallHandler.setDPMBase(this);
    domainHandler.setDPMBase(this);
    periodicBoundaryHandler.setDPMBase(this);
    //Initialise the handlers
    domainHandler.initialise();
    periodicBoundaryHandler.initialise();

    //setting contents equal to the other handlers!
    speciesHandler = other.speciesHandler;
    particleHandler = other.particleHandler;
    cgHandler = other.cgHandler;
    //cgHandler = other.cgHandler.copy(); //\todo
    //cgHandler.setDPMBase(this);
    wallHandler = other.wallHandler;
    boundaryHandler = other.boundaryHandler;
    interactionHandler = other.interactionHandler;
    vtkWriter_ = other.vtkWriter_;
    writeSuperquadricParticlesVTK_ = other.writeSuperquadricParticlesVTK_;
    writeParticlesVTK_ = other.writeParticlesVTK_;
    numberOfOMPThreads_ = other.numberOfOMPThreads_;
}

/*!
 * Constructor for the DPMBase class. Initialises a set of default parameters allowing
 * a simulation to be created 'off the shelf'. For full details of the parameters
 * initialised and their assigned values, see constructor()
 */
DPMBase::DPMBase() : wallVTKWriter_(wallHandler), interactionVTKWriter_(interactionHandler), boundaryVTKWriter_(boundaryHandler),
                     wallDetailsVTKWriter_(wallHandler)
{
constructor();
}

/*!
* \details Provides all the necessary default values for the DPMBase() constructor. When called, will initialise a two-dimensional simulation
* (<tt>setSystemDimensions(2), setParticleDimensions(2)</tt>)
* with "normal" vertical gravity
* (<tt>gravity_ = Vec3D(0.0, -9.8, 0.0);</tt>)
* as well as defining an arbitrary length (1s) and XBalls viewing domain (0.01 x 0.01) and other relevant viewing parameters (e.g. colourscheme, scale...).
* The first block of text creates the necessary handlers and sets their content according to the current ("this") instance of the DPMBase superclass.
*/
void DPMBase::constructor()
{
    //constructor();
    dataFile.getFstream().precision(10);
    fStatFile.getFstream().precision(10);
    eneFile.getFstream().precision(10);
    restartFile.getFstream().precision(
            std::numeric_limits<double>::digits10); //highly accurate, so the overlap is accurate
    statFile.getFstream().precision(10);
    statFile.getFstream().setf(std::ios::left);
    interactionFile.getFstream().precision(10);
    name_ = ""; // needs to be user-specified, otherwise checkSettings throws error
    //by default, the fileType of all files is ONE_FILE. However, by default we don't want an interaction file since it
    // is very large.
    interactionFile.setFileType(FileType::NO_FILE);

    runNumber_ = 0;

    //Decomposition direction for MPI
    numberOfDomains_ = {1, 1, 1};

    //Check if MPI is already initialised
    initialiseMPI();

    //This sets the maximum number of particles
    boundaryHandler.setDPMBase(this);
    periodicBoundaryHandler.setDPMBase(this);
    speciesHandler.setDPMBase(this);
    particleHandler.setDPMBase(this);
    cgHandler.setDPMBase(this);
    interactionHandler.setDPMBase(this);
    wallHandler.setDPMBase(this);
    interactionHandler.setDPMBase(this);
    domainHandler.setDPMBase(this);
    domainHandler.initialise();
    periodicBoundaryHandler.setDPMBase(this);
    periodicBoundaryHandler.initialise();

    //set defaults for DPMBase parameters
    setSystemDimensions(3);
    setParticleDimensions(3);
    setRestarted(false);
    setGravity(Vec3D(0, 0, 0));
    setBackgroundDrag(0.0);

    //This is the parameter of the numerical part
    setTime(0);
    numberOfTimeSteps_ = 0;
    setTimeMax(0);
    timeStep_ = 0; // needs to be user-specified, otherwise checkSettings throws error
    setSaveCount(20);

    //This sets the default xballs domain
    min_ = Vec3D(0, 0, 0);
    max_ = Vec3D(0, 0, 0); // needs to be user-specified, otherwise checkSettings throws error

    //sets the default write particles data in VTK format flag to false
    writeParticlesVTK_ = false;
    writeSuperquadricParticlesVTK_ = false;
    vtkWriter_ = nullptr;

    setName(""); // needs to be user-specified, otherwise checkSettings throws error

    //Default mode is energy with no scale of the vectors
    xBallsColourMode_ = 0;
    xBallsVectorScale_ = -1;
    xBallsScale_ = -1;
    xBallsAdditionalArguments_ = "";
    setAppend(false);

    //The default random seed is 0
    random.setRandomSeed(0);

    logger(DEBUG, "DPMBase problem constructor finished");

    readSpeciesFromDataFile_ = false;

    numberOfOMPThreads_ = 1;

    //Set number of elements to write to the screen if a user wants to output write information to the terminal
    nToWrite_ = 4;
}

/*!
 * \details A simple destructor for "DPMBase"-type objects, used to free-up memory when an object
 * is no longer necessary.
 */
DPMBase::~DPMBase()
{
delete vtkWriter_;
}

/*!
* \returns File& (A reference of object type File i.e. File& dataFile)
*/
File& DPMBase::getDataFile()
{
return dataFile;
}

/*!
* \returns File& (A reference of object type File i.e. File& eneFile)
*/
File& DPMBase::getEneFile()
{
return eneFile;
}

/*!
* \returns File& (A reference of object type File i.e. File& fStatFile)
*/
File& DPMBase::getFStatFile()
{
return fStatFile;
}

/*!
* \returns File& (A reference of object type File i.e. File& restartFile)
*/
File& DPMBase::getRestartFile()
{
return restartFile;
}

/*!
* \returns File& (A reference of object type File i.e. File& statFile)
*/
File& DPMBase::getStatFile()
{
return statFile;
}

/*!
* \return A reference of object type File i.e. File* interactionFile_
*/
File& DPMBase::getInteractionFile()
{
return interactionFile;
}


/*!
* \returns const File& (A const reference of object type File i.e. const File& dataFile)
*/
const File& DPMBase::getDataFile() const
{
return dataFile;
}

/*!
* \returns const File& (A const reference of object type File i.e. const File& eneFile)
*/
const File& DPMBase::getEneFile() const
{
return eneFile;
}

/*!
* \returns const File& (A const reference of object type File i.e. const File& fStatFile)
*/
const File& DPMBase::getFStatFile() const
{
return fStatFile;
}

/*!
* \returns const File& (A const reference of object type File i.e. const File& restartFile)
*/
const File& DPMBase::getRestartFile() const
{
return restartFile;
}

/*!
* \returns const File& (A const reference of object type File i.e. const File& statFile)
*/
const File& DPMBase::getStatFile() const
{
return statFile;
}
/*!
* \returns const File& (A const reference of object type std::string i.e. const std::string& name_)
*/
/// \bug The InteractionFile does not work across multifiles.
const File& DPMBase::getInteractionFile() const
{
return interactionFile;
}

const std::string& DPMBase::getName() const
{
return name_;
}

/*!
* \details sets the number of time steps skipped between each save for ALL data files, except for the interaction file.
* Note, that the interaction file is independent of time steps, and just writes when an interaction starts or ends.
*/
void DPMBase::setSaveCount(unsigned int saveCount)
{
dataFile.setSaveCount(saveCount);
fStatFile.setSaveCount(saveCount);
restartFile.setSaveCount(saveCount);
statFile.setSaveCount(saveCount);
eneFile.setSaveCount(saveCount);
for (auto cg : cgHandler)
cg->statFile.setSaveCount(saveCount);
}

/*!
* \param[in] name
*/
void DPMBase::setName(const std::string& name)
{
if (NUMBER_OF_PROCESSORS > 1)
{
name_ = name; // was before this->name_ = name
dataFile.setName(name_ + ".data" + std::to_string(PROCESSOR_ID));
fStatFile.setName(name_ + ".fstat" + std::to_string(PROCESSOR_ID));
restartFile.setName(name_ + ".restart" + std::to_string(PROCESSOR_ID));
statFile.setName(name_ + ".stat" + std::to_string(PROCESSOR_ID));
eneFile.setName(name_ + ".ene" + std::to_string(PROCESSOR_ID));
getInteractionFile().setName(name_ + ".interaction" + std::to_string(PROCESSOR_ID));
}
else
{
name_ = name; // was before this->name_ = name
dataFile.setName(name_ + ".data");
fStatFile.setName(name_ + ".fstat");
restartFile.setName(name_ + ".restart");
statFile.setName(name_ + ".stat");
eneFile.setName(name_ + ".ene");
interactionFile.setName(name_ + ".interaction");
}
}

/*!
* \param[in] name
*/
void DPMBase::setName(const char* name)
{
setName(std::string(name));
}

/*!
* \details Calls the setFileType() function from the File.h, which basically sets the File::fileType_. Note, this does
* not affect the interactionFile.
* \param[in] fileType (an object of enum class FileType)
*/
void DPMBase::setFileType(FileType fileType)
{
dataFile.setFileType(fileType);
fStatFile.setFileType(fileType);
restartFile.setFileType(fileType);
statFile.setFileType(fileType);
eneFile.setFileType(fileType);
}

/*!
* \details This implicitly calls the setCounter() function defined in File.h
*/
void DPMBase::resetFileCounter()
{
    dataFile.setCounter(0);
    fStatFile.setCounter(0);
    restartFile.setCounter(0);
    statFile.setCounter(0);
    eneFile.setCounter(0);
    interactionFile.setCounter(0);
    setLastSavedTimeStep(NEVER);
    if (vtkWriter_) vtkWriter_->setFileCounter(0);
    boundaryVTKWriter_.setFileCounter(0);
    interactionVTKWriter_.setFileCounter(0);
    wallVTKWriter_.setFileCounter(0);
    wallDetailsVTKWriter_.setFileCounter(0);
}

/*!
* \param[in] openmode
*/
void DPMBase::setOpenMode(std::fstream::openmode openMode)
{
dataFile.setOpenMode(openMode);
fStatFile.setOpenMode(openMode);
restartFile.setOpenMode(openMode);
statFile.setOpenMode(openMode);
eneFile.setOpenMode(openMode);
interactionFile.setOpenMode(openMode);
}

/*!
*
*/
void DPMBase::closeFiles()
{
dataFile.close();
fStatFile.close();
restartFile.close();
statFile.close();
eneFile.close();
interactionFile.close();
}

/*!
* \details Sets the time step when the files will next be saved, except for the interaction file.
* Note, that the interaction file is independent of time steps, and just writes when an interaction starts or ends.
* \param[in] nextSavedTimeStep
*
*/
void DPMBase::setLastSavedTimeStep(unsigned int nextSavedTimeStep)
{
dataFile.setLastSavedTimeStep(nextSavedTimeStep);
fStatFile.setLastSavedTimeStep(nextSavedTimeStep);
restartFile.setLastSavedTimeStep(nextSavedTimeStep);
//statFile.setLastSavedTimeStep(nextSavedTimeStep); //this one is not force-written
eneFile.setLastSavedTimeStep(nextSavedTimeStep);
}

/*!
* \details Using the three functions named above, the autoNumber() function acts to:
*
* 1) Use the \ref readRunNumberFromFile() function toead the current run number from the file COUNTER_DONOTDEL
* created by any script which utilises auto-numbering.
*
* 2) Set the  \ref runNumber_ counter to the value obtained from the above using the \ref setRunNumber() function.
*
* 3) Increment the value stored in the COUNTER_DONOTDEL file by one once the current value has been read
* using the \ref incrementRunNumberInFile() function.
*/
void DPMBase::autoNumber()
{
setRunNumber(readRunNumberFromFile());

if (!getRestarted())
{
incrementRunNumberInFile();
}
}

/*!
* \details Reads in the current counter in from the COUNTER_DONOTDEL file stored on the disk.
* If a COUNTER_DONOTDEL file does not already exist, creates one and initialises it with a value "1"
*/
int DPMBase::readRunNumberFromFile()
{
int counter;

FILE* counter_file;
//checking if there exists already a file named "COUNTER_DONOTDEL" which can be opened for
//input and output (returns "true" if no such file exists).
if ((counter_file = fopen("COUNTER_DONOTDEL", "r+")) == nullptr)
{
//if a file does not already exist, checks whether a new file can be created
//retutns "true" if a new file CANNOT be created
if ((counter_file = fopen("COUNTER_DONOTDEL", "w")) == nullptr)
{
    //If true, outputs an error message and ends the program
    fprintf(stderr, "\n\n\tERROR :: Counter File NOT found, please re-create\n\n");
    fclose(counter_file);
    exit(-1);
}
    //alternatively, if a new file CAN be created...
else
{
    //starts the new counter file, writing to it the value "1"
    fprintf(counter_file, "1");
    fprintf(stderr, "Counter File created\n");
    fclose(counter_file);
    return 1;
}
}
//alternatively, if a counter file DOES already exist...
else
{
//...checks if there exists only 1 value in the file (as would be expected from a COUNTER_DONOTDEL file...
if (fscanf(counter_file, "%d", &counter) != 1)
{
    //...and if not, returns an error.
    fprintf(stderr, "\n\n\tERROR :: Counter File found, but something went wrong with reading it\n\n");
    fclose(counter_file);
    exit(-1);
}
else
{
    //...otherwise, if all has been successful, returns the current value of the file!
    fclose(counter_file);
    return counter;
}
}

}

/*!
* A simple "set function" which allows the user to simply overwrite the current run number to any valid new value.
* \param[in] runNumber - the value to which we want to (re)set the internally stored run number parameter, \ref runNumber_
*/
void DPMBase::setRunNumber(int runNumber)
{
runNumber_ = runNumber;
}

/*!
* A simple "get function" which allows the user to retrieve the current value corresponding to the run number counter, \ref runNumber_
*
\returns \ref runNumber_ - the stored value of the current run number, i.e. the number of files corresponding to a given Mercury script that have been produced in
a given directory.
*/
int DPMBase::getRunNumber() const
{
return runNumber_;
}

/*!
* \details In order to increment the counter stored in COUNTER_DONOTDEL, we initialise two fstream objects counter_file, counter_file2 and
* an integer type temp_counter. First we open the file COUNTER_DONOTDEL, check if everything went fine with the opening. If yes, we extract the
* runNumber (counter) into the temp_counter. Increment the temp_counter and then write it into COUNTER_DONOTDEL. This is how we increment the
* counter in the file.
*/
void DPMBase::incrementRunNumberInFile()
{
//opening two filestreams - counter_file and counter_file2
std::fstream counter_file, counter_file2;
//declares an integer, temp_counter
int temp_counter;
//attempts to open the COUNTER_DONOTDEL text file
counter_file.open("COUNTER_DONOTDEL", std::ios::in);
//gives error message if file could not be successfully opened and ends the program
if (counter_file.fail())
{
fprintf(stderr, "\n\n\tERROR :: Counter File NOT found, please re-create\n\n");
counter_file.close();
exit(0);
}
// if opened successfully, reads in the counter corresponding to the current run number
//and stored it in the "temp_counter" variable
counter_file >> temp_counter;
counter_file.close();
//Increments the temp_counter
temp_counter++;
//opens an output stream to the COUNTER_DONOTDEL file
counter_file2.open("COUNTER_DONOTDEL", std::ios::out);
if (counter_file2.fail())
{
fprintf(stderr, "\n\n\tERROR :: Counter File NOT found, please re-create2\n\n");
counter_file2.close();
exit(0);
}
//writes the new valuer of the counter to COUNTER_DONOTDEL
counter_file2 << temp_counter;

counter_file2.close();
}

/*!
* \details Let's say sizeX = 5, counter stored in COUNTER_DONOTDEL = 1.
* Substituting these values into the algorithm below implies that studyNum = 0 or 1. Everytime the code is executed the
* counter gets incremented and the values of studyNum and i are updated, which is returned as std::vector<int>
* \param[in] sizeX The (integer) number of values to be tested in 1D parameter space.
* \returns std::vector<int> The current study numbers.
*/
std::vector<int> DPMBase::get1DParametersFromRunNumber(int sizeX) const
{
// Declare a vector of integers capable of storing 2 values
std::vector<int> temp(2);

// Declare and initialise for the current simulation run number
int counter = getRunNumber();

// Give studyNum value 0 if study is incomplete, otherwise value > 0
int studyNum = (counter-1)/sizeX;
counter = counter - sizeX*studyNum;

int i = ((counter - 1) % sizeX) + 1;
logger(INFO,"StudyNum: % \t Counter: % \t i: %", studyNum, counter, i);
temp[0] = studyNum;
temp[1] = i;

return temp;
}

/*!
* \details Let's say sizeX = 2 and sizeY = 5, counter stored in COUNTER_DONOTDEL =1. The studySize = 10.
* Substituting these values into the below algorithm implies that studyNum = 0 or 1, everytime the code is executed the counter gets incremented and hence determined
* the values of studyNum, i and j which is returned as a std::vector<int>
* \param[in] sizeX The (integer) number of values to be tested for one of the 2 parameters forming the 2D parameter space.
* \param[in] sizeY The (integer) number of values to be tested for the other of the 2 parameters forming the 2D parameter space.
* \returns std::vector<int>
*/
std::vector<int> DPMBase::get2DParametersFromRunNumber(int sizeX, int sizeY) const
{
//declares a vector of integers capable of storing 3 values,
std::vector<int> temp(3);
//declares and initialises an integer variable named "counter"
//with the current counter number, runNumber_
int counter = getRunNumber();
//calculates the total size of the study, i.e. the number of points
//in the 2D parameter space explored
int studySize = sizeX * sizeY;
//(counter - 1) / studySize gives a fraction comparing the number of runs conducted so far
//to the total size of the study, i.e. the total number of runs that need to be performed.
//since studyNum is an integer, will declare zero until an adequate number of runs has been performed,
//at which point it will equal 1
int studyNum = (counter - 1) / studySize;

counter = counter - studySize * studyNum;
int i = ((counter - 1) % sizeX) + 1;
int j = ((counter - i) / sizeX) + 1;
logger(INFO,"StudyNum: % \t Counter: % \t i: % \t j: %", studyNum, counter, i, j);

temp[0] = studyNum;
temp[1] = i;
temp[2] = j;

return (temp);
}

/*!
* \details Let's say sizeX = 2, sizeY = 5 and sizeZ = 3, counter stored in COUNTER_DONOTDEL =1. The studySize = 30.
* Substituting these values into the below algorithm implies that studyNum = 0 or 1, everytime the code is executed the counter gets incremented and hence determined
* the values of studyNum, i,j and k which is returned as a std::vector<int>
* \param[in] sizeX The (integer) number of values to be tested for one of the 3 parameters forming the 3D parameter space.
* \param[in] sizeY The (integer) number of values to be tested for one of the 3 parameters forming the 3D parameter space.
* \param[in] sizeZ The (integer) number of values to be tested for one of the 3 parameters forming the 3D parameter space.
* \returns std::vector<int>
*/
std::vector<int> DPMBase::get3DParametersFromRunNumber(int sizeX, int sizeY, int sizeZ) const
{
//declares a vector of integers capable of storing 4 values,
std::vector<int> temp(4);
//declares and initialises an integer variable named "counter"
//with the current counter number, runNumber_
int counter = getRunNumber();
//calculates the total size of the study, i.e. the number of points
//in the 3D parameter space explored
int studySize = sizeX * sizeY * sizeZ;
//(counter - 1) / studySize gives a fraction comparing the number of runs conducted so far
//to the total size of the study, i.e. the total number of runs that need to be performed.
//since studyNum is an integer, will declare zero until an adequate number of runs has been performed,
//at which point it will equal 1
int studyNum = (counter - 1) / studySize;

counter = counter - studySize * studyNum;
int i = ((counter-1) % sizeX) + 1;
int j = static_cast<int>(std::floor((counter-1)/sizeX)) % sizeY + 1;
int k = static_cast<int>(std::floor((counter-1)/(sizeX*sizeY))) % sizeZ + 1;
logger(INFO,"StudyNum: % \t Counter: % \t i: % \t j: % \t k: %", studyNum, counter, i, j, k);

temp[0] = studyNum;
temp[1] = i;
temp[2] = j;
temp[3] = k;

return (temp);
}

/*!
* \details Reads in the name of the command (code) to be launched.
* This name is then converted to a string stream and appended with " &" (such that command
* is run in the background), before being
* converted back to a C string and then fed to the system() command which will execute the named code
* from within the running Mercury program.
* \param[in] name The name of the code to be launched
* \param[in] quick
* \return int
*/
int DPMBase::launchNewRun(const char* name, bool quick UNUSED)
{
//defines an (empty) stringstream named "com"
std::stringstream com("");
//adds the name of the code to run (fed in as an argument)
//to the "com" string and appends the string with " &"
com << name << " &";
//converts the stringstream "com" to a standard string, and then
//converts this string to a C string
//the string is then fed to the "system" function, which will run the named command
return system(com.str().c_str());
}

/*!
* \details
* This method allows flags to be passed to Mercury from driver codes, such that variables can be
* altered without needing to alter the driver files - for example, if the user wishes to give
* specific commands for the manner in which the system will be displayed in xballs.
*
* After reading in the arguments provided, the normal 'solve()' routine is called. For full details
* see the documentation of the solve() function (linked).
* \param[in] argc
* \param[in] argv
*/
void DPMBase::solve(int argc, char* argv[])
{
readArguments(argc, argv);
solve();
}

/*!
* \return time_
*/
Mdouble DPMBase::getTime() const
{
return time_;
}

/*!
* \return time_
*/
Mdouble DPMBase::getNextTime() const
{
return time_ + timeStep_;
}

/*!
* \return numberOfTimeSteps_
*/
unsigned int DPMBase::getNumberOfTimeSteps() const
{
return numberOfTimeSteps_;
}

/*!
* \details This may be useful in codes where some initial set-up is required e.g. if a system of particles
* is first prepared and then exposed to excitation.
* In this situation, <TT>getNumberOfTimeSteps()</TT> may be used to reset the time to zero at the point at which excitation
* begins to be applied.
* \param[in] time
*/
void DPMBase::setTime(Mdouble time)
{
Mdouble diff = time_ - time;
time_ = time;
//this sets the interaction timestamp, so each interaction has the right time
for (auto i: interactionHandler)
{
i->setTimeStamp(i->getTimeStamp() - diff);
}
}


/*!
* \details sets nToWrite_. If a user wants to output e.g. particle information to the screen we limit the number of
* outputs to nToWrite elements to not get an overflow of information in the terminal.
* \param[in] nToWrite          Number of elements to write to the screen
*/
void DPMBase::setNToWrite(int nToWrite)
{
nToWrite_ = nToWrite;
}

/*!
* \details Gets nToWrite. If a user wants to output e.g. particle information to the screen we limit the number of
* outputs to nToWrite elements to not get an overflow of information in the terminal.
* \param[out] nToWrite_          Number of elements to write to the screen
*/

int DPMBase::getNToWrite() const
{
return nToWrite_;
}

/*!
* \details A sanity check is performed to ensure that the new maximum simulation duration is nonnegative.
* \param[in] newTMmax
*/
void DPMBase::setTimeMax(Mdouble newTMax)
{
if (newTMax >= 0)
{
timeMax_ = newTMax;
}
else
{
logger(ERROR, "Error in setTimeMax, new timeMax=% is not positive", newTMax);
}
}

/*!
* \return timeMax_
*/
Mdouble DPMBase::getTimeMax() const
{
return timeMax_;
}
/*!
* Uses the preprocessor directive <TT>ifdef</TT> to check if there exists a
* CONTACT_LIST_HGRID, before any code is compiled.
* If CONTACT_LIST_HGRID <B>does</B> exist, this function can be used to
* return the "possibleContactsList" - but not to alter it.
*/
#ifdef CONTACT_LIST_HGRID
PossibleContactList& DPMBase::getPossibleContactList()
{
return possibleContactList;
}
#endif

/*!
* \details
* The VTK file is used for visualisation in Paraview.
* \todo Move this (and the get) to WallHandler.
* \param[in] writeWallsVTK
*/
void DPMBase::setWallsWriteVTK(FileType writeWallsVTK)
{
    logger(WARN, "DPMBase.setWallsWriteVTK(FileType) is deprecated! Use wallHandler.setWriteVTK(FileType) instead.");
    wallHandler.setWriteVTK(writeWallsVTK);
}

/*!
* \details
* The VTK file is used for visualisation in Paraview.
* \todo Move this (and the get) to WallHandler.
* \param[in] writeWallsVTK
*/
void DPMBase::setWallsWriteVTK(bool writeVTK)
{
    FileType writeVTKFileType = writeVTK?FileType::MULTIPLE_FILES:FileType::NO_FILE;
    logger(WARN, "DPMBase.setWallsWriteVTK(bool) is deprecated! Use wallHandler.setWriteVTK(bool) instead.");
    wallHandler.setWriteVTK(writeVTKFileType);
}

void DPMBase::setInteractionsWriteVTK(bool writeVTK)
{
    logger(WARN, "DPMBase.setInteractionsWriteVTK(bool) is deprecated! Use interactionHandler.setWriteVTK(bool) instead.");
    interactionHandler.setWriteVTK(writeVTK);
}

/*!
* \details
* The VTK format is used for visualisation in Paraview.
* \todo Move this (and the get) to ParticleHandler.
* \param[in] writeParticlesVTK
*/
void DPMBase::setParticlesWriteVTK(bool writeParticlesVTK)
{
writeParticlesVTK_ = writeParticlesVTK;
if (writeParticlesVTK_)
{
writeSuperquadricParticlesVTK_ = false;
}
delete vtkWriter_;
vtkWriter_ = new SphericalParticleVtkWriter(particleHandler);
}

/*!
* \param[in] writeParticlesVTK
*/
void DPMBase::setSuperquadricParticlesWriteVTK(bool writeParticlesVTK)
{
writeSuperquadricParticlesVTK_ = writeParticlesVTK;
if (writeSuperquadricParticlesVTK_)
{
writeParticlesVTK_ = false;
}
delete vtkWriter_;
vtkWriter_ = new SuperQuadricParticleVtkWriter(particleHandler);
}

/*!
* \details
* The VTK file is used for visualisation in Paraview.
* \todo Move this (and the set) to WallHandler.
* \returns bool
*/
FileType DPMBase::getWallsWriteVTK() const
{
    logger(WARN, "DPMBase.getWallsWriteVTK() is deprecated! Use wallHandler.getWriteVTK() instead.");
    return wallHandler.getWriteVTK();
}

/*!
* \details
* The VTK format is used for visualisation in Paraview.
* \todo Move this (and the set) to ParticleHandler.
* \returns bool
*/
bool DPMBase::getParticlesWriteVTK() const
{
return writeParticlesVTK_;
}

/*!
* \returns bool
*/
bool DPMBase::getSuperquadricParticlesWriteVTK() const
{
return writeSuperquadricParticlesVTK_;
}

/*!
* \details An access function which allows the user to alter the value of the (private) x value of the "Vec3D"
* object "min_", thus setting the minimum x-value corresponding to the system, i.e. the lower limit of the domain in
* the x-direction.
* \n \n
* These bounds are used for display, plotting and statistics purposes. They do not affect the
* dynamics.
* \n \n
* The function also performs a sanity check to stop the user defining a lower bound that
* is higher than the corresponding upper bound (XMax), giving a (logged) warning if the user attempts to do so.
* \param[in] newXMin
*/
void DPMBase::setXMin(Mdouble newXMin)
{
if (newXMin <= getXMax())
{
min_.x() = newXMin;
}
else
{
logger(WARN, "Warning in setXMin(%): xMax=%", newXMin, getXMax());
}
}

/*!
* \details An access function which allows the user to alter the value of the (private) y value of the "Vec3D"
* object "min_", thus setting the minimum y-value corresponding to the system, i.e. the lower limit of the domain in
* the y-direction.
* \n \n
* These bounds are used for display, plotting and statistics purposes. They do not affect the
* dynamics.
* \n \n
* The function also performs a sanity check to stop the user defining a lower bound that
* is higher than the corresponding upper bound (YMax), giving a (logged) warning if the user attempts to do so.
* \param[in] newYMin
*/
void DPMBase::setYMin(Mdouble newYMin)
{
if (newYMin <= getYMax())
{
min_.y() = newYMin;
}
else
{
logger(WARN, "Warning in setYMin(%): yMax=%", newYMin, getYMax());
}
}

/*!
* An access function which allows the user to alter the value of the (private) z value of the "Vec3D"
* object "min_", thus setting the minimum z-value corresponding to the system, i.e. the lower limit of the domain in
* the z-direction.
* \n \n
* These bounds are used for display, plotting and statistics purposes. They do not affect the
* dynamics.
* \n \n
* The function also performs a sanity check to stop the user defining a lower bound that
* is higher than the corresponding upper bound (ZMax), giving a (logged) warning if the user attempts to do so.
* \param[in] newZMin
*/
void DPMBase::setZMin(Mdouble newZMin)
{

if (newZMin <= getZMax())
{
min_.z() = newZMin;
}
else
{
logger(WARN, "Warning in setZMin(%): zMax=%", newZMin, getZMax());
}

}

/*!
* \details
* This specifies one corner of the problem's cuboidal bounding box.
* These bounds are used for display, plotting and statistics purposes. They do not affect the
* dynamics.
* \n \n
* A sanity check is performed to verify that each of the maximum coordinates are greater than the
* corresponding minimum coordinates. It raises a (logged) warning if not.
* \param[in] newMax
*/
void DPMBase::setMax(const Vec3D& newMax)
{
if (min_.x() > newMax.x() ||
min_.y() > newMax.y() ||
min_.z() > newMax.z())
{
logger(WARN, "Warning in setMax: upper bound is smaller"
	     " than lower bound. (%,%,%) > (%,%,%)",
       min_.x(), min_.y(), min_.z(), newMax.x(), newMax.y(), newMax.z());
}
else
{
max_ = newMax;
}
}

void DPMBase::setDomain(const Vec3D& min, const Vec3D& max)
{

logger.assert_debug(min.X <= max.X, "lower x-bound (%) is larger than upper x-bound (%)", min.X, max.X);
logger.assert_debug(min.Y <= max.Y, "lower x-bound (%) is larger than upper x-bound (%)", min.Y, max.Y);
logger.assert_debug(min.Z <= max.Z, "lower x-bound (%) is larger than upper x-bound (%)", min.Z, max.Z);
min_ = min;
max_ = max;
}

/*!
* \details
* This specifies one corner of the problem's cuboidal bounding box.
* These bounds are used for display, plotting and statistics purposes. They do not affect the
* dynamics.
* \n \n
* A sanity check is performed to verify that each of the minimum coordinates are greater than the
* corresponding maximum coordinates. It raises a (logged) warning if not.
* \param[in] newMin
*/
void DPMBase::setMin(const Vec3D& newMin)
{
if (max_.x() < newMin.x() ||
max_.y() < newMin.y() ||
max_.z() < newMin.z())
{
logger(WARN, "Warning in setMin: lower bound is larger"
	     " than upper bound. (%,%,%) < (%,%,%)",
       max_.x(), max_.y(), max_.z(), newMin.x(), newMin.y(), newMin.z());
}
else
{
min_ = newMin;
}
}

/*!
* \details
* As in \ref setMin but the coordinates are passed as three Mdouble, not a Vec3D.
*/
void DPMBase::setMin(const Mdouble xMin, const Mdouble yMin, const Mdouble zMin)
{
setMin(Vec3D(xMin, yMin, zMin));
}

/*!
* \details
* As in \ref setMax but the coordinates are passed as three Mdouble, not a Vec3D.
*/
void DPMBase::setMax(const Mdouble xMax, const Mdouble yMax, const Mdouble zMax)
{
setMax(Vec3D(xMax, yMax, zMax));
}

/*!
* \details
* An access function which allows the user to alter the value of the (private) x value of the "Vec3D"
* object "max_", thus setting the maximum x-value corresponding to the system, i.e. the upper limit of the domain in
* the x-direction.
* \n \n
* These bounds are used for display, plotting and statistics purposes. They do not affect the
* dynamics.
* \n \n
* The function also performs a sanity check to stop the user defining an upper bound that
* is lower than the corresponding lower bound (XMin), giving a (logged) warning if the user attempts to do so.
* \param[in] newXMax
*/
void DPMBase::setXMax(Mdouble newXMax)
{

if (newXMax >= getXMin())
{
max_.x() = newXMax;
}
else
{
logger(WARN, "Warning in setXMax(%): xMax=%", newXMax, getXMin());
}

}

/*!
* An access function which allows the user to alter the value of the (private) y value of the "Vec3D"
* object "max_", thus setting the maximum y-value corresponding to the system, i.e. the upper limit of the domain in
* the y-direction.
* \n \n
* These bounds are used for display, plotting and statistics purposes. They do not affect the
* dynamics.
* \n \n
* The function also performs a sanity check to stop the user defining an upper bound that
* is lower than the corresponding lower bound (YMin), giving a (logged) warning if the user attempts to do so.
* \param[in] newYMax
*/
void DPMBase::setYMax(Mdouble newYMax)
{

if (newYMax >= getYMin())
{
max_.y() = newYMax;
}
else
{
logger(WARN, "Warning in setYMax(%): yMax=%", newYMax, getYMin());
}

}

/*!
* \details An access function which allows the user to alter the value of the (private) z value of the "Vec3D"
* object "max_", thus setting the maximum z-value corresponding to the system, i.e. the upper limit of the domain in
* the z-direction.
* \n \n
* These bounds are used for display, plotting and statistics purposes. They do not affect the
* dynamics.
* \n \n
* The function also performs a sanity check to stop the user defining an upper bound that
* is lower than the corresponding lower bound (ZMin), giving a (logged) warning if the user attempts to do so.
* \param[in] newZMax
*/
void DPMBase::setZMax(Mdouble newZMax)
{
if (newZMax >= getZMin())
{
max_.z() = newZMax;
}
else
{
logger(WARN, "Warning in setZMax(%): zMax=%", newZMax, getZMin());
}
}

/*!
* \details A sanity check is performed to ensure that the time step must be positive.
* \param[in] timeStep
* The (Mdouble) value of the desired new time step
*/
void DPMBase::setTimeStep(Mdouble timeStep)
{
if (timeStep > 0.0)
{
timeStep_ = timeStep;
}
else
{
logger(ERROR, "Error in setTimeStep: new timeStep % is not positive", timeStep);
}
}

/*!
* \return timeStep_
* The current (Mdouble) value of the simulation time step.
*/
Mdouble DPMBase::getTimeStep() const
{
return timeStep_;
}


/* Allows user to set the number of omp threads */
void DPMBase::setNumberOfOMPThreads(int numberOfOMPThreads)
{
logger.assert_always(numberOfOMPThreads > 0, "Number of OMP threads must be positive");
numberOfOMPThreads_ = numberOfOMPThreads;

#ifdef MERCURYDPM_USE_OMP
if(numberOfOMPThreads > omp_get_max_threads()) {
logger(INFO, "Number of omp threads set to the maximum number of threads allowed: %",
       omp_get_max_threads());
numberOfOMPThreads_ = numberOfOMPThreads = omp_get_max_threads();
}
#pragma omp parallel num_threads(getNumberOfOMPThreads())
{
if (omp_get_thread_num()==0)
    logger(INFO, "Using % of % omp threads; testing thread",
	   omp_get_num_threads(), omp_get_max_threads(), Flusher::NO_FLUSH);
}
#pragma omp parallel num_threads(getNumberOfOMPThreads())
{
        logger(INFO, " %", std::to_string(omp_get_thread_num()), Flusher::NO_FLUSH);
    }
    logger(INFO,"");

#else
    logger(WARN, "You are setting the number of omp threads to %, but OMP is not turned on", getNumberOfOMPThreads());
#endif
}

/* Returns the number of omp threads */
int DPMBase::getNumberOfOMPThreads() const
{
    //logger.assert_debug(numberOfOMPThreads_,"You need to set the number of OMP threads");
    return numberOfOMPThreads_;
}

/*!
 * \details Allows the user to change the default "cmode" (colour mode) variable in MercuryDPM's built-in visualiser, "XBalls".
 * cmode takes an integer value between 1 and 27 (1 and 14 for 2D problems); to each integer is assigned
 * a different colour mode, which can be used to assign colours to individual particles based on a parameter such as the particles radius, velocity,
 * rotational energy etc. etc.
 * For further details, refer to the \ref xballs
 * \param[in] newCMode The numerical value corresponding to the colour mode (see \ref xballs) you want to use when visually reconstructing Mercury simulations.
 */
void DPMBase::setXBallsColourMode(int newCMode)
{
    xBallsColourMode_ = newCMode;
}

/*!
 * Returns the integer value corresponding to the colour scheme used by the XBalls visualisation software.
 * See also \ref setXBallsColourMode and the \ref xballs
 * \return int xBallsColourMode_ The integer value corresponding to the colour scheme used by the XBalls visualisation software.
 */
int DPMBase::getXBallsColourMode() const
{
    return xBallsColourMode_;
}

/*!
 * Allows the user to choose the default vector scaling, i.e. the length of the vectors representing particle velocities in the XBalls visualisation software.
 * Further details may be found in the \ref xballs
 * \param[in] newVScale The value of the desired vector length - a value of 100 sets the length to 1 particle radius, 1000 sets it to 10 particle radii etc.
 */
void DPMBase::setXBallsVectorScale(double newVScale)
{
    xBallsVectorScale_ = newVScale;
}

/*!
 * Returns the length of the vectors which represent particle velocities in XBalls visualisations (see also \ref setXBallsVectorScale and the \ref xballs).
 * \return double xBallsVectorScale_ The value of the vector length used in XBalls visualisations.
 * A value of 100 sets the length to 1 particle radius, 1000 sets it to 10 particle radii etc.
 */
double DPMBase::getXBallsVectorScale() const
{
    return xBallsVectorScale_;
}

/*!
 * \details Used to set, from the driver code itself, arguments to control the visualisation
 * of the simulation run using xballs.
 * <B>All arguments</B> can be passed as a <B>single string</B>, for example:
 *
 * setXBallsAdditionalArguments("-cmode 8 -solidf");
 *
 * will set the colour mode (cmode) to 8 (colour dependent on species)
 *  and draw particles with solid lines (solidf)
 *
 * \param[in] newXBArgs
 */
void DPMBase::setXBallsAdditionalArguments(std::string xBallsAdditionalArguments)
{
    xBallsAdditionalArguments_ = xBallsAdditionalArguments;
}

/*!
 * \return xBallsAdditionalArguments_
 */
std::string DPMBase::getXBallsAdditionalArguments() const
{
    return xBallsAdditionalArguments_;
}

/*!
 * \param[in] newScale The desired new scaling or "zoom". Values > 1 act to "zoom out", values < 1 act to "zoom in".
 */
void DPMBase::setXBallsScale(Mdouble newScale)
{
    xBallsScale_ = newScale;
}

/*!
 * \return double xBallsScale_ The scaling or "zoom" - corresponds tol the XBalls "-s" flag. Values > 1 mean a "zoomed out" view,
 * values < 1 give a "zoomed in" view.
 */
double DPMBase::getXBallsScale() const
{
    return xBallsScale_;
}

/*!
 * Allows the user to set the value of the gravitational acceleration (g) to which their simulated system is exposed.
 * The gravity is passed to the function as a (Vec3D) vector value, thus providing all three components of g
 * and hence allowing both its direction and strength to be altered.
 * \param[in] newGravity The desired new value of the gravitational acceleration as a Vec3D vector.
 */
void DPMBase::setGravity(Vec3D newGravity)
{
    gravity_ = newGravity;
}

/*!
 * \return Vec3D gravity_ The desired new value of the gravitational acceleration as a Vec3D vector.
 */
Vec3D DPMBase::getGravity() const
{
    return gravity_;
}

/*!
 * \details <B>Note:</B> In MercuryDPM, it is possible to simulate, for example, 3D particles in a 2D system by setting
 * "setSystemDimensions(newDim)" and "setParticleDimensions(newDim)" individually.
 * \n \n
 * The relevant sanity checks exist within the "setSystemDimensions(newDim)" and "setParticleDimensions(newDim)" functions.
 * \param[in] newDim The desired dimensionality of the system and the particles: 1 &rarr;  1D, 2 &rarr;  2D, 3 &rarr;  3D.
 */
void DPMBase::setDimension(unsigned int newDim)
{
    setSystemDimensions(newDim);
    setParticleDimensions(newDim);
}

/*!
 * \details Sets the dimensionality of the <B>system</B>.
 * (Note that <b>particles</p> may possess a different dimensionality.)
 * \n\n
 * A sanity check is performed to ensure that the dimensionality is 1, 2 or 3. If not, an error
 * message is logged and the simulation terminates.
 * \param[in] newDim The desired dimensionality of the system: 1 &rarr;  1D, 2 &rarr;  2D, 3 &rarr;  3D.
 */
void DPMBase::setSystemDimensions(unsigned int newDim)
{
    if (newDim >= 1 && newDim <= 3)
        systemDimensions_ = newDim;
    else
    {
        logger(ERROR, "Error in setSystemDimensions; newDim % is not 1, 2 or 3", newDim);
    }
}

/*!
 * \return systemDimensions_ The dimensionality of the <B>system</B>. (Note that <B>particles</B> may possess a different dimensionality.)
 */
unsigned int DPMBase::getSystemDimensions() const
{
    return systemDimensions_;
}

/*!
 * \details Sets the dimensionality of the <B>particles</B>.
 * This affects whether particles are treated as rods, discs or spheres. This is used e.g. for
 * calculating masses and moments of inertia.
 * \n \n
 * Note that the <i>system</i> may possess a different dimensionality.
 * \n \n
 * The particles can be chosen to be one- two- or three-dimensional by passing, respectively, arguments equal to 1, 2 or 3.
 * The code also performs a sanity check to ensure that the dimensionality cannot be higher than three or lower than one.
 * If an unphysical dimensionality is chosen, the particle dimension is set to be equal to the system dimension. Note,
 * that we cannot throw an error here, since in some old restart files the "particleDimensions_" are 0.
 * \param[in] particleDimensions An integer value representing the desired dimensionality of the particles: 1 &rarr;  1D, 2 &rarr;  2D, 3 &rarr;  3D.
 */
void DPMBase::setParticleDimensions(unsigned int particleDimensions)
{
    if (particleDimensions >= 1 && particleDimensions <= 3)
    {
        particleDimensions_ = particleDimensions;
        particleHandler.computeAllMasses();
    }
    else
    {
        logger(WARN, "Error in setParticleDimensions; particleDimensions % is not 1, 2 or 3, setting it to "
                     "systemDimension now (%)", particleDimensions, systemDimensions_);
        particleDimensions_ = systemDimensions_;
    }
}

/*!
 * \return particleDimensions_ The dimensionality of the <B>particles</B>. (Note that the <B>system</B> may possess a different dimensionality).
 */

unsigned int DPMBase::getParticleDimensions() const
{
    return particleDimensions_;
}

/*!
 * \return restartVersion_
 */
std::string DPMBase::getRestartVersion() const
{
    return restartVersion_;
}

/*!
 * \param[in] newRV
 */

void DPMBase::setRestartVersion(std::string newRV)
{
    restartVersion_ = newRV;
}

/*!
 * \return restarted_
 */

bool DPMBase::getRestarted() const
{
    return restarted_;
}

/*!
 * \param[in] newRestartedFlag
 */
void DPMBase::setRestarted(bool newRestartedFlag)
{
    restarted_ = newRestartedFlag;
    //setAppend(new_);
}

/*!
 * \return <tt>true</tt> if the "append" option is on; <tt>false</tt> if the "append" option is off.
 */
bool DPMBase::getAppend() const
{
    return append_;
}

/*!
 * Used when a code is restarted. If set to <tt>true</tt>, existing files (other than
 * restart files, which are always overwritten) will be appended, rather than
 * new files being created and/or old files being overwritten.
 * If <tt>false</tt>, files will simply be overwritten.
* \param[in] newAppendFlag
 */
void DPMBase::setAppend(bool newAppendFlag)
{
    append_ = newAppendFlag;
}

/*!
 * \return elasticEnergy The total elastic energy of all current particle interactions.
 */
Mdouble DPMBase::getElasticEnergy() const
{
    Mdouble elasticEnergy = 0.0;
    // JMFT: Note that we do count the elastic energy of fixed particles here.
    for (const BaseInteraction* c : interactionHandler)
    {
        elasticEnergy += c->getElasticEnergy();
    }
    return elasticEnergy;
}

/*!
 * \return kineticEnergy The total kinetic energy of all particles.
 */
Mdouble DPMBase::getKineticEnergy() const
{
    Mdouble kineticEnergy = 0;
    for (const BaseParticle* const p : particleHandler)
    {
        if (!(p->isFixed()))
        {
            kineticEnergy += .5 * p->getMass() * p->getVelocity().getLengthSquared();
        }
    }
    return kineticEnergy;
}

/*!
 * \return gravitationalEnergy The total gravitational potential energy of all particles (relative
 * to the origin).
 */
Mdouble DPMBase::getGravitationalEnergy() const
{
    Mdouble gravitationalEnergy = 0;
    for (const BaseParticle* const p : particleHandler)
    {
        // Don't consider fixed particles. 'Fixed' particles aren't necessarily
        // stationary; it just means their position is prescribed.
        if (!(p->isFixed()))
        {
            gravitationalEnergy += p->getMass() * Vec3D::dot((-getGravity()), p->getPosition());
        }
    }
    return gravitationalEnergy;
}

Mdouble DPMBase::getRotationalEnergy() const
{
    Mdouble kineticEnergy = 0;
    for (const BaseParticle* const p : particleHandler)
    {
        if (!(p->isFixed()))
        {
            kineticEnergy += .5 * Vec3D::dot(p->getAngularVelocity(), p->getInertia() * p->getAngularVelocity());
        }
    }
    return kineticEnergy;
}

Mdouble DPMBase::getTotalEnergy() const {
    return getElasticEnergy()+getKineticEnergy()+getGravitationalEnergy()+getRotationalEnergy();
}

/*!
 * \return double
 */
Mdouble DPMBase::getTotalMass() const
{
    /*
    double mass_sum = 0;
    for (std::vector<BaseParticle*>::const_iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
        if (!(*it)->isFixed())
            mass_sum += (*it)->getMass();
    return mass_sum;
    */
    return particleHandler.getMass();
}

/*!
 * \details
 * Returns the centre of mass of particles, excluding fixed particles.
 */
Vec3D DPMBase::getCentreOfMass() const
{
    return particleHandler.getCentreOfMass();
}

/*!
 * \details
 * Returns the total momentum in the system, excluding fixed particles (which will usually, but not always, have velocity 0)
 * \return Vec3D
 */
Vec3D DPMBase::getTotalMomentum() const
{
    return particleHandler.getMomentum();
    /*
    Vec3D total_momentum = Vec3D(0,0,0);
    for (std::vector<BaseParticle*>::const_iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
        if (!(*it)->isFixed())
            total_momentum += (*it)->getMass() * (*it)->getVelocity();
    return total_momentum;
    */
}

/*!
 * \return double
 */
/* JMFT: This information is placed on the last entry on each line of the .data
 * file. That space is, in general, reserved for 'additional' information.
 */
Mdouble DPMBase::getInfo(const BaseParticle& p) const
{
//    return p.getSpecies()->getId(); // was getIndex()
    return p.getInfo();
}

/*!
 * \details Determines whether two particles are distinct and in contact by comparing the separation of their centres to their (interaction) radii.
 * \n \n
 * Firstly, checks if the two particles are different (if pI == pJ, the result is "false").
 * Secondly, if the two particles are distinct, finds the distance between the two particles' centres
 * (<TT>getDistanceSquared(pI->getPosition(), pJ->getPosition()))</TT>)
 * and tests whether the separation of the particles is less than the sum of their radii
 * (<TT>pI->getInteractionRadius() + pJ->getInteractionRadius()</TT>).
 * If so, the bool returns  "true", i.e. the particles are in contact.
 * \param[in] pI A pointer to a particle
 * \param[in] pJ A pointer to a second particle
 * \return bool (True or False) - lets the user know whether two particles are in contact
 */
bool DPMBase::areInContact(const BaseParticle* pI, const BaseParticle* pJ)
{
    return (pI != pJ && pI->isInContactWith(pJ));
}

/*!
 * \details no implementation but can be overriden in its derived classes.
 */
void DPMBase::actionsBeforeTimeLoop()
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::hGridActionsBeforeTimeLoop()
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::actionsOnRestart()
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::hGridActionsBeforeTimeStep()
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::hGridInsertParticle(BaseParticle* obj UNUSED)
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::hGridUpdateParticle(BaseParticle* obj UNUSED)
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::hGridRemoveParticle(BaseParticle* obj UNUSED)
{
}

/*!
 * \return bool (True or False)
 */
bool DPMBase::getHGridUpdateEachTimeStep() const
{
    return true;
}

/*!
 * \brief Function that checks if the mpi particle should really be inserted by the current domain
 * \details When adding a particle, all domains "add" the particle to enable communication between processors
 * However not very domain should _add_ the particle, only the domain that actually contains the particle
 * There is one exception, if an MPI Particle is added (which is not physically in the current domain),
 * this has already been approved by the domain and hence it should return true.
 * \param[in] P Pointer to a baseParticle that requires an insertion check
 * \return Returns if the baseParticle should be inserted or not
 */
bool DPMBase::mpiInsertParticleCheck(BaseParticle* P)
{
#ifdef MERCURYDPM_USE_MPI
    //If only one core is used (i.e. domainHandler is empty) then the result is always true
    if (domainHandler.getSize() == 0)
    {
        return true;
    }
    //Get the current domain
    Domain* domain = domainHandler.getCurrentDomain();

    //Check if the particle is in the current domain
    if(domain->containsParticle(P))
    {
        //When adding a particle inside the domain, this should always be true
        return true;
    }
    else
    {
        //MPI particles that are inserted in the communication zone should still be inserted
        return (P->isMPIParticle());
    }
#else
    return false;
#endif
}

/*!
 * \brief Checks if the position of the particle is in an mpi communication zone or not
 * \param[in] particle Pointer to a base particle
 * \return Returns if the particle is in the communication zone (true) or not (false)
 */
bool DPMBase::mpiIsInCommunicationZone(BaseParticle* particle)
{

    bool insideCommunicationZone = false;
#ifdef MERCURYDPM_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();

    //Check for the current domain if the particle is within the communication domain
    int val = domainHandler.getCurrentDomain()->isInCommunicationZone(particle);

    //The root gathers all results
    int *list = nullptr;
    if (PROCESSOR_ID == 0)
    {
        list = new int [NUMBER_OF_PROCESSORS];
    }
    communicator.gather(val,list);

    //Compute the global value
    //if on any processor the val is true, we have to do the communcation step
    /// \todo MX: In future this might be replaced with an allReduce command and the appropriate operation.
    int result = 0;
    if (PROCESSOR_ID == 0)
    {
        for (int i = 0; i< NUMBER_OF_PROCESSORS; i++)
        {
            if (list[i] == 1)
            {
                result = 1;
                break;
            }
        }
    }

    //The root now tells the other processors what the global value for the interaction is
    communicator.broadcast(result);

    //Convert the result back to bool
    insideCommunicationZone = result;
#endif
    return insideCommunicationZone;
}

/*!
 * \brief This function inserts a particle in the mpi communication boundaries
 * \param[in] particle Pointer to a base particle that needs to be inserted in the communication boundaries
 */
void DPMBase::insertGhostParticle(BaseParticle* particle)
{
#ifdef MERCURYDPM_USE_MPI
    //mpi particles only exist when there is more than one domain
    if (domainHandler.getSize() > 0)
    {
        //Add the particle to the mpi domain
        domainHandler.getCurrentDomain()->addParticle(particle);
    }

    //If periodic boundaries are present..
    if (periodicBoundaryHandler.getSize() > 0)
    {
        periodicBoundaryHandler.addNewParticle(particle);
    }
#endif
}

/*!
 * \brief Checks if the domain and periodicBoundaryHandler need to update the interaction distance.
 * If this is the case it will update the ghost particles as well.
 * \details When adding a new particle with a larger interaction radius than all previous particles,
 *  the domain needs to update the interactionDistance and initialise all particles that are now also
 *  included in the communication zones. All other domains are updated as well.
 * \param[in] P Pointer to a baseParticle that recently has been added to the simulation
 */
void DPMBase::updateGhostGrid(BaseParticle* P)
{
#ifdef MERCURYDPM_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1) { return; }

    //Check if the interactionRadius of the BaseParticle is larger than given in the domain
    Domain* domain = domainHandler.getCurrentDomain();
    if(2.0*P->getMaxInteractionRadius() > domainHandler.getInteractionDistance())
    {
        logger(VERBOSE,"Processor % | Updating mpi grid. Old interactionDistance: %, new interactionDistance %.",
            PROCESSOR_ID,domainHandler.getInteractionDistance(),2.0*P->getMaxInteractionRadius());

        //Update the interactionDistance in the domain and periodicBoundaryHandler
        domainHandler.setInteractionDistance(2.0*P->getMaxInteractionRadius());
        periodicBoundaryHandler.setInteractionDistance(2.0*P->getMaxInteractionRadius());

        //Find new ghost particless
        domainHandler.addNewParticles();
        periodicBoundaryHandler.addNewParticles();
    }
#endif
}


/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::actionsBeforeTimeStep()
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::actionsAfterSolve()
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::actionsAfterTimeStep()
{
}

/*!
 * \details This function is overridden by StatisticsVector
 */
void DPMBase::initialiseStatistics()
{
    cgHandler.initialise();
}

/*!
 * \details This function is overridden by StatisticsVector
 */
void DPMBase::outputStatistics()
{
    //cgHandler.evaluate();
}

void DPMBase::gatherContactStatistics()
{
    for (BaseInteraction* c : interactionHandler)
    {
        c->gatherContactStatistics();
    }
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::gatherContactStatistics(unsigned int index1 UNUSED, int index2 UNUSED, Vec3D Contact UNUSED,
                                      Mdouble delta UNUSED, Mdouble ctheta UNUSED, Mdouble fdotn UNUSED,
                                      Mdouble fdott UNUSED, Vec3D P1_P2_normal_ UNUSED, Vec3D P1_P2_tangential UNUSED)
{
}

/*!
 * \details This function is overridden by StatisticsVector
 */
void DPMBase::processStatistics(bool usethese UNUSED)
{
}

/*!
 * \details This function is overridden by StatisticsVector
 */
void DPMBase::finishStatistics()
{
    cgHandler.finish();
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::hGridUpdateMove(BaseParticle*, Mdouble)
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::hGridActionsBeforeIntegration()
{
}

/*!
 * \details no implementation but can be overidden in its derived classes.
 */
void DPMBase::hGridActionsAfterIntegration()
{
}

/*!
 * \details That is to say that the first <TT>n</TT> particles in the handler will have their masses
 * and inertiae set to infinity and their velocities set to zero, making them effectively immobile and immovable
 * (see also \ref  BaseParticle::fixParticle() ).
 * If <TT>n</TT> exceeds the number of particles in the handler, then all particles will be "fixed".
 * \param[in] n The number of particles to "fix".
 * \todo  Question: can this function be called during the run of the program? i.e. can particles
	be made to suddenly "become fixed"?
 */
void DPMBase::setFixedParticles(unsigned int n)
{
    for (unsigned int i = 0; i < std::min(particleHandler.getSize(), n); i++)
        particleHandler.getObject(i)->fixParticle();
}

/*!
 * \details Gets and prints the current simulation time (<tt>getTime()</tt>) and the currently set maximum simulation time
 * (<tt>getTimeMax()</tt>) .
 */
void DPMBase::printTime() const
{
#ifdef MERCURYDPM_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    if (communicator.getProcessorID() != 0)
    {
        return;
    }

#endif
  logger(INFO, "t=%3.6, tmax=%3.6", getTime(), getTimeMax());
}

/*!
 * \details Used within the main loop of the 'solve()' routine to let the code know whether or
 * not the time step should continue to be advanced, i.e. whether the simulation should be continued.
 * By default this is always <tt>true</tt> but the user may redefine it to return <tt>false</tt>
 * under certain desired circumstances.
 * \return bool (True or False)
 */
bool DPMBase::continueSolve() const
{
    return true;
}

/*!
 * \brief A virtual function with no implementation but can be overriden
 */
void DPMBase::setupInitialConditions()
{
}

/*!
 * \details If the "append" flag is off (false) - i.e. the file has not been restarted - creates a header for the output
 * ".ene" file. The headers simply give titles to each of the columns output to this file.
 * \n \n
 * If the "append" flag is on (true), i.e. if the file has restarted, simply ends the function without writing
 * a header - otherwise restarted files will have a random header at the point of restart, making
 * data processing more difficult...
 * \n \n
 * The function takes as an argument an output stream, "os", which tells the compiler where to output the headers,
 * if they are required.
* \param[in] os The output stream to the file in which the header should be written.
*/
void DPMBase::writeEneHeader(std::ostream& os) const
{
    //only write if we don't restart
    if (getAppend())
        return;

    /// \todo JMFT: Get rid of aligned columns. They make things too wide. (changed back)
    /// \todo {Why is there a +6 here? TW: to get the numbers and title aligned}
    /// \todo Add number of particles to this file (change from Jonny to be added later)
    long width = os.precision() + 6;
    os << std::setw(width)
       << "time " << std::setw(width)
       << "gravitEnergy " << std::setw(width) //gravitational potential energy
       << "traKineticEnergy " << std::setw(width) //translational kinetic energy
       << "rotKineticEnergy " << std::setw(width) //rotational kE
       << "elasticEnergy " << std::setw(width)
       << "centerOfMassX " << std::setw(width)
       << "centerOfMassY " << std::setw(width)
       << "centerOfMassZ\n";
}

/*!
 * \details The function takes as an argument an output stream and, to the corresponding file,
 * outputs the relevant header <B>and data</B> for the ".fstat"-type output files.
 * For more information regarding the fstat file format, please refer to the user guide on the
 * <A HREF="http://mercurydpm.org/assets/downloads/MercuryLesson/MercuryDPMLessonSlides.pdf">Mercury website</A>
 * \param[in] os The output stream (e.g. a file stream or cout)
 */
void DPMBase::writeFstatHeader(std::ostream& os) const
{

    // line #1: time, volume fraction
    // line #2: wall box: wx0, wy0, wz0, wx1, wy1, wz1
    // line #3: radii-min-max & moments: rad_min, rad_max, r1, r2, r3, r4
    os << "#"
       << " " << getTime()
       << " " << 1 //marker that these are fstat files with contact point instead of center point
       << '\n';
    os << "#"
       << " " << getXMin()
       << " " << getYMin()
       << " " << getZMin()
       << " " << getXMax()
       << " " << getYMax()
       << " " << getZMax()
       << '\n';
    os << "#"
       << " ";

    if (!(particleHandler.getSmallestParticleLocal() == nullptr))
    {
        os << particleHandler.getSmallestParticleLocal()->getRadius();
    }
    else
    {
        os << std::numeric_limits<double>::quiet_NaN();
    }
    os << " ";
    if (!(particleHandler.getLargestParticleLocal() == nullptr))
    {
        os << particleHandler.getLargestParticleLocal()->getRadius();
    }
    else
    {
        os << std::numeric_limits<double>::quiet_NaN();
    }

    os << " " << 0
       << " " << 0
       << " " << 0
       << " " << 0
       << '\n';
    //B: write data
    for (BaseInteraction* c : interactionHandler)
    {
        c->writeToFStat(os, getTime());
    }
    //os << std::flush;
}

/*!
 * \details The function cycles over all particles within the system (or rather, the particleHandler), creating sums of the relevant energies and
 * "mass lengths" (m.x, m.y, m.z) from which the system's centre of mass can also be calculated. The summed energy values and calculated centre of
 * mass values are then output to the file corresponding to "os" alongside the current time step.
 * \n \n
 * A check is performed - <TT>if (!p->isFixed())</TT> - to ensure that calculations are not performed on fixed particles, as these are assigned an effectively infinite mass
 * and would hence cause compiler issues.
 * \param[in] os The output stream to which the data will be written
 */
void DPMBase::writeEneTimeStep(std::ostream& os) const
{
    if (eneFile.getCounter() == 1 || eneFile.getFileType() == FileType::MULTIPLE_FILES ||
        eneFile.getFileType() == FileType::MULTIPLE_FILES_PADDED)
        writeEneHeader(os);

    const Mdouble m = particleHandler.getMass();
    const Vec3D com = particleHandler.getMassTimesPosition();
    //Ensure the numbers fit into a constant width column: for this we need the precision given by the operating system,
    //plus a few extra characters for characters like a minus and scientific notation.
    const static int width = os.precision() + 6;
    os << std::setw(width) << getTime()
       << " " << std::setw(width) << -Vec3D::dot(getGravity(), com)
       << " " << std::setw(width) << particleHandler.getKineticEnergy()
       << " " << std::setw(width) << particleHandler.getRotationalEnergy()
       << " " << std::setw(width) << getElasticEnergy()
       // we need to write x, y and z coordinates separately, otherwise the width of the columns is incorrect
       << " " << std::setw(width)
       << (m == 0 ? constants::NaN : com.X / m) //set to nan because 0/0 implementation in gcc and clang differs
       << " " << std::setw(width) << (m == 0 ? constants::NaN : com.Y / m)
       << " " << std::setw(width) << (m == 0 ? constants::NaN : com.Z / m)
       << std::endl;
}

void DPMBase::writeVTKFiles() const
{
    if (wallHandler.getWriteVTK() == FileType::ONE_FILE && wallVTKWriter_.getFileCounter() == 0)
    {
        wallVTKWriter_.writeVTK();
    } else if (wallHandler.getWriteVTK() == FileType::MULTIPLE_FILES
        || wallHandler.getWriteVTK() == FileType::MULTIPLE_FILES_PADDED) {
        wallVTKWriter_.writeVTK();
    } // else do nothing

    if (getParticlesWriteVTK() || getSuperquadricParticlesWriteVTK())
    {
        vtkWriter_->writeVTK();
    }  // else do nothing

    if (interactionHandler.getWriteVTK() != FileType::NO_FILE)
    {
        interactionVTKWriter_.writeVTK();
    }

    if (boundaryHandler.getWriteVTK())
    {
        boundaryVTKWriter_.writeVTK();
    }

    if (wallHandler.getWriteDetailsVTKAny())
    {
        wallDetailsVTKWriter_.writeVTK();
    }

    //only write once
    bool writePython = getParticlesWriteVTK() || wallHandler.getWriteVTK() != FileType::NO_FILE ||
                       interactionHandler.getWriteVTK() != FileType::NO_FILE ||
                       wallHandler.getWriteDetailsVTKAny();

    writePython &= forceWritePythonFileForVTKVisualisation_;

    if (writePython && getTime() == 0)
    {
        writePythonFileForVTKVisualisation();
    }
}

void DPMBase::writePythonFileForVTKVisualisation() const
{
#ifdef MERCURYDPM_USE_MPI
    if (PROCESSOR_ID == 0)
    {
#endif
    logger(INFO, "Writing python script for paraview visualisation");

    std::string script = "#script to visualise the output of data2pvd of MercuryDPM in paraview.\n"
                         "#usage: change the path below to your own path, open paraview\n"
                         "#Tools->Python Shell->Run Script->VisualisationScript.py\n"
                         "#or run paraview --script=VisualisationScript.py \n"
                         "\n"
                         "from paraview.simple import *\n"
                         "import os\n"
                         "import glob\n"
                         "import re # for natural sorting\n\n"
                         "# Path to directory containing vtu files\n"
                         "# Hardcoded (for windows you might need to change start of path, for example: from /mnt/c/MyFolder to c:/MyFolder)\n"
                         "dirPath = '" + helpers::getPath() + "'\n"
                         "# Or directory containing this script\n"
                         "#dirPath = os.path.dirname(os.path.realpath(__file__))\n"
                         "# Change to directory\n"
                         "os.chdir(dirPath)\n\n";

    script += "# Some functions to achieve natural sorting for the numbering of the files (regardless of them being padded or not etc.)\n"
              "# See https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside for more info\n"
              "def atoi(text):\n"
              "\treturn int(text) if text.isdigit() else text\n"
              "\n"
              "def natural_keys(text):\n"
              "\treturn [ atoi(c) for c in re.split(r'(\\d+)', text) ]\n\n";

    script += "simName = '" + getName() + "' # name of the simulation\n\n";

    if (getParticlesWriteVTK())
    {
#ifdef MERCURYDPM_USE_MPI
        if (NUMBER_OF_PROCESSORS > 1)
        {
            script += "# PARTICLES ##########################################################\n"
                      "for processorNumber in range(0, " + std::to_string(NUMBER_OF_PROCESSORS) + "):\n"
                      "\t# Load data in any order and sort it\n"
                      "\tDataParticles = glob.glob('./' + simName + 'Processor_' + str(processorNumber) + '_Particle_*.vtu')\n"
                      "\tDataParticles.sort(key = natural_keys)"
                      "\t# Load the data and visualise it in paraview\n"
                      "\tparticles = XMLUnstructuredGridReader(FileName = DataParticles, registrationName = simName + 'Processor_' + str(processorNumber) + '_Particle_*')\n"
                      "\tglyphP = Glyph(particles, registrationName = 'Glyph Particle')\n" // From here on no differences with (non-)MPI
                      "\tglyphP.GlyphType = 'Sphere'\n"
                      "\tglyphP.ScaleArray = 'Radius'\n"
                      "\tglyphP.ScaleFactor = 2\n"
                      "\tglyphP.GlyphMode = 'All Points'\n"
                      "\tparticlesDisplay = Show(glyphP)\n"
                      "\t#ColorBy(particlesDisplay, ('POINTS', 'Velocity', 'Magnitude'))\n"
                      "\t#particlesDisplay.Opacity = 0.5\n"
                      "\n\n";
        }
        else
        {
#endif
        script += "# PARTICLES ##########################################################\n"
                  "# Load data in any order and sort it\n"
                  "DataParticles = glob.glob('./' + simName + 'Particle_*.vtu')\n"
                  "DataParticles.sort(key = natural_keys)"
                  "# Load the data and visualise it in paraview\n"
                  "particles = XMLUnstructuredGridReader(FileName = DataParticles, registrationName = simName + 'Particle_*')\n"
                  "glyphP = Glyph(particles, registrationName = 'Glyph Particle')\n" // From here on no differences with (non-)MPI
                  "glyphP.GlyphType = 'Sphere'\n"
                  "glyphP.ScaleArray = 'Radius'\n"
                  "glyphP.ScaleFactor = 2\n"
                  "glyphP.GlyphMode = 'All Points'\n"
                  "particlesDisplay = Show(glyphP)\n"
                  "#ColorBy(particlesDisplay, ('POINTS', 'Velocity', 'Magnitude'))\n"
                  "#particlesDisplay.Opacity = 0.5\n"
                  "\n\n";
#ifdef MERCURYDPM_USE_MPI
        }
#endif
    }

    if (wallHandler.getWriteVTK() != FileType::NO_FILE)
    {
        script += "# WALLS ##############################################################\n"
                  "# Load data in any order and sort it\n"
                  "DataWalls = glob.glob('./' + simName + 'Wall_*.vtu')\n"
                  "DataWalls.sort(key = natural_keys)"
                  "# Load the data and visualise it in paraview\n"
                  "walls = XMLUnstructuredGridReader(FileName = DataWalls, registrationName = simName + 'Wall_*')\n"
                  "wallsDisplay = Show(walls)\n"
                  "#wallsDisplay.Opacity = 0.5\n"
                  "\n\n";
    }

    if (interactionHandler.getWriteVTK() != FileType::NO_FILE)
    {
#ifdef MERCURYDPM_USE_MPI
        if (NUMBER_OF_PROCESSORS > 1)
        {
            script += "# INTERACTIONS #######################################################\n"
                      "for processorNumber in range(0, " + std::to_string(NUMBER_OF_PROCESSORS) + "):\n"
                      "\t# Load data in any order and sort it\n"
                      "\tDataInteractions = glob.glob('./' + simName + 'Processor_' + str(processorNumber) + '_Interaction_*.vtu')\n"
                      "\tDataInteractions.sort(key = natural_keys)"
                      "\t# Load the data and visualise it in paraview\n"
                      "\tinteractions = XMLUnstructuredGridReader(FileName = DataInteractions, registrationName = simName + 'Processor_' + str(processorNumber) + '_Interaction_*')\n"
                      "\tglyphI = Glyph(interactions, registrationName = 'Glyph Interaction')\n" // From here on no differences with (non-)MPI
                      "\tglyphI.GlyphType = 'Arrow'\n"
                      "\tglyphI.OrientationArray = 'Normal'\n"
                      "\tglyphI.ScaleArray = 'Force'\n"
                      "\tglyphI.ScaleFactor = 0.01\n"
                      "\tglyphI.GlyphMode = 'All Points'\n"
                      "\tinteractionsDisplay = Show(glyphI)\n"
                      "\t#ColorBy(interactionsDisplay, ('POINTS', 'Force'))\n"
                      "\t#interactionsDisplay.Opacity = 0.5\n"
                      "\n\n";
        }
        else
        {
#endif
        script += "# INTERACTIONS #######################################################\n"
                  "# Load data in any order and sort it\n"
                  "DataInteractions = glob.glob('./' + simName + 'Interaction_*.vtu')\n"
                  "DataInteractions.sort(key = natural_keys)"
                  "# Load the data and visualise it in paraview\n"
                  "interactions = XMLUnstructuredGridReader(FileName = DataInteractions, registrationName = simName + 'Interaction_*')\n"
                  "glyphI = Glyph(interactions, registrationName = 'Glyph Interaction')\n" // From here on no differences with (non-)MPI
                  "glyphI.GlyphType = 'Arrow'\n"
                  "glyphI.OrientationArray = 'Normal'\n"
                  "glyphI.ScaleArray = 'Force'\n"
                  "glyphI.ScaleFactor = 0.01\n"
                  "glyphI.GlyphMode = 'All Points'\n"
                  "interactionsDisplay = Show(glyphI)\n"
                  "#ColorBy(interactionsDisplay, ('POINTS', 'Force'))\n"
                  "#interactionsDisplay.Opacity = 0.5\n"
                  "\n\n";
#ifdef MERCURYDPM_USE_MPI
        }
#endif
    }

    if (wallHandler.getWriteDetailsVTK(WallHandler::DetailsVTKOptions::BOUNDINGBOX) != FileType::NO_FILE)
    {
        script += "# WALL DETAILS - BOUNDING BOX ########################################\n"
                  "# Load data in any order and sort it\n"
                  "DataWDBoundingBox = glob.glob('./' + simName + 'WallDetailsBoundingBox_*.vtu')\n"
                  "DataWDBoundingBox.sort(key = natural_keys)"
                  "# Load the data and visualise it in paraview\n"
                  "wdBoundingBox = XMLUnstructuredGridReader(FileName = DataWDBoundingBox, registrationName = simName + 'WallDetailsBoundingBox_*')\n"
                  "wdBoundingBoxDisplay = Show(wdBoundingBox)\n"
                  "wdBoundingBoxDisplay.Representation = 'Wireframe'\n"
                  "\n\n";
    }

    if (wallHandler.getWriteDetailsVTK(WallHandler::DetailsVTKOptions::NURBSWALL) != FileType::NO_FILE)
    {
        script += "# WALL DETAILS - NURBS WALL ##########################################\n"
                  "# Load data in any order and sort it\n"
                  "DataWDNurbsWall = glob.glob('./' + simName + 'WallDetailsNurbsWall_*.vtu')\n"
                  "DataWDNurbsWall.sort(key = natural_keys)"
                  "# Load the data and visualise it in paraview\n"
                  "wdNurbsWall = XMLUnstructuredGridReader(FileName = DataWDNurbsWall, registrationName = simName + 'WallDetailsNurbsWall_*')\n"
                  "wdNurbsWallDisplay = Show(wdNurbsWall)\n"
                  "wdNurbsWallDisplay.Representation = 'Wireframe'\n"
                  "ColorBy(wdNurbsWallDisplay, ('POINTS', 'ID'))\n"
                  "glyphWDNW = Glyph(wdNurbsWall, registrationName='Glyph NURBS Control Point')\n"
                  "glyphWDNW.GlyphType = 'Sphere'\n"
                  "glyphWDNW.ScaleArray = 'Weight'\n"
                  "glyphWDNW.ScaleFactor = 0.05\n"
                  "glyphWDNW.GlyphMode = 'All Points'\n"
                  "glyphWDNWDisplay = Show(glyphWDNW)\n"
                  "ColorBy(glyphWDNWDisplay, ('POINTS', 'ID'))\n"
                  "\n\n";
    }

    script += "GetAnimationScene().PlayMode = 'Snap To TimeSteps'\n"
              "Render()\n"
              "ResetCamera()\n";

    helpers::writeToFile(getName() + ".py", script);
#ifdef MERCURYDPM_USE_MPI
    } // end of communicator is root statement
#endif
}


/*!
 * \param[in] os
 */
void DPMBase::outputXBallsData(std::ostream& os) const
{


    //Set the correct formation based of dimension if the formation is not specified by the user

    unsigned int format;
    switch (getSystemDimensions())
    {
        case 2:
            format = 8;
            break;
        case 3:
            format = 14;
            break;
        default:
            logger(ERROR, "Unknown system dimension");
    }

    unsigned int numberOfParticles = particleHandler.getNumberOfRealObjectsLocal();

    // This outputs the location of walls and how many particles there are to file this is required by the xballs plotting
    if (format != 14) // dim = 1 or 2
    {
        os << numberOfParticles
           << " " << getTime()
           << " " << getXMin()
           << " " << getYMin()
           << " " << getXMax()
           << " " << getYMax()
           << " " << std::endl;
    }
    else
    {
        //dim==3
        os << numberOfParticles
           << " " << getTime()
           << " " << getXMin()
           << " " << getYMin()
           << " " << getZMin()
           << " " << getXMax()
           << " " << getYMax()
           << " " << getZMax()
           << " " << std::endl;
    }

    // This outputs the particle data
    for (unsigned int i = 0; i < particleHandler.getSize(); i++)
    {
#ifdef MERCURYDPM_USE_MPI
        if (!particleHandler.getObject(i)->isPeriodicGhostParticle() && !particleHandler.getObject(i)->isMPIParticle())
        {
            outputXBallsDataParticle(i, format, os);
        }
#else
        outputXBallsDataParticle(i, format, os);
#endif
    }
#ifdef DEBUG_OUTPUT
    logger(DEBUG, "Have output the properties of the problem to disk ");
#endif
}

/*!
 * \details This function reads a .data file,
 *          which contains info about each particle's
 *          position, velocity, angular velocity, radius ...info.
 *          See also MD::readRestartFile
 *          For XBalls:
 *          Can read in format_ 14 - 8 or format_ 7 data format.
 *          This code saves in format_ 8 for 2D and format_ 14 for 3D.
 *          So if no extra parameters are specified it will assume many parameters,
 *          like density cannot be set using the data file.
 * use of string instead of string& b/c this function is often used with a string literal
 * \param[in] fileName
 * \param[in] format (format for specifying if its for 2D or 3D data)
 * \return bool (True or False)
 */
bool DPMBase::readDataFile(std::string fileName, unsigned int format)
{
    //default value: dataFile.getFullName()
    if (!fileName.compare(""))
        fileName = dataFile.getFullName();

    std::string oldFileName = dataFile.getName();
    unsigned oldCounter = dataFile.getCounter();
    //Updates the name of the data file to the user-input from the argument.
    dataFile.setName(fileName);
    //opens a filestream of the input type
    dataFile.open(std::fstream::in);
    //Checks if the file has been successfully opened...
    if (!dataFile.getFstream().is_open() || dataFile.getFstream().bad())
    {
        //...and if not, ends the function and returns "false"
        logger(WARN, "Loading data file % failed.", fileName);
        return false;
    }

    //retrieves and saves the "FileType" of the file
    FileType fileTypeData = dataFile.getFileType();
    dataFile.setFileType(FileType::ONE_FILE);
    readNextDataFile(format);
    dataFile.setFileType(fileTypeData);
    dataFile.close();
    dataFile.setName(oldFileName);
    dataFile.setCounter(oldCounter);
    return true;
}

/*!
 * use of string instead of string& b/c this function is often used with a string literal
 * \param[in] fileName
 * \return bool (True or False)
 */
bool DPMBase::readParAndIniFiles(const std::string fileName)
{
    //Opens the par.ini file
    std::fstream file;
    file.open(fileName, std::fstream::in);
    if (!file.is_open() || file.bad())
    {
        //std::cout << "Loading par.ini file " << filename << " failed" << std::endl;
        return false;
    }

    Mdouble doubleValue;
    int integerValue;

    // inputfile par.ini
    // line 1 =============================================================
    // Example: 1 1 0
    //   1: integer (0|1) switches from non-periodic to periodic
    //      integer (5|6) does 2D integration only (y-coordinates fixed)
    //                    and switches from non-periodic to periodic
    //      integer (11) uses a quarter system with circular b.c.
    file >> integerValue;
    //~ std::cout << "11" << integerValue << std::endl;
    if (integerValue == 0)
    {
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, -1, 0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 0, 1), Vec3D(0, 0, getZMax()));
        wallHandler.copyAndAddObject(w0);
    }
    else if (integerValue == 1)
    {
        PeriodicBoundary b0;
        b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0, 0, 1), getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(b0);
    }
    else if (integerValue == 5)
    {
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1, 0, 0), Vec3D(-getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, -1, 0), Vec3D(0, -getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);

    }
    else if (integerValue == 6)
    {
        PeriodicBoundary b0;
        b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0, 0, 1), getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(b0);
    }
    else
    {
        logger(ERROR, "Error in par.ini: line 1, value 1 is %", integerValue);
    }

    //   2: integer (0|1) dont use | use the search pattern for linked cells
    file >> integerValue; //ignore

    //   3: real - gravity in z direction: positive points upwards
    file >> doubleValue;
    setGravity(Vec3D(0.0, 0.0, doubleValue));

    // line 2 =============================================================
    // Example: -10000 .5e-2
    //   1: time end of simulation - (negative resets start time to zero
    //                                and uses -time as end-time)
    file >> doubleValue;
    if (doubleValue < 0)
        setTime(0);
    setTimeMax(fabs(doubleValue));

    //   2: time-step of simulation
    file >> doubleValue;
    setTimeStep(doubleValue);

    // line 3 =============================================================
    // Example: 1e-1 100
    file >> doubleValue;
    if (doubleValue >= 0)
    {
        //   1: time-step for output on time-series protocoll file  -> "ene"
        unsigned int savecount = static_cast<unsigned int>(round(doubleValue / getTimeStep()));
        setSaveCount(savecount);

        //   2: time-step for output on film (coordinate) file      -> "c3d"
        //      (fstat-output is coupled to c3d-output time-step)
        file >> doubleValue;
        savecount = static_cast<unsigned int>(round(doubleValue / getTimeStep()));
        dataFile.setSaveCount(savecount);
        fStatFile.setSaveCount(savecount);
    }
    else
    {
        //  or: ---------------------------------------------------------------
        //   1: negative number is multiplied to the previous log-output-time
        //   2: requires initial log-output time
        //   3: negative number is multiplied to the previous film-output-time
        //   4: requires initial film-output time
        logger(ERROR, "Error in par.ini: line 3, value 1 is %", doubleValue);
    }

    // line 4 =============================================================
    // Example: 2000 1e5 1e3 1e2
    //   1: particle density (mass=4/3*constants::pi*density*rad^3)
    file >> doubleValue;

    //clear species handler
    speciesHandler.clear();
    auto S = new LinearViscoelasticSlidingFrictionSpecies;
    speciesHandler.addObject(S);

    setParticleDimensions(3);
    S->setDensity(doubleValue);

    //   2: linear spring constant
    file >> doubleValue;
    S->setStiffness(doubleValue);

    //   3: linear dashpot constant
    file >> doubleValue;
    S->setDissipation(doubleValue);

    //   4: background damping dashpot constant
    file >> doubleValue;
    if (doubleValue != 0.0)
        logger(WARN, "Warning in par.ini: ignored background damping %", doubleValue);

    // line 5 =============================================================
    // Example: 0 0
    //   1: growth rate:  d(radius) = xgrow * dt
    file >> doubleValue;
    if (doubleValue != 0.0)
        logger(WARN, "Warning in par.ini: ignored growth rate %", doubleValue);

    //   2: target volume_fraction
    file >> doubleValue;
    if (doubleValue != 0.0)
        logger(WARN, "Warning in par.ini: ignored target volume_fraction %", doubleValue);

    file.close();
    //std::cout << "Loaded par.ini file " << filename << std::endl;
    return true;
}

/*!
 * \details
 * First, checks to see if the file type is \ref MULTIPLE_FILES or \ref MULTIPLE_FILES_PADDED and, if so, whether the file contains data (ending
 * the function if not).
 * Then, checks if the time corresponding to the current file exceeds the minimum value entered (tMin). If not, keeps looking through subsequent data files.
 * When a data file that satisfies t > tMin is found and successfully opened, the function returns true.
 * \n \n
 * Useful when fileType is chosen as \ref MULTIPLE_FILES or \ref MULTIPLE_FILES_PADDED, which write
 * data corresponding to each time step as a separate, consecutively numbered file (see \ref FileType).
 * \param[in] tMin Compared with the t value belonging to the file being checked to see if it is viable.
 * \param[in] verbose Allows the function to give output to the screen if desired.
 * \return bool - true if the next file is found, false if not.
 */
bool DPMBase::findNextExistingDataFile(Mdouble tMin, bool verbose)
{
    if (dataFile.getFileType() == FileType::MULTIPLE_FILES || dataFile.getFileType() == FileType::MULTIPLE_FILES_PADDED)
    {
        while (true)// This true corresponds to the if s
        {
            dataFile.open();
            //check if file exists and contains data
            int N;
            dataFile.getFstream() >> N;
            if (dataFile.getFstream().eof() || dataFile.getFstream().peek() == -1)
            {
                logger(WARN, "file % not found", dataFile.getName());
                return false;
            }
            //check if tmin condition is satisfied
            Mdouble t;
            dataFile.getFstream() >> t;
            if (t > tMin)
            {
                //set_file_counter(get_file_counter()-1);
                return true;
            }
            if (verbose)
                logger(VERBOSE, "Jumping file counter: %", dataFile.getCounter());
        }
    }
    return true;
}

/*!
 * \details Used to load data from existing .data files (particle positions, velocities, sizes etc.)
 * when restarting a simulation.
 *  The 'format' to choose depends on the formatting of the data files to be read in which varies,
 * for instance, if the code is 2D (format 8) or 3D (format 14).
 * \param[in] format
 */
bool DPMBase::readNextDataFile(unsigned int format)
{
    dataFile.open(std::fstream::in);
    //logger(INFO,"Reading %",dataFile.getFullName());
    //fStatFile.open();
    //Set the correct format based of dimension if the format is not explicitly specified by the user
    if (format == 0)
    {
        //checking the dimensionality of the system
        /// \todo make systemDimensions enum (2 or 3)
        switch (getSystemDimensions())
        {
            case 1:
                //if 2D, sets format 8 (data file has 8 columns for 2D)
            case 2:
                format = 8;
                break;
            case 3:
                //if 3D, sets format 14 (data file has 14 columns for 3D)
                format = 14;
                break;
        }
        //end case
    }
    //end if

    // read in the particle number (as a double)
    double doubleN = -1;
    dataFile.getFstream() >> doubleN;

    // If N cannot be read, we reached the end of the file
    if (doubleN == -1) return false;

    // converting N to an integer; skipping the line if there is a problem (this happens when there is a corrupt data file)
    unsigned N = doubleN;
    while (doubleN != N) {
        std::string dummy;
        getline(dataFile.getFstream(),dummy,'\n');
        logger(WARN,"Skipping bad line in data file: % %",doubleN, dummy);
        dataFile.getFstream() >> doubleN;
        N = doubleN;
    }
    //store the parameters you want to preserve:
    const size_t nHistory = std::min(N,particleHandler.getSize());
    std::vector<const ParticleSpecies*> species(nHistory);
    std::vector<bool> fix(nHistory);
    std::vector<Vec3D> positionHistory(nHistory);
    for (size_t i=0; i<nHistory; ++i) {
        const BaseParticle *p = particleHandler.getObject(i);
        species[i] = p->getSpecies();
        fix[i] = p->isFixed();
        positionHistory[i] = p->getPosition();
        //store from reading the restart file which particles are fixed and which species each particle has
    }

    BaseParticle* p;
    if (particleHandler.getSize() == 0) {
        logger.assert_always(speciesHandler.getSize()>0,"readData: species needs to be set first");
        p = new SphericalParticle(speciesHandler.getObject(0));
        p->unfix();
    } else {
        p = particleHandler.getObject(0)->copy();
    }

    //empty the particle handler
    particleHandler.clear();
    particleHandler.setStorageCapacity(N);

    //now fill it again with the read-in particles
    std::stringstream line;
    helpers::getLineFromStringStream(dataFile.getFstream(), line);
    Mdouble radius;
    Vec3D position, velocity, angle, angularVelocity;
    size_t indSpecies;
    //read all other data available for the time step
    if (format==7) {
        line >> time_ >> min_.x() >> min_.y() >> min_.z() >> max_.x() >> max_.y() >> max_.z();
        for (size_t i = 0; i < N; ++i) {
            helpers::getLineFromStringStream(dataFile.getFstream(), line);
            line >> position.X >> position.Z >> position.Y >> velocity.X >> velocity.Z
                 >> velocity.Y >> radius
                 >> indSpecies;
            p->setPosition(position);
            p->setVelocity(velocity);
            p->setOrientation({1, 0, 0, 0});
            p->setAngularVelocity({0.0, 0.0, 0.0});
            p->setRadius(radius);
            if (readSpeciesFromDataFile_)
                p->setSpecies(speciesHandler.getObject(indSpecies));
            else if (i < nHistory)
                p->setSpecies(species[i]);
            if (i < nHistory && fix[i])
                p->fixParticle();
            if (i < nHistory)
                p->setPreviousPosition(positionHistory[i]);
            else
                p->setPreviousPosition(p->getPosition());
            particleHandler.copyAndAddObject(*p);
            p->unfix();
        }
    } else if (format==8) {
        line >> time_ >> min_.x() >> min_.y() >> max_.x() >> max_.y();
        min_.z() = 0.0; max_.z() = 0.0;//For 2d functions we define Z to be of no interest
        for (size_t i = 0; i < N; ++i) {
            helpers::getLineFromStringStream(dataFile.getFstream(), line);
            line >> position.X >> position.Y >> velocity.X >> velocity.Y >> radius >> angle.Z >> angularVelocity.Z >> indSpecies;
            Quaternion q;
            q.setAngleZ(angle.Z);
            p->setPosition(position);
            p->setVelocity(velocity);
            p->setOrientation(q);
            p->setAngularVelocity(-angularVelocity);
            p->setRadius(radius);
            if (readSpeciesFromDataFile_) p->setSpecies(speciesHandler.getObject(indSpecies));
            else if (i < nHistory) p->setSpecies(species[i]);
            if (i < nHistory && fix[i]) p->fixParticle();
            if (i < nHistory) p->setPreviousPosition(positionHistory[i]);
            else p->setPreviousPosition(p->getPosition());
            particleHandler.copyAndAddObject(*p);
            p->unfix();
        } //end for all particles
    } else if (format==14) {
        //This is a 3D format_
        // \todo : Check bounds or get rid of this function
        line >> time_ >> min_.x() >> min_.y() >> min_.z() >> max_.x() >> max_.y() >> max_.z();
        for (size_t i = 0; i < N; ++i) {
            helpers::getLineFromStringStream(dataFile.getFstream(), line);
            line >> position >> velocity >> radius >> angle >> angularVelocity >> indSpecies;
#ifdef MERCURYDPM_USE_MPI
            //This is required for the CG tool. When reading the data file it is neseccary to know if a particle is an MPIParticle or not
            bool isMPIParticle = false;
            bool isPeriodicGhostParticle = false;
            if (NUMBER_OF_PROCESSORS > 1)
            {
                ///\todo TW should be checked it doesn't break marnix' code
                line >> isMPIParticle >> isPeriodicGhostParticle;
            }
#endif
            p->setPosition(position);
            p->setVelocity(velocity);
            p->setOrientationViaEuler(angle);
            p->setAngularVelocity(angularVelocity);
            p->setRadius(radius);
            if (readSpeciesFromDataFile_) {
                if (indSpecies<speciesHandler.getSize()) {
                    p->setSpecies(speciesHandler.getObject(indSpecies));
                } else {
                    logger(WARN, "Read in bad species data; species is not set");
                }
            } else if (i < nHistory)
                p->setSpecies(species[i]);
            if (i < nHistory && fix[i])
                p->fixParticle();
            if (i < nHistory) p->setPreviousPosition(positionHistory[i]);
            else p->setPreviousPosition(p->getPosition());
#ifdef MERCURYDPM_USE_MPI
            if (NUMBER_OF_PROCESSORS)
            {
                p->setMPIParticle(isMPIParticle);
                p->setPeriodicGhostParticle(isPeriodicGhostParticle);
            }
#endif
            particleHandler.copyAndAddObject(*p);
            p->unfix();
        } //end read into existing particles logger(INFO, "read % particles", particleHandler.getNumberOfObjects());
    } else if (format==15) {
        line >> time_ >> min_.x() >> min_.y() >> min_.z() >> max_.z() >> max_.y() >> max_.z();
        for (size_t i = 0; i < N; ++i) {
            helpers::getLineFromStringStream(dataFile.getFstream(), line);
            line >> position >> velocity >> radius >> angle >> angularVelocity >> indSpecies >> indSpecies;
            Quaternion q;
            q.setEuler(angle);
            p->setPosition(position);
            p->setVelocity(velocity);
            p->setOrientation(q);
            p->setAngularVelocity(angularVelocity);
            p->setRadius(radius);
            if (readSpeciesFromDataFile_) p->setSpecies(speciesHandler.getObject(indSpecies));
            else if (i < nHistory) p->setSpecies(species[i]);
            if (i < nHistory && fix[i]) p->fixParticle();
            if (i < nHistory) p->setPreviousPosition(positionHistory[i]);
            else p->setPreviousPosition(p->getPosition());
            particleHandler.copyAndAddObject(*p);
        } //end for all particles
    } //end if format

    particleHandler.computeAllMasses();
    return true;
}

void DPMBase::readNextFStatFile()
{
    fStatFile.open(std::fstream::in);
    std::string line;
    std::fstream& in = fStatFile.getFstream();
    // read the first three lines
    getline(in, line);
    getline(in, line);
    getline(in, line);
    Mdouble time;
    unsigned int indexP;
    int indexI; //could be negative
    interactionHandler.clear();
    while ((in.peek() != -1) && (in.peek() != '#'))
    {
        /* # 1: time
         # 2: particle Number i
         # 3: contact partner j (particles >= 0, walls < 0)
         # 4: x-position \
		 # 5: y-position  > of the contact point (I hope)
         # 6: z-position /
         # 7: delta = overlap at the contact
         # 8: ctheta = length of the tangential spring
         # 9: P1_P2_normal force |f^n|
         # 10: remaining (tangential) force |f^t|=|f-f^n|
         # 11-13: P1_P2_normal unit vector nx, ny, nz
         # 14-16: tangential unit vector tx, ty, tz
         */
        in >> time >> indexP >> indexI;
        BaseParticle* P = particleHandler.getObject(indexP);
        BaseInteraction* C;
        if (indexI >= 0)
        {
            BaseParticle* I = particleHandler.getObject(static_cast<const unsigned int>(indexI));
            //read only one of the two fstat lines reported
            if (indexI >= indexP or I->isFixed())
            {
                // particle pair contact
                C = interactionHandler.addInteraction(P, I, getNumberOfTimeSteps() + 1);
                C->setFStatData(in, P, I);
                // skip next line
                //in.ignore(256, '\n');
            }
        }
        else
        {
            // wall-particle contact
            while (wallHandler.getNumberOfObjects() <= -indexI - 1)
            {
                wallHandler.copyAndAddObject(InfiniteWall(speciesHandler.getLastObject()));
                logger(WARN, "Added new wall because .fstat file indicates contact with wall % that doesn't exist",
                       -indexI - 1);
            }
            BaseWall* I = wallHandler.getObject(static_cast<const unsigned int>(-indexI - 1));
            C = interactionHandler.addInteraction(P, I, getNumberOfTimeSteps() + 1);
            C->setFStatData(in, P, I);
        }
        //skip the rest of the line, to allow additional output and to ignore the second time a particle-particle contact is printed
        in.ignore(256, '\n');
    }
    //logger(INFO,"read % contacts at t = % (N=%)",counter,timeStamp,interactionHandler.getNumberOfObjects());
    //interactionHandler.write(std::cout);
    //logger(INFO,"normal % %",interactionHandler.getObject(0)->getNormal(),interactionHandler.getObject(0)->getContactPoint());
    //logger(INFO,"normal % %",interactionHandler.getLastObject()->getNormal(),interactionHandler.getLastObject()->getContactPoint());
}

/*!
 * \details Calls the \ref write() function in order to output all relevant data
 * (particle positions and velocities, system dimensions, positions of walls and boundaries...)
 * to a restart file.
 * \n \n
 * See also \ref readRestartFile
 */
void DPMBase::writeRestartFile()
{
    if (restartFile.openWriteNoAppend(getNumberOfTimeSteps()))
    {
        //logger(DEBUG, "Writing restart file %th time step",getNumberOfTimeSteps());
        write(restartFile.getFstream());
        restartFile.close();
    }
}

void DPMBase::writeDataFile()
{
    if (dataFile.openWrite(getNumberOfTimeSteps()))
    {
        outputXBallsData(dataFile.getFstream());
        dataFile.close();
    }
}

void DPMBase::writeEneFile()
{
    if (eneFile.openWrite(getNumberOfTimeSteps()))
    {
        //If the file type is "multiple files, writes a header for each individual files. If not, only writes for the first time step
        writeEneTimeStep(eneFile.getFstream());
        eneFile.close();
    }
}

void DPMBase::writeFStatFile()
{
    if (fStatFile.openWrite(getNumberOfTimeSteps()))
    {
        writeFstatHeader(fStatFile.getFstream());
        //fStatFile.getFstream().ignore(2000,'\t');
        fStatFile.close();
    }
}

/**
 * Inserts particles in the whole domain.
 * THis is useful if you want to check whether the wall visualisation or wall computation is correct:
 * First insert the walls, then the particles, then check in paraview if the walls and particles overlap
 * @param N
 */
void DPMBase::fillDomainWithParticles(unsigned N) {
    logger.assert_always(speciesHandler.getSize()>0,"There needs to be at least one species");
    ParticleSpecies* s = speciesHandler.getLastObject();
    SphericalParticle p(s);
    CubeInsertionBoundary b;
    b.set(p, 100, getMin(), getMax(), {0, 0, 0}, {0, 0, 0});
    b.insertParticles(this);
    logger(INFO,"Inserted % particles",particleHandler.getSize());
}


/*!
 * \details Opens a file input stream corresponding to the restart file to be opened. If the file can be successfully opened,
 * uses the \ref read() function to extract all data from the restart file.
 * Once the read-in has been successfully completed, sets the "restarted_" flag to true such that, where necessary, other functions will
 * know that the file has been restarted.
 *
 * \return int
 */
bool DPMBase::readRestartFile(ReadOptions opt)
{
    //Assuming a filename corresponding to "restartFile" has already been established,
    //opens an input filestream to the relevant file
    if (restartFile.open(std::fstream::in))
    {
        //reads the input stream line-by-line
        read(restartFile.getFstream(), opt);
        logger(INFO, "Loaded restart file %", restartFile.getFullName());
        restartFile.close();
        //sets the flag such that other functions can tell that the file has been restarted
        //e.g. does not run "setUpInitialConditions" or add headers to the .ene files etc.
        setRestarted(true);
        return true;
    }
    else /* if the file could not be opened */
    {
        logger(INFO, "% could not be loaded.", restartFile.getFullName());
        return false;
    }
}

/*!
 * \details Reads in the name of a (.restart) file and then opens and reads in the data corresponding to this file using the
 * argument-less \ref readRestartFile() function.
 * Note that this function should be called before setupInitialConditions().
 * \param[in] fileName The name of the (.restart) file to be read in.
 * \return int
 */
int DPMBase::readRestartFile(std::string fileName, ReadOptions opt)
{
    //add ".restart" if necessary
    if (fileName.find(".restart") == std::string::npos)
    {
        fileName.append(".restart");
    }

    //If padded or numbered files are used, we need to extract the counter and remove it from the filename
    //First find the last point in the restart name
    unsigned int pos = fileName.find('.');
    while (fileName.find('.', pos + 1) != std::string::npos)
    {
        pos = fileName.find('.', pos + 1);
    }
    //If the next char after the last . is a digit we are using numbered files
    std::string counter;
    if (isdigit(fileName[pos + 1]))
    {
        for (int i = pos + 1; i < fileName.length(); i++)
        {
            counter.push_back(fileName[i]);
        }
        //Set counter in restart file
        restartFile.setCounter(std::stoi(counter));
        logger(INFO, "Counter: %", std::stoi(counter));
    }

#ifdef MERCURYDPM_USE_MPI
    //Correct for the processor number
    if (NUMBER_OF_PROCESSORS > 1 && !helpers::fileExists(fileName))
    {
        //Modify file name
        const unsigned int length = fileName.length();
        if (isdigit(fileName[pos + 1]))
        {
            for (int i = pos + 1; i < length + 1; i++)
            {
                fileName.pop_back();
            }
        }
        fileName.append(std::to_string(PROCESSOR_ID));
        if (counter.size() > 0)
        {
            fileName.append(".");
            fileName.append(counter);
        }
    }
#endif

    restartFile.setName(fileName);

    logger(INFO, "Restarting from %", fileName);
    return readRestartFile(opt);
}

/*!
 * \details
 * Firstly, checks the types of particles involved in order to ensure that only
 * viable interactions are counted.
 *
 * Secondly, if the particle combination is viable, checks if the particles are interacting.
 *
 * Finally, if the particles are found to be interacting, calculates the relevant forces (as well as torques, if the "rotation"
 * flag is turned "on") acting
 * between the particles, and applies them to each particle.
 * \param[in] P1
 * \param[in] P2
 */
void DPMBase::computeInternalForce(BaseParticle* const P1, BaseParticle* const P2)
{
    //Does not compute forces if particles are fixed
    //this is necessary because the rough bottom allows overlapping fixed particles
    if (P1->isFixed() && P2->isFixed())
    {
        return;
    }
//Ensures that interactions between the "ghost" particles used to implement periodic behaviour
    //are not included in calculations
    //i.e. ends the function if both particles are "ghosts".
    if ((P1->getPeriodicFromParticle() != nullptr) && (P2->getPeriodicFromParticle() != nullptr))
    {
        return;
    }
//if statement below ensures that the PI has the lower id than PJ
    BaseParticle* PI, * PJ;
    if (P1->getId() > P2->getId())
    {
        PI = P2;
        PJ = P1;
    }
    else
    {
        PI = P1;
        PJ = P2;
    }
    //checks if the two particles are interacting
    //("getInteractionWith" returns the relevant pointer if PI and PJ are interacting,
    //zero if not)
    //if statement above ensures that the PI has the lower id than PJ
    BaseInteraction* i = PJ->getInteractionWith(PI, getNumberOfTimeSteps() + 1,
                                                &interactionHandler);
    if (i!= nullptr) {
        //calculates the force corresponding to the interaction
        i->computeForce();

        //Applies the relevant calculated forces to PI and PJ
        PI->addForce(i->getForce());
        PJ->addForce(-i->getForce());

        //checks if particle rotation is turned on...
        if (getRotation()) {
            //...and, if so, performs equivalent calculations for the torque as were
            //performed for the force.
            PI->addTorque(i->getTorque() - Vec3D::cross(PI->getPosition() - i->getContactPoint(), i->getForce()));
            PJ->addTorque(-i->getTorque() + Vec3D::cross(PJ->getPosition() - i->getContactPoint(), i->getForce()));
        }
    }
}

/*!
 * By default we have gravity and background drag as the external forces
 * \param[in] CI The BaseParticle object to which the relevant external forces are applied.
 */
void DPMBase::computeExternalForces(BaseParticle* CI)
{
    //Checks that the current particle is not "fixed"
    //and hence infinitely massive!
    if ((!CI->isFixed()) && (!CI->isPebble()))
    {
        // Applying the force due to gravity (F = m.g) Do not apply gravity to slave particles
        CI->addForce(getGravity() * CI->getMass()-getBackgroundDrag()*CI->getVelocity());
        // Still calls this in compute External Forces.
        // computeForcesDueToWalls(CI);
    }
}

/*!
 * Checks if a particle pI is currently in contact - i.e. interacting - with any of the
 * walls within the system using the \ref BaseParticle::getInteractionWith() function.
 * If an interaction <B>is</B> detected, computes the force acting between particle and wall and applies the relevant
 * torques and forces to both particle and wall(s).
 * \param[in] pI The BaseParticle object to which the wall forces are applied.
 */
void DPMBase::computeForcesDueToWalls(BaseParticle* pI, BaseWall* w)
{

    //No need to compute interactions between periodic particle images and walls
    if (pI->getPeriodicFromParticle() != nullptr)
        return;

    //Checks if the particle is interacting with the current wall
    BaseInteraction* i = w->getInteractionWith(pI, getNumberOfTimeSteps() + 1,
                                               &interactionHandler);
    if (i!=nullptr) {
        //...calculates the forces between the two objects...
        i->computeForce();

        //...and applies them to each of the two objects (wall and particle).
        pI->addForce(i->getForce());
        w->addForce(-i->getForce());

        //If the rotation flag is on, also applies the relevant torques
        //(getRotation() returns a boolean).
        if (getRotation()) // getRotation() returns a boolean.
        {
            pI->addTorque(i->getTorque() - Vec3D::cross(pI->getPosition() - i->getContactPoint(), i->getForce()));
            ///\todo TW: I think this torque has the wrong sign
            w->addTorque(-i->getTorque() + Vec3D::cross(w->getPosition() - i->getContactPoint(), i->getForce()));
        }
    }
}

/*!
 * \details
 * Performs integration - i.e. updating particle's positions, velocities and accelerations - for all
 * particles and walls within the system (i.e. in the particleHandler and wallHandler). Integration is performed
 * using the \ref BaseParticle::integrateBeforeForceComputation() function.
 * \n \n
 * The velocity Verlet algorithm requires us to integrate twice each time step: both before and after
 * the force computation. This method is therefore used in conjunction with \ref
 * DPMBase::integrateAfterForceComputation().
 * See http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet for details.
 */
void DPMBase::integrateBeforeForceComputation()
{
    //cycling through all particles, p, in the particleHandler
    //for_each(particleHandler.begin(), particleHandler.end(), [this](BaseParticle* p)
    //for (BaseParticle* p : particleHandler) {

    #pragma omp parallel for num_threads(getNumberOfOMPThreads()) //schedule(dynamic)
    for (int k = 0; k < particleHandler.getSize(); ++k) {
        BaseParticle *p = particleHandler.getObject(k);
#ifdef MERCURYDPM_USE_MPI
        //MPI particles are not integrated, they are purely ghost particles and get their new velocity and position from an MPI update
        if (!(p->isMPIParticle() || p->isPeriodicGhostParticle()))
        {
            p->integrateBeforeForceComputation(getTime(), getTimeStep());
        }
#else
        //using the particle p's internal "integrateBeforeForceComputation" function
        //to update the relevant parameters concerning the particle's position and motion
        p->integrateBeforeForceComputation(getTime(), getTimeStep());
#endif
    }
    //});
    //cycling through all walls, w, in the wallHandler
    //for_each(wallHandler.begin(), wallHandler.end(), [this](BaseWall* w)
    //for (BaseWall* w : wallHandler) {
    #pragma omp parallel for num_threads(getNumberOfOMPThreads()) //schedule(dynamic)
    for (int k = 0; k < wallHandler.getSize(); k++) {
        BaseWall *w = wallHandler.getObject(k);
        //using the wall's internal "integrateBeforeForceComputation" function
        //to update the relevant parameters concerning its position and motion
        w->integrateBeforeForceComputation(getTime(), getTimeStep());
    }
    //});
}

/*!
 * \details For each boundary, checks whether each particle in the system has "passed"
 * it and performs an action according to the type of boundary involved.
 *
 * For instance, if the boundary is a periodic boundary, the periodic boundary version of "checkBoundaryAfterParticleMoved"
 * will be called ( \ref PeriodicBoundary::checkBoundaryAfterParticleMoved()) and in turn apply the \ref shiftPosition() function
 * to the particle. If the boundary is a deletion boundary ( \ref DeletionBoundary::checkBoundaryAfterParticleMoved	()), any particle
 * passing the boundary will be deleted. Further details can be seen in the in-code comments below.
 */
void DPMBase::checkInteractionWithBoundaries()
{

    //Cycling over all boundaries within the system...
    for (BaseBoundary* b : boundaryHandler)
    {
        //check all boundaries...
        b->checkBoundaryAfterParticlesMove(particleHandler);


#ifdef MERCURYDPM_USE_MPI
        //When ghost particles are deleted by deletion boundaries they need to be removed
        //from their communication lists to avoid segfaults
        if (NUMBER_OF_PROCESSORS > 1)
        {
            //Flush deleted particles from mpi communication zones
            getCurrentDomain()->flushParticles(b->getParticlesToBeDeleted());
            getCurrentDomain()->cleanCommunicationLists();
            periodicBoundaryHandler.flushParticles(b->getParticlesToBeDeleted());
            periodicBoundaryHandler.cleanCommunicationLists();
        }

        //Delete particles that were in communication zone
        for (auto p_it = b->getParticlesToBeDeleted().begin(); p_it != b->getParticlesToBeDeleted().end(); p_it++)
        {
            particleHandler.removeGhostObject((*p_it)->getIndex());
        }
#endif
    }
}

/*!
 * \details
 * Performs integration - i.e. updating particle's positions, velocities and accelerations - for all
 * particles and walls within the system (i.e. in the particleHandler and wallHandler). Integration is performed
 * using the \ref BaseParticle::integrateBeforeForceComputation() function.
 * \n \n
 * The velocity Verlet algorithm requires us to integrate twice each time step: both before and after
 * the force computation. This method is therefore used in conjunction with \ref
 * DPMBase::integrateAfterForceComputation().
 * See http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet for details.
 */
void DPMBase::integrateAfterForceComputation()
{
    //cycling through all particles, p, in the particleHandler
    //for_each(particleHandler.begin(), particleHandler.end(), [this](BaseParticle* p){
    #pragma omp parallel for num_threads(getNumberOfOMPThreads()) //schedule(dynamic)
    for (int k = 0; k < particleHandler.getSize(); ++k) {
        BaseParticle *p = particleHandler.getObject(k);
#ifdef MERCURYDPM_USE_MPI
        //MPI particles do not require integration - they are updated by the communication step
        if (!(p->isMPIParticle() || p->isPeriodicGhostParticle()))
        {
            p->integrateAfterForceComputation(getTime(), getTimeStep());
        }
#else
        //using the particle p's internal "integrateAfterForceComputation" function
        //to update the relevant parameters concerning the particle's position and motion
        p->integrateAfterForceComputation(getTime(), getTimeStep());
#endif
    }
    //});
    //cycling through all walls, w, in the wallHandler
    //for_each(wallHandler.begin(), wallHandler.end(), [this](BaseWall* w){
    #pragma omp parallel for num_threads(getNumberOfOMPThreads()) //schedule(dynamic)
    for (int k = 0; k < wallHandler.getSize(); k++) {
        BaseWall *w = wallHandler.getObject(k);
        //using the wall's internal "integrateAfterForceComputation" function
        //to update the relevant parameters concerning its position and motion
        w->integrateAfterForceComputation(getTime(), getTimeStep());
    }
    //});
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////statisticsFromRestartData
///////////////////////////////////////////////////////////////////////////////////////////////////////
//void DPMBase::statisticsFromRestartData(const char *name)
//{
//    /// \todo {Check this whole function}
//    //This function loads all MD data
//    readRestartFile();
//
//    //This creates the file statistics will be saved to
//    std::stringstream ss("");
//    ss << name << ".stat";
//    statFile.setName(ss.str());
//    statFile.setOpenMode(std::fstream::out);
//    statFile.open();
//
//    // Sets up the initial conditions for the simulation
//    // setupInitialConditions();
//    // Setup the previous position arrays and mass of each particle.
//    computeParticleMasses();
//    // Other routines required to jump-start the simulation
//    actionsBeforeTimeLoop();
//    initialiseStatistics();
//    hGridActionsBeforeTimeLoop();
//    writeEneHeader(eneFile.getFstream());
//
//    while (readDataFile(dataFile.getName().c_str()))
//    {
//        hGridActionsBeforeTimeLoop();
//        actionsBeforeTimeStep();
//        checkAndDuplicatePeriodicParticles();
//        hGridActionsBeforeTimeStep();
////        dataFile.setSaveCurrentTimeStep(true);
////        eneFile.setSaveCurrentTimeStep(true);
////        statFile.setSaveCurrentTimeStep(true);
////        fStatFile.setSaveCurrentTimeStep(true);
//        computeAllForces();
//        removeDuplicatePeriodicParticles();
//        actionsAfterTimeStep();
//        writeEneTimeStep(eneFile.getFstream());
//        std::cout << std::setprecision(6) << getTime() << std::endl;
//    }
//
//    dataFile.close();
//    statFile.close();
//}
/*!
 * Initially, resets all forces to zero for all particles and all walls.
 * For each particle in turn, the function searches for particle interactions, and computes
 * the relevant internal forces, followed by the relevant external forces (e.g. gravity).
 */
void DPMBase::computeAllForces()
{
    //Resetting all forces on both particles and walls to zero
    #pragma omp parallel num_threads(getNumberOfOMPThreads())
    {
        #pragma omp for
        for (int k = 0; k < particleHandler.getSize(); ++k) {
            particleHandler.getObject(k)->resetForceTorque(getNumberOfOMPThreads());
        }
        #pragma omp for
        for (int k = 0; k < wallHandler.getSize(); k++) {
            wallHandler.getObject(k)->resetForceTorque(getNumberOfOMPThreads());
        }
    }
    logger(DEBUG,"All forces set to zero");

    // for omp simulations, reset the newObjects_ variable (used for reduction)
    interactionHandler.resetNewObjectsOMP();

    // compute all internal and external forces; for omp simulations, this can be done in parallel
    #pragma omp parallel num_threads(getNumberOfOMPThreads())
    {
        //logger(INFO, "Number of omp threads = %", getNumberOfOMPThreads());
        ///Now loop over all particles contacts computing force contributions
        #pragma omp for schedule(dynamic)
        for (int k = 0; k < particleHandler.getSize(); ++k) {
            BaseParticle *p = particleHandler.getObject(k);
            //computing both internal forces (e.g. due to collisions)
            //and external forces (e.g. gravity)
            //(compute internal forces compares the current particle p
            //with all others in the handler!)
            computeInternalForces(p);
            // body forces
            computeExternalForces(p);
        }

        // wall-forces
        #pragma omp for schedule(dynamic)
        for (int k = 0; k < wallHandler.getSize(); k++) {
            BaseWall *w = wallHandler.getObject(k);
            computeWallForces(w);
            w->computeWear();
        }

    }

    #ifdef MERCURYDPM_TRIANGLE_WALL_CORRECTION
    // This statement deals with interactions of a particle with a group of walls:
    // If there are multiple contact forces between the particle and the group of walls, we assume that these contacts share a concave edge or a vertex. Thus, we multiply each force by a weight $w_{\alpha\gamma}=|\vec{f}_{\alpha\gamma}^{\text{w}}|/\sum_\gamma|\vec{f}_{\alpha\gamma}^{\text{w}}|$. Note that these weights add up to 1, thus, the total contact force between the particle and the group of walls will be a weighted average of the individual contact forces.
    // A second modification is done in #BaseWall::getInteractionWith
    ///\todo add OMP here
    for (BaseParticle *p : particleHandler) {
        // find all current interactions between this particle and a TriangleWall
        std::vector<BaseInteraction *> wallInters;
        double overlapSum = 0.0;
        for (const auto i : p->getInteractions()) {
            if (i->getI()->getName() == "TriangleWall" && i->getTimeStamp() > getNumberOfTimeSteps()) {
                wallInters.push_back(i);
                overlapSum += i->getOverlap();
            }
        }
        // if more than one interactions are detected, we have a concave wall interaction problem.
        if (wallInters.size() > 1) {
            /// Stefan's weighting approach
            // get the multiple contact points to define the new contact plane, normal, and overlap
            for (BaseInteraction *i : wallInters) {
                auto q = i->getI();
                // undo the force/torque computation from multiple TriangleWalls
                p->addForce(-i->getForce());
                q->addForce(i->getForce());
                if (getRotation()) {
                    p->addTorque(-i->getTorque() + Vec3D::cross(p->getPosition() - i->getContactPoint(), i->getForce()));
                    q->addTorque( i->getTorque() - Vec3D::cross(q->getPosition() - i->getContactPoint(), i->getForce()));
                }
                // reset the forces
                double w = i->getOverlap() / overlapSum;
                i->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * i->getOverlap()) * w * i->getNormal());
                i->setForce(i->getForce() * w);
                i->setTorque(i->getTorque() * w);
                // applies the reset forces to the particle and the wall
                p->addForce(i->getForce());
                q->addForce(-i->getForce());
                if (getRotation()) {
                    p->addTorque( i->getTorque() - Vec3D::cross(p->getPosition() - i->getContactPoint(), i->getForce()));
                    q->addTorque(-i->getTorque() + Vec3D::cross(q->getPosition() - i->getContactPoint(), i->getForce()));
                }
            }
        }
    }
    #endif

#ifdef CONTACT_LIST_HGRID
    PossibleContact* Curr=possibleContactList.getFirstPossibleContact();
    while(Curr)
    {
        computeInternalForces(Curr->getP1(),Curr->getP2());
        Curr=Curr->getNext();
    }
#endif

    // Check wall forces
    #pragma omp for schedule(dynamic)
    for (int k = 0; k < wallHandler.getNumberOfObjects(); k++) {
        BaseWall *w = wallHandler.getObject(k);
        w->checkInteractions(&interactionHandler, getNumberOfTimeSteps() + 1);
    }

    computeAdditionalForces();

    // for omp simulations, sum up all forces and add all newObjects_ (needed since both are using reduction)
    #ifdef MERCURYDPM_USE_OMP
    if (getNumberOfOMPThreads()>1) {
        interactionHandler.addNewObjectsOMP();
    }
    //Resetting all forces on both particles and walls to zero
    #pragma omp parallel num_threads(getNumberOfOMPThreads())
    {
        #pragma omp for
        for (int k = 0; k < particleHandler.getSize(); k++) {
            particleHandler.getObject(k)->sumForceTorqueOMP();
        }
        #pragma omp for
        for (int k = 0; k < wallHandler.getSize(); k++) {
            wallHandler.getObject(k)->sumForceTorqueOMP();
        } //end reset forces loop
    }
    #endif

    //end outer loop over contacts.
}

/*!
 * \details Taking a single BaseParticle object as an argument, passes it to the \ref broadPhase() function
 * which then loops over all other particles in the particleHandler and computes the relevant forces
 * for any particle pairing found to be in contact.
 * \param[in] i A BaseParticle object for which we want to calculate the internal forces.
 */
void DPMBase::computeInternalForces(BaseParticle* i)
{
    for (auto it = particleHandler.begin(); (*it) != i; ++it)
    {
        computeInternalForce(i, *it);
    }
}


/*!
 * \details Writes out all relevant information - e.g. system dimensions, run duration, particle information -
 * for a .restart file to the chosen output stream, <TT>os</TT>. More detailed comments may be seen in the body of the code.
 * \param[in] writeAllParticles A boolean which decides whether or not all particle information is written
 * to the output file. (Otherwise, only a small number of particles are printed.)
 * \param[in] nToWrite
 */
void DPMBase::write(std::ostream& os, bool writeAllParticles) const
{
    os << "MercuryDPM " << getVersion();
    //which outputs basic information regarding the various files  (.data, .fstat etc. etc.)
    //only writes the run number if it is different from 0
    if (runNumber_ != 0)
        os << " runNumber " << runNumber_;
    os << " name " << name_;
    os << " revision " << getRevision();
    os << " repository " << getRepositoryURL() << '\n';
    os << "dataFile    " << dataFile << '\n';
    os << "fStatFile   " << fStatFile << '\n';
    os << "eneFile     " << eneFile << '\n';
    os << "restartFile " << restartFile << '\n';
    os << "statFile    " << statFile << '\n';
    os << "interactionFile " << interactionFile << '\n';
    //Outputs the "domain" corresponding to the system for
    //use with XBalls, as well as other information regarding the system as a whole
    os << "xMin " << getXMin()
       << " xMax " << getXMax()
       << " yMin " << getYMin()
       << " yMax " << getYMax()
       << " zMin " << getZMin()
       << " zMax " << getZMax() << '\n'
       << "timeStep " << getTimeStep()
       << " time " << getTime()
       << " ntimeSteps " << numberOfTimeSteps_
       << " timeMax " << getTimeMax() << '\n'
       << "systemDimensions " << getSystemDimensions()
       << " particleDimensions " << getParticleDimensions()
       << " gravity " << getGravity()
       << " backgroundDrag " <<getBackgroundDrag();
    os << " writeVTK " << writeParticlesVTK_
        << " " << wallHandler.getWriteVTK()
        << " " << interactionHandler.getWriteVTK()
        << " " << boundaryHandler.getWriteVTK()
        << " " << (vtkWriter_?vtkWriter_->getFileCounter():0)
        << " " << wallVTKWriter_.getFileCounter()
        << " " << interactionVTKWriter_.getFileCounter()
        << " " << boundaryVTKWriter_.getFileCounter();

    // Outputting the wall details vtk enum values and file types.
    // For example, we have three wall details options with enum values 4, 12 and 25, and three corresponding filetypes
    // NO_FILE, ONE_FILE and MULTIPLE_FILES. Let's also say the file counter is at 101.
    // The output will be: writeWallDetailsVTK 3 4 NO_FILE 12 ONE_FILE 25 MULTIPLE_FILES 101
    // Only write when any of them have something other than NO_FILE. In this way the selftests restart files all keep working.
    if (wallHandler.getWriteDetailsVTKAny())
    {
        std::unordered_map<WallHandler::DetailsVTKOptions, FileType> writeWallDetailsVTK = wallHandler.getWriteWallDetailsVTKAll();
        os << " writeWallDetailsVTK " << writeWallDetailsVTK.size();
        for (const auto& p : writeWallDetailsVTK)
            os << " " << static_cast<int>(p.first) << " " << p.second;
        os << " " << wallDetailsVTKWriter_.getFileCounter();
    }

    os << " random ";
    random.write(os);
#ifdef MERCURYDPM_USE_OMP
    //Write number of OMP threads
    if(getNumberOfOMPThreads() > 1) {
        os << " numberOfOMPThreads " << getNumberOfOMPThreads();
    }
#endif
#ifdef MERCURYDPM_USE_MPI
    //Check if we are dealing with multiple cores
    if (NUMBER_OF_PROCESSORS > 1 )
    {
        os << " numberOfProcessors " << NUMBER_OF_PROCESSORS
            << " numberOfDomains " << numberOfDomains_[Direction::XAXIS] << " " << numberOfDomains_[Direction::YAXIS] << " " << numberOfDomains_[Direction::ZAXIS];
    }
#endif
    //only write xBallsArguments if they are nonzero
    if (getXBallsAdditionalArguments().compare(""))
        os << " xBallsArguments " << getXBallsAdditionalArguments();
    os << '\n';
    //writes all species (including mixed species) to an output stream

    speciesHandler.write(os);

    //outputs the number of walls in the system
    os << "Walls " << wallHandler.getNumberOfObjects() << std::endl;
    if (writeAllParticles || wallHandler.getSize() < 9) {
        for (BaseWall* w : wallHandler)
            os << (*w) << std::endl;
    } else {
        for (int i=0; i<2; ++i)
            os << *wallHandler.getObject(i) << std::endl;
        os << "...\n";
    }

    //outputs the number of boundaries in the system
    os << "Boundaries " << boundaryHandler.getNumberOfObjects() << std::endl;
    if (writeAllParticles || boundaryHandler.getSize() < 9)
    {
        for (BaseBoundary* b: boundaryHandler)
            os << (*b) << std::endl;
    }
    else
    {
        for (int i = 0; i < 2; ++i)
            os << *boundaryHandler.getObject(i) << std::endl;
        os << "...\n";
    }

    if (writeAllParticles || particleHandler.getSize() < getNToWrite())
    {
        //if the "writeAllParticles" bool == true, or there are fewer than 4 particles
        //calls the particleHandler version of the "write" function and also
        //outputs to file all relevant particle information for all particles in the system
        particleHandler.write(os);
    }
    else
    {
        //otherwise, only prints out limited information
        os << "Particles " << particleHandler.getSize() << '\n';
        for (unsigned int i = 0; i < getNToWrite(); i++)
            os << *(particleHandler.getObject(i)) << '\n';
        os << "..." << '\n';
    }
    // Similarly, print out interaction details (all of them, or up to nToWrite of them)
    if (writeAllParticles || interactionHandler.getNumberOfObjects() < getNToWrite())
    {
        interactionHandler.write(os);
    }
    else
    {
        os << "Interactions " << interactionHandler.getNumberOfObjects() << '\n';
        for (unsigned int i = 0; i < getNToWrite(); i++)
            os << *(interactionHandler.getObject(i)) << '\n';
        os << "..." << '\n';
    }
}

/*!
 * Reads in an <B>existing</B> .restart file line-by-line and passes all relevant parameters to the current instance of DPMBase.
 * The data stream corresponding to the desired input file is passed as an argument.
 * \param[in] is The data stream from which the particle data will be read.
 */
void DPMBase::read(std::istream& is, ReadOptions opt)
{
#ifdef MERCURYDPM_USE_MPI
    int previousNumberOfProcessors;
#endif
    //Declares...
    std::string dummy;
    //...and reads in a dummy variable from the start of the stream "is"
    is >> dummy;
    //compare the string read in to the phrase "restart_version" to see if the stream corresponds
    //to a restart file (all restart files begin with this phrase)
    //if both strings match, strcmp(dummy.c_str(), "restart_version") returns 0 (here read as "false")
    if (dummy != "restart_version" && dummy != "MercuryDPM")
    {
        //If the strings do not match, if statement is fulfilled  and the error logged
        //Note: only very old files did not have a restart_version
        logger(FATAL, "Error in DPMBase::read(is): this is not a valid restart file");
    }
    else
    {
        //reads in the restart version (earlier versions of Mercury possess different file formats!)
        is >> restartVersion_;
        //checking which version the current data file corresponds to, and reads the data in
        //accordingly
        if (restartVersion_ == "1.0" || restartVersion_ == "0.14")
        {
            //reads in and saves the relevant values from the data file to the current instance of DPMBase
            std::stringstream line;

            // Store path (if restart file is nonlocal)
            auto slash = restartFile.getName().rfind('/');
            std::string path;
            if (slash != std::string::npos)
            {
                path = restartFile.getName().substr(0, slash + 1);
            }
            if (!path.empty())
            {
                logger(INFO, "Adding path information (%) to file names", path);
            }

            //line 1
            helpers::getLineFromStringStream(is, line);
            //discards the whitespace (ws) at the start of the stream
            line >> std::ws;
            //uses the "peek" function to access the stream's first
            //non-whitespace character, and check if it is an "r"
            if (line.peek() == 'r')
                //if so, reads in the current run number
                line >> dummy >> runNumber_;
            //In either case, then uses the "Files" version of the read function
            //to read in the rest of the relevant information.
            line >> dummy >> name_;
            setName(name_);

            //Read line 2-7 (definition of i/o files)
            helpers::getLineFromStringStream(is, line);
            line >> dummy >> dataFile;
            helpers::getLineFromStringStream(is, line);
            line >> dummy >> fStatFile;
            helpers::getLineFromStringStream(is, line);
            line >> dummy >> eneFile;
            helpers::getLineFromStringStream(is, line);
            line >> dummy >> restartFile;
            helpers::getLineFromStringStream(is, line);
            line >> dummy >> statFile;

            // Add the file path from the restart file to the file names
            dataFile.setName(path + dataFile.getName());
            fStatFile.setName(path + fStatFile.getName());
            eneFile.setName(path + eneFile.getName());
            restartFile.setName(path + restartFile.getName());
            statFile.setName(path + statFile.getName());

            // Get current position
            //check if the next line starts with 'interactionFile'; otherwise, skip interaction
            if (helpers::compare(is, "interactionFile"))
            {
                helpers::getLineFromStringStream(is, line);
                line >> interactionFile;
                interactionFile.setName(path + interactionFile.getName());
            }

            helpers::getLineFromStringStream(is, line);
            line >> dummy >> min_.x()  /// \todo : Bound checking
                 >> dummy >> max_.x()  /// \todo : Same order as other file format, please?
                 >> dummy >> min_.y()
                 >> dummy >> max_.y()
                 >> dummy >> min_.z()
                 >> dummy >> max_.z();

            helpers::getLineFromStringStream(is, line);
            line >> dummy >> timeStep_
                 >> dummy >> time_
                 >> dummy >> numberOfTimeSteps_
                 >> dummy >> timeMax_;

            helpers::getLineFromStringStream(is, line);
            line >> dummy >> systemDimensions_
                 >> dummy >> particleDimensions_
                 >> dummy >> gravity_
                 >> dummy >> backgroundDrag_;
            line >> dummy;
            if (!dummy.compare("writeVTK"))
            {
                FileType writeInteractionsVTK = FileType::NO_FILE, writeWallsVTK;
                unsigned particlesCounter, wallCounter, interactionCounter, boundaryCounter;
                bool writeBoundaryVTK, writeParticlesVTK;
                line >> writeParticlesVTK >> writeWallsVTK >> writeInteractionsVTK >> writeBoundaryVTK >> particlesCounter >> wallCounter >> interactionCounter >> boundaryCounter;
                line.clear();//because the number of arguments  in writeVTK has changed
                line >> dummy;
                setParticlesWriteVTK(writeParticlesVTK);
                wallHandler.setWriteVTK(writeWallsVTK);
                interactionHandler.setWriteVTK(writeInteractionsVTK);
                boundaryHandler.setWriteVTK(writeBoundaryVTK);
                vtkWriter_->setFileCounter(particlesCounter);
                wallVTKWriter_.setFileCounter(wallCounter);
                interactionVTKWriter_.setFileCounter(interactionCounter);
                boundaryVTKWriter_.setFileCounter(boundaryCounter);
            }

            if (!dummy.compare("writeWallDetailsVTK"))
            {
                // Example input: writeWallDetailsVTK 3 4 NO_FILE 12 ONE_FILE 25 MULTIPLE_FILES 101
                // Number of options is 3, each option having a enum value and file type (4 NO_FILE), (12 ONE_FILE) and
                // (25 MULTIPLE_FILES). At the end the file counter (101).
                unsigned numberOfOptions, fileCounter, enumType;
                FileType fileType;
                line >> numberOfOptions;
                for (unsigned  i = 0; i < numberOfOptions; i++)
                {
                    line >> enumType >> fileType;
                    wallHandler.setWriteDetailsVTK(static_cast<WallHandler::DetailsVTKOptions>(enumType), fileType);
                }
                line >> fileCounter;
                wallDetailsVTKWriter_.setFileCounter(fileCounter);
                line.clear();
                line >> dummy;
            }

            if (!dummy.compare("random"))
            {
                random.read(line);
                line >> dummy;
            }

#ifdef MERCURYDPM_USE_OMP
            //Read the number of OMP threads
            if (!dummy.compare("numberOfOMPThreads")) {
                int numberOfOMPThreads;
                line >> numberOfOMPThreads;
                setNumberOfOMPThreads(numberOfOMPThreads);
                //logger(INFO," Check the number of OMP threads = % ", getNumberOfOMPThreads());
            }
#endif
#ifdef MERCURYDPM_USE_MPI
            if (!dummy.compare("numberOfProcessors"))
            {
                line  >> previousNumberOfProcessors
                    >> dummy >> numberOfDomains_[Direction::XAXIS]
                    >> numberOfDomains_[Direction::YAXIS]
                    >> numberOfDomains_[Direction::ZAXIS];
            }
            else
            {
                logger(INFO,"Reading a serial restart file");
                //numberOfDomains_ = {1,1,1};
            }
#endif
            if (!dummy.compare("xBallsArguments")) {
                helpers::getLineFromStringStream(line, line);
                setXBallsAdditionalArguments(line.str());
            }

            speciesHandler.read(is);

#ifdef MERCURYDPM_USE_MPI
            //Initialise MPI structures and perform domain decomposition
            decompose();
#endif

            //reading in the various relevant handlers
            unsigned int N;
            is >> dummy >> N;
            if (dummy.compare("Walls"))
                logger(ERROR, "DPMBase::read(is): Error during restart: 'Walls' argument could not be found.");
            wallHandler.clear();
            wallHandler.setStorageCapacity(N);
            helpers::getLineFromStringStream(is, line);
            for (unsigned int i = 0; i < N; i++)
            {
                helpers::getLineFromStringStream(is, line);
                wallHandler.readAndAddObject(line);
            }

            is >> dummy >> N;
            boundaryHandler.clear();
            boundaryHandler.setStorageCapacity(N);
            if (dummy.compare("Boundaries"))
                logger(ERROR, "DPMBase::read(is): Error during restart: 'Boundaries' argument could not be found.");
            helpers::getLineFromStringStream(is, line);
            for (unsigned int i = 0; i < N; i++)
            {
                helpers::getLineFromStringStream(is, line);
                boundaryHandler.readAndAddObject(line);
            }

            if (opt==ReadOptions::ReadNoParticlesAndInteractions) return;

            is >> dummy >> N;
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            //display a message if a large amount o fparticles is read
            if (N>2.5e5) logger(INFO, "Reading % particles (may take a while)",N);
            logger.assert_always(dummy.compare("Particles")==0, "DPMBase::read(is): Error during restart: 'Particles' argument could not be found. %",dummy);
            particleHandler.clear();
            particleHandler.setStorageCapacity(N);
            for (unsigned int i = 0; i < N; i++)
            {
                //ParticleHandler::readAndCreateObject reads line-by-line
                particleHandler.readAndAddObject(is);
                //skip the remaining data in line
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                /// \todo {Do we want to calculate the mass?}
                //particleHandler.getLastObject()->computeMass();
            }
#ifdef MERCURYDPM_USE_MPI
            //Interaction distances of the domainHandler and periodicBoundaryHandler need to be set
            Mdouble interactionRadius = particleHandler.getLargestInteractionRadius();
            domainHandler.setInteractionDistance(2.0*interactionRadius);
            periodicBoundaryHandler.setInteractionDistance(2.0*interactionRadius);

            if (NUMBER_OF_PROCESSORS > 1)
            {
                //Create ghost particles
                domainHandler.addNewParticles();
                periodicBoundaryHandler.addNewParticles();
            }
#endif
            //Add interactions to particles and ghost particles
            if (opt==ReadOptions::ReadNoInteractions) return;
            interactionHandler.read(is);
        }
            //reading in for older versions of the Mercury restart file.
        else if (!restartVersion_.compare("3"))
        {
            logger(INFO, "DPMBase::read(is): restarting from an old restart file (restart_version %).",
                   restartVersion_);
            readOld(is);
        }
            //returning an error if there is no restart file to read in due to the use of outdated files.
        else
        {
            //only very old files did not have a restart_version
            logger(FATAL,
                   "Error in DPMBase::read(is): restart_version % cannot be read; use an older version of Mercury to upgrade the file",
                   restartVersion_);
        }
    }
}

/*!
 * \param[in] is
 */
void DPMBase::readOld(std::istream& is)
{
    std::string dummy;
    is >> dummy >> dummy;
    setName(dummy);

    unsigned int saveCountData, saveCountEne, saveCountStat, saveCountFStat;
    unsigned int fileTypeFstat, fileTypeData, fileTypeEne, fileTypeRestart;
    is >> dummy >> min_.x()
       >> dummy >> max_.x()
       >> dummy >> min_.y()
       >> dummy >> max_.y()
       >> dummy >> min_.z()
       >> dummy >> max_.z()
       >> dummy >> timeStep_
       >> dummy >> time_
       >> dummy >> timeMax_
       >> dummy >> saveCountData
       >> dummy >> saveCountEne
       >> dummy >> saveCountStat
       >> dummy >> saveCountFStat
       >> dummy >> systemDimensions_
       >> dummy >> gravity_
       >> dummy >> fileTypeFstat
       >> dummy >> fileTypeData
       >> dummy >> fileTypeEne;
    dataFile.setSaveCount(saveCountData);
    eneFile.setSaveCount(saveCountEne);
    statFile.setSaveCount(saveCountStat);
    fStatFile.setSaveCount(saveCountFStat);

    fStatFile.setFileType(static_cast<FileType>(fileTypeFstat));
    dataFile.setFileType(static_cast<FileType>(fileTypeData));
    eneFile.setFileType(static_cast<FileType>(fileTypeEne));

    //this is optional to allow restart files with and without restartFile.getFileType()
    is >> dummy;
    if (!strcmp(dummy.c_str(), "options_restart"))
    {
        is >> fileTypeRestart;
        restartFile.setFileType(static_cast<FileType>(fileTypeRestart));
    }

    speciesHandler.read(is);
    wallHandler.read(is);
    boundaryHandler.read(is);
    particleHandler.read(is);
    setRestarted(true);
}

/*!
 *\details
 * 1. Checks if at least one species exists in the SpeciesHandler.
 * 2. Checks if the time step is set or not.
 * \n \n
 * If any of the above checks fail, gives an error message to the user and terminates the program.
 */
void DPMBase::checkSettings()
{
    //check if name is set
    logger.assert_always(getName() != "",
                         "File name not set: use setName()");
    //check if time step is set
    logger.assert_always(getTimeStep() != 0,
                         "Time step undefined: use setTimeStep()");
    //check if domain is set
    logger.assert_always(getXMax() > getXMin(),
                         "Domain size not set: use setXMin() and setXMax()");
    logger.assert_always(getYMax() > getYMin(),
                         "Domain size not set: use setYMin() and setYMax()");
    logger.assert_always(systemDimensions_ == 3 ? (getZMax() > getZMin()) : (getZMax() >= getZMin()),
                         "Domain size not set: use setZMin() and setZMax()", systemDimensions_);

    //check for species parameters
    logger.assert_always(speciesHandler.getNumberOfObjects() > 0,
                         "No species defined: use speciesHandler.copyAndAddObject()");
    for (BaseParticle* p : particleHandler)
    {
        logger.assert_always(p->getSpecies() != nullptr, "particle % has no species", p->getId());
    }
    for (BaseWall* w : wallHandler)
    {
        logger.assert_always(w->getSpecies() != nullptr, "% with index % has no species", w->getName(), w->getId());
    }
}

void DPMBase::forceWriteOutputFiles()
{
    setLastSavedTimeStep(NEVER);
    writeOutputFiles();
}

/*!
 *\details Writes headers and all relevant information to the relevant output files. Note that the \ref writeFstatHeader()
 * actually contains within it the functionality to write the full fstat data, whereas for .ene files
 * the functions to write the headers and main data are separate. Note that the interaction file is not written here: it
 * is written with the start and end of each interaction.
 *
 * The function [X].saveCurrentTimeStep(numberOfTimeSteps_) returns <tt>true</tt> if:
 *
 * a) The current time step is greater than or equal to the time step at which the next write or read operation is supposed to happen.
 *
 * b) The \ref FileType is <B>not</B> "NO_FILE".
 *
 * c) The file is open.
 *
 */
void DPMBase::writeOutputFiles()
{
    //writing fstat data if .saveCurrentTimeStep(numberOfTimeSteps_) is true
    if (fStatFile.saveCurrentTimeStep(numberOfTimeSteps_))
    {
        writeFStatFile();
    }

    //writing .ene data if .saveCurrentTimeStep(numberOfTimeSteps_) is true
    if (eneFile.saveCurrentTimeStep(numberOfTimeSteps_))
    {
        writeEneFile();
    }
    //writing .data data if .saveCurrentTimeStep(numberOfTimeSteps_) is true
    if (dataFile.saveCurrentTimeStepNoFileTypeCheck(numberOfTimeSteps_))
    {
        if (dataFile.getFileType() != FileType::NO_FILE) {
            if (getRestarted() || dataFile.getCounter() == 0)
                writeXBallsScript();
            writeDataFile();
        } else {
            dataFile.setLastSavedTimeStep(numberOfTimeSteps_);
        }
        printTime();
        writeVTKFiles();
    }
    cgHandler.evaluate();


    //write restart file last, otherwise the output counters are wrong
    if (restartFile.saveCurrentTimeStep(numberOfTimeSteps_))
    {
        writeRestartFile();
    }
}

/*!
 * \details This function takes the simulation domain boundaries and decomposes it into sub domains ready for parallel computations
 */
void DPMBase::decompose()
{
#ifdef MERCURYDPM_USE_MPI

    //If running in parallel build, but just running with one core - no domain decomposition required
    int numberOfRequiredProcessors = numberOfDomains_[Direction::XAXIS]*
                                     numberOfDomains_[Direction::YAXIS]*
                                     numberOfDomains_[Direction::ZAXIS];
    if (NUMBER_OF_PROCESSORS != numberOfRequiredProcessors)
    {
        logger(ERROR,"The domain decompositions expects % processors, but only % are requested.\n"
                     "Either run your process using \"mpirun -np % [executable]\", "
                     "or change the domain decomposition to e.g. setNumberOfDomains({%,1,1}).", numberOfRequiredProcessors, NUMBER_OF_PROCESSORS, numberOfRequiredProcessors, NUMBER_OF_PROCESSORS);
    }

    if (NUMBER_OF_PROCESSORS == 1) {return;}

    //Check if the simulation domain has been set
    logger.assert_always(getXMax() - getXMin() > 0,"Please set your simulation domain (setXMax(),setXmin()) before calling solve()");
    logger.assert_always(getYMax() - getYMin() > 0,"Please set your simulation domain (setYMax(),setYmin()) before calling solve()");
    logger.assert_always(getZMax() - getZMin() > 0,"Please set your simulation domain (setZMax(),setZmin()) before calling solve()");

    //Grab simulation domains
    std::vector<Mdouble> simulationMin{getXMin(), getYMin(), getZMin()};
    std::vector<Mdouble> simulationMax{getXMax(), getYMax(), getZMax()};

    //Check if the user input decomposition is correct
    logger.assert_always(numberOfDomains_[Direction::XAXIS] > 0,"Number of domain in x-direction incorrect: %",numberOfDomains_[Direction::XAXIS]);
    logger.assert_always(numberOfDomains_[Direction::YAXIS] > 0,"Number of domain in y-direction incorrect: %",numberOfDomains_[Direction::YAXIS]);
    logger.assert_always(numberOfDomains_[Direction::ZAXIS] > 0,"Number of domain in z-direction incorrect: %",numberOfDomains_[Direction::ZAXIS]);

    //Open domain decomposition, closed is not implemented
    bool open = true;

    //Check if the number of domains is equal to the number of processors
    logger.assert_always(numberOfDomains_[Direction::XAXIS]*numberOfDomains_[Direction::YAXIS]*numberOfDomains_[Direction::ZAXIS] == NUMBER_OF_PROCESSORS,
             "Number of Processors is not equal to number of domains.    Processors %, domains, %",
             NUMBER_OF_PROCESSORS,
             numberOfDomains_[Direction::XAXIS]*numberOfDomains_[Direction::YAXIS]*numberOfDomains_[Direction::ZAXIS]);

    //Create all processor domains

    domainHandler.setDPMBase(this); /// \todo Find out why this has to be done explicit once again
    domainHandler.createMesh(simulationMin, simulationMax, numberOfDomains_, open);
    logger(VERBOSE,"Number of domains: % | Number of processors: %",domainHandler.getNumberOfObjects(), NUMBER_OF_PROCESSORS);
    //logger.assert_always(domainHandler.getNumberOfObjects() == numberOfProcessors, "Invalid decomposition: Number of domains and processors are different");

    //Tell the current processor to which domain it belongs
    for (Domain* domain : domainHandler)
    {
        if (domain->getRank() == PROCESSOR_ID)
        {
            logger(VERBOSE,"processor: %, domain index: %",PROCESSOR_ID, domain->getIndex());
            domainHandler.setCurrentDomainIndex(domain->getIndex());
        }
    }

    //Define the mpi transfer types, which requires a definition of the species already
    /// \todo update this function for Jonny
    logger.assert_always(speciesHandler.getNumberOfObjects() > 0, "Please create a particle species before calling solve()");
    MPIContainer::Instance().initialiseMercuryMPITypes(speciesHandler);

    //Make sure all processors are done with decomposition before proceeding
    logger(VERBOSE,"processor %: #real particles: %, #total particles: %", PROCESSOR_ID, particleHandler.getNumberOfRealObjects(), particleHandler.getSize());
    MPIContainer::Instance().sync();
#endif
}


/*!
 * \details  - Initialises the time, sets up the initial conditions for the simulation by
 *             calling the setupInitialConditions() and resets the counter using
 *             setNExtSavedTimeStep().
 *          -  HGrid operations which is the contact detection algorithm.
 *          -  Checks if the basic essentials are set for carrying out the
 *             simulations using checkSettings()
 *          -  And many more vital operations.
 *
 * Further details are included in the body of the code, below.
 * \todo Is it necessary to reset initial conditions here and in setTimeStepByParticle
 *       (i.e. should it be in constructor) Thomas: I agree, setTimeStepByParticle should be
 *       rewritten to work without calling setupInitialConditions
 */
void DPMBase::initialiseSolve() {
    logger(DEBUG, "Entered solve");
#ifdef CONTACT_LIST_HGRID
    logger(INFO,"Using CONTACT_LIST_HGRID");
#endif

    /// Initialise the time and
    /// sets up the initial conditions for the simulation
    /// \todo Is it necessary to reset initial conditions here and in setTimeStepByParticle (i.e. should it be in constructor)?
    ///Thomas: I agree, setTimeStepByParticle should be rewritten to work without calling setupInitialConditions
    if (!getRestarted()) {
        // If the simulation is "new" (i.e. not restarted):
        // - set time, nTimeSteps to zero
        // - reset the file counter etc.
        // - decompose the domain based on XMin, XMax, ....
        // - run user-defined setupInitialConditions
        numberOfTimeSteps_ = 0;
        setTime(0.0);
        resetFileCounter();
        decompose();
        // \todo tw there was a function combining the next two lines, why is it back to the old version?
        //setLastSavedTimeStep(NEVER); //reset the counter
        //this is to ensure that the interaction time stamps agree with the resetting of the time value
        for (auto &i: interactionHandler)
            i->setTimeStamp(0);
        setupInitialConditions();
        logger(DEBUG, "Have created the particles initial conditions");
    } else {
        // If the simulation is "restarted" (i.e. not restarted):

        // - run wall-defined actionsOnRestart
        for (auto w: wallHandler) {
            w->actionsOnRestart();
        }

        // - run user-defined actionsOnRestart
        actionsOnRestart();
    }

    // Check that the code has been correctly set up,
    // i.e. system dimensions, particles and time steps are sensibly implemented
    checkSettings();

    // If the simulation is "new" and the runNumber is used, append the run number to the problem name
    if (getRunNumber() > 0 && !getRestarted()) {
        std::stringstream name;
        name << getName() << "." << getRunNumber();
        setName(name.str());
    }

    //If append is true, files are appended, not overwritten
    if (getAppend()) {
        setOpenMode(std::fstream::out | std::fstream::app);
        //Restart files should always be overwritten.
        restartFile.setOpenMode(std::fstream::out);
    } else {
        setOpenMode(std::fstream::out);
    }

    //sets the hgrid, writes headers to the .stat output file
    initialiseStatistics();

    if (getInteractionFile().getFileType() == FileType::ONE_FILE) {
        logger(WARN, "Warning: interaction file will take up a lot of disk space!");
        getInteractionFile().open();
    }

    // Sets the mass of all particle.
    /// \todo MX: Why does the mass get computed here? if a particle is assigned a radius, it automatically also computes its mass.
    /// IFCD: commenting out this line does not make any test fail on my system.
    particleHandler.computeAllMasses();

    // Other initialisations
    //max_radius = getLargestParticle()->getRadius();

    actionsBeforeTimeLoop();
    boundaryHandler.boundaryActionsBeforeTimeLoop();
    hGridActionsBeforeTimeLoop();

    // Performs a first force computation
    checkAndDuplicatePeriodicParticles();

#ifdef MERCURYDPM_USE_MPI
    if (NUMBER_OF_PROCESSORS > 1)
    {
        //Find new mpi particles
        domainHandler.addNewParticles();
        //Periodic particles in parallel
        periodicBoundaryHandler.addNewParticles();
    }
#endif

    hGridActionsBeforeTimeStep();
    computeAllForces();
    removeDuplicatePeriodicParticles();
    interactionHandler.actionsAfterTimeStep();
    logger(DEBUG, "Have computed the initial values for the forces ");
}

void DPMBase::solve()
{
    initialiseSolve();
    setSoftStop();

    // Can be used to measure simulation time
    clock_.tic();
    // This is the main loop over advancing time
    while (getTime() < getTimeMax() && continueSolve() && continueFlag_)
    {
        computeOneTimeStep();
    }
    // Can be used to measure simulation time
    clock_.toc();

    finaliseSolve();
}

void DPMBase::finaliseSolve() {

    //force writing of the last time step
    forceWriteOutputFiles();

    //end loop over interaction count
    actionsAfterSolve();

    //To make sure getTime gets the correct time for outputting statistics
    finishStatistics();

    closeFiles();
}

/*!
 * \details Performs one time step in the time loop, including updating the time. It is made public, since this makes
 * coupling multiple DPM simulations easier in the future.
 */
void DPMBase::computeOneTimeStep()
{
    logger(DEBUG, "starting computeOneTimeStep()");

    logger(DEBUG, "about to call writeOutputFiles()");
    writeOutputFiles(); //everything is written at the beginning of the time step!

    logger(DEBUG, "about to call hGridActionsBeforeIntegration()");
    hGridActionsBeforeIntegration();

    //Computes the half-time step velocity and full time step position and updates the particles accordingly
    logger(DEBUG, "about to call integrateBeforeForceComputation()");

    integrateBeforeForceComputation();
    //New positions require the MPI and parallel periodic boundaries to do things
    logger(DEBUG, "about to call performGhostParticleUpdate()");
    performGhostParticleUpdate();
    // Some walls need to be aware of the new positions
    wallHandler.actionsAfterParticleGhostUpdate();

    /// \todo MX: this is not true anymore. all boundaries are handled here.
    /// particles have received a position update, so here the deletion boundary deletes particles
    /// \todo add particles need a periodic check

    logger(DEBUG, "about to call checkInteractionWithBoundaries()");
    checkInteractionWithBoundaries(); // INSERTION boundaries handled

    logger(DEBUG, "about to call hGridActionsAfterIntegration()");
    hGridActionsAfterIntegration();

    // Compute forces
    ///\bug{In chute particles are added in actions_before_time_set(), however they are not written to the xballs data yet, but can have a collision and be written to the fstat data}
    // INSERTION/DELETION boundary flag change
    for (BaseBoundary* b : boundaryHandler)
    {
        b->checkBoundaryBeforeTimeStep(this);
    }

    logger(DEBUG, "about to call actionsBeforeTimeStep()");
    actionsBeforeTimeStep();
    particleHandler.actionsBeforeTimeStep();
    //interactionHandler.actionsBeforeTimeStep();

    logger(DEBUG, "about to call checkAndDuplicatePeriodicParticles()");
    checkAndDuplicatePeriodicParticles();

    logger(DEBUG, "about to call hGridActionsBeforeTimeStep()");
    hGridActionsBeforeTimeStep();

    //Creates and updates interactions and computes forces based on these
    logger(DEBUG, "about to call computeAllForces()");
    computeAllForces();

    logger(DEBUG, "about to call removeDuplicatePeriodicParticles()");
    removeDuplicatePeriodicParticles();

    logger(DEBUG, "about to call actionsAfterTimeStep()");
    actionsAfterTimeStep();

    //Computes new velocities and updates the particles accordingly
    logger(DEBUG, "about to call integrateAfterForceComputation()");
    integrateAfterForceComputation();

    //erase interactions that have not been used during the last time step
    //logger(DEBUG, "about to call interactionHandler.eraseOldInteractions(getNumberOfTimeSteps())");
    interactionHandler.eraseOldInteractions(getNumberOfTimeSteps());
    logger(DEBUG, "about to call interactionHandler.actionsAfterTimeStep()");
    interactionHandler.actionsAfterTimeStep();
    particleHandler.actionsAfterTimeStep();

    time_ += timeStep_;
    numberOfTimeSteps_++;

    logger(DEBUG, "finished computeOneTimeStep()");
}

/*!
 * Interprets commands <B>passed in the command line</B> (e.g. -tmin 0 -tmax 100 ...).
 *
 * <TT>argc</TT> gives the number of commands passed, while <TT>argv</TT> stores the commands themselves
 * (as strings).
 *
 * Outputs the name and value of each flag passed, then calls the \ref readNextArgument() function to actually
 * interpret and implement the relevant arguments. Will raise an error if an unknown flag is passed.
 *
 * \param[in] argc "Argument count" - number of individual elements that argv will possess
 * \param[in] *argv[] An array of length argc - specifically an array of strings (or, in C terminology, a character array)
 */
bool DPMBase::readArguments(int argc, char* argv[])
{
    bool isRead = true;
    // Cycles over every second element. (Most flags will contain both name labels and actual data.
    // Those that don't will have to do i--; some examples in readNextArgument.)
    for (int i = 1; i < argc; i += 2)
    {
        logger(INFO, "interpreting input argument %\n", argv[i], Flusher::NO_FLUSH);
        for (int j = i + 1; j < argc; j++)
        {
            //looks for the next string that starts with a minus sign
            //i.e. the next flag, as each flag may take 0, 1 , 2, 3... arguments
            //and we need to make sure all are read in!
            if (argv[j][0] == '-')
                break;
            logger(INFO, " %\n", argv[j], Flusher::NO_FLUSH);
        }
        logger(INFO, "");
        //if "isRead"is true and "readNextArgument" is also true...
        //(i.e. checking if any argument is false)
        isRead &= readNextArgument(i, argc, argv);

        // If the read was unsuccessful, raise an error and quit. (JMFT: Used to just be a warning.)
        if (!isRead)
        {
            logger(ERROR, "Warning: not all arguments read correctly!");
        }
    }
    return isRead;
}

void DPMBase::removeOldFiles() const
{
    //logger(INFO,"ID %",PROCESSOR_ID);
    //if (PROCESSOR_ID!=0) return;

    logger(INFO,"Removing old files named %.*",getName());
    std::ostringstream filename;

    // add processor id to file extension for mpi jobs
    std::string p = (NUMBER_OF_PROCESSORS > 1)?std::to_string(PROCESSOR_ID):"";
    // all the file extensions that should be deleted
    std::vector<std::string> ext{".restart"+p, ".stat"+p, ".fstat"+p, ".data"+p, ".ene"+p, ".xballs"};
    for (const auto& j : ext)
    {
        // remove files with given extension for FileType::ONE_FILE
        filename.str("");
        filename << getName() << j;
        if (!remove(filename.str().c_str()))
        {
            logger(INFO,"  File % successfully deleted",filename.str());
        }
        // remove files with given extension for FileType::MULTIPLE_FILES
        unsigned k = 0;
        filename.str("");
        filename << getName() << j << '.' << k;
        while (!remove(filename.str().c_str()))
        {
            if (k<3) logger(INFO,"  File % successfully deleted",filename.str());
            filename.clear();
            filename.str("");
            filename << getName() << j << '.' << ++k;
        }
        // remove files with given extension for FileType::MULTIPLE_FILES_PADDED
        k = 0;
        filename.str("");
        filename << getName() << j << '.' << to_string_padded(k);
        while (!remove(filename.str().c_str()))
        {
            if (k<3) logger(INFO,"  File % successfully deleted",filename.str());
            filename.clear();
            filename.str("");
            filename << getName() << j << '.' << to_string_padded(++k);
        }
    }
    // remove vtk files
    // add processor id to file extension for mpi jobs
    std::string q = (NUMBER_OF_PROCESSORS > 1)?("Processor_"+std::to_string(PROCESSOR_ID)+"_"):"";
    // all the file extensions that should be deleted
    ext = {"Wall_", q+"Particle_", q+"Interaction_", q+"Boundary_"};
    for (const auto& j : ext)
    {
        // remove files with given extension for FileType::ONE_FILE
        filename.str("");
        filename << getName() << j << ".vtu";
        if (!remove(filename.str().c_str()))
        {
            logger(INFO,"  File % successfully deleted",filename.str());
        }
        // remove files with given extension for FileType::MULTIPLE_FILES
        unsigned k = 0;
        filename.str("");
        filename << getName() << j << k << ".vtu";
        while (!remove(filename.str().c_str()))
        {
            if (k<3) logger(INFO,"  File % successfully deleted",filename.str());
            filename.str("");
            filename << getName() << j << ++k << ".vtu";
        }
        //std::cout << "  File " << filename.str() << " not found" << std::endl;
    }
}

/*!
 * \brief Reads, recognises and applies all valid flags passed when starting or restarting a Mercury simulation.
 * \details For all of the N = <TT>argc</TT> (argument count) command line arguments passed when starting/restarting a
 * code (e.g. -tmax, -tmin ...), compares them to the "known" arguments understood by Mercury (note that further
 * recognised arguments can be added in derived classes).
 * If a match is found, the relevant parameter is set to the corresponding value(s) following the flag and <tt>true</tt>
 * is returned. Otherwise, <tt>false</tt> is returned.
 *
 * For instance, if the flag <TT>-xmin 0</TT>
 * is passed, the code's second if statement will recognise the flag, convert the subsequent string in <TT>argv</TT> to
 * a double, and then call the \ref setXMin() function to implement the new value (0) of XMin.
 *
 * For developers: note the use of strcmp here. This cannot be replaced with a simpler ==, as we are comparing c-style
 * strings (char*), instead of std::string. Thus, == would return equality of the pointers instead of the contents
 * of the string. strcmp returns 0 if the strings are the same, and another number if they are different. This is then
 * implicitly cast to a bool, where 0->false and other numbers will give true. Finally, the !-operator makes sure that
 * the expression in the if-statements are true if the strings are the same, and false otherwise.
 *
 * \param[in] i the position of the element that will be read, note that the count starts at 1, as element 0 is the name
 * of the executable
 * \param[in] argc number of arguments the user has given
 * \param[in] *argv[] the command-line arguments the user has given when calling the executable
 * \return <tt>true</tt> if the argument is successfully read, and <tt>false</tt> otherwise.
 */
bool DPMBase::readNextArgument(int& i, int argc, char* argv[])
{
    // The argument argv[i] identifies the label of the flag, and subsequent arguments (usually 1)
    // contain the content.
    //
    // For example...
    // Checks if the "-name" flag has been passed
    // The strcmp returns 0 if "argv[i]" is "-name" (i.e. !strcmp(argv[i], "-name") --> 1)
    // In this case, the setName function is run with the relevant input (i.e. the value which
    // immediately follows the "-name" flag
    if (!strcmp(argv[i], "-name"))
    {
        setName(argv[i + 1]);
    }
        // The above process is repeated for all viable flags.
    else if (!strcmp(argv[i], "-xmin"))
    {
        setXMin(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-ymin"))
    {
        setYMin(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-zmin"))
    {
        setZMin(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-xmax"))
    {
        setXMax(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-ymax"))
    {
        setYMax(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-zmax"))
    {
        setZMax(atof(argv[i + 1]));
        //} else if (!strcmp(argv[i],"-svn")) {
        //	std::cout << "svn version " << SVN_VERSION << std::endl;
        //	i--;
    }
    else if (!strcmp(argv[i], "-dt"))
    {
        Mdouble old = getTimeStep();
        setTimeStep(atof(argv[i + 1]));
        logger(INFO, "  reset dt from % to %", old, getTimeStep());
    }
//    else if (!strcmp(argv[i], "-Hertz"))
//    {
//        speciesHandler.getObject(0)->setForceType(ForceType::HERTZ);
//        i--;
//    }
    else if (!strcmp(argv[i], "-tmax"))
    {
        Mdouble old = getTimeMax();
        setTimeMax(atof(argv[i + 1]));
        logger(INFO, "  reset timeMax from % to %", old, getTimeMax());
    }
    else if (!strcmp(argv[i], "-saveCount"))
    {
        Mdouble old = dataFile.getSaveCount();
        setSaveCount(static_cast<unsigned int>(atoi(argv[i + 1])));
        logger(INFO, "  reset saveCount from & to %", old, dataFile.getSaveCount());
    }
    else if (!strcmp(argv[i], "-saveCountData"))
    {
        dataFile.setSaveCount(static_cast<unsigned int>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-saveCountFStat"))
    {
        fStatFile.setSaveCount(static_cast<unsigned int>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-saveCountStat"))
    {
        statFile.setSaveCount(static_cast<unsigned int>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-saveCountEne"))
    {
        eneFile.setSaveCount(static_cast<unsigned int>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-saveCountRestart"))
    {
        restartFile.setSaveCount(static_cast<unsigned int>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-dim"))
    {
        setSystemDimensions(static_cast<unsigned int>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-gravity"))
    {
        ///-gravity_ requires three arguments
        setGravity(Vec3D(atof(argv[i + 1]), atof(argv[i + 2]), atof(argv[i + 3])));
        i += 2;
    }
    else if (!strcmp(argv[i], "-fileType"))
    { //uses int input
        setFileType(static_cast<FileType>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-fileTypeFStat"))
    { //uses int input
        fStatFile.setFileType(static_cast<FileType>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-fileTypeRestart"))
    {
        restartFile.setFileType(static_cast<FileType>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-fileTypeData"))
    {
        dataFile.setFileType(static_cast<FileType>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-fileTypeStat"))
    {
        statFile.setFileType(static_cast<FileType>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-fileTypeEne"))
    {
        eneFile.setFileType(static_cast<FileType>(atoi(argv[i + 1])));
    }
    else if (!strcmp(argv[i], "-auto_number"))
    {
        autoNumber();
        i--;
    }
//    else if (!strcmp(argv[i], "-number_of_saves"))
//    {
//        set_number_of_saves_all(atof(argv[i + 1]));
//    }
    else if (!strcmp(argv[i], "-restart") || !strcmp(argv[i], "-r"))
    {
        ///-restart or -r loads a restart file.
        ///By default, it loads <name>.restart.
        ///If an argument "arg" is given it loads the file "arg", or "arg".restart (if the ending is not given).
        std::string filename;

        //use default filename if no argument is given
        if (i + 1 >= argc || argv[i + 1][0] == '-')
        {
            i--;
            filename = getName();
            logger(INFO, "%", getName());
        }
        else
        {
            filename = argv[i + 1];
        }

        //add ".restart" if necessary
        if (filename.find(".restart") == std::string::npos)
        {
            filename = filename + ".restart";
        }

        readRestartFile(filename);
    }
    else if (!strcmp(argv[i], "-clean") || !strcmp(argv[i], "-c"))
    {
        logger(INFO, "Remove old %.* files", getName());
        removeOldFiles();
        i--;
    }
    else if (!strcmp(argv[i], "-data"))
    {
        std::string filename = argv[i + 1];
        readDataFile(filename);
    }
    else if (!strcmp(argv[i], "-readSpeciesFromDataFile"))
    {
        readSpeciesFromDataFile_ = true;
        i--;
        logger(INFO, "Last column of data file will be interpreted as species index");
    }
//    else if (!strcmp(argv[i], "-k"))
//    {
//        Mdouble old = getSpecies(0)->getStiffness();
//        getSpecies(0)->setStiffness(atof(argv[i + 1]));
//        std::cout << "  reset k from " << old << " to " << getSpecies(0)->getStiffness() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-dissipation") || !strcmp(argv[i], "-disp"))
//    {
//        Mdouble old = getSpecies(0)->getDissipation();
//        getSpecies(0)->setDissipation(atof(argv[i + 1]));
//        std::cout << "  reset getDissipation() from " << old << " to " << getSpecies(0)->getDissipation() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-kt"))
//    {
//        Mdouble old = getSpecies(0)->getSlidingStiffness();
//        getSpecies(0)->setSlidingStiffness(atof(argv[i + 1]));
//        std::cout << "  reset kt from " << old << " to " << getSpecies(0)->getSlidingStiffness() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-dispt"))
//    {
//        Mdouble old = getSpecies(0)->getSlidingDissipation();
//        getSpecies(0)->setSlidingDissipation(atof(argv[i + 1]));
//        std::cout << "  reset dispt from " << old << " to " << getSpecies(0)->getSlidingDissipation() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-krolling"))
//    {
//        Mdouble old = getSpecies(0)->getRollingStiffness();
//        getSpecies(0)->setRollingStiffness(atof(argv[i + 1]));
//        std::cout << "  reset krolling from " << old << " to " << getSpecies(0)->getRollingStiffness() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-disprolling"))
//    {
//        Mdouble old = getSpecies(0)->getRollingDissipation();
//        getSpecies(0)->setRollingDissipation(atof(argv[i + 1]));
//        std::cout << "  reset disprolling from " << old << " to " << getSpecies(0)->getRollingDissipation() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-mu"))
//    {
//        Mdouble old = getSpecies(0)->getSlidingFrictionCoefficient();
//        getSpecies(0)->setSlidingFrictionCoefficient(atof(argv[i + 1]));
//        std::cout << "  reset mu from " << old << " to " << getSpecies(0)->getSlidingFrictionCoefficient() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-murolling"))
//    {
//        Mdouble old = getSpecies(0)->getRollingFrictionCoefficient();
//        getSpecies(0)->setRollingFrictionCoefficient(atof(argv[i + 1]));
//        std::cout << "  reset murolling from " << old << " to " << getSpecies(0)->getRollingFrictionCoefficient() << std::endl;
//    }
    else if (!strcmp(argv[i], "-randomise") || !strcmp(argv[i], "-randomize"))
    {
        random.randomise();
        i--;
    }
//    else if (!strcmp(argv[i], "-k0"))
//    {
//        Mdouble old = speciesHandler.getObject(0)->getAdhesionStiffness();
//        speciesHandler.getObject(0)->setAdhesionStiffness(atof(argv[i + 1]));
//        std::cout << "  reset k0 from " << old << " to " << speciesHandler.getObject(0)->getAdhesionStiffness() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-f0"))
//    {
//        Mdouble old = speciesHandler.getObject(0)->getBondForceMax();
//        speciesHandler.getObject(0)->setBondForceMax(atof(argv[i + 1]));
//        std::cout << "  reset f0 from " << old << " to " << speciesHandler.getObject(0)->getBondForceMax() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-AdhesionForceType"))
//    {
//        AdhesionForceType old = speciesHandler.getObject(0)->getAdhesionForceType();
//        speciesHandler.getObject(0)->setAdhesionForceType(argv[i + 1]);
//        std::cout << "  reset AdhesionForceType from "
//                << static_cast<signed char>(old) << " to "
//                << static_cast<signed char>(speciesHandler.getObject(0)->getAdhesionForceType()) << std::endl;
//    }
    else if (!strcmp(argv[i], "-append"))
    {
        setAppend(true);
        i--;
    }
    else if (!strcmp(argv[i], "-fixedParticles"))
    {
        setFixedParticles(static_cast<unsigned int>(atoi(argv[i + 1])));
    }
//    else if (!strcmp(argv[i], "-rho"))
//    {
//        Mdouble old = speciesHandler.getObject(0)->getDensity();
//        speciesHandler.getObject(0)->setDensity(atof(argv[i + 1]));
//        std::cout << "  reset rho from " << old << " to " << speciesHandler.getObject(0)->getDensity() << std::endl;
//    }
//    else if (!strcmp(argv[i], "-dim_particle"))
//    {
//        setParticleDimensions(atoi(argv[i + 1]));
//    }
    else if (!strcmp(argv[i], "-counter"))
    {
        setRunNumber(atoi(argv[i + 1]));
    }
    else
    {
        //returns false if the flag passed does not match any of the currently used flags.
        return false;
    }
    //returns true if argv is found
    return true;
}

/*!
 * \details
 * A very useful feature. For example, when one wants to have an initial condition with particles
 * free of interactions with other particles or walls, one could use this to see if a particle about
 * to be inserted would have interactions. If yes, then the particle would not be considered for
 * insertion.
 * \n\n
 * However can prove expensive if the number of particles is large.
 *
 * \param[in] p The particle for which one wants to detect collisions (or the lack thereof).
 * \return <tt>true</tt> if and only if there are no interactions with other particles or walls.
 */
bool DPMBase::checkParticleForInteraction(const BaseParticle& p)
{
#ifdef MERCURYDPM_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
        return checkParticleForInteractionLocalPeriodic(p);
    }

    int localInteraction = checkParticleForInteractionLocal(p);
    //The root gathers all values and computes the global value
    int *interactionList = nullptr;
    if (PROCESSOR_ID == 0)
    {
        interactionList = new int [NUMBER_OF_PROCESSORS];
    }

    //Gather all local values
    MPIContainer::Instance().gather(localInteraction,interactionList);

    //Compute the global value
    int globalInteraction = 1;
    if (PROCESSOR_ID == 0)
    {
        for (int i = 0; i < NUMBER_OF_PROCESSORS; i++)
        {
            if (interactionList[i] == 0)
            {
                globalInteraction = 0;
                break;
            }
        }
    }
    //The root now tells the other processors what the global value for the interaction is
    MPIContainer::Instance().broadcast(globalInteraction);

    //Convert the result back to bool
    bool interaction = globalInteraction;
#else
    bool interaction = checkParticleForInteractionLocalPeriodic(p);
#endif
    return interaction;
}

/*!
 * Extends the capability of detecting intersecting particles to periodic systems
 * \todo TW the implementation of this function is not very efficient and should be improved
 * @param p
 * @return
 */
bool DPMBase::checkParticleForInteractionLocalPeriodic(const BaseParticle& p)
{
    //A vector of ghost particles of the particle that is to be inserted (empty if no periodic boundaries are present)
    std::vector<Vec3D> pPeriodic;
    for (BaseBoundary* b : boundaryHandler)
    {
        PeriodicBoundary* pb = dynamic_cast<PeriodicBoundary*>(b);
        if (pb && particleHandler.getNumberOfObjects() > 0 )
        {
            const Mdouble maxDistance = p.getMaxInteractionRadius() + particleHandler.getLargestParticle()->getMaxInteractionRadius();
            for (int i = pPeriodic.size() - 1; i >= 0; --i)
            {
                if (pb->getDistance(pPeriodic[i]) < maxDistance)
                {
                    pPeriodic.push_back(pPeriodic[i]);
                    pb->shiftPosition(pPeriodic.back());
                }
            }
            if (pb->getDistance(p) < maxDistance)
            {
                pPeriodic.push_back(p.getPosition());
                pb->shiftPosition(pPeriodic.back());
            }
        }
    }
    //check the particle AND the ghost particles for intersection problems
    bool insertable = checkParticleForInteractionLocal(p);
    if (!pPeriodic.empty()) {
        BaseParticle* q = p.copy();
        for (const Vec3D& pos : pPeriodic) {
            q->setPosition(pos);
            insertable &= checkParticleForInteractionLocal(*q);
        }
        delete q;
    }
    return insertable;
}

/*!
 * \details A very useful feature. For example, when one wants to have an initial condition
 *          with particles free of interactions with other particles or walls, one could use this
 *          method and whether particles are interacting. If yes, then it would not
 *          consider this particle for insertion and continue onto the next particle.
 *          However can prove expensive if the number of particles is large.
 *
 * Returns true if and only if there are no interactions with other particles in the local domain <B>or</B> walls.
 *
 * \param[in] p The particle for which one wants to detect collisions (or the lack thereof).
 * \return bool - true if particle P has no interactions, false if P has one or more interactions with other particles or walls.
 */
bool DPMBase::checkParticleForInteractionLocal(const BaseParticle& p)
{
    Mdouble distance;
    Vec3D normal;

    //Check if it has no collision with walls
    for (BaseWall* w : wallHandler)
    {
        //returns false if the function getDistanceAndNormal returns true,
        //i.e. if there exists an interaction between wall and particle
        // \todo TW getDistanceAndNormal(p,distance,normal) should ideally be replaced by a inContact(p) function, as it doesn't require distance and normal for anything (and walls now can have multiple contacts, soon particles can have it too.)
        if (w->getDistanceAndNormal(p, distance, normal))
        {
            //std::cout<<"failure: Collision with wall: "<<**it<<std::endl;
            return false;
        }
        else
        {
            //std::cout<<"No collision with wall: "<<**it<<std::endl;
        }
    }

    //Check if it has no collision with other particles
    for (BaseParticle* q : particleHandler)
    {
        //returns false if the particle separation is less than the relevant sum of interaction radii
        //(i.e. another particle is in contact with P)
        if (Vec3D::getDistanceSquared(q->getPosition(), p.getPosition())
            < mathsFunc::square(p.getSumOfInteractionRadii(q)))
        {
            //std::cout<<"failure: Collision with particle "<<**it<<std::endl;
            return false;
        }
        else
        {
            //std::cout<<"No collision with particle "<<**it<<std::endl;
        }
    }
    return true;
    /// \todo tw check against periodic copies (see ShearCell3DInitialConditions.cpp)
}

/*!
 * \details Copies particles, interactions assigning species from a local simulation to a global one; useful for the creation
 *              of a cluster.
 *
 * \param[in] particleH: the particle handler from wich particles are copied,
 * \param[in] particleH: the interaction handler from wich interactions are copied,
 * \param[in] species: the species that will be assigned to the particle.
 */
void DPMBase::importParticlesAs(ParticleHandler& particleH, InteractionHandler& interactionH, const ParticleSpecies* species )
{
    size_t nParticlesPreviouslyIn = particleHandler.getSize();
    for (auto k = particleH.begin(); k != particleH.end(); ++k) {
        auto p = particleHandler.copyAndAddObject( *k );
        p->setSpecies(species);
    }

    for (std::vector<BaseInteraction*>::const_iterator i = interactionH.begin(); i != interactionH.end(); ++i) {
        if ( (*i)->getP()->getInvMass() != 0.0 && (*i)->getI()->getInvMass() != 0.0 ) {
            auto j = interactionHandler.copyAndAddObject(*i);
            j->importP(particleHandler.getObject(nParticlesPreviouslyIn + j->getP()->getIndex()));
            j->importI(particleHandler.getObject(nParticlesPreviouslyIn + j->getI()->getIndex()));
            j->setTimeStamp(getNumberOfTimeSteps());
        }
    }
}

/*!
 *
 * \details Removes particles created by \ref checkAndDuplicatePeriodicParticles().
 *
 * Loops (from back to front) through the particle handler, removing particles marked as
 * 'periodic' (i.e. BaseParticle::getPeriodicFromParticle() returns a real value,
 * as opposed to a null pointer).
 *
 * Since the duplicated particles will be stored at the 'back'  of the particle handler,
 * the function ceases when the first non-duplicated particle is found.
 *
 * Note that between these two functions it is not allowed to create additional functions.
 * \image html Walls/periodicBoundary.png
 */
void DPMBase::removeDuplicatePeriodicParticles()
{
#ifdef MERCURYDPM_USE_MPI
    /// \note: when there is a truely parallel computation (cpu > 2) then periodic particles are done in a different manner (see computeOneTimeStep)
    if (NUMBER_OF_PROCESSORS == 1)
    {
#endif
    //Looping from the *back* of the particle handler, where all the periodic or "ghost" particles are stored.
    //The loop ends when there exist no more periodic particles (i.e. "getPeriodicFromParticle()" returns a null pointer)
    for (unsigned int i = particleHandler.getSize();
         i >= 1 && particleHandler.getObject(i - 1)->getPeriodicFromParticle() != nullptr; i--)
    {
        while (!particleHandler.getObject(i - 1)->getInteractions().empty())
        {
            interactionHandler.removeObjectKeepingPeriodics(
                    particleHandler.getObject(i - 1)->getInteractions().front()->getIndex());
        }
        particleHandler.removeObject(i - 1);
    }

    // OMP parallelism
    /*#pragma omp parallel for num_threads(getNumberOfOMPThreads()) //schedule(dynamic)
    for (unsigned int i = particleHandler.getSize(); i >= 1 ; i--)
    {
        if (particleHandler.getObject(i - 1)->getPeriodicFromParticle() != nullptr) {
            while (!particleHandler.getObject(i - 1)->getInteractions().empty()) {
                interactionHandler.removeObjectKeepingPeriodics(
                        particleHandler.getObject(i - 1)->getInteractions().front()->getIndex());
            }
            particleHandler.removeObject(i - 1);
        }
    }*/

#ifdef MERCURYDPM_USE_MPI
    }
#endif
}

/*!
 * \details For all particles in the system, checks their proximity to all periodic boundaries. If a particle is found
 * to be near a periodic boundary, creates and adds a periodic ("ghost") particle.
 * \image html Walls/periodicBoundary.png
 */
void DPMBase::checkAndDuplicatePeriodicParticles()
{
    //Looping over all boundaries in the boundaryHandler
    for (BaseBoundary* boundary : boundaryHandler)
    {
        //Calls the createPeriodicParticles() function which checks if a particle is adequately
        //close to a periodic particle that a periodic (ghost) particle should be created and,
        //if so, creates one and adds it to the system (hence the necessity to keep "N" variable).
        //
        // (The loop is over all boundaries, but if a boundary is not a PeriodicBoundary, then
        // this does nothing.)
        boundary->createPeriodicParticles(particleHandler);
    }

    // OMP parallelism
    /*#pragma omp parallel for num_threads(getNumberOfOMPThreads()) //schedule(dynamic)
    for (int k = 0; k < boundaryHandler.getNumberOfObjects(); k++)
    {
        //Calls the createPeriodicParticles() function which checks if a particle is adequately
        //close to a periodic particle that a periodic (ghost) particle should be created and,
        //if so, creates one and adds it to the system (hence the necessity to keep "N" variable).
        //
        // (The loop is over all boundaries, but if a boundary is not a PeriodicBoundary, then
        // this does nothing.)

        BaseBoundary* boundary = boundaryHandler.getObject(k);
        #pragma omp critical
        boundary->createPeriodicParticles(particleHandler);
    }*/
}

/*!
 * \todo MX: Under construction
 */
void DPMBase::performGhostParticleUpdate()
{
#ifdef MERCURYDPM_USE_MPI
    //MPIContainer& communicator = MPIContainer::Instance();
    if (NUMBER_OF_PROCESSORS == 1) {return;}

    //Update the postion and velocity data of ghosts and perform some bookkeeping
    std::set<BaseParticle*> particlesToBeDeleted;
    domainHandler.updateStatus(particlesToBeDeleted);
    periodicBoundaryHandler.updateStatus(particlesToBeDeleted);

    //Delete particles
    deleteGhostParticles(particlesToBeDeleted);

    //Add new particles
    domainHandler.addNewParticles();
    periodicBoundaryHandler.addNewParticles();
#endif
}

/*!
 * \todo: doc
 */
void DPMBase::deleteGhostParticles(std::set<BaseParticle*>& particlesToBeDeleted)
{
    //Flush mixed particles from lists (MP particles are located in both structures)
    if (periodicBoundaryHandler.getSize() > 0)
    {
        //Flush particles from boundaries
        domainHandler.getCurrentDomain()->flushParticles(particlesToBeDeleted);
        periodicBoundaryHandler.flushParticles(particlesToBeDeleted);
    }

    //Clean communication lists
    domainHandler.getCurrentDomain()->cleanCommunicationLists();
    periodicBoundaryHandler.cleanCommunicationLists();

    //Delete the particles
    for (auto particle_it : particlesToBeDeleted)
    {
        particleHandler.removeGhostObject(particle_it->getIndex());
    }
}


//This function takes a base particle and then copies all the particle information from the root particle to the other particles
// on neighbouring domains.
void DPMBase::synchroniseParticle(BaseParticle* p, unsigned fromProcessor)
{
#ifdef MERCURYDPM_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();

    //The processor that contains the particle that needs to be copied needs to identify the target, and communicate this
    MPIParticle pInfo;
    if (communicator.getProcessorID() == fromProcessor)
    {
        pInfo.copyDataFromParticleToMPIParticle(p);
    }

    //Broadcast from processor i
    communicator.broadcast(&pInfo,MercuryMPIType::PARTICLE,fromProcessor);
    copyDataFromMPIParticleToParticle(&pInfo, p, &particleHandler);
#endif
}

void DPMBase::performGhostVelocityUpdate()
{
#ifdef MERCURYDPM_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1) {return;}
    // \todo If required, I can implement this for periodic particles, first discuss with Thomas if it is actually requiredf
    //periodicDomainHandler.updateVelocity()
    //domainHandler.updateVelocity();
#endif
}

/*!
* \brief signal handler function.
*/
void DPMBase::signalHandler(int signal)
{
    switch (signal)
    {
        case SIGINT:
            if (!continueFlag_)
            {
                logger(INFO, "SIGINT has been captured for the second time!\nMercuryDPM will forcefully exit!");
                std::exit(SIGINT);
            }

            logger(INFO, "SIGINT has been captured!\nMercuryDPM will finish the current time step, then it will stop!\nHit Ctrl+C again to forcefully exit.");
            continueFlag_ = false;
            return;

        case SIGTERM:
            logger(INFO, "\nSIGTERM has been captured!\nMercuryDPM will finish the current time step, then it will stop!");
            continueFlag_ = false;
            return;

        default:
            logger(INFO, "No Signal to Capture!");
    }
}


/*!
* \brief function for setting sigaction constructor.
*/
void DPMBase::setSoftStop()
{
    if (disableSoftStop_)
    {
        logger(INFO, "Soft stop disabled.");
        return;
    }

    struct sigaction act{};
    memset(&act, 0, sizeof(act));

    act.sa_handler = &signalHandler;

    sigaction(SIGINT, &act, nullptr);
    sigaction(SIGTERM, &act, nullptr);

    logger(VERBOSE, "Soft stop enabled.");
}

/*!
 * \details Skims through all the object pointers of type BaseInteraction in the interaction handler. Outputs the type of
 *          interaction between two particles P and I.
 */
void DPMBase::outputInteractionDetails() const
{
    logger(INFO, "Interactions currently in the handler:\n", Flusher::NO_FLUSH);
    //looping over all individual objects in the interactionHandler
    for (BaseInteraction* p : interactionHandler)
    {
        p->write(std::cout);
        logger(INFO, "\nInteraction % % between % and %",
               p->getName(), p->getId(), p->getP()->getId(), p->getI()->getId());
    }
}

/*!
 * \details
 * Returns true if and only if the "time" argument passed to the function is equal to the current simulation time
 * i.e. if "time" is either exactly equal to the current simulation time (\ref getTime() ) or at least lies between this time step and the next increment
 * (this nicely avoids rounding errors!)
 * \return <tt>true</tt> if "time" and \ref getTime() are equal, otherwise <tt>false</tt>.
 *
 */
bool DPMBase::isTimeEqualTo(Mdouble time) const
{
    return getTime() <= time && getTime() + getTimeStep() > time;
}

/*!
 * \brief Sets the number of domains in DPMbase
 * \details The parallel code decomposes the domain according to these values
 * \param[in] numberOfDomains
 */
void DPMBase::setNumberOfDomains(std::vector<unsigned> numberOfDomains)
{
#ifdef MERCURYDPM_USE_MPI
    numberOfDomains_ = numberOfDomains;
    logger(INFO, "Split domain into a %x%x% grid",numberOfDomains[0],numberOfDomains[1],numberOfDomains[2]);
#else
    logger(WARN, "Setting number of domains, but code is not compiled with MPI on");
#endif
}

void DPMBase::splitDomain(DomainSplit domainSplit) {
    //one-d problems
    if (domainSplit == DomainSplit::X) {
        setNumberOfDomains({NUMBER_OF_PROCESSORS,1,1});
        return;
    } else if (domainSplit == DomainSplit::Y) {
        setNumberOfDomains({1,NUMBER_OF_PROCESSORS,1});
        return;
    } else if (domainSplit == DomainSplit::Z) {
        setNumberOfDomains({1,1,NUMBER_OF_PROCESSORS});
        return;
    }
    //two-d problems
    // split into axb grid with a the largest integer that divides NUMBER_OF_PROCESSORS and is smaller than sqrt(NUMBER_OF_PROCESSORS)
    unsigned a;
    for (unsigned n = floor(sqrt(NUMBER_OF_PROCESSORS));n>0; n--) {
        if (NUMBER_OF_PROCESSORS % n == 0) {
            a = n;
            break;
        }
    }
    if (domainSplit == DomainSplit::XY) {
        setNumberOfDomains({NUMBER_OF_PROCESSORS/a,a,1});
        return;
    } else if (domainSplit == DomainSplit::XZ) {
        setNumberOfDomains({NUMBER_OF_PROCESSORS/a,1,a});
        return;
    } else if (domainSplit == DomainSplit::YZ) {
        setNumberOfDomains({1,NUMBER_OF_PROCESSORS/a,a});
        return;
    }
    //three-d problems
    // split into axbxc grid with
    //  - a the largest integer that divides NUMBER_OF_PROCESSORS and is smaller than cbrt(NUMBER_OF_PROCESSORS)
    //  - b the largest integer that divides NUMBER_OF_PROCESSORS/a and is smaller than sqrt(NUMBER_OF_PROCESSORS/a)
    unsigned b;
    for (unsigned n = floor(cbrt(NUMBER_OF_PROCESSORS));n>0; n--) {
        if (NUMBER_OF_PROCESSORS % n == 0) {
            a = n;
            break;
        }
    }
    for (unsigned n = floor(sqrt(NUMBER_OF_PROCESSORS/a));n>0; n--) {
        if (NUMBER_OF_PROCESSORS % (n*a) == 0) {
            b = n;
            break;
        }
    }
    setNumberOfDomains({NUMBER_OF_PROCESSORS/a/b,a,b});
}


/*!
 * \brief returns the number of domains
 * \details number of domains in parallel code in terms of domains in x,y,z direction
 * \return Returns the number of domains in cartesian cooridates of the parallel mesh
 */
std::vector<unsigned> DPMBase::getNumberOfDomains()
{
    return numberOfDomains_;
}

/*!
 * \brief Function that returns a pointer to the domain correseponding to the processor
 * \return Pointer to a domain corresponding to the processor
 */
Domain* DPMBase::getCurrentDomain()
{
    return domainHandler.getCurrentDomain();
}

ParticleVtkWriter* DPMBase::getVtkWriter() const
{
    return vtkWriter_;
}

/*!
 * \details The function first generates random velocities for every particle in the system and then
 * injects the desired kinetic energy and sets the desired mean velocity in the system.
 * \param[in] V_mean_goal The mean velocity you want to set after injecting energy
 * \param[in] Ek_goal  The kinetic energy you want to inject into the system
 */
void DPMBase::setMeanVelocity(Vec3D V_mean_goal)
{
    Vec3D meanVelocity = getTotalMomentum() / getTotalMass();

    //correct the mean velocity to zero
    for (auto& p : particleHandler)
    {
        p->addVelocity(-meanVelocity);
    }
}

/*!
 * \details The function first generates random velocities for every particle in the system and then
 * injects the desired kinetic energy and sets the desired mean velocity in the system.
 * \param[in] V_mean_goal The mean velocity you want to set after injecting energy
 * \param[in] Ek_goal  The kinetic energy you want to inject into the system
 */
void DPMBase::setMeanVelocityAndKineticEnergy(Vec3D V_mean_goal, Mdouble Ek_goal)
{
    Vec3D V_mean;
    Mdouble Ek_mean_goal = 0, Ek_fluct_factor = 0;
    RNG rng;

    //assign random velocity to each particle
    for (auto& p : particleHandler)
    {
        p->setVelocity(Vec3D(rng.getRandomNumber(-1, 1), rng.getRandomNumber(-1, 1), rng.getRandomNumber(-1, 1)));
    }

    //calculate the mean velocity in the system now
    Ek_mean_goal = 0.5 * getTotalMass() * V_mean_goal.getLengthSquared();
    V_mean = getTotalMomentum() / getTotalMass();
    //check if the user input mean kinetic energy is larger than the total kinetic energy input, then return error
    logger.assert_always(0.5 * getTotalMass() * V_mean_goal.getLengthSquared() < Ek_goal,
                         "Too large mean velocity input, Kinetic energy from mean velocity part is larger than the "
                                 "total kinetic energy you want to set");

    //correct the mean velocity to zero
    for (auto& p : particleHandler)
    {
        p->addVelocity(-V_mean);
    }

    //set the new fluctuating velocity based on the goal fluctuating kinetic energy
    Ek_fluct_factor = std::sqrt((Ek_goal - Ek_mean_goal) / getKineticEnergy());
    for (auto& p : particleHandler)
    {
        p->setVelocity(Ek_fluct_factor * p->getVelocity());
    }

    //correct the mean velocity finally to the user set values
    V_mean = getTotalMomentum() / getTotalMass();
    for (auto& p : particleHandler)
    {
        p->addVelocity(V_mean_goal - V_mean);
    }

    //check the final mean velocity and kinetic energy
    logger(INFO, "In DPMBase::setMeanVelocityAndKineticEnergy,\nV_mean_final %\n Ek_final %",
           getTotalMomentum() / getTotalMass(), getKineticEnergy());
}

/*!
 * \return The total volume of the domain.
 */
Mdouble DPMBase::getTotalVolume() const
{
    return (getXMax() - getXMin()) * (getYMax() - getYMin()) * (getZMax() - getZMin());
}

/*!
 * \details The function calculate the kinetic stress tensor based on particle fluctuation velocity.
 * \return The kinetic stress of the whole system (all particles).
 */
Matrix3D DPMBase::getKineticStress() const
{
    Matrix3D F; //set the kinetic energy tensor, this is in terms of Sum(m*v^2)
    Vec3D J; //set the momentum tensor

    //calculate stress for kinetic part
    for (const auto& p : particleHandler)
    {
        F += Matrix3D::dyadic(p->getVelocity(), p->getVelocity()) * p->getMass();
        J += p->getVelocity() * p->getMass();
    }

    Matrix3D stressKinetic = F - Matrix3D::dyadic(J, J) / getTotalMass();
    stressKinetic /= getTotalVolume();
    return stressKinetic;
}

/*!
 * \details The function calculate the static stress tensor based on particle contact force and
 * contact normal branch vector.
 * \return The static stress of the whole system (all interactions).
 */
Matrix3D DPMBase::getStaticStress() const
{
    //stress components calculation variables
    Matrix3D stressStatic;

    //calculate the static stress tensor based on all the interactions
    for (const auto i : interactionHandler)
    {
        stressStatic += Matrix3D::dyadic(i->getForce(), i->getNormal()) * i->getDistance();
    }

    stressStatic /= getTotalVolume();
    return stressStatic;
}

/*!
 * \details The function calculate the total stress tensor which is
 * the sum of kinetic and static stress tensors.
 * \return The total stress of the whole system (all particles and all interactions).
 */
Matrix3D DPMBase::getTotalStress() const
{
    return getKineticStress() + getStaticStress();
}


void DPMBase::computeWallForces(BaseWall* const w)
{
    //compute forces for all particles that are neither fixed or ghosts
    for (auto p : particleHandler)
    {
        if (!p->isFixed() && p->getPeriodicFromParticle() == nullptr)
        {
            //w->computeForces(p);
            computeForcesDueToWalls(p, w);
        }
    }
}

/*!
 * \details sets the number of time steps skipped between each save for ALL data files, except for the interaction file.
 * And the increament of this number of time steps is based on a user input logarithmicSaveCountBase, e.g. if you put 10
 * as your input, the saving point will be at 10, 10^1, 10^2, 10^3, which results equal distance on a log scale.
 * Note, that the interaction file is independent of time steps, and just writes when an interaction starts or ends.
 *  \param[in] logarithmicSaveCountBase
 */
void DPMBase::setLogarithmicSaveCount(const Mdouble logarithmicSaveCountBase)
{
    dataFile.setlogarithmicSaveCount(logarithmicSaveCountBase);
    fStatFile.setlogarithmicSaveCount(logarithmicSaveCountBase);
    restartFile.setlogarithmicSaveCount(logarithmicSaveCountBase);
    statFile.setlogarithmicSaveCount(logarithmicSaveCountBase);
    eneFile.setlogarithmicSaveCount(logarithmicSaveCountBase);
}

/*!
 * \details This function is called by ParticleHandler::removeObject and
 * ParticleHandler::removeGhostObject to broadcast the removal of a particle
 * from the handler. It passes the information on to the walls in the wallHandler.
 *  \param[in] id The id of the removed particle.
 */
void DPMBase::handleParticleRemoval(unsigned int id)
{
    for (auto w: wallHandler)
    {
        w->handleParticleRemoval(id);
    }
}

/*!
 * \details This function is called by ParticleHandler::addObject and
 * ParticleHandler::addGhostObject to broadcast the addition of a particle
 * from the handler. It passes the information on to the walls in the wallHandler.
 *  \param[in] id The id of the removed particle.
 *  \param[in] p A pointer to the particle.
 */
void DPMBase::handleParticleAddition(unsigned int id, BaseParticle* p)
{
    for (auto w: wallHandler)
    {
        w->handleParticleAddition(id, p);
    }
}

// initialise static member variables
volatile sig_atomic_t DPMBase::continueFlag_ = true;

/// \todo When restarting the indexMax should be reset

/*!
* \details
* Enables/disables the writePythonFileForVTKVisualisation() function.
* \param[in] forceWritePythonFileForVTKVisualisation
*/
void DPMBase::setWritePythonFileForVTKVisualisation(bool forceWritePythonFileForVTKVisualisation)
{
    forceWritePythonFileForVTKVisualisation_ = forceWritePythonFileForVTKVisualisation;
}

bool DPMBase::getWritePythonFileForVTKVisualisation() const
{
    return forceWritePythonFileForVTKVisualisation_;
}

WallVTKWriter & DPMBase::getWallVTKWriter()
{
    return wallVTKWriter_;
}
