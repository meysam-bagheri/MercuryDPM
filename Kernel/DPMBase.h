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

#ifndef DPMBase_H
#define DPMBase_H

//so that the user doesn't have to include string/io manipulations:
#include <string>
#include <iomanip>
// include file for signal handler types
#include <csignal>
//The vector class contains a 3D vector class.
#include "Math/Vector.h"
//This class defines the particle handler
#include "ParticleHandler.h"
//This class defines the base particle (such that not every Driver has to include it)
#include "Particles/BaseParticle.h"
#include "Particles/SphericalParticle.h"
//This class defines the wall handler
#include "WallHandler.h"
//This class defines the boundary handler
#include "BoundaryHandler.h"
#include "PeriodicBoundaryHandler.h"
#include "DomainHandler.h"
//This class defines the interaction handler
#include "InteractionHandler.h"
//This class defines the Species handler
#include "SpeciesHandler.h"
//This class defines the cg handler
#include "CG/CGHandler.h"
//This class defines the possibleContact lists
#ifdef CONTACT_LIST_HGRID
#include "PossibleContactList.h"
#endif
//This class defines the random number generator
#include "Math/RNG.h"
#include "Domain.h"
#include "VTKWriter/ParticleVtkWriter.h"
#include "VTKWriter/WallVTKWriter.h"
#include "VTKWriter/InteractionVTKWriter.h"
#include "VTKWriter/BoundaryVTKWriter.h"
#include "VTKWriter/WallDetailsVTKWriter.h"
#include "MercuryTime.h"


/*!
 * \class DPMBase
 * \brief The DPMBase header includes quite a few header files, defining all the
 *        handlers, which are essential. Moreover, it defines and solves a DPM problem.
 *        It is inherited from FilesAndRunNumber (public). 
 * \bug When restarting the first time step is not saved, therefore there is a missing
 *      time step after a restart
 */
class DPMBase
{
    /*!
     *  A simple vector of vectors for collecting and ordering interactions in the OpenMP parallel environment
     */

public:
    
    /*!
     * \brief A function which initialises the member variables to default values, so
     *        that the problem can be solved off the shelf; sets up a basic two dimensional
     *        problem which can be solved off the shelf. It is called in the constructor DPMBase().
     */
    void constructor();
    
    /*!
     * \brief Constructor that calls the "void constructor()"
     */
    DPMBase();
    
    /*!
     * \brief Copy constructor type-2
     */
    DPMBase(const DPMBase& other);
    
    /*!
     * \brief virtual destructor
     */
    virtual ~DPMBase();
    
    /*!
     * \brief Increment the run Number (counter value) stored in the file_counter (COUNTER_DONOTDEL) by 1 and store the new value
     * in the counter file
     */
    static void incrementRunNumberInFile();
    
    /*!
     * \brief Read the run number or the counter from the counter file (COUNTER_DONOTDEL)
     */
    static int readRunNumberFromFile();
    
    /*!
     * \brief The autoNumber() function calls three functions: setRunNumber(), readRunNumberFromFile() and
     * incrementRunNumberInFile().
     */
    void autoNumber();

    /*!
     * \brief This turns a counter into 1 index, which is a useful feature for performing 1D parameter study.
     * The index run from 1:size_x, while the study number starts at 0 (initially the counter=1 in COUNTER_DONOTDEL)
     */
    std::vector<int> get1DParametersFromRunNumber(int size_x) const;

    /*!
     * \brief This turns a counter into 2 indices which is a very useful feature for performing a 2D study.
     * The indices run from 1:size_x and 1:size_y, while the study number starts at 0 ( initially the counter=1 in COUNTER_DONOTDEL)
     */
    std::vector<int> get2DParametersFromRunNumber(int size_x, int size_y) const;

    /*!
     * \brief This turns a counter into 3 indices, which is a useful feature for performing a 3D parameter study.
     * The indices run from 1:size_x, 1:size_y and 1:size_z, while the study number starts at 0 ( initially the counter=1 in COUNTER_DONOTDEL)
     */
    std::vector<int> get3DParametersFromRunNumber(int size_x, int size_y, int size_z) const;
    
    /*!
     * \brief This launches a code from within this code. Please pass the name of the code to run.
     */
    int launchNewRun(const char* name, bool quick = false);
    
    //setters and getters
    
    /*!
     * \brief This sets the counter/Run number, overriding the defaults
     */
    void setRunNumber(int runNumber);
    
    /*!
     * \brief This returns the current value of the counter (runNumber_)
     */
    int getRunNumber() const;
    
    /*!
     * \brief Sends particles from processorId to the root processor
     */
//    void sendParticlesToRoot(unsigned int &numberOfLocalParticles, unsigned int processorId, ParticleHandler &tempParticleHandler) const;
    
    /*!
     * \brief Decomposes the simulation domain in a structured cube mesh of domains for parallel processing
     */
    virtual void decompose();
    
    /*!
     * \brief The work horse of the code
     * \todo IFCD AT, TW: Consider moving some things before the time loop to actionsBeforeTimeLoop
     */
    void solve();

    /*!
     * \brief Beginning of the solve routine, before time stepping
     */
    void initialiseSolve();

    /*!
     * \brief End of the solve routine, after time stepping
     */
    void finaliseSolve();

    /*!
     * \brief Performs everything needed for one time step, used in the time-loop of solve().
     * \todo IFCD @AT, TW: please check what should be in here, and whether it should be virtual or not.
     */
    virtual void computeOneTimeStep();
    
    /*!
     * \brief Checks if the essentials are set properly to go ahead with solving the problem.
     *
     */
    void checkSettings();
    
    /*!
     * \brief Writes output files immediately, even if the current time step was
     * not meant to be written. Also resets the last saved time step. */
    void forceWriteOutputFiles();
    
    /*!
     * \brief Writes simulation data to all the main Mercury files: .data, .ene, .fstat, .xballs and .restart (see the
     * <A HREF="http://mercurydpm.org/assets/downloads/MercuryLesson/MercuryDPMLessonSlides.pdf">Mercury website</A> for more details regarding these files).
     */
    virtual void writeOutputFiles();
    
    /*!
     * \brief The work horse of the code. Can handle flags from the command line.
     */
    void solve(int argc, char* argv[]);
    
    /*!
     * \brief This function allows to set the initial conditions for our problem to be solved, 
     *        by default particle locations are randomly set. Remember particle properties must
     *        also be defined here.
     * \todo I (Anthony) wants to change this to be an external function. This has a lot of advantages
     *       especially when using copy-constructors. This is a major change and will break other codes,
     *       so therefore has to be done carefully.
     * \details This sets up the particles initial conditions it is as you expect the user to override
     *          this. By default the particles are randomly distributed
     */
    virtual void setupInitialConditions();
    
    /*!
     * \brief This writes a script which can be used to load the xballs problem to display the data just generated
     * \todo Implement or make pure virtual
     */
    virtual void writeXBallsScript() const;
    
    /*!
     * \brief A virtual function that returns some user-specified information about a particle.
     */
    virtual Mdouble getInfo(const BaseParticle& P) const;
    
    ParticleVtkWriter* getVtkWriter() const;
    
    /*!
     * \brief Stores all the particle data for current save time step to a "restart" file, which is a file
     * simply intended to store all the information necessary to "restart" a simulation from a given time step
     * (see also <A HREF=http://mercurydpm.org/assets/downloads/MercuryLesson/MercuryDPMLessonSlides.pdf> MercuryDPM.org</A> for
     * more information on restart files).
     *
     */
    virtual void writeRestartFile();
    
    void writeDataFile();
    
    void writeEneFile();
    
    void writeFStatFile();

    void fillDomainWithParticles(unsigned N=50);

    enum class ReadOptions : int {
        ReadAll,
        ReadNoInteractions,
        ReadNoParticlesAndInteractions
    };

    /*!
     * \brief Reads all the particle data corresponding to a given, existing . restart file (for more details regarding restart files,
     * refer to the <A HREF="http://mercurydpm.org/assets/downloads/MercuryLesson/MercuryDPMLessonSlides.pdf">training materials</A> on the
     *        MercuryDPM website).Returns true if it is successful, false otherwise.
     */
    bool readRestartFile(ReadOptions opt = ReadOptions::ReadAll);

    /*!
     * \brief The same as readRestartFile(bool), but also reads all the particle data corresponding to the current saved
     * time step.
     */
    int readRestartFile(std::string fileName, ReadOptions opt = ReadOptions::ReadAll);

//    /*!
//     * \brief Loads all MD data and plots statistics for all time steps in the .data file
//     */
//    void statisticsFromRestartData(const char *name);
///\todo what to do with statisticsFromRestartData?

    /*!
     * \brief Writes all data into a restart file
     */
    virtual void write(std::ostream& os, bool writeAllParticles = true) const;

    /*!
     * \brief Reads all data from a restart file, e.g. domain data and particle data
     * \todo warning: hides non-virtual function from the class 'Files'.
     */
    virtual void read(std::istream& is, ReadOptions opt = ReadOptions::ReadAll);

    /*!
     * \brief Allows you to read in a wall defined in a Driver directory; see USER/Luca/ScrewFiller
     */
    virtual BaseWall* readUserDefinedWall(const std::string& type) const
    { return nullptr; }

    /*!
     * \brief Reads all data from a restart file, e.g. domain data and particle data; old version
     * \deprecated Use read or the new magic driver instead
     */
    virtual void readOld(std::istream& is);

    /*!
     * \brief This allows particle data to be reloaded from data files
     */
    bool readDataFile(std::string fileName = "", unsigned int format = 0);

    /*!
     * \brief Allows the user to read par.ini files (useful to read files produced by the MDCLR simulation code - external to MercuryDPM)
     */
    bool readParAndIniFiles(std::string fileName);

    /*!
     * \brief Reads the next data file with default format=0. However, one can
     *        modify the format based on whether the particle data corresponds to
     *        3D or 2D data- see \ref xballs.
     */
    bool readNextDataFile(unsigned int format = 0);

    /*!
     * \brief Reads the next fstat file.
     */
    void readNextFStatFile();

    /*!
     * \brief Finds and opens the next data file, if such a file exists.
     */
    bool findNextExistingDataFile(Mdouble tMin, bool verbose = true);

    /*!
     * \brief Can interpret main function input arguments that are passed by the
     *        driver codes
     */
    bool readArguments(int argc, char* argv[]);

    /*!
     * \brief Interprets the i^th command-line argument
     */
    virtual bool readNextArgument(int& i, int argc, char* argv[]);

    /*!
     * \brief Checks whether a particle P has any interaction with walls or other particles.
     */
    virtual bool checkParticleForInteraction(const BaseParticle& P);

    /*!
     * \brief Checks if a particle P has any interaction with walls or other particles in the local domain.
     */
    virtual bool checkParticleForInteractionLocal(const BaseParticle& P);

    bool checkParticleForInteractionLocalPeriodic(const BaseParticle& P);

    void readSpeciesFromDataFile(bool read = true){readSpeciesFromDataFile_=read;}

    /*!
     * \brief Copies particles, interactions assigning species from a local simulation to a global one.
     *        Useful for the creation of a cluster.
     */
    void importParticlesAs(ParticleHandler& particleHandler, InteractionHandler& interactionHandler, const ParticleSpecies* species );

    //getters and setters

    /*!
 * \brief The non const version. Allows one to edit the File::dataFile
 * \deprecated dataFile is now protected, so it can be used by all applications.
 * Please don't use getDataFile() anymore.
 */
    MERCURYDPM_DEPRECATED
    File& getDataFile();

    /*!
     * \brief The non const version. Allows to edit the File::eneFile
     * \deprecated eneFile is now protected, so it can be used by all applications.
     * Please don't use getEneFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    File& getEneFile();

    /*!
     * \brief The non const version. Allows to edit the File::fStatFile
     * \deprecated fStatFile is now protected, so it can be used by all applications.
     * Please don't use getFStatFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    File& getFStatFile();

    /*!
     * \brief The non const version. Allows to edit the File::restartFile
     * \deprecated restartFile is now protected, so it can be used by all applications.
     * Please don't use getRestartFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    File& getRestartFile();

    /*!
     * \brief The non const version. Allows to edit the File::statFile
     * \deprecated statFile is now protected, so it can be used by all applications.
     * Please don't use getStatFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    File& getStatFile();

    /*!
     * \brief Return a reference to the file InteractionsFile
     */
    File& getInteractionFile();

    /*!
     * \brief The const version. Does not allow for any editing of the File::dataFile
     * \deprecated dataFile is now protected, so it can be used by all applications.
     * Please don't use getDataFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    const File& getDataFile() const;

    /*!
     * \brief The const version. Does not allow for any editing of the File::eneFile
     * \deprecated eneFile is now protected, so it can be used by all applications.
     * Please don't use getEneFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    const File& getEneFile() const;

    /*!
     * \brief The const version. Does not allow for any editing of the File::fStatFile
     * \deprecated fStatFile is now protected, so it can be used by all applications.
     * Please don't use getFStatFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    const File& getFStatFile() const;

    /*!
     * \brief The const version. Does not allow for any editing of the File::restartFile
     * \deprecated restartFile is now protected, so it can be used by all applications.
     * Please don't use getRestartFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    const File& getRestartFile() const;

    /*!
     * \brief The const version. Does not allow for any editing of the File::statFile
     * \deprecated statFile is now protected, so it can be used by all applications.
     * Please don't use getStatFile() anymore.
     */
    MERCURYDPM_DEPRECATED
    const File& getStatFile() const;

    /*!
     * \brief Returns a constant reference to an Interactions file
     */
    const File& getInteractionFile() const;

    /*!
     * \brief Returns the name of the file. Does not allow to change it though.
     */
    const std::string& getName() const;

    /*!
     * \brief Allows to set the name of all the files (ene, data, fstat, restart, stat)
     */
    void setName(const std::string& name);

    /*!
     * \brief Calls setName(std::string)
     */
    void setName(const char* name);

    /*!
     * \brief Sets File::saveCount_ for all files (ene, data, fstat, restart, stat)
     */
    void setSaveCount(unsigned int saveCount);

    /*!
     * \brief Sets File::fileType_ for all files (ene, data, fstat, restart, stat)
     */
    void setFileType(FileType fileType);

    /*!
     * \brief Sets File::openMode_ for all files (ene, data, fstat, restart, stat)
     */
    void setOpenMode(std::fstream::openmode openMode);


    //other member functions

    /*!
     * \brief Resets the file counter for each file i.e. for ene, data, fstat, restart, stat)
     */
    void resetFileCounter();

    /*!
     * \brief Closes all files (ene, data, fstat, restart, stat) that were opened to read or write.
     */
    void closeFiles();

    /*!
     * \brief Sets the next time step for all the files (ene, data, fstat, restart, stat) at which the data is to be written or saved.
     */
    void setLastSavedTimeStep(unsigned int nextSavedTimeStep);

    /*!
     * \brief Returns the current simulation time.
     */
    Mdouble getTime() const;

    /*!
     * \brief Returns the current simulation time.
     */
    Mdouble getNextTime() const;

    /*!
     * \brief Returns the current counter of time-steps, i.e. the number of time-steps that the
     * simulation has undergone so far.
     */
    unsigned int getNumberOfTimeSteps() const;

    /*!
     * \brief Sets a new value for the current simulation time.
     */
    void setTime(Mdouble time);

    /*!
     * \brief Sets a new value for the maximum simulation duration.
     */
    void setTimeMax(Mdouble newTMax);

    /*!
     * \brief Returns the maximum simulation duration.
     */
    Mdouble getTimeMax() const;

    /*!
     * \brief Sets File::logarithmicSaveCount_ for all files (ene, data, fstat, restart, stat)
     */
    void setLogarithmicSaveCount(Mdouble logarithmicSaveCountBase);

    /*!
     * \brief set the number of elements to write to the screen
     */
    void setNToWrite(int nToWrite);

    /*!
     * \brief get the number of elements to write to the
     */
    int getNToWrite() const;

#ifdef CONTACT_LIST_HGRID
    /*!
     * \brief Returns the linked list of all possible contacts between particles
     */
        PossibleContactList& getPossibleContactList();
#endif

    /*!
     * \todo{these functions should also update the mixed species}
     */

    /*!
     * \brief Sets whether particle rotation is enabled or disabled.
     * \details * Passing <tt>true</tt> will enable particle rotation.
     * Passing <tt>false</tt> will disable particle rotation.
     * \param[in] newRotFlag
     */
    void setRotation(bool rotation)
    { rotation_ = rotation; }

    /*!
     * \brief Indicates whether particle rotation is enabled or disabled.
     * \returns <tt>true</tt> if particle rotation is enabled; <tt>false</tt> if particle rotation is disabled.
     */
    bool getRotation() const
    { return rotation_; }

    MERCURYDPM_DEPRECATED
    /*!
     * \brief Sets whether walls are written into a VTK file.
     * \deprecated Use wallHandler.setWriteVTK(FileType) instead (since 21-07-2021)
     */
    void setWallsWriteVTK(FileType writeWallsVTK);

    MERCURYDPM_DEPRECATED
    /*!
     * \brief Sets whether walls are written into a VTK file.
     * \deprecated Use wallHandler.setWriteVTK(bool) instead (since 21-07-2021)
     */
    void setWallsWriteVTK(bool);

    MERCURYDPM_DEPRECATED
    /*!
     * \brief Sets whether interactions are written into a VTK file.
     * \deprecated Use interactionHandler.setWriteVTK(bool) instead (since 21-07-2021)
     */
    void setInteractionsWriteVTK(bool);

    /*!
     * \brief Sets whether particles are written in a VTK file.
     */
    void setParticlesWriteVTK(bool writeParticlesVTK);

    void setSuperquadricParticlesWriteVTK(bool writeSuperquadricParticlesVTK);

    MERCURYDPM_DEPRECATED
    /*!
     * \brief Returns whether walls are written in a VTK file.
     * \deprecated Use wallHandler.getWriteVTK() instead (since 21-07-2021)
     */
    FileType getWallsWriteVTK() const;

    /*!
     * \brief Returns whether particles are written in a VTK file.
     */
    bool getParticlesWriteVTK() const;

    bool getSuperquadricParticlesWriteVTK() const;

    /*!
     * \brief If the length of the problem domain in x-direction is XMax - XMin,
     *        then getXMin() returns XMin
     */
    Mdouble getXMin() const
    { return min_.x(); }

    /*!
     * \brief If the length of the problem domain in x-direction is XMax - XMin,
             then getXMax() returns XMax
     */
    Mdouble getXMax() const
    { return max_.x(); }

    /*!
     * \brief If the length of the problem domain in y-direction is YMax - YMin, then getYMin() returns YMin
     */
    Mdouble getYMin() const
    { return min_.y(); }

    /*!
     * \brief If the length of the problem domain in y-direction is YMax - YMin, then getYMax() returns XMax
     */
    Mdouble getYMax() const
    { return max_.y(); }

    /*!
     * \brief If the length of the problem domain in z-direction is ZMax - ZMin, then getZMin() returns ZMin
     */
    Mdouble getZMin() const
    { return min_.z(); }

    /*!
     * \brief If the length of the problem domain in z-direction is ZMax - ZMin, then getZMax() returns ZMax
     */
    Mdouble getZMax() const
    { return max_.z(); }

    Mdouble getXCenter() const
    { return 0.5*(max_.X+min_.X); }

    Mdouble getYCenter() const
    { return 0.5*(max_.Y+min_.Y); }

    Mdouble getZCenter() const
    { return 0.5*(max_.Z+min_.Z); }
    /*!
     * \brief Returns the minimum coordinates of the problem domain.
     */
    Vec3D getMin() const
    { return min_; }

    /*!
     * \brief Returns the maximum coordinates of the problem domain.
     */
    Vec3D getMax() const
    { return max_; }

    /*!
     * \deprecated
     * \brief Sets the value of XMin, the lower bound of the problem domain in the x-direction
     */
    void setXMin(Mdouble newXMin);

    /*!
     * \deprecated
     * \brief Sets the value of YMin, the lower bound of the problem domain in the y-direction
     */
    void setYMin(Mdouble newYMin);

    /*!
     * \deprecated
     * \brief Sets the value of ZMin, the lower bound of the problem domain in the z-direction
     */
    void setZMin(Mdouble newZMin);

    /*!
     * \deprecated
     * \brief Sets the value of XMax, the upper bound of the problem domain in the x-direction
     */
    void setXMax(Mdouble newXMax);

    /*!
     * \deprecated
     * \brief Sets the value of YMax, the upper bound of the problem domain in the y-direction
     */
    void setYMax(Mdouble newYMax);

    /*!
     * \deprecated
     * \brief Sets the value of ZMax, the upper bound of the problem domain in the z-direction
     */
    void setZMax(Mdouble newZMax);

    /*!
     * \brief Sets the maximum coordinates of the problem domain.
     */
    void setMax(const Vec3D& max);

    /*!
     * \brief Sets the maximum coordinates of the problem domain.
     */
    void setMax(Mdouble, Mdouble, Mdouble);

    /*!
     * \brief Sets the minimum coordinates of the problem domain.
     */
    void setDomain(const Vec3D& min, const Vec3D& max);

    /*!
     * \brief Sets the minimum coordinates of the problem domain.
     */
    void setMin(const Vec3D& min);

    /*!
     * \brief Sets the minimum coordinates of the problem domain.
     */
    void setMin(Mdouble, Mdouble, Mdouble);


    /*!
     * \brief Sets a new value for the simulation time step.
     */
    void setTimeStep(Mdouble newDt);

    /*!
     * \brief Returns the simulation time step.
     */
    Mdouble getTimeStep() const;

    /*! \brief Sets the number of omp threads */
    void setNumberOfOMPThreads(int numberOfOMPThreads);

    /*! \brief Returns the number of omp threads */
    int getNumberOfOMPThreads() const;

    /*!
     * \brief Set the xballs output mode.
     */
    void setXBallsColourMode(int newCMode);

    /*!
     * \brief Get the xballs colour mode (CMode).
     */
    int getXBallsColourMode() const;

    /*!
     * \brief Set the scale of vectors in xballs.
     */
    void setXBallsVectorScale(double newVScale);

    /*!
     * \brief Returns the scale of vectors used in xballs.
     */
    double getXBallsVectorScale() const;

    /*!
     * \brief Set the additional arguments for xballs
     */
    void setXBallsAdditionalArguments(std::string newXBArgs);

    /*!
     * \brief Returns the additional arguments for xballs
     */
    std::string getXBallsAdditionalArguments() const;

    /*!
     * \brief Sets the scale of the view (either normal, zoom in or zoom out) to
     *        display in xballs. The default is fit to screen
     */
    void setXBallsScale(Mdouble newScale);

    /*!
     * \brief Returns the scale of the view in xballs.
     */
    double getXBallsScale() const;

    /*!
     * \brief Sets a new value for the gravitational acceleration.
     */
    void setGravity(Vec3D newGravity);

    /*!
     * \brief Returns the gravitational acceleration.
     */
    Vec3D getGravity() const;

    /*!
     * \brief Simple access function to turn on a background drag.
     * The force of particleVelocity*drag is applied (note, it allowd to be negaitve i.e. create energy)
     */
     void setBackgroundDrag(Mdouble backgroundDrag){backgroundDrag_=backgroundDrag;}

     /*!
      * \brief Return the background drag
      * @return double which is the background drag. Positive number removes energy, negative gains.
      */
     const Mdouble getBackgroundDrag() const {return backgroundDrag_;}

    /*!
     * \brief Sets both the system dimensions and the particle dimensionality.
     */
    void setDimension(unsigned int newDim);

    /*!
     * \brief Sets the system dimensionality.
     */
    void setSystemDimensions(unsigned int newDim);

    /*!
     * \brief Returns the system dimensionality.
     */
    unsigned int getSystemDimensions() const;

    /*!
    * \brief Sets the particle dimensionality.
    */
    void setParticleDimensions(unsigned int particleDimensions);

    /*!
     * \brief Returns the particle dimensionality.
     */
    unsigned int getParticleDimensions() const;

    /*!
     * \brief This is to take into account for different Mercury versions.
     *        Returns the version of the restart file.
     */
    std::string getRestartVersion() const;

    /*!
     * \brief Sets restart_version
     */
    void setRestartVersion(std::string newRV);

    /*!
     * \brief Returns the flag denoting if the simulation was restarted or not.
     */
    bool getRestarted() const;

    /*!
     * \brief Allows to set the flag stating if the simulation is to be restarted or not.
     */
    void setRestarted(bool newRestartedFlag);

    /*!
     * \brief Returns whether the "append" option is on or off.
     */
    bool getAppend() const;

    /*!
     * \brief Sets whether the "append" option is on or off.
     */
    void setAppend(bool newAppendFlag);

    /*!
     * \brief Returns the global elastic energy within the system.
     */
    Mdouble getElasticEnergy() const;

    /*!
     * \brief Returns the global kinetic energy stored in the system.
     */
    Mdouble getKineticEnergy() const;

    /*!
     * \brief Returns the global gravitational potential energy stored in the system.
     */
    Mdouble getGravitationalEnergy() const;

    /*!
     * \brief JMFT Returns the global rotational energy stored in the system.
     */
    Mdouble getRotationalEnergy() const;

    Mdouble getTotalEnergy() const;

    /*!
     * \brief JMFT: Return the total mass of the system, excluding fixed particles.
     */
    Mdouble getTotalMass() const;

    /*!
     * \brief JMFT: Return the centre of mass of the system, excluding fixed particles.
     */
    Vec3D getCentreOfMass() const;

    /*!
     * \brief JMFT: Return the total momentum of the system, excluding fixed particles.
     */
    Vec3D getTotalMomentum() const;

    /*!
     * \brief
     */
    double getCPUTime() { return clock_.getCPUTime(); }

    /*!
     * \brief
     */
    double getWallTime() { return clock_.getWallTime(); }


    /*!
     * \brief Checks if two particle are in contact or is there any positive overlap
     */
    static bool areInContact(const BaseParticle* pI, const BaseParticle* pJ);

    /// \bug Why are the hGRID actions public, this seems wrong. Someone please comment [Ant].
    /*!
     * \brief
     */
    virtual void hGridInsertParticle(BaseParticle* obj UNUSED);

    /*!
     * \brief
     */
    virtual void hGridUpdateParticle(BaseParticle* obj UNUSED);

    /*!
     * \brief
     */
    virtual void hGridRemoveParticle(BaseParticle* obj UNUSED);

    /*!
     * \brief
     */
    virtual void hGridUpdateMove(BaseParticle*, Mdouble);

    /*!
     * \brief Checks if the position of the particle is in an mpi communication zone or not
     */
    bool mpiIsInCommunicationZone(BaseParticle* particle);

    /*!
     * \brief Function that checks if the mpi particle should really be inserted by the current domain
     */
    bool mpiInsertParticleCheck(BaseParticle* P);

    /*!
     * \brief This function inserts a particle in the mpi communication boundaries
     */
    void insertGhostParticle(BaseParticle* P);

    /*!
     * \brief Checks if the Domain/periodic interaction distance needs to be updated and updates it accordingly.
     */
    void updateGhostGrid(BaseParticle* P);

    /*!
     * \brief
     * //Not unsigned index because of possible wall collisions.
     */
    virtual void gatherContactStatistics(unsigned int index1, int index2, Vec3D Contact, Mdouble delta, Mdouble ctheta,
                                         Mdouble fdotn, Mdouble fdott, Vec3D P1_P2_normal_, Vec3D P1_P2_tangential);

    /*!
     * \brief Sets the number of domains in x-,y- and z-direction. Required for parallel computations
     */
    void setNumberOfDomains(std::vector<unsigned> direction);

    enum class DomainSplit {X, Y, Z, XY, XZ, YZ, XYZ};

    /*!
     * Splits domain as neatly as possible.
     * e.g. splitDomain(XY) splits domain into a 6x5x1 grid for 30 processors, a 6x6x1 grid for 36 processors
     * \todo the function needs improvement: for non-cubic domains (with different domain length in x,y,z), having a equal number of grid cell is not the best choice.
     */
    void splitDomain(DomainSplit domainSplit);

    /*!
     * \brief returns the number of domains
     */
    std::vector<unsigned> getNumberOfDomains();

    /*!
     * \brief Function that returns a pointer to the domain corresponding to the processor
     */
    Domain* getCurrentDomain();

    void removeOldFiles() const;

    /*!
     * \brief Creates a list of neighbour particles obtained from the hgrid
     */
    virtual void hGridGetInteractingParticleList(BaseParticle* obj, std::vector<BaseParticle*>& list)
    {};

    virtual void computeWallForces(BaseWall* w);

    /*!
     * \brief
     */
    virtual bool getHGridUpdateEachTimeStep() const;

    /// \brief This function will help you set a fixed kinetic energy and mean velocity in your system.
    void setMeanVelocity(Vec3D V_mean_goal);

    /// \brief This function will help you set a fixed kinetic energy and mean velocity in your system.
    void setMeanVelocityAndKineticEnergy(Vec3D V_mean_goal, Mdouble Ek_goal);

    /// \brief Get the total volume of the cuboid system.
    Mdouble getTotalVolume() const;

    /// \brief Calculate the kinetic stress tensor in the system averaged over the whole volume.
    Matrix3D getKineticStress() const;

    /// \brief Calculate the static stress tensor in the system averaged over the whole volume.
    Matrix3D getStaticStress() const;

    /// \brief Calculate the total stress tensor in the system averaged over the whole volume.
    Matrix3D getTotalStress() const;

    /*!
     * \brief Handles the removal of particles from the particleHandler
     */
    virtual void handleParticleRemoval(unsigned int id);

    /*!
     * \brief Handles the addition of particles to the particleHandler
     */
    virtual void handleParticleAddition(unsigned int id, BaseParticle* p);

    //functions that should only be used in the class definitions
protected:

    /*!
     * \brief Computes all the forces acting on the particles using the \ref BaseInteractable::setForce() and
     * \ref BaseInteractable::setTorque()
     */
    virtual void computeAllForces();

    /*!
     * \brief Computes the internal forces on particle i (internal in the sense that the
     *        sum over all these forces is zero i.e. fully modelled forces)
     */
    virtual void computeInternalForces(BaseParticle*);

    /*!
     * \brief Computes the forces between two particles (internal in the sense that
     *        the sum over all these forces is zero i.e. fully modelled forces)
     */
    virtual void computeInternalForce(BaseParticle*, BaseParticle*);

    /*!
     * \brief Computes the external forces, such as gravity, acting on particles.
     */
    virtual void computeExternalForces(BaseParticle*);

    /*!
     * \brief Computes the forces on the particles due to the walls (normals are outward normals)
     */
    virtual void computeForcesDueToWalls(BaseParticle*, BaseWall*);

    /*!
     * \brief A virtual function where the users can add extra code which is executed
     *  only when the code is restarted.
     */
    virtual void actionsOnRestart();

    /*!
     * \brief A virtual function. Allows one to carry out any operations before the start
     *        of the time loop.
     */
    virtual void actionsBeforeTimeLoop();

    /*!
     * \brief A virtual function that allows one to carry out hGrid operations before the
     *       start of the time loop.
     */
    virtual void hGridActionsBeforeTimeLoop();

    /*!
     * \brief A virtual function that allows one to set or execute hGrid parameters or
     *        operations before every simulation time step.
     */
    virtual void hGridActionsBeforeTimeStep();

    /*!
     * \brief A virtual function which allows to define operations to be executed before
     *        the new time step.
     */
    virtual void actionsBeforeTimeStep();

    /*!
     * \brief A virtual function which allows to define operations to be executed
     *       prior to the OMP force collect.
     * \details no implementation but can be overidden in its derived classes.
     * This function is called by DPMBase::computeAllForces. It may be used to add
     * additional user defined forces to any MercuryObject. Within this function,
     * OMP parallelization may be used, as it is called prior to sumForceTorqueOMP.
     */
    virtual void computeAdditionalForces() {}

    /*!
     * \brief A virtual function which allows to define operations to be executed after
     *       the solve().
     */
    virtual void actionsAfterSolve();

    /*!
     * \brief A virtual function which allows to define operations to be executed after
     *       time step.
     */
    virtual void actionsAfterTimeStep();

    void writeVTKFiles() const;

    /*!
     * \brief This function writes the location of the walls and particles in a format the
     *        XBalls program can read. For more information on the XBalls program, see \ref xballs.
     */
    virtual void outputXBallsData(std::ostream& os) const;

    /*!
     * \brief This function writes out the particle locations into an output stream in a
     *        format the XBalls program can read. For more information on the XBalls program, see \ref xballs.
     */
    virtual void outputXBallsDataParticle(unsigned int i, unsigned int format, std::ostream& os) const;

    /*!
     * \brief Writes a header with a certain format for ENE file.
     */
    virtual void writeEneHeader(std::ostream& os) const;

    /*!
     * \brief Writes a header with a certain format for FStat file.
     */
    virtual void writeFstatHeader(std::ostream& os) const;

    /*!
     * \brief Write the global kinetic, potential energy, etc. in the system.
     */
    virtual void writeEneTimeStep(std::ostream& os) const;

    // Functions for statistics
    /*!
     *
     */
    virtual void initialiseStatistics();

    /*!
     * \brief
     */
    virtual void outputStatistics();

    /*!
     *
     */
    void gatherContactStatistics();

    /*!
     * \brief
     */
    virtual void processStatistics(bool);

    /*!
     * \brief
     */
    virtual void finishStatistics();

    /*!
     * \brief Update particles' and walls' positions and velocities before force computation.
     * \details This is where the integration is done, at the moment it is
     *        velocity Verlet integration and is done before the forces are
     *        computed. See http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
     */
    virtual void integrateBeforeForceComputation();

    /*!
     * \brief Update particles' and walls' positions and velocities after force computation.
     */
    virtual void integrateAfterForceComputation();

    /*!
     * \brief There are a range of boundaries one could implement depending on ones' problem.
     *        This methods checks for interactions between particles and such range of boundaries.
     *        See BaseBoundary.h and all the boundaries in the Boundaries folder
     */
    virtual void checkInteractionWithBoundaries();

    /*!
     * \brief This function has to be called before integrateBeforeForceComputation.
     */
    virtual void hGridActionsBeforeIntegration();

    /*!
     * \brief This function has to be called after integrateBeforeForceComputation.
     */
    virtual void hGridActionsAfterIntegration();

    /*!
     * \brief Sets a number, <TT>n</TT>, of particles in the particleHandler as "fixed particles".
     */
    void setFixedParticles(unsigned int n);

    /*!
     * \brief Displays the current simulation time and the maximum simulation duration
     */
    virtual void printTime() const;

    /*!
     * \brief A virtual function for deciding whether to continue the simulation, based on a user-specified criterion.
     * \details Alternatively, the discontinueSolve function can be used as well, which can be called directly from any place in the code.
     */
    virtual bool continueSolve() const;

    /*!
     * \brief Displays the interaction details corresponding to the pointer objects
     *        in the interaction handler.
     */
    void outputInteractionDetails() const;

    /*!
     * \brief Checks whether the input variable "time" is the current time in the simulation
     */
    bool isTimeEqualTo(Mdouble time) const;

    /*!
     * \brief Removes periodic duplicate Particles
     */
    void removeDuplicatePeriodicParticles();

    /*!
     * \brief For simulations using periodic boundaries, checks and adds particles
     *        when necessary into the particle handler. See DPMBase.cc and PeriodicBoundary.cc
     *        for more details.
     */
    void checkAndDuplicatePeriodicParticles();

    /*!
     * \brief When the Verlet scheme updates the positions and velocities of particles,
     * ghost particles will need an update as wel. Their status will also be updated
     * accordingly.
     */
    void performGhostParticleUpdate();

    void deleteGhostParticles(std::set<BaseParticle*>& particlesToBeDeleted);

    void synchroniseParticle(BaseParticle*, unsigned fromProcessor = 0);

    /*!
     * \brief updates the final time-step velocity of the ghost particles
     */
    void performGhostVelocityUpdate();

    /*!
     * \brief This prevents the initialisation of the soft stop feature.
     */
    void disableSoftStop() { disableSoftStop_ = true; }

    /*!
     * \brief This discontinues the solve loop after the current time step is finished.
     * \details Use this to gracefully exit the simulation early, i.e. write the output files etc., by calling it from
     * any place in the code. It is an alternative to overriding the continueSolve function.
     */
    void discontinueSolve() { continueFlag_ = false; }

private:

    /*!
    * \brief signal handler function.
    */
    static void signalHandler(int signal);

    /*!
    * \brief function for setting sigaction constructor.
    */
    void setSoftStop();

    /**
     * Stores whether code should continue solving, see #discontinueSolve.
     */
    static volatile sig_atomic_t continueFlag_;

    /**
     * Stores whether the soft stop should not be initialised.
     */
    bool disableSoftStop_ = false;

    /*The number of openmp (symmetric multiprocessing threads)*/
    int numberOfOMPThreads_;

    /*!
     * \brief The dimensions of the simulation i.e. 2D or 3D
     */
    unsigned int systemDimensions_;

    /*!
     * \brief determines if 2D or 3D particle volume is used for mass calculations
     */
    unsigned int particleDimensions_;

    /*!
     * \brief Gravity vector
     */
    Vec3D gravity_;

    /*!
     * \brief This is a global background drag to help retard particles, or to do a simple fluid drag.
     * The force particleVeclocity*backgroundDrag is applied.
     */
     Mdouble backgroundDrag_;

    /*!
     * \brief Vector containing the number of domains in x-,y- and z-direction, required for parallel computations
     */
    std::vector<unsigned> numberOfDomains_;

    /*!
     * \brief These vectors are used for the XBalls domain, and occasionally people use it to add walls.
     */
    Vec3D min_;
    Vec3D max_;

    /*!
     * \brief Stores the current simulation time
     */
    Mdouble time_;

    /*!
     * \brief Stores the number of time steps
     */
    unsigned int numberOfTimeSteps_;

    /*!
     * \brief Stores the simulation time step
     */
    Mdouble timeStep_;

    /*!
     * \brief Stores the duration of the simulation
     */
    Mdouble timeMax_;

    /*!
     * \brief Previous versions of MercuryDPM had a different restart file format,
     *        the below member variable allows one to specify the version in order
     *        to choose between the available version support.
     */
    std::string restartVersion_;

    /*!
     * \brief A bool to check if the simulation was restarted or not, ie. if setupInitialConditionsShould be run and the fileCounters reset.
     */
    bool restarted_;

    /*!
     * \brief A flag to determine if the file has to be appended or not. See DPMBase::Solve()
     *        for example.
     */
    bool append_;

    /*!
     * \brief A flag to turn on/off particle rotation.
     * <tt>true</tt> will enable particle rotation.
     * <tt>false</tt> will disable particle rotation.
     */
    bool rotation_;

    /*!
     * \brief A flag to turn on/off the vtk writer for particles.
     */
    bool writeParticlesVTK_;

    bool writeSuperquadricParticlesVTK_;

    ParticleVtkWriter* vtkWriter_;

    WallVTKWriter wallVTKWriter_;

    InteractionVTKWriter interactionVTKWriter_;

    BoundaryVTKWriter boundaryVTKWriter_;

    WallDetailsVTKWriter wallDetailsVTKWriter_;

    //This is the private data that is only used by the xballs output

    /*!
     * \brief XBalls is a package to view the particle data. As an alternative MercuryDPM also supports ParaView.
     * The below variable is used to set the argument cmode in xballs script (see XBalls/xballs.txt)
     */
    int xBallsColourMode_;

    /*!
     * \brief sets the xballs argument vscale (see XBalls/xballs.txt)
     */
    Mdouble xBallsVectorScale_;

    /*!
     * \brief sets the xballs argument scale (see XBalls/xballs.txt)
     */
    Mdouble xBallsScale_;

    /*!
     * \brief A string of additional arguments for xballs can be specified (see XBalls/xballs.txt). e.g. "-solidf -v0"
     */
    std::string xBallsAdditionalArguments_;

    /*!
     * \brief This stores the run number for saving
     */
    int runNumber_;

    /*!
     * \brief the name of the problem, used, e.g., for the output files
     */
    std::string name_;

    // defines a Macro for creating an instance of class PossibleContactList. See PossibleContactList.h
#ifdef CONTACT_LIST_HGRID
    PossibleContactList possibleContactList;
#endif

    /*!
     * \brief Determines if the last column of the data file is interpreted as the info parameter during restart
     */
    bool readSpeciesFromDataFile_;

    /*!
     * \brief number of elements to write to a screen
     */
    int nToWrite_;

public:
    /*!
     * \brief A handler to that stores the species type i.e. LinearViscoelasticSpecies, etc.
     */
    SpeciesHandler speciesHandler;

    /*!
     * \brief This is a random generator, often used for setting up the initial conditions etc...
     */
    RNG random;

    /*!
     * \brief An object of the class ParticleHandler, contains the pointers to all the particles created.
     */
    ParticleHandler particleHandler;

    /*!
     * \brief Fake particleHandler created by Paolo needed temporary by just Paolo.
     */
    ParticleHandler paoloParticleHandler;

    /*!
     * \brief An object of the class WallHandler. Contains pointers to all the walls created.
     */
    WallHandler wallHandler;

    /*!
     * \brief An object of the class BoundaryHandler which concerns insertion and deletion of particles into or from regions.
     */
    BoundaryHandler boundaryHandler;

    /*!
     * \brief Internal handler that deals with periodic boundaries, especially in a parallel build
     */
    PeriodicBoundaryHandler periodicBoundaryHandler;

    /*!
     * \brief An object of the class DomainHandler which deals with parallel code
     */
    DomainHandler domainHandler;

    /*!
     * \brief An object of the class InteractionHandler
     */
    InteractionHandler interactionHandler;


    /*!
     * \brief Object of the class cgHandler
     */
    CGHandler cgHandler;

    /*!
     * \brief An instance of class File to handle in- and output into a .data file
     */
    File dataFile;

    /*!
     * \brief An instance of class File to handle in- and output into a .fstat file
     */
    File fStatFile;

    /*!
     * \brief An instance of class File to handle in- and output into a .ene file
     */
    File eneFile;

    /*!
     * \brief An instance of class File to handle in- and output into a .restart file
     */
    File restartFile;

    /*!
     * \brief An instance of class File to handle in- and output into a .stat file
     */
    File statFile;

    /*!
     * \brief File class to handle in- and output into .interactions file. This file hold
     * information about interactions.
     */
    File interactionFile;

    /*!
     * \brief record when the simulation started
     */
    Time clock_;

    void writePythonFileForVTKVisualisation() const;

    void setWritePythonFileForVTKVisualisation(bool forceWritePythonFileForVTKVisualisation);
    bool getWritePythonFileForVTKVisualisation() const;

    WallVTKWriter &getWallVTKWriter();

private:
    bool forceWritePythonFileForVTKVisualisation_ { true };
};

#endif
