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


#include "SpeciesHandler.h"
#include "DPMBase.h"

#include "Species/LinearViscoelasticSpecies.h"
#include "Species/LinearPlasticViscoelasticSpecies.h"
#include "Species/SinterSpecies.h"
#include "Species/SinterReversibleAdhesiveSpecies.h"
#include <Species/SinterFrictionSpecies.h>
#include "Species/SinterFrictionReversibleAdhesiveSpecies.h"
#include "Species/HertzianSinterSpecies.h"
#include "Species/HertzianSinterFrictionSpecies.h"
#include "Species/HertzianSinterSlidingFrictionSpecies.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Species/LinearPlasticViscoelasticFrictionSpecies.h"
#include "Species/LinearPlasticViscoelasticFrictionBondedSpecies.h"
#include "Species/LinearPlasticViscoelasticSlidingFrictionSpecies.h"
#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include "Species/HertzianViscoelasticMindlinRollingTorsionSpecies.h"
#include "Species/HertzianViscoelasticMindlinLiquidBridgeClassicalWilletSpecies.h"
#include "Species/HertzianViscoelasticMindlinLiquidBridgeBagheriSpecies.h"
#include "Species/SPHSpecies.h"

#include "Species/LinearViscoelasticBondedSpecies.h"
#include "Species/LinearViscoelasticSlidingFrictionBondedSpecies.h"
#include "Species/LinearViscoelasticIrreversibleAdhesiveSpecies.h"
#include "Species/LinearPlasticViscoelasticIrreversibleAdhesiveSpecies.h"
#include "Species/LinearViscoelasticFrictionIrreversibleAdhesiveSpecies.h"
#include "Species/LinearViscoelasticSlidingFrictionIrreversibleAdhesiveSpecies.h"
#include "Species/LinearPlasticViscoelasticFrictionIrreversibleAdhesiveSpecies.h"
#include "Species/LinearPlasticViscoelasticSlidingFrictionIrreversibleAdhesiveSpecies.h"
#include "Species/HertzianViscoelasticFrictionChargedBondedSpecies.h"
#include "Species/HertzianViscoelasticMindlinLiquidMigrationWilletSpecies.h"

#include "Species/LinearViscoelasticFrictionChargedBondedSpecies.h"
#include "Species/LinearPlasticViscoelasticSlidingFrictionLiquidMigrationWilletSpecies.h"
#include "Species/LinearPlasticViscoelasticFrictionLiquidBridgeWilletSpecies.h"

#include "Species/LinearViscoelasticReversibleAdhesiveSpecies.h"
#include "Species/LinearPlasticViscoelasticReversibleAdhesiveSpecies.h"
#include "Species/LinearViscoelasticFrictionReversibleAdhesiveSpecies.h"

#include "Species/LinearViscoelasticFrictionJKRAdhesiveSpecies.h"
#include "Species/MeltableSpecies.h"
#include "Species/MeltableFrictionSpecies.h"

#include "Species/LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies.h"
#include "Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h"
#include "Species/LinearPlasticViscoelasticSlidingFrictionReversibleAdhesiveSpecies.h"


#include "Species/LinearViscoelasticFrictionLiquidBridgeWilletSpecies.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
#include "Species/LinearPlasticViscoelasticFrictionLiquidMigrationWilletSpecies.h"
#include "Species/LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationLSSpecies.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletViscousSpecies.h"
#include "Species/HertzianViscoelasticSlidingFrictionParhamiMcMeekingSinterSpecies.h"
#include "Species/NormalForceSpecies/ThermalSpecies.h"
#include "Species/HeatFluidCoupledLinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
#include "Species/HeatFluidCoupledLinearViscoelasticFrictionSolidifyingLiquidMigrationWilletSpecies.h"

#include "Species/NormalForceSpecies/SinterLinNormalSpecies.h"
#include "Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h"

#include "Species/HertzianBSHPViscoelasticSpecies.h"
#include "Species/HertzianBSHPViscoelasticMindlinSpecies.h"
/*!
 * \details Constructor of the SpeciesHandler class. It creates an empty SpeciesHandler.
 */
SpeciesHandler::SpeciesHandler()
{
    logger(DEBUG, "SpeciesHandler::SpeciesHandler() finished");
}

/*!
 * \param[in] other The SpeciesHandler that has to be copied.
 * \details This is not a copy constructor! This constructor copies only all
 *          BaseSpecies and MixedSpecies and copies the pointer to the DPMBase.
 *          It sets all other data members to 0 or nullptr.
 */
SpeciesHandler::SpeciesHandler(const SpeciesHandler& other)
{
    clear();
    setDPMBase(other.getDPMBase());
    copyContentsFromOtherHandler(other);
    for (BaseSpecies* mixSpec : other.mixedObjects_)
    {
        mixedObjects_.push_back(mixSpec->copy());
        mixedObjects_.back()->setHandler(this);
    }
    logger(DEBUG, "SpeciesHandler::SpeciesHandler(const SpeciesHandler &other) finished");
}

/*!
 * \param[in] rhs The BoundaryHandler on the right hand side of the assignment.
 * \return The SpeciesHandler that is a copy of the input SpeciesHandler rhs.
 * \details This is not a copy assignment operator! It copies only all
 *          BaseSpecies and MixedSpecies and copies the pointer to the DPMBase.
 *          It sets all other data members to 0 or nullptr.
 */
SpeciesHandler& SpeciesHandler::operator=(const SpeciesHandler& rhs)
{
    if (this != &rhs)
    {
        clear();
        copyContentsFromOtherHandler(rhs);
        for (BaseSpecies* mixSpec : mixedObjects_)
        {
            delete mixSpec;
        }
        mixedObjects_.clear();
        for (BaseSpecies* mixSpec : rhs.mixedObjects_)
        {
            mixedObjects_.push_back(mixSpec->copy());
            mixedObjects_.back()->setHandler(this);
        }
    }

    logger(DEBUG, "SpeciesHandler SpeciesHandler::operator =(const SpeciesHandler& rhs)");
    return *this;
}

/*!
 * \details Destructor: first destroys the objects of the BaseHandler, then destroys the mixedObjects
 * \todo TW Note: deleting the species does not delete the particles and walls of this species.
 */
SpeciesHandler::~SpeciesHandler()
{
    clear(); //this deletes everything that comes from the BaseHandler.
    /*for (BaseSpecies* o : mixedObjects_)
    {
        delete o;
    }*/
    for (BaseSpecies* mixedObject : mixedObjects_)
    {
        delete mixedObject;
    }
    mixedObjects_.clear();
    logger(DEBUG, "SpeciesHandler::~SpeciesHandler() finished");
}

/*!
 * \param[in] is The input stream from which the information is read.
 * \details First determine the type of the object we want to read, then read
 * the actual object. After that, clear the mixed objects and read the mixed objects.
 */
void SpeciesHandler::readAndAddObject(std::istream& is)
{
    std::string type;
    is >> type;
    logger(DEBUG, "SpeciesHandler::readAndAddObject(is): reading type %.", type);
    if (type == "LinearViscoelasticSpecies")
    {
        LinearViscoelasticSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "SPHSpecies")
    {
        SPHSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticSpecies")
    {
        LinearPlasticViscoelasticSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "SinterSpecies")
    {
        SinterSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "SinterReversibleAdhesiveSpecies")
    {
        SinterReversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "SinterFrictionReversibleAdhesiveSpecies")
    {
        SinterFrictionReversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "SinterFrictionSpecies")
    {
        SinterFrictionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianSinterSpecies")
    {
        HertzianSinterSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianSinterFrictionSpecies")
    {
        HertzianSinterFrictionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianSinterSlidingFrictionSpecies")
    {
        HertzianSinterSlidingFrictionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticSlidingFrictionSpecies")
    {
        LinearViscoelasticSlidingFrictionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticSlidingFrictionSpecies")
    {
        LinearPlasticViscoelasticSlidingFrictionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionSpecies")
    {
        LinearViscoelasticFrictionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticFrictionSpecies")
    {
        LinearPlasticViscoelasticFrictionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticFrictionBondedSpecies")
    {
        LinearPlasticViscoelasticFrictionBondedSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticIrreversibleAdhesiveSpecies")
    {
        LinearViscoelasticIrreversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticBondedSpecies")
    {
        LinearViscoelasticBondedSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticIrreversibleAdhesiveSpecies")
    {
        LinearPlasticViscoelasticIrreversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticSlidingFrictionIrreversibleAdhesiveSpecies")
    {
        LinearViscoelasticSlidingFrictionIrreversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticSlidingFrictionIrreversibleAdhesiveSpecies")
    {
        LinearPlasticViscoelasticSlidingFrictionIrreversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionIrreversibleAdhesiveSpecies")
    {
        LinearViscoelasticFrictionIrreversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticFrictionIrreversibleAdhesiveSpecies")
    {
        LinearPlasticViscoelasticFrictionIrreversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticReversibleAdhesiveSpecies")
    {
        LinearViscoelasticReversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticReversibleAdhesiveSpecies")
    {
        LinearPlasticViscoelasticReversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies")
    {
        LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticSlidingFrictionReversibleAdhesiveSpecies")
    {
        LinearPlasticViscoelasticSlidingFrictionReversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionReversibleAdhesiveSpecies")
    {
        LinearViscoelasticFrictionReversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionJKRAdhesiveSpecies")
    {
        LinearViscoelasticFrictionJKRAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "MeltableFrictionSpecies")
    {
        MeltableFrictionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "MeltableSpecies")
    {
        MeltableSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies")
    {
        LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionLiquidBridgeWilletSpecies")
    {
        LinearViscoelasticFrictionLiquidBridgeWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionLiquidMigrationWilletSpecies")
    {
        LinearViscoelasticFrictionLiquidMigrationWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianViscoelasticMindlinLiquidMigrationWilletSpecies")
    {
        HertzianViscoelasticMindlinLiquidMigrationWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianViscoelasticMindlinLiquidBridgeClassicalWilletSpecies")
    {
        HertzianViscoelasticMindlinLiquidBridgeClassicalWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianViscoelasticMindlinLiquidBridgeBagheriSpecies")
    {
        HertzianViscoelasticMindlinLiquidBridgeBagheriSpecies species;
        is >> species;
        copyAndAddObject(species);
    }    
    else if (type == "LinearPlasticViscoelasticFrictionLiquidMigrationWilletSpecies")
    {
        LinearPlasticViscoelasticFrictionLiquidMigrationWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies")
    {
        LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionLiquidMigrationLSSpecies")
    {
        LinearViscoelasticFrictionLiquidMigrationLSSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionLiquidMigrationWilletViscousSpecies")
    {
        LinearViscoelasticFrictionLiquidMigrationWilletViscousSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticSlidingFrictionLiquidMigrationWilletSpecies")
    {
        LinearPlasticViscoelasticSlidingFrictionLiquidMigrationWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HeatFluidCoupledLinearViscoelasticFrictionLiquidMigrationWilletSpecies")
    {
        HeatFluidCoupledLinearViscoelasticFrictionLiquidMigrationWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HeatFluidCoupledLinearViscoelasticFrictionSolidifyingLiquidMigrationWilletSpecies")
    {
        HeatFluidCoupledLinearViscoelasticFrictionSolidifyingLiquidMigrationWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianViscoelasticMindlinSpecies")
    {
        HertzianViscoelasticMindlinSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianViscoelasticMindlinRollingTorsionSpecies" ||
             type == "HertzianViscoelasticFrictionSpecies")
    {
        HertzianViscoelasticMindlinRollingTorsionSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianViscoelasticFrictionChargedBondedSpecies")
    {
        HertzianViscoelasticFrictionChargedBondedSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianBSHPViscoelasticMindlinSpecies")
    {
        HertzianBSHPViscoelasticMindlinSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "HertzianBSHPViscoelasticSpecies")
    {
        HertzianBSHPViscoelasticSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearViscoelasticFrictionChargedBondedSpecies")
    {
        LinearViscoelasticFrictionChargedBondedSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "ThermalSinterSlidingFrictionSpecies")
    {
        Species<ThermalSpecies<SinterNormalSpecies>, SlidingFrictionSpecies> species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "ThermalSinterFrictionSpecies")
    {
        Species<ThermalSpecies<SinterNormalSpecies>, FrictionSpecies> species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "LinearPlasticViscoelasticFrictionLiquidBridgeWilletSpecies")
    {
        LinearPlasticViscoelasticFrictionLiquidBridgeWilletSpecies species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "ThermalSinterLinFrictionReversibleAdhesiveSpecies")
    {
        Species<ThermalSpecies<SinterLinNormalSpecies>, FrictionSpecies, ReversibleAdhesiveSpecies> species;
        is >> species;
        copyAndAddObject(species);
    }
    else if (type == "k") //for backwards compatibility
    {
        addObject(readOldObject(is));
    }
    else
    {
        logger(WARN,
               "Species type % not understood in restart file: You need to add this species to SpeciesHandler::readObject.",
               type);
        std::stringstream line;
        helpers::getLineFromStringStream(is, line);
        copyAndAddObject(LinearViscoelasticSpecies());
    }

    //remove the default mixed species
    for (unsigned int i = 0; i + 1 < getNumberOfObjects(); i++)
    {
        ///\todo TW why does deleting these objects create a segmentation fault
        ///How do you create the segmentation fault? \author weinhartt
        ///\todo IFCD how does the numbering of mixedSpecies_ work?
        ///the numbering of mixed species is 01, 02, 12, 03, 13, 23, 04, 14, 24, 34,
        ///i.e. if you add the n-th ParticleSpecies, then you have to add n-1 MixedSpecies.
        ///So here I remove the last n-1 MixedSpecies and add n-1 new ones. \author weinhartt
        delete mixedObjects_.back();
        mixedObjects_.pop_back();
    }

    //Read the mixed species.
    for (unsigned int i = 0; i + 1 < getNumberOfObjects(); i++)
    {
        is >> type;
        if (type == "LinearViscoelasticMixedSpecies")
        {
            LinearViscoelasticMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticMixedSpecies")
        {
            LinearPlasticViscoelasticMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "SinterMixedSpecies")
        {
            SinterMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "SinterReversibleAdhesiveMixedSpecies")
        {
            SinterReversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "SinterFrictionReversibleAdhesiveMixedSpecies")
        {
            SinterFrictionReversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "SinterFrictionMixedSpecies")
        {
            SinterFrictionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianSinterMixedSpecies")
        {
            HertzianSinterMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianSinterFrictionMixedSpecies")
        {
            HertzianSinterFrictionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianSinterSlidingFrictionMixedSpecies")
        {
            HertzianSinterSlidingFrictionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticSlidingFrictionMixedSpecies")
        {
            LinearViscoelasticSlidingFrictionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticSlidingFrictionMixedSpecies")
        {
            LinearPlasticViscoelasticSlidingFrictionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionMixedSpecies")
        {
            LinearViscoelasticFrictionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticFrictionMixedSpecies")
        {
            LinearPlasticViscoelasticFrictionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticFrictionBondedMixedSpecies")
        {
            LinearPlasticViscoelasticFrictionBondedMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticIrreversibleAdhesiveMixedSpecies")
        {
            LinearViscoelasticIrreversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticBondedMixedSpecies")
        {
            LinearViscoelasticBondedMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticIrreversibleAdhesiveMixedSpecies")
        {
            LinearPlasticViscoelasticIrreversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticSlidingFrictionIrreversibleAdhesiveMixedSpecies")
        {
            LinearViscoelasticSlidingFrictionIrreversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticSlidingFrictionIrreversibleAdhesiveMixedSpecies")
        {
            LinearPlasticViscoelasticSlidingFrictionIrreversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionIrreversibleAdhesiveMixedSpecies")
        {
            LinearViscoelasticFrictionIrreversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticFrictionIrreversibleAdhesiveMixedSpecies")
        {
            LinearPlasticViscoelasticFrictionIrreversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticSlidingFrictionLiquidMigrationWilletMixedSpecies")
        {
            LinearPlasticViscoelasticSlidingFrictionLiquidMigrationWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticFrictionLiquidBridgeWilletMixedSpecies")
        {
            LinearPlasticViscoelasticFrictionLiquidBridgeWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianViscoelasticMindlinLiquidMigrationWilletMixedSpecies")
        {
            HertzianViscoelasticMindlinLiquidMigrationWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianViscoelasticMindlinLiquidBridgeClassicalWilletMixedSpecies")
        {
            HertzianViscoelasticMindlinLiquidBridgeClassicalWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }  
        else if (type == "HertzianViscoelasticMindlinLiquidBridgeBagheriMixedSpecies")
        {
            HertzianViscoelasticMindlinLiquidBridgeBagheriMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }               
        else if (type == "HeatFluidCoupledLinearViscoelasticFrictionLiquidMigrationWilletMixedSpecies")
        {
            HeatFluidCoupledLinearViscoelasticFrictionLiquidMigrationWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HeatFluidCoupledLinearViscoelasticFrictionSolidifyingLiquidMigrationWilletMixedSpecies")
        {
            HeatFluidCoupledLinearViscoelasticFrictionSolidifyingLiquidMigrationWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticReversibleAdhesiveMixedSpecies")
        {
            LinearViscoelasticReversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticReversibleAdhesiveMixedSpecies")
        {
            LinearPlasticViscoelasticReversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticSlidingFrictionReversibleAdhesiveMixedSpecies")
        {
            LinearViscoelasticSlidingFrictionReversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticSlidingFrictionReversibleAdhesiveMixedSpecies")
        {
            LinearPlasticViscoelasticSlidingFrictionReversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionReversibleAdhesiveMixedSpecies")
        {
            LinearViscoelasticFrictionReversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticFrictionReversibleAdhesiveMixedSpecies")
        {
            LinearPlasticViscoelasticFrictionReversibleAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionJKRAdhesiveMixedSpecies")
        {
            LinearViscoelasticFrictionJKRAdhesiveMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "MeltableFrictionMixedSpecies")
        {
            MeltableFrictionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionLiquidBridgeWilletMixedSpecies")
        {
            LinearViscoelasticFrictionLiquidBridgeWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionLiquidMigrationWilletMixedSpecies")
        {
            LinearViscoelasticFrictionLiquidMigrationWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearPlasticViscoelasticFrictionLiquidMigrationWilletMixedSpecies")
        {
            LinearPlasticViscoelasticFrictionLiquidMigrationWilletMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticSlidingFrictionLiquidMigrationLSMixedSpecies")
        {
            LinearViscoelasticSlidingFrictionLiquidMigrationLSMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionLiquidMigrationLSMixedSpecies")
        {
            LinearViscoelasticFrictionLiquidMigrationLSMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionLiquidMigrationWilletViscousMixedSpecies")
        {
            LinearViscoelasticFrictionLiquidMigrationWilletViscousMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianViscoelasticMindlinMixedSpecies")
        {
            HertzianViscoelasticMindlinMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianViscoelasticMindlinRollingTorsionMixedSpecies" ||
                 type == "HertzianViscoelasticFrictionMixedSpecies")
        {
            HertzianViscoelasticMindlinRollingTorsionMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianViscoelasticFrictionChargedBondedMixedSpecies")
        {
            HertzianViscoelasticFrictionChargedBondedMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianBSHPViscoelasticMindlinMixedSpecies")
        {
            HertzianBSHPViscoelasticMindlinMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "HertzianBSHPViscoelasticMixedSpecies")
        {
            HertzianBSHPViscoelasticMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "LinearViscoelasticFrictionChargedBondedMixedSpecies")
        {
            LinearViscoelasticFrictionChargedBondedMixedSpecies species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "ThermalSinterSlidingFrictionMixedSpecies")
        {
            MixedSpecies<ThermalSpecies<SinterNormalSpecies>, SlidingFrictionSpecies> species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "ThermalSinterFrictionMixedSpecies")
        {
            MixedSpecies<ThermalSpecies<SinterNormalSpecies>, FrictionSpecies> species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "ThermalSinterLinFrictionReversibleAdhesiveMixedSpecies")
        {
            MixedSpecies<ThermalSpecies<SinterLinNormalSpecies>, FrictionSpecies, ReversibleAdhesiveSpecies> species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else if (type == "ThermalSinterFrictionMixedSpecies")
        {
            MixedSpecies<ThermalSpecies<SinterNormalSpecies>, FrictionSpecies> species;
            is >> species;
            mixedObjects_.push_back(species.copy());
        }
        else
        {
            logger(WARN,
                   "Species type % not understood in restart file: You need to add this species to SpeciesHandler::readMixedObject.",
                   type);
            std::stringstream line;
            helpers::getLineFromStringStream(is, line);
            LinearViscoelasticMixedSpecies species;
            mixedObjects_.push_back(species.copy());
        }
    }
}

/*!
 * \param[in] is The input stream from which the information is read.
 * \return A pointer to the ParticleSpecies that has just been read.
 * \details To read the old object, we first make a stringstream of the line that
 * describes this ParticleSpecies. After that, we read the properties one by one,
 * first the stiffness and after that the other properties. We stop when we either
 * reach the end of the file(eof) or if a string is not recognized as a property.
 */
ParticleSpecies* SpeciesHandler::readOldObject(std::istream& is)
{
    //read in next line
    std::stringstream line;
    helpers::getLineFromStringStream(is, line);

    //read each property
    std::string property;
    unsigned int particleDimension = 0;
    Mdouble density = 0.0, stiffness = 0.0, dissipation = 0.0, slidingFrictionCoefficient = 0.0, slidingFrictionCoefficientStatic = 0.0, slidingStiffness = 0.0, slidingDissipation = 0.0;
    line >> stiffness;
    while (true)
    {
        line >> property;
        if (property == "disp")
            line >> dissipation;
        else if (property == "rho")
            line >> density;
        else if (property == "kt")
            line >> slidingStiffness;
        else if (property == "dispt")
            line >> slidingDissipation;
        else if (property == "mu")
            line >> slidingFrictionCoefficient;
        else if (property == "mus")
            line >> slidingFrictionCoefficientStatic;
        else if (property == "dim_particle")
        {
            line >> particleDimension;
            getDPMBase()->setParticleDimensions(particleDimension);
        }
        else if (property == "(mixed)")
        {
            density = 0;
        }
        else
        {
            logger(WARN, "Warning: % is not a species property", property);
            break;
        }
        if (line.eof())
            break;
    }

    //create the correct species
    if (slidingFrictionCoefficient == 0.0)
    {
        LinearViscoelasticSpecies* species = new LinearViscoelasticSpecies;
        species->setDensity(density);
        species->setStiffness(stiffness);
        species->setDissipation(dissipation);
        return species;
    }
    else
    {
        LinearViscoelasticSlidingFrictionSpecies* species = new LinearViscoelasticSlidingFrictionSpecies;
        species->setDensity(density);
        species->setStiffness(stiffness);
        species->setDissipation(dissipation);
        species->setSlidingStiffness(slidingStiffness);
        species->setSlidingDissipation(slidingDissipation);
        species->setSlidingFrictionCoefficient(slidingFrictionCoefficient);
        if (slidingFrictionCoefficientStatic == 0.0)
            slidingFrictionCoefficientStatic = slidingFrictionCoefficient;
        species->setSlidingFrictionCoefficientStatic(slidingFrictionCoefficientStatic);
        return species;
    }
}

/*!
 * \param[in] id1 Id of the first species.
 * \param[in] id2 Id of the second species.
 * \return An unsigned integer that denotes the Id of the mixed species.
 * \details The numbering of the mixed species is 0-1,  0-2, 1-2,  0-3, 1-3, 2-3,  0-4, 1-4, 2-4, 3-4, ...,
 * where each pair of numbers a and b denotes the mixed species between ParticleSpecies a and b.
 * Thus, first compute which id has a higher value, the id of the mixed species
 * is then given by (maxId*(maxId-1))/2 + minId.
 */
unsigned int SpeciesHandler::getMixedId(const unsigned int id1, const unsigned int id2) const
{
    unsigned int maxId = std::max(id1, id2);
    return (maxId * (maxId - 1)) / 2 + std::min(id1, id2);
}

/*!
 * \return A reference to the vector of pointers of all the mixedObjects.
 */
const std::vector<BaseSpecies*>& SpeciesHandler::getMixedObjects() const
{
    return mixedObjects_;
}

/*!
 * \param[in] id1 Id of the first BaseSpecies.
 * \param[in] id2 Id of the second BaseSpecies.
 * \return A pointer to an object that is a MixedSpecies of both input Species.
 * \todo This function should probably be made private. The user should use the function
 * SpeciesHandler::getMixedObject(const U* S, const U* T), which deals with pointers.
 */
BaseSpecies* SpeciesHandler::getMixedObject(const unsigned int id1, const unsigned int id2)
{
    if (id1 == id2)
    {
        return getObject(id1);
    }
    else
    {
        const unsigned int mixedId = getMixedId(id1, id2);
        if (mixedId>=mixedObjects_.size())
        {
            logger(WARN,
                   "In: Object* SpeciesHandler::getMixedObject(const unsigned int id) const. No Object exist with index %, number of objects is %",
                   std::max(id1, id2), getNumberOfObjects());
            return nullptr;
        }
        else
        {
            return mixedObjects_[mixedId];
        }
    }
}

/*!
 * \details Updates mixed species by calling the function mixAll() with the
 * corresponding original species.
 * \sa mixAll(Object1, Object2) which creates sensible default values for mixed species based on the two species.
 * Remember mixedSpecies properties can be set independently if required.
 *
 * This may (should) be called, if a value is changed
 * in the original species and the mixed species need to be aware of this change and update by the default rule.
 * Note this only make sense if you are using the default mixed species and not a custom implementation.
 */
void SpeciesHandler::updateMixedObjects()
{
    for (unsigned int index1 = 0; index1 < getNumberOfObjects(); ++index1)
    {
        for (unsigned int index2 = index1+1; index2 < getNumberOfObjects(); ++index2)
        {
            mixedObjects_[getMixedId(index1, index2)]->mixAll(getObject(index1), getObject(index2));
        }
    }
}

/*!
 * \param[in] S A pointer to the ParticleSpecies that has to be added.
 * \details First, add the ParticleSpecies to the vector of ParticleSpecies (object_),
 * then construct all MixedSpecies. Tell the ParticleSpecies that this is its
 * handler and compute all masses and whether it should use angular degrees of freedom.
 *
 * Note: The MixedSpecies objects are initialized with
 * averaged values from both species: e.g., the mixedSpecies between Species A
 * and B will have a stiffness \$fk=(1/k_a+1/k_b)^{-1}\$f, you have to change
 * the MixedSpecies properties if you don't like these defaults.
 */
void SpeciesHandler::addObject(ParticleSpecies* const S)
{
    BaseHandler<ParticleSpecies>::addObject(S);
    //logger(INFO, "Part / Mix: % / %", objects_.size(), mixedObjects_.size());
    ///\todo TW don't put logger messages that only make sense for one application!
    S->setHandler(this);
    for (unsigned int id = 0; id + 1 < getNumberOfObjects(); ++id)
    {
        mixedObjects_.push_back(S->copyMixed());
        mixedObjects_.back()->setIndex(id);
        mixedObjects_.back()->setId(getNumberOfObjects() - 1);
        mixedObjects_.back()->mixAll(S, getObject(id));
    }
    getDPMBase()->particleHandler.computeAllMasses(S->getIndex());
    getDPMBase()->setRotation(useAngularDOFs());
    S->setMaxInteractionDistance(S->getInteractionDistance());
}

/*!
 * \param[in] index The index of which ParticleSpecies has to be removed from this
 *                  ParticleHandler.
 * \details         The ParticleSpecies with index is removed and the last
 *                  ParticleSpecies in the vector is moved to its position.
 *                  It also removes all mixed species for this ParticleSpecies.
 */
void SpeciesHandler::removeObject(unsigned const int index)
{
    BaseHandler<ParticleSpecies>::removeObject(index);
    for (unsigned int index2 = 0; index2 < getNumberOfObjects(); ++index2)
    {
        mixedObjects_.erase(mixedObjects_.begin() + getMixedId(index, index2));
    }
    getDPMBase()->setRotation(useAngularDOFs());
}

/*!
 * \details First write "Species" and the amount of species in this handler,
 * then write all ParticleSpecies and MixedSpecies.
 */
void SpeciesHandler::write(std::ostream& os) const
{
    os << "Species " << getNumberOfObjects() << "\n";
    unsigned idMixed = 0;
    for (const ParticleSpecies* species : objects_)
    {
        os << *species << "\n";
        for (unsigned int id2 = 0; id2 < species->getIndex(); id2++)
        {
            os << *mixedObjects_[idMixed] << "\n";
            idMixed++;
        }
    }
}

/*!
 * \return The string "SpeciesHandler"
 */
std::string SpeciesHandler::getName() const
{
    return "SpeciesHandler";
}

/*!
 * \return The boolean which says whether or not AnuglarDOFs must be used in this handler.
 */
bool SpeciesHandler::useAngularDOFs()
{
    for (unsigned int i = 0; i < getSize(); i++)
    {
        if (getObject(i)->getUseAngularDOFs())
            return true;
        for (unsigned int j = 0; j + 1 < i; j++)
        {
            if (getMixedObject(i, j)->getUseAngularDOFs())
                return true;
        }
    }
    return false;
}
