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

#ifndef MECURYDPM_PSD_H
#define MECURYDPM_PSD_H

#include <fstream>
#include <Logger.h>
#include <vector>
#include <Tools/csvReader.h>
#include <random>
#include <utility>
#include <algorithm>
#include "iostream"
#include "ExtendedMath.h"
#include "RNG.h"
#include "Math/DistributionElements.h"

using mathsFunc::square;
using mathsFunc::cubic;

/*!
 * \brief Contains a vector with radii and probabilities of a user defined particle size distribution (PSD)
 *
 * \details Stores radii and probabilities of a particle size distribution (PSD) in a vector of type DistributionElements and
 * converts them to a PSD which can be used by insertionBoundaries which insert particles into a simulation:
 *
 * Cumulative distribution function (CDF): gives percentage of particles p_i whose radius r is less than
 * a certain radius r_i. It requires p_0 = 0 and p_end = 1. The CDF's probabilities can be number, length, area or
 * volume based. The cumulative number distribution function (CUMULATIVE_NUMBER_DISTRIBUTION) is used as PSD in MercuryDPM, as it can be
 * interpreted as the probability that a particle's radius is less than r_i:
 *     CDF(r<r_i) = p_i.
 *
 * Probability density function (PDF): p_i is the percentage of particles whose radius is between r_i-1 and
 * r_i. It requires p_0=0 and sum(p_i)=1. The PDF's probabilities can also be number, length, area or volume based.
 * PDF's are utilized to convert any type of PDF to a probability number density function (PROBABILITYDENSITY_NUMBER_DISTRIBUTION) which in a next step
 * are converted to the default CUMULATIVE_NUMBER_DISTRIBUTION.
 *
 * Sieve data: p_i is the percentage of particles whose radius is between r_i and r_i+1. It requires sum(p_i)=1.
 * relation to PDF: p_i = p(PDF)_i-1. Sieve data is not yet used in this class.
 *
 * Default distribution is the CUMULATIVE_NUMBER_DISTRIBUTION.
 */
class PSD
{

public:
    /*!
     * \brief Enum class which stores the possible types of CDFs and PDFs.
     * Particle size distributions can be represented by different probabilities based on number of particles, length
     * of particles, surface area of particles or volume of particles.
     */
    enum class TYPE
    {
        CUMULATIVE_NUMBER_DISTRIBUTION,
        CUMULATIVE_LENGTH_DISTRIBUTION,
        CUMULATIVE_VOLUME_DISTRIBUTION,
        CUMULATIVE_AREA_DISTRIBUTION,
        PROBABILITYDENSITY_NUMBER_DISTRIBUTION,
        PROBABILITYDENSITY_LENGTH_DISTRIBUTION,
        PROBABILITYDENSITY_AREA_DISTRIBUTION,
        PROBABILITYDENSITY_VOLUME_DISTRIBUTION
    };
    
    /*!
     * \brief Constructor; sets everything to 0 or default.
     */
    PSD();

    /*!
     * \brief Copy constructor with deep copy.
     */
    PSD(const PSD& other);

    /*!
     * \brief Destructor; default destructor.
     */
    ~PSD();

    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    PSD* copy() const;

    /*!
     * \brief Prints radii and probabilities of the PSD vector.
     */
    void printPSD() const;

    /*!
     * \brief Draw a sample radius from a CUMULATIVE_NUMBER_DISTRIBUTION.
     */
    Mdouble drawSample();

    /*!
     * \brief Draw sample radius manually per size class and check the volumeAllowed of each size class to insert the
     * PSD as accurate as possible.
     */
    Mdouble insertManuallyByVolume(Mdouble volume);

    /*!
     * \brief Validates if a CDF starts with zero and adds up to unity.
    */
    void validateCumulativeDistribution();

    /*!
     * \brief Validates if the integral of the PDF equals to unity.
     */
    void validateProbabilityDensityDistribution();

    /*!
     * \brief Sets the PSD from a vector with DistributionElements.
     */
    void setPSDFromVector(std::vector<DistributionElements> psd, TYPE PSDType);

    /*!
     * \brief read in the PSD vector with probabilities and radii saved in a .csv file.
     */
    void setPSDFromCSV(const std::string& fileName, TYPE PSDType, bool headings = false, Mdouble
    unitScalingFactorRadii = 1.0);

    /*!
     * \brief create a PSD vector for a uniform distribution.
     */
    void setDistributionUniform(Mdouble radMin, Mdouble radMax, int numberOfBins);

    /*!
     * \brief create a PSD vector for a normal distribution.
     */
    void setDistributionNormal(Mdouble mean, Mdouble standardDeviation, int numberOfBins);
    
    static PSD getDistributionNormal(Mdouble mean, Mdouble standardDeviation, int numberOfBins) {
        PSD psd;
        psd.setDistributionNormal(mean, standardDeviation, numberOfBins);
        return psd;
    }

    /*!
     * \brief create a PSD vector for a normal distribution.
     */
    void setDistributionLogNormal(Mdouble mean, Mdouble standardDeviation, int numberOfBins);

    static PSD getDistributionLogNormal(Mdouble mean, Mdouble standardDeviation, int numberOfBins) {
        PSD psd;
        psd.setDistributionLogNormal(mean, standardDeviation, numberOfBins);
        return psd;
    }
    
     /*!
     * \brief create a PSD vector for a normal distribution in Phi Units that has the demanded D50.
     */
    void setDistributionPhiNormal(Mdouble D50, Mdouble standardDeviationinPhi, int numberOfBins);
    
    /*!
     * \brief Converts a PDF to a CDF by integration.
     */
    void convertProbabilityDensityToCumulative();

    /*!
     * \brief Converts a CDF to a PDF by derivation.
     */
    void convertCumulativeToProbabilityDensity();

    /*!
     * \brief Convert any PDF to a PROBABILITYDENSITY_NUMBER_DISTRIBUTION.
     */
    void convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE PDFType);

    /*!
     * \brief Convert a PROBABILITYDENSITY_NUMBER_DISTRIBUTION to any PDF.
     */
    void convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE PDFType);

    /*!
     * \brief convert a PROBABILITYDENSITY_NUMBER_DISTRIBUTION to a PROBABILITYDENSITY_VOLUME_DISTRIBUTION.
     * \deprecated Use convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION) instead.
     */
    MERCURYDPM_DEPRECATED
    void convertProbabilityDensityNumberDistributionToProbabilityDensityVolumeDistribution();

    /*!
     * \brief Convert any CDF to a CUMULATIVE_NUMBER_DISTRIBUTION.
     */
    void convertCumulativeToCumulativeNumberDistribution(TYPE CDFType);

    /*!
     * \brief Convert a CUMULATIVE_NUMBER_DISTRIBUTION to any CDF.
     */
    void convertCumulativeNumberDistributionToCumulative(TYPE CDFType);

    /*!
     * \brief cutoff the PSD at given quantiles.
     * \deprecated Instead use cutoffCumulativeDistribution(PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION, quantileMin, quantileMax, minPolydispersity)
     */
    MERCURYDPM_DEPRECATED
    void cutoffCumulativeNumber(Mdouble quantileMin, Mdouble quantileMax, Mdouble minPolydispersity = 0.1)
     {
         cutoffCumulativeDistribution(TYPE::CUMULATIVE_NUMBER_DISTRIBUTION, quantileMin, quantileMax, minPolydispersity);
     }

    /*!
     * \brief Cut off the PSD at the given quantiles, with the PSD in the provided CDF form.
     */
    void cutoffCumulativeDistribution(TYPE CDFType, Mdouble quantileMin, Mdouble quantileMax, Mdouble minPolydispersity = 0.1);

    /*!
     * \brief cutoff the PSD at given quantiles and make it less polydisperse by squeezing it.
     * \deprecated Instead use cutoffAndSqueezeCumulativeDistribution(TYPE::CUMULATIVE_NUMBER_DISTRIBUTION, quantileMin, quantileMax, squeeze, minPolydispersity)
     */
    MERCURYDPM_DEPRECATED
    void cutoffAndSqueezeCumulative(Mdouble quantileMin, Mdouble quantileMax, Mdouble squeeze,
                                    Mdouble minPolydispersity = 0.1)
    {
        cutoffAndSqueezeCumulativeDistribution(TYPE::CUMULATIVE_NUMBER_DISTRIBUTION, quantileMin, quantileMax, squeeze, minPolydispersity);
    }

    /*!
     * \brief Cut off the PSD at the given quantiles, with the PSD in the provided CDF form, and make it less polydisperese by squeezing it.
     */
    void cutoffAndSqueezeCumulativeDistribution(TYPE CDFType, Mdouble quantileMin, Mdouble quantileMax, Mdouble squeeze, Mdouble minPolydispersity = 0.1);

    /*!
     * \brief Get smallest radius of the PSD.
     */
    Mdouble getMinRadius() const;

    /*!
     * \brief Get largest radius of the PSD.
     */
    Mdouble getMaxRadius() const;
    
    /*!
     * \brief Get the PSD vector.
     */
    std::vector<DistributionElements> getParticleSizeDistribution() const;

    /*!
     * \brief Gets the PSD vector, converted to the preferred type.
     */
    std::vector<DistributionElements> getParticleSizeDistributionByType(TYPE psdType, Mdouble scalingFactor = 1.0) const;
    
    /*!
     * \brief set the PSD by a suitable vector.
     * \deprecated Use setPSDFromVector() instead.
     */
    MERCURYDPM_DEPRECATED
    void setParticleSizeDistribution(std::vector<DistributionElements>);

    /*!
     * \brief Scales all particle sizes by a factor
     */
    void scaleParticleSize(double scale);

    /*!
     * \brief Scales all particle sizes, such that the total volume of N particles approximately equals the target volume.
     * @param numberOfParticles The number of particles N.
     * @param targetVolume The volume to match.
     * @param allowScaleDown Whether or not to scale down the particle sizes when the scale factor is smaller than 1.
     * @return The scale factor used.
     */
    double scaleParticleSizeAuto(int numberOfParticles, double targetVolume, bool allowScaleDown = false);

    /*!
     * \brief Get the number of particles already inserted into the simulation.
     */
    int getInsertedParticleNumber() const;
    
    /*!
    * \brief Decrement nParticlesPerClass_ counter.
    */
    void decrementNParticlesPerClass();
    
    /*!
     * \brief Decrement volumePerClass_ counter.
     */
    void decrementVolumePerClass(Mdouble volume);

/*!
     * \brief get a radius at a certain index of the particleSizeDistribution_
     */
    Mdouble getRadius(int index) const;
    
    /*!
     * \brief Calculate a certain diameter (e.g. D10, D50, D90, etc.) from a percentile x of the number based PSD.
     */
    Mdouble getNumberDx(Mdouble x) const;

    /*!
     * \brief Calculate a certain diameter (e.g. D10, D50, D90, etc.) from a percentile x of the length based PSD.
     */
    Mdouble getLengthDx(Mdouble x) const;

    /*!
     * \brief Calculate a certain diameter (e.g. D10, D50, D90, etc.) from a percentile x of the area based PSD.
     */
    Mdouble getAreaDx(Mdouble x) const;
    
    /*!
     * \brief Calculate a certain diameter (e.g. D10, D50, D90, etc.) from a percentile x of the volume based PSD.
     */
    Mdouble getVolumeDx(Mdouble x) const;

    /*!
     * \brief Calculate the quantile of the PSD.
     */
    Mdouble getRadiusByQuantile(Mdouble quantile) const;

    /*!
     * \brief Calculates the quantile corresponding to a certain radius.
     */
    Mdouble getQuantileByRadius(Mdouble radius) const;

    /*!
     * \brief get a volumetric mean radius of the PSD.
     */
    Mdouble getVolumetricMeanRadius() const;
    
    /*!
     * \brief get the size ratio (width) of the PSD.
     */
    Mdouble getSizeRatio() const;
    
    /*!
     * \brief Check if the size ratio is too high and cut it
     */
    void cutHighSizeRatio();

//    /*!
//     * \brief compute raw momenta of the user defined PSD.
//     */
//    void computeRawMomenta();
//
//    /*!
//     * \brief compute central momenta of the user defined PSD.
//     */
//    void computeCentralMomenta();
//
//    /*!
//     * \brief compute standardised momenta of the user defined PSD.
//     */
//    void computeStandardisedMomenta();
//
//    /*!
//     * \brief get momenta of the user defined PSD.
//     */
//    std::array<Mdouble, 6> getMomenta() const;
    

    /*!
     * \brief convert the probabilityDensity of diameter in Phi to the probabilityDensity of diameter in Meter.
     */
    std::vector< std::pair <Mdouble,Mdouble> > convertPdfPhiToPdfMeter(Mdouble meaninPhi, Mdouble standardDeviationinPhi, int numberOfBins);

    /*!
     * \brief compute the mass-based D_50 of the paired diameter and probabilityDensity.
     */
    Mdouble computeD50(std::vector< std::pair <Mdouble,Mdouble> > vectorOfPair);
    
    /*!
     * \brief create a vector of linearly spaced values.
     */
    static std::vector<Mdouble> linspace(Mdouble Min, Mdouble Max, int numberOfBins);
    
    /*!
     * \brief set a fixed seed for the random number generator; this is used for the PSDSelfTest to reproduce results
     */
    void setFixedSeed(int seed)
    {
//        random_.setLinearCongruentialGeneratorParmeters(0,0,1);
        random_.setRandomSeed(seed);
    }
    
    /*!
     * \brief determines if a certain value of the PSD vector is lower than another one. Used for std::lower_bound()
     */
    friend bool operator<(const DistributionElements& l, const DistributionElements& r);
    
    /*!
     * \brief determines if a certain value of the PSD vector is lower than a double.
     */
    friend bool operator<(const DistributionElements& l, Mdouble r);
    
    /*!
     * \brief Writes to output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, DistributionElements& p);
    
    /*!
     * \brief Reads from input stream.
     */
    friend std::istream& operator>>(std::istream& is, DistributionElements& p);
    
    /*!
     * \brief Determines if a certain value of the PSD vector is equal to a double.
     */
    friend Mdouble operator==(DistributionElements l, Mdouble r);

private:
    /*!
     * Vector of the DistributionElements class which stores radii as internalVariable and probabilities of a PSD.
     */
    std::vector<DistributionElements> particleSizeDistribution_;

    /*!
     * Vector of integers which represents the number of inserted particles in each size class.
     * The classes in this vector are defined to contain the particles between size r_i and r_i-1. (e.g. size class
     * 12 consists of particles between size class 12 and 11 of the PDF)
     */
    std::vector<int> nParticlesPerClass_;
    
    /*!
     * Vector of doubles which stores the volume of inserted particles for each size class. This vector is used in
     * the insertManuallyByVolume() function to check if the volumeAllowed per class is exceeded and thus no further
     * particles should be added to a certain class.
     * The classes in this vector are defined to contain the volume of particles between size r_i and r_i-1. (e.g. size
     * class 12 consists of the particles' volume between size class 12 and 11 of the PDF)
     */
    std::vector<Mdouble> volumePerClass_;
    
    /*!
     * Integer which determines the class in which a particle has to be inserted for the manual insertion routine.
     */
    int index_;
    
    /*!
     * Mercury random number generator object used to draw random numbers from a random initial seed
     */
    RNG random_;
};


#endif //MECURYDPM_PSD_H
