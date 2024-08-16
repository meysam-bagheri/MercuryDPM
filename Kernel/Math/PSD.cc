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

#include "PSD.h"

/*!
 * \details Default constructor; sets every data member to 0 or default.
 */
PSD::PSD()
{
//    for (auto& momenta: momenta_)
//        momenta = 0;
    index_ = 0;
    // default seed is zero for reproducibility
    random_.setRandomSeed(0);
}

/*!
 * \details Copy constructor
 */
PSD::PSD(const PSD& other)
{
    particleSizeDistribution_ = other.particleSizeDistribution_;
//    momenta_ = other.momenta_;
    index_ = other.index_;
    random_ = other.random_;
}

/*!
 * \details Destructor. Since there are no pointers in this class, there is no
 *          need for any actions here.
 */
PSD::~PSD()
= default;


/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 * \return pointer to the copy on the heap
 */
PSD* PSD::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "PSD::copy() const finished" << std::endl;
#endif
    return new PSD(*this);
}

/*!
 * \details Prints the radii [m] and probabilities [%] of the psd vector. It currently supports nanometers, micrometers,
 *          milimeters and meters as size units.
 */
void PSD::printPSD() const
{
    if (particleSizeDistribution_.front().internalVariable > 1e-1)
    {
        for (const auto p: particleSizeDistribution_)
        {
            logger(INFO, "%m\t %\%", p.internalVariable, p.probability * 100);
        }
    }
    else if (particleSizeDistribution_.front().internalVariable > 1e-4)
    {
        for (const auto p: particleSizeDistribution_)
        {
            logger(INFO, "%mm\t %\%", p.internalVariable * 1000, p.probability * 100);
        }
    }
    else if (particleSizeDistribution_.front().internalVariable > 1e-7)
    {
        for (const auto p: particleSizeDistribution_)
        {
            logger(INFO, "%um\t %\%", p.internalVariable * 1e6, p.probability * 100);
        }
    }
    else
    {
        for (const auto p: particleSizeDistribution_)
        {
            logger(INFO, "%nm\t %\%", p.internalVariable * 1e9, p.probability * 100);
        }
    }
}

/*!
 * \details Draws a sample probability of a real uniform distribution. A random number is
 * generated in range [0,1] and by linear interpolation a suitable radius is returned. This function is only valid
 * for CumulativeNumberDistributions as particles are drawn and not volumes, areas or lengths.
 * \return A Radius for a randomly assigned probability.
 */
Mdouble PSD::drawSample()
{
    // draw a number between 0 and 1, uniformly distributed
    Mdouble prob = random_.getRandomNumber(0, 1);
    // find the interval [low,high] of the psd in which the number sits
    auto high = std::lower_bound(particleSizeDistribution_.begin(), particleSizeDistribution_.end(), prob);
    auto low = std::max(particleSizeDistribution_.begin(), high - 1);
    Mdouble rMin = low->internalVariable;
    Mdouble rMax = high->internalVariable;
    Mdouble pMin = low->probability;
    Mdouble pMax = high->probability;
    // interpolate linearly between [low.radius,high.radius] (assumption: CDF is piecewise linear)
    Mdouble a = (prob - pMin) / (pMax - pMin);
    return a * rMin + (1 - a) * rMax;
}

/*!
 * \details Draws a sample probability of a real uniform distribution within a given size class. A random number is
 * generated in the size class range [r_i,r_i+1] and by linear interpolation a suitable radius is returned. from a
 * CUMULATIVE_VOLUME_DISTRIBUTION the volumeAllowed of each class is checked to insert the PSD as accurate as possible.
 * This function is only valid for cumulativeNumberDistributions as particles are drawn and not volumes, areas or
 * lengths. Furthermore this insertion routine is most accurate for non-continuous particle insertion.
 * \param[in] volume                         volume of the geometry to be filled.
 * \return A Radius for a randomly assigned probability of a specific size class.
 */
Mdouble PSD::insertManuallyByVolume(Mdouble volume)
{
    // initialize the particleNumberPerClass vector if empty.
    if (nParticlesPerClass_.empty())
        nParticlesPerClass_.resize(particleSizeDistribution_.size());
    // initialize the particleNumberPerClass vector if empty.
    if (volumePerClass_.empty())
        volumePerClass_.resize(particleSizeDistribution_.size());

    for (auto it = particleSizeDistribution_.end() - 1; it != particleSizeDistribution_.begin(); --it)
    {
        static std::mt19937 gen(0);
        static std::uniform_real_distribution<Mdouble> dist(0, 1);
        Mdouble prob = dist(gen);
        // map the probability to the current class interval
        prob = (it - 1)->probability + (it->probability - (it - 1)->probability) * prob;
        // find the interval [low,high] of the psd in which the number sits
        auto high = std::lower_bound(particleSizeDistribution_.begin(), particleSizeDistribution_.end(), prob);
        auto low = std::max(particleSizeDistribution_.begin(), high - 1);
        Mdouble rMin = low->internalVariable;
        Mdouble rMax = high->internalVariable;
        Mdouble pMin = low->probability;
        Mdouble pMax = high->probability;
        // interpolate linearly between [low.internalVariable,high.internalVariable] (assumption: CDF is piecewise linear)
        index_ = std::distance(particleSizeDistribution_.begin(), high);
        Mdouble a = (prob - pMin) / (pMax - pMin);
        Mdouble rad = a * rMin + (1 - a) * rMax;
        //\todo generalize this to work for non-spherical particles
        // Compare inserted volume against volumeAllowed (only works for spherical at the moment)
        Mdouble radVolume = 4.0 / 3.0 * constants::pi * rad * rad * rad;
        convertCumulativeToProbabilityDensity();
        convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
        Mdouble volumeAllowed = particleSizeDistribution_[index_].probability * volume;
        convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
        convertProbabilityDensityToCumulative();
        if (volumePerClass_[index_] > volumeAllowed)
        {
            continue;
        }
        else if (radVolume + volumePerClass_[index_] > volumeAllowed)
        {
            Mdouble differenceVolumeLow = -(volumePerClass_[index_] - volumeAllowed);
            Mdouble differenceVolumeHigh = radVolume + volumePerClass_[index_] - volumeAllowed;
            Mdouble volumeRatio = differenceVolumeLow / differenceVolumeHigh;
            // If volumeRatio > 1 it will insert anyways, because it should be no problem for the distribution
            prob = dist(gen);
            if (prob <= volumeRatio)
            {
                ++nParticlesPerClass_[index_];
                volumePerClass_[index_] += radVolume;
                return rad;
            }
            else
            {
                continue;
            }
        }
        else
        {
            ++nParticlesPerClass_[index_];
            volumePerClass_[index_] += radVolume;
            return rad;
        }
        
    }
    return 0;
}

/*!
 * \details Validates if a distribution is cumulative by first checking if the psd vector is empty and that the
 * scaling of probabilities is in the range [0,1]. Also it checks if each consecutive value is higher than the latter
 * value (i.e. assuming piecewise linear CDF: p_i > p_i-1). Further it replaces probabilities if the CDF does not
 * start with zero or does not end with unity (i.e. p_0=0 and p_end=1)
 */
void PSD::validateCumulativeDistribution()
{
    // check whether the distribution is cumulative
    for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
    {
        logger.assert_always(it->probability >= (it - 1)->probability,
                             "psd is not cumulative: radius % %, probabilities % %", (it - 1)->internalVariable,
                             it->internalVariable, (it - 1)->probability, it->probability);
    }
    // set negative probabilities to 0
    for (auto& rp : particleSizeDistribution_)
    {
        if (rp.probability < 0)
        {
            logger(INFO, "negative probability % has been set to 0", rp.probability);
            rp.probability = 0;
        }
    }
    // ensure interval of probabilities to be [0,1]
    for (auto p = 0; p < particleSizeDistribution_.size(); ++p)
        particleSizeDistribution_[p].probability /= particleSizeDistribution_.back().probability;
    // cdf needs to start with a probability of zero
    if (particleSizeDistribution_[0].probability != 0)
    {
        logger(INFO, "adding a zero at the beginning of the psd");
        particleSizeDistribution_.insert(particleSizeDistribution_.begin(),
                                         {std::max(1e-3 * particleSizeDistribution_[0].internalVariable,
                                                   2 * particleSizeDistribution_[0].internalVariable -
                                                   particleSizeDistribution_[1].internalVariable), 0});
    }
    // cdf needs to end with a probability of one
    if (particleSizeDistribution_.back().probability < 1)
    {
        logger(INFO, "adding a one at the end of the psd");
        particleSizeDistribution_.push_back({particleSizeDistribution_.back().internalVariable, 1});
    }
    // cut off equal subsequent values on the end of the cdf to remove zero size classes.
    for (auto i = particleSizeDistribution_.end() - 1; i != particleSizeDistribution_.begin(); i--)
    {
        if (i->probability == (i - 1)->probability)
        {
            particleSizeDistribution_.erase(particleSizeDistribution_.end() + std::distance(particleSizeDistribution_
                                                                                                    .end(), i));
        }
        else
        {
            break;
        }
    }
}

/*!
 * \details Validates if a PDF is valid by first checking if the PDF starts with zero (i.e. p_0=0). Further it checks
 * if the integral is equal to unity (i.e. assuming piecewise constant PDF: sum(p_i)=1).
 */
void PSD::validateProbabilityDensityDistribution()
{
    // set negative probabilities to 0
    for (auto& rp : particleSizeDistribution_)
    {
        if (rp.probability < 0)
        {
            logger(INFO, "negative probability % has been set to 0", rp.probability);
            rp.probability = 0;
        }
    }
    Mdouble sum = 0;
    // check if the psd starts with zero probability (i.e. p_0=0)
    if (particleSizeDistribution_[0].probability != 0)
    {
        logger(INFO, "adding a zero at the beginning of PDF");
        particleSizeDistribution_.insert(particleSizeDistribution_.begin(),
                                         {std::max(1e-3 * particleSizeDistribution_[0].internalVariable,
                                                   2 * particleSizeDistribution_[0].internalVariable -
                                                   particleSizeDistribution_[1].internalVariable), 0});
    }
    // sum up probabilities to check if it equals to unity (sum(p_i)=1)
    for (auto& it : particleSizeDistribution_)
    {
        sum += it.probability;
    }
    // ensure interval of probabilities to be [0,1] by normalization.
    for (auto& it : particleSizeDistribution_)
        it.probability /= sum;
    // if the sum of probabilities is not equal to unity, sum up the normalized probabilities again
    if (sum - 1 > 1e-6)
    {
        sum = 0;
        for (auto& it: particleSizeDistribution_)
        {
            sum += it.probability;
        }
    }
    logger.assert_debug(sum - 1 < 1e-6, "PDF is not valid: Integral of PDF is not equal to unity");
}

/*!
 * \details sets the particle size distribution to a discretised uniform (linear) cumulative number distribution
 * \param[in] radMin            Double representing The smallest particle radius of the particle size distribution
 * \param[in] radMax            Double representing the biggest particle radius of the particle size distribution
 * \param[in] numberOfBins      Integer determining the number of bins (aka. particle size classes or resolution) of the
 *                              particle size distribution
 */
void PSD::setDistributionUniform(Mdouble radMin, Mdouble radMax, int numberOfBins)
{
    if (!particleSizeDistribution_.empty())
    {
        particleSizeDistribution_.clear();
    }
    std::vector<Mdouble> probabilities = helpers::linspace(0.0, 1.0, numberOfBins);
    std::vector<Mdouble> radii = helpers::linspace(radMin, radMax, numberOfBins);
    // combine radii and probabilities
    for (int i = 0; i < radii.size(); ++i)
    {
        particleSizeDistribution_.push_back({radii[i], probabilities[i]});
    }
    validateCumulativeDistribution();
}

/*!
 * \details sets the particle size distribution to a discretised normal (gaussian) cumulative number distribution,
 * which covers the range of 3 * standardDeviation (99,73% of all values covered).
 * \param[in] mean                  Double representing the mean of the particle size distribution
 * \param[in] standardDeviation     Double representing the standard deviation of the particle size distribution
 * \param[in] numberOfBins          Integer determining the number of bins (aka. particle size classes or resolution)
 *                                  of the particle size distribution
 */
void PSD::setDistributionNormal(Mdouble mean, Mdouble standardDeviation, int numberOfBins)
{
    logger.assert_always(mean > 3 * standardDeviation,
                         "Reduce standardDeviation of your normal distribution to avoid negative radii; The mean "
                         "should be greater than 3*standardDeviation");
    if (!particleSizeDistribution_.empty())
    {
        particleSizeDistribution_.clear();
    }
    Mdouble radMin = mean - 3 * standardDeviation;
    Mdouble radMax = mean + 3 * standardDeviation;
    std::vector<Mdouble> radii = helpers::linspace(radMin, radMax, numberOfBins);
    std::vector<Mdouble> probabilities;
    for (int i = 0; i < radii.size(); i++)
    {
        Mdouble probability = 0.5 * (1 + erf((radii[i] - mean) / (sqrt(2) * standardDeviation)));
        probabilities.push_back(probability);
    }
    for (int j = 0; j < radii.size(); j++)
    {
        particleSizeDistribution_.push_back({radii[j], probabilities[j]});
    }
    validateCumulativeDistribution();
}

/*!
 * \details sets the particle size distribution to a discretised lognormal cumulative number distribution,
 * which generates random values of radius following a lognormal distribution and iterates
 * until the mass-based D50 is close to two times the user-defined mean radius.
 * \param[in] mean                  Double representing the mean of the particle size distribution in meters (radius)
 * \param[in] standardDeviation     Double representing the standard deviation of the particle size distribution in meters (radius)
 * \param[in] numberOfBins          Integer determining the number of bins (aka. particle size classes or resolution)
 *                                  of the particle size distribution
 * See https://en.cppreference.com/w/cpp/numeric/random/lognormal_distribution
 */
void PSD::setDistributionLogNormal(Mdouble mean, Mdouble standardDeviation, int numberOfBins)
{
    //setDistributionNormal(mean, standardDeviation, numberOfBins);
    if (!particleSizeDistribution_.empty())
    {
        particleSizeDistribution_.clear();
    }
    Mdouble radMin = mean - 3 * standardDeviation;
    Mdouble radMax = mean + 3 * standardDeviation;
    std::vector<Mdouble> radii = helpers::linspace(radMin, radMax, numberOfBins);
    std::vector<Mdouble> probabilities;
    for (int i = 0; i < radii.size(); i++)
    {
        Mdouble probability = 0.5 * (1 + erf((radii[i] - mean) / (sqrt(2) * standardDeviation)));
        probabilities.push_back(probability);
    }
    for (int j = 0; j < radii.size(); j++)
    {
        particleSizeDistribution_.push_back({radii[j], probabilities[j]});
    }
    validateCumulativeDistribution();
    for (auto& p: particleSizeDistribution_)
    {
        p.internalVariable = exp(p.internalVariable);
    }
}

/*!
 * \details sets the particle size distribution to a discretised cumulative number distribution,
 * which is a normal distribution in phi units and has the demanded D50. D50 is the mass-based median grain size, while
 * the Phi scale is a logarithmic scale computed by: Phi = -log2(1000*D). It can be used for visualizing sediment grain
 * size distributions over a wide range of grain sizes. The combination of D50 and standardDeviation in Phi units is a
 * representative measure of grain size distribution in the field of Sedimentology.
 * Source: Donoghue, J.F. (2016). Phi Scale. In: Kennish, M.J. (eds) Encyclopedia of Estuaries. Encyclopedia of Earth Sciences Series.
 * Springer, Dordrecht. https://doi.org/10.1007/978-94-017-8801-4_277
 * \param[in] D50                       Double representing the mass-based D50 of the particle size distribution in meters
 * \param[in] standardDeviationInPhi    Double representing the standard deviation of the particle size distribution (diameter) in phi units
 * \param[in] numberOfBins              Integer determining the number of bins (aka. particle size classes or resolution)
 *                                      of the particle size distribution
 */
void PSD::setDistributionPhiNormal(Mdouble D50, Mdouble standardDeviationInPhi, int numberOfBins){
    Mdouble PhiMean = -log2(1000*D50);
    std::vector<std::pair<Mdouble,Mdouble>> vect = convertPdfPhiToPdfMeter(PhiMean, standardDeviationInPhi, numberOfBins);
    std::sort(vect.begin(), vect.end());
    Mdouble D50Test = computeD50(vect);
    Mdouble res = D50Test - D50;
    Mdouble tol = D50*1e-5;
    int itter = 0;
    Mdouble dPhiMean, slope, res_old;
    while (res > tol && itter < 100) {
        if(itter==0) dPhiMean = 0.1;
        else {
            slope = (res - res_old)/dPhiMean;
            dPhiMean = -res/slope;
        }
        PhiMean = PhiMean+dPhiMean;
        vect = convertPdfPhiToPdfMeter(PhiMean, standardDeviationInPhi, numberOfBins);
        std::sort(vect.begin(), vect.end());
        D50Test = computeD50(vect);
        res_old = res;
        res = D50Test - D50;
        itter++;
    }
    logger(INFO, "D50 %, Final D50Test %", D50, D50Test);

    Mdouble probTot=0;
    for(int i = 0; i < vect.size(); i++){
        probTot += vect[i].second;
    }
    for (int i = 0; i < vect.size(); i++)
    {
        vect[i].second = vect[i].second/probTot;
    }

    if (!particleSizeDistribution_.empty())
    {
        particleSizeDistribution_.clear();
    }
    for (int j = 0; j < vect.size(); j++)
    {
        particleSizeDistribution_.push_back({vect[j].first*0.5, vect[j].second});
    }
    convertProbabilityDensityToCumulative();
    validateCumulativeDistribution();
}

/*!
 * \details generates a normal particle size distribution in Phi units and converts it to the distribution in meters.
 * \param[in] meanInPhi                  Double representing the mean diameter of the particle size distribution in phi units
 * \param[in] standardDeviationInPhi     Double representing the standard deviation of the particle size distribution (diameter) in phi units
 * \param[in] numberOfBins               Integer determining the number of bins (aka. particle size classes or resolution)
 *                                       of the particle size distribution
 */
std::vector< std::pair <Mdouble,Mdouble> > PSD::convertPdfPhiToPdfMeter(Mdouble meanInPhi, Mdouble standardDeviationInPhi, int numberOfBins) {
    Mdouble dMin = meanInPhi - 6 * standardDeviationInPhi;
    Mdouble dMax = meanInPhi + 6 * standardDeviationInPhi;
    std::vector<Mdouble> diameterPhi = helpers::linspace(dMin, dMax, numberOfBins);
    std::vector<Mdouble> pdfPhi;
    for (int i = 0; i < diameterPhi.size(); i++) {
        Mdouble probability = 1 / (standardDeviationInPhi * sqrt(2.0 * constants::pi)) *
                              exp(-0.5 * pow((diameterPhi[i] - meanInPhi) / standardDeviationInPhi, 2));
        pdfPhi.push_back(probability);
    }
    std::vector<Mdouble> diameter;
    std::vector<Mdouble> pdf(pdfPhi.size());
    for(int i = 0; i < diameterPhi.size(); i++){
        Mdouble D = pow(2, -diameterPhi[i])/1000.0;
        diameter.push_back(D);
    }
    pdf[0] = 0;
    for(int i= 0; i < diameterPhi.size()-1; i++){
        Mdouble area = (diameterPhi[i+1]-diameterPhi[i]) * 0.5*(pdfPhi[i+1]+pdfPhi[i]);
        Mdouble dD = fabs(diameter[i+1]-diameter[i]);
        pdf[i+1] = 2.0*area/dD - pdf[i];
    }
    std::vector< std::pair <Mdouble,Mdouble> > vect;
    for (int i=0; i<diameter.size(); i++)
        vect.push_back(std::make_pair(diameter[i],pdf[i]));
    return vect;
}

/*!
 * \details computes the mass-based D50 of the paired diameter and probabilityDensity
 * \param[in]  vectorOfPair         vector of pair containing paired data of diameter and probabilityDensity
 */
Mdouble PSD::computeD50(std::vector< std::pair <Mdouble,Mdouble> > vectorOfPair)
    {
        Mdouble D_50 = 0.0;
        int numberOfBins = vectorOfPair.size();
        std::vector<Mdouble> Vcum(numberOfBins, 0);
        Mdouble Vint=0;
        for (int i = 0; i < Vcum.size()-1; i++)
        {
                Mdouble dD = fabs(vectorOfPair[i+1].first-vectorOfPair[i].first);
                Mdouble V = 0.5*(pow(vectorOfPair[i].first, 3.0) * vectorOfPair[i].second +
                        pow(vectorOfPair[i+1].first, 3.0) * vectorOfPair[i+1].second) * dD;
                Vint = Vint + V;
                Vcum[i+1] = Vint;
        }

        for(int j=0; j<Vcum.size()-1; ++j){
            if(Vcum[j]/Vint*100<=50 && Vcum[j+1]/Vint*100>50){
                D_50 = vectorOfPair[j].first + fabs(vectorOfPair[j+1].first-vectorOfPair[j].first)/(Vcum[j+1]/Vint*100-Vcum[j]/Vint*100)
                        *(50-Vcum[j]/Vint*100);
            }
        }
        return D_50;
    }

/*!
 * \details creates the psd vector from radii and probabilities filled in by hand. The Type of PSD will be converted
 * to the default cumulative number distribution function (CUMULATIVE_NUMBER_DISTRIBUTION) for further processing.
 * \param[in] psdVector                         Vector containing radii and probabilities ([0,1]).
 * \param[in] PSDType                           Type of the PSD: CUMULATIVE_VOLUME_DISTRIBUTION,
 *                                              CUMULATIVE_NUMBER_DISTRIBUTION, CUMULATIVE_LENGTH_DISTRIBUTION,
 *                                              CUMULATIVE_AREA_DISTRIBUTION, PROBABILITYDENSITY_VOLUME_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_NUMBER_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_LENGTH_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_AREA_DISTRIBUTION.
 */
void PSD::setPSDFromVector(std::vector<DistributionElements> psdVector, TYPE PSDType)
{
    particleSizeDistribution_ = psdVector;

    // To prevent errors when the PSD is empty, throw a warning and return.
    ///\todo JWB: Is a PSD ever allowed to be empty? drawSample() will cause a segfault in that case. I'm not throwing
    /// an error here, as that causes the UnitTests/FullRestartSelfTest to fail, due to the InsertionBoundary being
    /// allowed to use (and defaults to using) an empty PSD.
    if (particleSizeDistribution_.empty())
    {
        logger(WARN, "Setting an empty PSD!");
        return;
    }

    switch (PSDType)
    {
        // default case is a CUMULATIVE_NUMBER_DISTRIBUTION
        default:
            validateCumulativeDistribution();
            break;
        case TYPE::CUMULATIVE_LENGTH_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(
                    TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::CUMULATIVE_AREA_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::CUMULATIVE_VOLUME_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(
                    TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_NUMBER_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
    }
}

/*!
 * \details creates the PSD vector from probabilities and radii saved in a .csv file. Radii should be located at
 *          the first column and probabilities at the second column of the file. The Type of PSD will be converted to the
 *          default cumulative number distribution function (CUMULATIVE_NUMBER_DISTRIBUTION) for further processing.
 * \param[in] fileName                          Name of the .csv file containing the radii and probabilities of the PSD.
 * \param[in] PSDType                           Type of the PSD: CUMULATIVE_VOLUME_DISTRIBUTION,
 *                                              CUMULATIVE_NUMBER_DISTRIBUTION, CUMULATIVE_LENGTH_DISTRIBUTION,
 *                                              CUMULATIVE_AREA_DISTRIBUTION, PROBABILITYDENSITY_VOLUME_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_NUMBER_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_LENGTH_DISTRIBUTION PROBABILITYDENSITY_AREA_DISTRIBUTION.
 * \param[in] headings                          If TRUE the file is assumed to have headings and the first row will
 *                                              be skipped. If FALSE the file has no headings and the file will be
 *                                              read in as is. Default is FALSE.
 * \param[in] unitScalingFactorRadii            Scaling factor of radii to match SI-units.
 */
void PSD::setPSDFromCSV(const std::string& fileName, TYPE PSDType, bool headings, Mdouble unitScalingFactorRadii)
{
    // read in csv file using the csvReader class
    csvReader csv;
    csv.setHeader(headings);
    csv.read(fileName);
    std::vector<Mdouble> radii = csv.getFirstColumn(unitScalingFactorRadii);
    std::vector<Mdouble> probabilities = csv.getSecondColumn(1.0);
    logger.assert_always(radii.size() == probabilities.size(), "The radii and probabilities vector have to be the "
                                                               "same size");
    // combine radii and probabilities
    std::vector<DistributionElements> psd;
    for (int i = 0; i < radii.size(); ++i)
    {
        psd.push_back({radii[i], probabilities[i]});
    }
    logger.assert_always(!psd.empty(), "PSD cannot be empty");
    setPSDFromVector(psd, PSDType);
    logger(INFO, "A PSD with size ratio % is now ready to be set", getSizeRatio());
}

/*!
 * \details Convert any type of PDF to a CDF. Probabilities are integrated for each radius by cumulatively summing
 * them up (i.e. p_i = p_i + p_i-1).
 */
void PSD::convertProbabilityDensityToCumulative()
{
    // add up probabilities
    for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
    {
        it->probability += (it - 1)->probability;
    }
    // check whether cumulative distribution is valid
    validateCumulativeDistribution();
}

/*!
 * \details Convert any type of CDF to a PDF. Probabilities are derivated for each radius by substracting the CDF
 * probabilities (i.e. pPDF_i = pCDF_i - pCDF_i-1).
 */
void PSD::convertCumulativeToProbabilityDensity()
{
    // subtract cumulative probabilities
    Mdouble probabilityOld = 0, probabilityCDF;
    for (auto& it : particleSizeDistribution_)
    {
        probabilityCDF = it.probability;
        it.probability -= probabilityOld;
        probabilityOld = probabilityCDF;
    }
    // check whether probability density distribution is valid
    validateProbabilityDensityDistribution();
}

/*!
 * \details Converts any of the PDFTypes to a PROBABILITYDENSITY_NUMBER_DISTRIBUTION based on their TYPE. This is a helper function
 * which enables the convertProbabilityDensityToCumulative function to convert the psd into the default TYPE (CUMULATIVE_NUMBER_DISTRIBUTION).
 * \param[in] PDFType Type of the PDF: PROBABILITYDENSITY_LENGTH_DISTRIBUTION, PROBABILITYDENSITY_AREA_DISTRIBUTION or PROBABILITYDENSITY_VOLUME_DISTRIBUTION.
 */
void PSD::convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE PDFType)
{
    Mdouble sum = 0;
    switch  (PDFType)
    {
        default:
            logger(ERROR, "Wrong PDFType");
            break;
        case TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability /= 0.5 * (it->internalVariable + (it - 1)->internalVariable);
                // old conversion
                //p.probability /= p.internalVariable;
                // sum up probabilities
                sum += it->probability;
            }
            break;
        case TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability /=
                        1.0 / 3.0 * (square(it->internalVariable) + it->internalVariable * (it - 1)->internalVariable +
                                     square((it - 1)->internalVariable));
                // old conversion
                //p.probability /= square(p.internalVariable);
                // sum up probabilities
                sum += it->probability;
            }
            break;
        case TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability /=
                        0.25 * (square(it->internalVariable) + square((it - 1)->internalVariable)) *
                        (it->internalVariable + (it - 1)->internalVariable);
                // old conversion
                //p.probability /= cubic(p.internalVariable);
                // sum up probabilities
                sum += it->probability;
            }
            break;
    }

    // normalize
    for (auto& p : particleSizeDistribution_)
    {
        p.probability /= sum;
    }
}

/*!
 * \details This is a helper function to have full flexibility in converting from a
 * PROBABILITYDENSITY_NUMBER_DISTRIBUTION type to any other PDF type. This can be useful e.g. when wanting to output the
 * PSD with a certain type.
 * @param PDFType The PDF type to convert to.
 */
void PSD::convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE PDFType)
{
    Mdouble sum = 0;
    switch  (PDFType)
    {
        default:
            logger(ERROR, "Wrong PDFType");
            break;
        case TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability *= 0.5 * (it->internalVariable + (it - 1)->internalVariable);
                // old conversion
                //p.probability *= p.internalVariable;
                // sum up probabilities
                sum += it->probability;
            }
            break;
        case TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability *=
                        1.0 / 3.0 * (square(it->internalVariable) + it->internalVariable * (it - 1)->internalVariable +
                                     square((it - 1)->internalVariable));
                // old conversion
                //p.probability *= square(p.internalVariable);
                // sum up probabilities
                sum += it->probability;
            }
            break;
        case TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability *=
                        0.25 * (square(it->internalVariable) + square((it - 1)->internalVariable)) *
                        (it->internalVariable + (it - 1)->internalVariable);
                // old conversion
                //p.probability *= cubic(p.internalVariable);
                // sum up probabilities
                sum += it->probability;
            }
            break;
    }

    // normalize
    for (auto& p : particleSizeDistribution_)
    {
        p.probability /= sum;
    }
}

/*!
 * \details converts a PROBABILITYDENSITY_NUMBER_DISTRIBUTION to a PROBABILITYDENSITY_VOLUME_DISTRIBUTION. Used for
 * getVolumetricMeanRadius() and insertManuallyByVolume().
 */
MERCURYDPM_DEPRECATED
void PSD::convertProbabilityDensityNumberDistributionToProbabilityDensityVolumeDistribution()
{
    convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
}

/*!
 * \details Converts any of the CDFTypes to a the default CUMULATIVE_NUMBER_DISTRIBUTION based on their TYPE.
 * \param[in] CDFType Type of the PDF: CUMULATIVE_LENGTH_DISTRIBUTION, CUMULATIVE_AREA_DISTRIBUTION or CUMULATIVE_VOLUME_DISTRIBUTION.
 */
void PSD::convertCumulativeToCumulativeNumberDistribution(TYPE CDFType)
{
    // Ignore the default case.
    if (CDFType == TYPE::CUMULATIVE_NUMBER_DISTRIBUTION)
        return;

    convertCumulativeToProbabilityDensity();

    if (CDFType == TYPE::CUMULATIVE_LENGTH_DISTRIBUTION)
        convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION);
    else if (CDFType == TYPE::CUMULATIVE_AREA_DISTRIBUTION)
        convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION);
    else if (CDFType == TYPE::CUMULATIVE_VOLUME_DISTRIBUTION)
        convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
    else
        logger(ERROR, "Wrong CDFType");

    convertProbabilityDensityToCumulative();
}

/*!
 * @param CDFType The CDF type to convert to.
 */
void PSD::convertCumulativeNumberDistributionToCumulative(TYPE CDFType)
{
    // Ignore the default case.
    if (CDFType == TYPE::CUMULATIVE_NUMBER_DISTRIBUTION)
        return;

    convertCumulativeToProbabilityDensity();

    if (CDFType == TYPE::CUMULATIVE_LENGTH_DISTRIBUTION)
        convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION);
    else if (CDFType == TYPE::CUMULATIVE_AREA_DISTRIBUTION)
        convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION);
    else if (CDFType == TYPE::CUMULATIVE_VOLUME_DISTRIBUTION)
        convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
    else
        logger(ERROR, "Wrong CDFType");

    convertProbabilityDensityToCumulative();
}

/*!
 * \details Cuts off the CDF at given minimum and maximum quantiles and applies a minimum polydispersity at the base.
 * \param[in] CDFType                           The CDF form the PSD should have when applying the cut off.
 * \param[in] quantileMin                       Undersize quantile to cut off the lower part of the CDF.
 * \param[in] quantileMax                       Oversize quantile to cut off the upper part of the CDF.
 * \param[in] minPolydispersity                 Applies a minimum polydispersity ([0,1]) at the base of the CDF.
 */
void PSD::cutoffCumulativeDistribution(PSD::TYPE CDFType, Mdouble quantileMin, Mdouble quantileMax, Mdouble minPolydispersity)
{
    // Convert the default CUMULATIVE_NUMBER_DISTRIBUTION to the provided CDF type.
    convertCumulativeNumberDistributionToCumulative(CDFType);

    Mdouble radiusMin = getRadiusByQuantile(quantileMin);
    Mdouble radiusMax = getRadiusByQuantile(quantileMax);
    // to get a minimum polydispersity at the base
    Mdouble radiusMinCut = std::min(radiusMin * (1 + minPolydispersity), radiusMax);
    // cut off min
    while (particleSizeDistribution_.front().internalVariable <= radiusMinCut)
    {
        particleSizeDistribution_.erase(particleSizeDistribution_.begin());
    }
    particleSizeDistribution_.insert(particleSizeDistribution_.begin(), {radiusMinCut, quantileMin});
    particleSizeDistribution_.insert(particleSizeDistribution_.begin(), {radiusMin, 0});
    // cut off max
    while (particleSizeDistribution_.back().internalVariable >= radiusMax)
    {
        particleSizeDistribution_.pop_back();
    }
    Mdouble radiusMaxCut = std::max(radiusMax - (1 - minPolydispersity) * (radiusMax - particleSizeDistribution_.back()
            .internalVariable), radiusMin);
    particleSizeDistribution_.push_back({radiusMaxCut, quantileMax});
    particleSizeDistribution_.push_back({radiusMax, 1});

    // Convert the provided CDF type back to the default CUMULATIVE_NUMBER_DISTRIBUTION.
    convertCumulativeToCumulativeNumberDistribution(CDFType);

    logger(INFO, "PSD was cut successfully and now has a size ratio of %", getSizeRatio());
}

/*!
 * \details Cuts off the CDF at given minimum and maximum quantiles, applies a minimum polydispersity at the base and
 * squeezes the distribution to make it less polydisperse.
 * \param[in] CDFType                           The CDF form the PSD should have when applying the cut off and squeeze.
 * \param[in] quantileMin                       Undersize quantile to cut off the lower part of the CDF.
 * \param[in] quantileMax                       Oversize quantile to cut off the upper part of the CDF.
 * \param[in] squeeze                           Applies a squeezing factor ([0,1]) which determines the degree the
 *                                              PDF gets squeezed.
 * \param[in] minPolydispersity                 Applies a minimum polydispersity ([0,1]) at the base of the CDF.
 */
void PSD::cutoffAndSqueezeCumulativeDistribution(PSD::TYPE CDFType, Mdouble quantileMin, Mdouble quantileMax,
                                                 Mdouble squeeze, Mdouble minPolydispersity)
{
    Mdouble r50;
    if (CDFType == TYPE::CUMULATIVE_NUMBER_DISTRIBUTION)
        r50 = 0.5 * PSD::getNumberDx(50);
    else if (CDFType == TYPE::CUMULATIVE_LENGTH_DISTRIBUTION)
        r50 = 0.5 * PSD::getLengthDx(50);
    else if (CDFType == TYPE::CUMULATIVE_AREA_DISTRIBUTION)
        r50 = 0.5 * PSD::getAreaDx(50);
    else if (CDFType == TYPE::CUMULATIVE_VOLUME_DISTRIBUTION)
        r50 = 0.5 * PSD::getVolumeDx(50);
    else
        logger(ERROR, "Wrong CDFType");

    // cut off
    cutoffCumulativeDistribution(CDFType, quantileMin, quantileMax, minPolydispersity);
    // squeeze psd
    for (auto& p: particleSizeDistribution_)
    {
        p.internalVariable = r50 + (p.internalVariable - r50) * squeeze;
    }
}

/*!
 * \details Gets the diameter from a certain percentile of the number based PSD.
 * \return A double which is the diameter corresponding to a certain percentile.
 * \param[in] x double [0, 100] which determines the obtained diameter as a percentile of the PSD.
 */
Mdouble PSD::getNumberDx(Mdouble x) const
{
    logger.assert_always(x >= 0 && x <= 100, "Percentile % is not between 0 and 100.", x);
    return 2.0 * getRadiusByQuantile(x / 100);
}

/*!
 * \details Gets the diameter from a certain percentile of the length based PSD.
 * \return A double which is the diameter corresponding to a certain percentile.
 * \param[in] x double [0, 100] which determines the obtained diameter as a percentile of the PSD.
 */
Mdouble PSD::getLengthDx(Mdouble x) const
{
    logger.assert_always(x >= 0 && x <= 100, "Percentile % is not between 0 and 100.", x);
    PSD psd = *this;
    psd.convertCumulativeNumberDistributionToCumulative(TYPE::CUMULATIVE_LENGTH_DISTRIBUTION);
    return 2.0 * psd.getRadiusByQuantile(x / 100);
}

/*!
 * \details Gets the diameter from a certain percentile of the area based PSD.
 * \return A double which is the diameter corresponding to a certain percentile.
 * \param[in] x double [0, 100] which determines the obtained diameter as a percentile of the PSD.
 */
Mdouble PSD::getAreaDx(Mdouble x) const
{
    logger.assert_always(x >= 0 && x <= 100, "Percentile % is not between 0 and 100.", x);
    PSD psd = *this;
    psd.convertCumulativeNumberDistributionToCumulative(TYPE::CUMULATIVE_AREA_DISTRIBUTION);
    return 2.0 * psd.getRadiusByQuantile(x / 100);
}

/*!
 * \details Gets the diameter from a certain percentile of the volume based PSD.
 * \return A double which is the diameter corresponding to a certain percentile.
 * \param[in] x double [0, 100] which determines the obtained diameter as a percentile of the PSD.
 */
Mdouble PSD::getVolumeDx(Mdouble x) const
{
    logger.assert_always(x >= 0 && x <= 100, "Percentile % is not between 0 and 100.", x);
    PSD psd = *this;
    psd.convertCumulativeNumberDistributionToCumulative(TYPE::CUMULATIVE_VOLUME_DISTRIBUTION);
    return 2.0 * psd.getRadiusByQuantile(x / 100);
}

/*!
 * \details gets the radius from a certain quantile of the PSD
 * \return A double which is the radius corresponding to a certain quantile of the PSD.
 * \param[in] quantile                 double which determines the returned radius as a quantile of the PSD.
 */
Mdouble PSD::getRadiusByQuantile(Mdouble quantile) const
{
    logger.assert_always(quantile <= 1 && quantile >= 0, "quantile is not between 0 and 1");
    // find the quantile corresponding to the PSD
    auto high = std::lower_bound(particleSizeDistribution_.begin(), particleSizeDistribution_.end(), quantile);
    auto low = std::max(particleSizeDistribution_.begin(), high - 1);
    if (high->probability == low->probability)
        return high->internalVariable;
    else
        return low->internalVariable +
               (high->internalVariable - low->internalVariable) * (quantile - low->probability) /
               (high->probability - low->probability);
}

/*!
 * \details Calculates the quantile corresponding to a certain radius.
 * @param radius The radius to find the quantile for.
 * @return The quantile found.
 */
Mdouble PSD::getQuantileByRadius(Mdouble radius) const
{
    auto high = std::lower_bound(particleSizeDistribution_.begin(), particleSizeDistribution_.end(),
                                 radius, [](const DistributionElements& rp, Mdouble value){ return rp.internalVariable < value; });
    auto low = std::max(particleSizeDistribution_.begin(), high - 1);
    if (high->internalVariable == low->internalVariable)
        return high->probability;
    else
        return low->probability + (high->probability - low->probability) * (radius - low->internalVariable) / (high->internalVariable - low->internalVariable);
}

/*!
 * \details Gets a radius such that a monodisperse system has the same number of particles as a polydisperse system.
 * (i.e. mean += p_i * 0.5*(r_i^3 + r_i-1^3)
 * \return A double which corresponds to the volumetric mean radius.
 */
Mdouble PSD::getVolumetricMeanRadius() const
{
    PSD psd = *this;
    psd.convertCumulativeToProbabilityDensity();
    psd.convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
    Mdouble mean = 0;
    for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
        mean += it->probability * 0.5 * (cubic(it->internalVariable) + cubic((it - 1)->internalVariable));
    mean = pow(mean, 1. / 3.);
    return mean;
}

/*!
 * \details Decrement the nParticlesPerClass_ counter of the last inserted class by 1. This is utilized when a
 * particle is generated but failed to be placed into an insertionBoundary.
 */
void PSD::decrementNParticlesPerClass()
{
    --nParticlesPerClass_[index_];
}

/*!
 * \details Decrement the volumePerClass_ counter of the last inserted class by the volume of the particle. This is
 * utilized when a particle is generated but failed to be placed into an insertionBoundary.
 * \param[in] volume        Volume of the particle which failed to be inserted.
 */
void PSD::decrementVolumePerClass(Mdouble volume)
{
    volumePerClass_[index_]-= volume;
}

/*!
 * \details Gets the size ratio (width) of the PSD defined by the ratio of minimum to maximum particle radius.
 * \return A double which corresponds to the size ratio of the PSD.
 */
Mdouble PSD::getSizeRatio() const
{
    Mdouble rMin = getMaxRadius();
    Mdouble rMax = getMinRadius();
    return rMin / rMax;
}

/*!
 * \details Checks if the Size ratio of the PSD is too high and cuts the PSD at head and tail by 10 percent to avoid
 * inaccurate results.
 */
void PSD::cutHighSizeRatio()
{
    Mdouble SR = getSizeRatio();
    if (SR > 100)
    {
        logger(WARN, "Size ratio > 100; Cutting the PSD to avoid inaccurate results");
        cutoffCumulativeDistribution(TYPE::CUMULATIVE_NUMBER_DISTRIBUTION, 0.1, 0.9, 0.5);
    }
}

/*!
 * \details Gets the minimal radius of the PSD.
 * \return A double which corresponds to the minimal radius of the PSD.
 */
Mdouble PSD::getMinRadius() const
{
    return particleSizeDistribution_[0].internalVariable;
}

/*!
 * \details Gets the maximum radius of the PSD.
 * \return A double which corresponds to the maximum radius of the PSD.
 */
Mdouble PSD::getMaxRadius() const
{
    return particleSizeDistribution_.back().internalVariable;
}

/*!
 * \details Gets the vector containing radii and probabilities of the PSD.
 * \return A vector containing radii and probabilities of the PSD.
 */
std::vector<DistributionElements> PSD::getParticleSizeDistribution() const
{
    return particleSizeDistribution_;
}

/*!
 * \details This returns the PSD vector converted to the specified type, with an optional scale factor for the internal
 * variable. This can be used e.g. to output the PSD with a certain type.
 * @param PSDType The type the PSD should be.
 * @param scalingFactor The scale factor by which the the internal variable is multiplied.
 * @return The converted PSD vector.
 */
std::vector<DistributionElements> PSD::getParticleSizeDistributionByType(TYPE PSDType, Mdouble scalingFactor) const
{
    PSD psd = *this;

    switch (PSDType)
    {
        // default case is a CUMULATIVE_NUMBER_DISTRIBUTION
        default:
            break;
        case TYPE::CUMULATIVE_LENGTH_DISTRIBUTION:
            psd.convertCumulativeToProbabilityDensity();
            psd.convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION);
            psd.convertProbabilityDensityToCumulative();
            break;
        case TYPE::CUMULATIVE_AREA_DISTRIBUTION:
            psd.convertCumulativeToProbabilityDensity();
            psd.convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION);
            psd.convertProbabilityDensityToCumulative();
            break;
        case TYPE::CUMULATIVE_VOLUME_DISTRIBUTION:
            psd.convertCumulativeToProbabilityDensity();
            psd.convertProbabilityDensityNumberDistributionToProbabilityDensity(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
            psd.convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_NUMBER_DISTRIBUTION:
            psd.convertCumulativeToProbabilityDensity();
            break;
        case TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION:
            psd.convertCumulativeToProbabilityDensity();
            psd.convertProbabilityDensityNumberDistributionToProbabilityDensity(PSDType);
            break;
        case TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION:
            psd.convertCumulativeToProbabilityDensity();
            psd.convertProbabilityDensityNumberDistributionToProbabilityDensity(PSDType);
            break;
        case TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION:
            psd.convertCumulativeToProbabilityDensity();
            psd.convertProbabilityDensityNumberDistributionToProbabilityDensity(PSDType);
            break;
    }

    auto psdVector = psd.getParticleSizeDistribution();
    for (auto& de : psdVector)
        de.internalVariable *= scalingFactor;

    return psdVector;
}

/*!
 * \details sets the PSD from a vector of the DistributionElements class; used to read in PSDs from a restart file.
 * \param[in] particleSizeDistribution      vector containing the radii and probabilities specifying the PSD.
 */
void PSD::setParticleSizeDistribution(std::vector<DistributionElements> particleSizeDistribution)
{
    particleSizeDistribution_ = particleSizeDistribution;
}

void PSD::scaleParticleSize(double scale)
{
    for (auto& p : particleSizeDistribution_) {
        p.internalVariable *= scale;
    }
}

double PSD::scaleParticleSizeAuto(int numberOfParticles, double targetVolume, bool allowScaleDown)
{
    // Calculate the mean particle volume, using the Mean Value Theorem for Integrals. I.e. the mean particle volume
    // over an infinite number of particles.
    double meanParticleVolume = 0.0;
    for (auto it = particleSizeDistribution_.begin(); it < particleSizeDistribution_.end() - 1; it++)
    {
        double radiusLeft = it->internalVariable;
        double radiusRight = (it+1)->internalVariable;
        double probabilityLeft = it->probability;
        double probabilityRight = (it+1)->probability;

        double meanVolume;
        if (radiusLeft != radiusRight)
        {
            // The CDF is assumed piecewise linear, so the radius is linearly interpolated.
            // Integrate the volume from radius left to right, i.e. integral of 4/3 pi r^3 gives 1/3 pi (r_r^4 - r_l^4).
            double volumeIntegral = (std::pow(radiusRight, 4) - std::pow(radiusLeft, 4)) * constants::pi / 3.0;
            // Calculate the mean volume that represents this bin.
            meanVolume = volumeIntegral / (radiusRight - radiusLeft);
        }
        else
        {
            meanVolume = std::pow(radiusLeft, 3) * constants::pi * 4.0 / 3.0;
        }

        // Calculate the volume that this bin accounts for, and add to total.
        meanParticleVolume += meanVolume * (probabilityRight - probabilityLeft);
    }

    // The volume that the number of particles currently fills, and the scale factor to make it match the target volume.
    double totalParticleVolume = numberOfParticles * meanParticleVolume;
    double scale = std::cbrt(targetVolume / totalParticleVolume);

    // The number of (unscaled) particles needed to fill the target volume.
    unsigned long long numberOfParticlesExpected = std::pow(scale, 3) * numberOfParticles;
    logger(INFO, "Rescaling the particle size based on number of particles %, and target volume %.", numberOfParticles, targetVolume);
    logger(INFO, "Expected number of particles %, scale factor %.", numberOfParticlesExpected, scale);

    // Always scale up, only scale down when requested.
    if (scale > 1.0 || allowScaleDown)
    {
        scaleParticleSize(scale);
        logger(INFO, "Rescaling applied.");
    }
    else
    {
        logger(INFO, "Rescaling not applied, since scale factor is less than or equal to 1.");
        scale = 1.0; // Reassign so that the returned value is correct.
    }

    return scale;
}

/*!
 * \details Gets the number of particles already inserted into the simulation by summing up the particles inserted in
 * each class.
 * \return An integer containing the number of particles already inserted into the simulation.
 */
int PSD::getInsertedParticleNumber() const
{
    int sum = 0;
    for (auto& it: nParticlesPerClass_)
        sum += it;
    return sum;
}

// /*!
// * \details Compute the raw momenta of the inserted PSD by converting it to a PDF, calculating the moments m_i according
// * to \f$ m_i = \sum_{j=0}^n \f$ scaling it by the number
// * of particles inserted into the simulation (moments are computed from a PROBABILITYDENSITY_NUMBER_DISTRIBUTION) and
// * converting it back to a CDF.
// * See <a href="https://en.wikipedia.org/wiki/Moment_(mathematics)#Significance_of_the_moments">Wikipedia</a> for
// * details.
// */
//void PSD::computeRawMomenta()
//{
//    convertCumulativeToProbabilityDensity();
//    for (size_t im = 0; im < momenta_.size(); ++im)
//    {
//        // prevent summing up of moments for each time step of the simulation
//        momenta_[im] = 0;
//        for (auto& it : particleSizeDistribution_)
//        {
//            momenta_[im] += std::pow(it.internalVariable, im) * it.probability;
//        }
//    }
//    // zeroth-moment equals particle number for a PROBABILITYDENSITY_NUMBER_DISTRIBUTION
//    momenta_[0] *= getInsertedParticleNumber();
//    convertProbabilityDensityToCumulative();
//}
//
// /*!
// * \details Compute the central momenta of the PSD from their respective raw momenta.
// */
//void PSD::computeCentralMomenta()
//{
//    computeRawMomenta();
//    Mdouble mean = momenta_[1];
//    momenta_[5] += -5 * mean * momenta_[4] + 10 * mean * mean * momenta_[3]
//                   - 10 * mean * mean * mean * momenta_[2] + 4 * mean * mean * mean * mean * mean;
//    momenta_[4] += -4 * mean * momenta_[3] + 6 * mean * mean * momenta_[2] - 3 * mean * mean * mean * mean;
//    momenta_[3] += -3 * mean * momenta_[2] + 2 * mean * mean * mean;
//    momenta_[2] += -mean * mean;
//}

// /*!
// * \details Compute the standardised momenta of the PSD from their respective central momenta.
// */
//void PSD::computeStandardisedMomenta()
//{
//    computeCentralMomenta();
//    Mdouble std = std::sqrt(momenta_[2]);
//    momenta_[3] /= std * std * std;
//    momenta_[4] /= std * std * std * std;
//    momenta_[5] /= std * std * std * std * std;
//}

/*!
 * \details Gets the radius of a particle from the PSD.
 * \param[in] i     An integer which is the index of the vector containing the radius.
 * \return A double which corresponds to the radius of the particle.
 */
Mdouble PSD::getRadius(int index) const
{
    return particleSizeDistribution_[index].internalVariable;
}

/*!
 * \details Required to use std::lower_bound for finding when the probability is higher than a certain value.
 * \return TRUE if probability from a vector of type DistributionElements is higher than a certain value from a vector of type
 * DistributionElements and FALSE in the opposite case.
 */
bool operator<(const DistributionElements& l, const DistributionElements& r)
{
    return l.probability < r.probability;
}

/*!
 * \details required to use std::lower_bound for finding when the probability provided as a double is higher than a
 * certain value.
 * \return TRUE if probability as double is higher than a certain value from a DistributionElements vector and FALSE in the opposite
 * case.
 */
bool operator<(const DistributionElements& l, const Mdouble prob)
{
    return l.probability < prob;
}

/*!
 * \details Required to use std::distance to find the index of the PSD size class in which a particle has to be inserted
 * \return A double which determines the size class (radius) a particle will be inserted to.
 */
Mdouble operator==(DistributionElements l, const Mdouble r)
{
    return l.internalVariable == r;
}

/*!
 * \details reads from input stream. This function is used for restart files.
 * \return a reference to an input stream.
 */
std::istream& operator>>(std::istream& is, DistributionElements& p)
{
    is >> p.internalVariable;
    is >> p.probability;
    return is;
}

/*!
 * \details Writes to output stream. This function is used for restart files.
 * \return a reference to an output stream.
 */
std::ostream& operator<<(std::ostream& os, DistributionElements& p)
{
    os << p.internalVariable << ' ' << p.probability << ' ';
    return os;
}
