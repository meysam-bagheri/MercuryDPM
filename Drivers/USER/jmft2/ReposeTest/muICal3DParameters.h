#ifndef MERCURYDPM_MUICAL3DPARAMETERS_H
#define MERCURYDPM_MUICAL3DPARAMETERS_H
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <fstream>

#define DEGREES (M_PI / 180.)

class muICal3DParameters {
public:
    double timeStep, timeMax;
    int saveEvery;
    double length, width, height;
    double theta;

    double rho, g;
    double collisionTime, restitutionCoefficient;
    double particleRadius, dispersity;
    double betaslide, betaroll, betators;
    double baseRadius, baseDispersity, baseConc;

    static muICal3DParameters *read(std::string paramsFileName) {
        auto params = muICal3DParameters::presets();

        std::ifstream file(paramsFileName);
        std::string paramName;
        double paramValue;
        while (file >> paramName >> paramValue) {
            /* Parse the line. The first word will be a parameter paramName and the
             * second will be a parameter value. The type of parameter value is a
             * string; depending on which parameter it is, we may have to use atof
             * or atoi to convert it into a number. */
            if (paramName == "timeStep")
                params->timeStep = paramValue;
            else if (paramName == "timeMax")
                params->timeMax = paramValue;
            else if (paramName == "saveEvery")
                params->saveEvery = paramValue;
            else if (paramName == "length")
                params->length = paramValue;
            else if (paramName == "width")
                params->width = paramValue;
            else if (paramName == "height")
                params->height = paramValue;
            else if (paramName == "theta")
                params->theta = paramValue * DEGREES;

            else if (paramName == "rho")
                params->rho = paramValue;
            else if (paramName == "g")
                params->g = paramValue;
            else if (paramName == "collisionTime")
                params->collisionTime = paramValue;
            else if (paramName == "restitutionCoefficient")
                params->restitutionCoefficient = paramValue;

            else if (paramName == "particleRadius")
                params->particleRadius = paramValue;
            else if (paramName == "dispersity")
                params->dispersity = paramValue;
            else if (paramName == "betaslide")
                params->betaslide = paramValue * DEGREES;
            else if (paramName == "betaroll")
                params->betaroll = paramValue * DEGREES;
            else if (paramName == "betators")
                params->betators = paramValue * DEGREES;

            else if (paramName == "baseRadius")
                params->baseRadius = paramValue;
            else if (paramName == "baseDispersity")
                params->baseDispersity = paramValue;
            else if (paramName == "baseConc")
                params->baseConc = paramValue;

            else
                std::cerr << "Warning: Ignored parameter " << paramName << ", value " << paramValue << std::endl;
        }

        fprintf(stderr, "Read the following parameters:\n");
        params->write(stderr);

        file.close();

        return params;
    }

    static muICal3DParameters *presets() {
        auto params = new muICal3DParameters();
        params->timeStep = 1e-5;
        params->timeMax = 10;
        params->saveEvery = 1000;

        params->length = 1;
        params->width = 0.2;
        params->height = 0.25;
        params->theta = 18. * DEGREES;

        params->rho = 1;
        params->g = 1;
        params->collisionTime = 1e-4;
        params->restitutionCoefficient = 0.66;

        params->particleRadius = 0.008;
        params->dispersity = 0.1;
        params->betaslide = 21.8 * DEGREES;
        params->betaroll = 0.0 * DEGREES;
        params->betators = 0.0 * DEGREES;
        params->baseRadius = 0.005;
        params->baseDispersity = 0.1;
        params->baseConc = 1.0;

        return params;
    }

    void write(FILE *f) {
        fprintf(f, "timeStep %g\n", this->timeStep);
        fprintf(f, "timeMax %g\n", this->timeMax);
        fprintf(f, "saveEvery %d\n", this->saveEvery);

        fprintf(f, "length %lf\n", this->length);
        fprintf(f, "width %lf\n", this->width);
        fprintf(f, "height %lf\n", this->height);
        fprintf(f, "theta %lf\n", this->theta * 180 / M_PI);

        fprintf(f, "rho %f\n", this->rho);
        fprintf(f, "g %f\n", this->g);
        fprintf(f, "collisionTime %lf\n", this->collisionTime);
        fprintf(f, "restitutionCoefficient %lf\n",
                this->restitutionCoefficient);

        fprintf(f, "particleRadius %g\n", this->particleRadius);
        fprintf(f, "dispersity %lf\n", this->dispersity);
        fprintf(f, "betaslide %lf\n", this->betaslide * 180 / M_PI);
        fprintf(f, "betaroll %lf\n", this->betaroll * 180 / M_PI);
        fprintf(f, "betators %lf\n", this->betators * 180 / M_PI);

        fprintf(f, "baseRadius %g\n", this->baseRadius);
        fprintf(f, "baseDispersity %lf\n", this->baseDispersity);
        fprintf(f, "baseConc %lf\n", this->baseConc);
    }
};

#endif
