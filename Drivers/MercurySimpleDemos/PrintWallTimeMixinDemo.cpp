#include <thread>

#include "Mercury3D.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Mixins/PrintWallTimeMixin.h"


class PrintWallTimeMixinDemo : public Mercury3D, public PrintWallTimeMixin {
public:
    PrintWallTimeMixinDemo(long delayMs) {
        this->delayMs_ = delayMs;
    }

    void actionsAfterTimeStep() override {
        std::this_thread::sleep_for(std::chrono::milliseconds(delayMs_));
    }

private:
    long delayMs_;
};


int main(const int argc, char *argv[]) {
    long delayMs = argc > 1 ? std::stol(argv[1]) : 10;
    auto sim = new PrintWallTimeMixinDemo(delayMs);
    sim->setName("PrintWallTimeMixinDemo");
    sim->setMin(0, 0, 0);
    sim->setMax(1, 1, 1);
    sim->setTimeMax(1);
    sim->setTimeStep(0.01);
    sim->setSaveCount(10);

    sim->speciesHandler.copyAndAddObject(new LinearViscoelasticSpecies());
    if (argc > 1) {
        argv[1] = argv[0];
        sim->solve(argc - 1, argv + 1);
    } else
        sim->solve();

    //deallocate memory
    delete sim;

    return 0;
}
