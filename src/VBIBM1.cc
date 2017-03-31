#include "VBIBM1.h"

VBIBM1::VBIBM1() : IBM1()
{
}

VBIBM1::~VBIBM1()
{
}

void VBIBM1::TrainModel(ParallelCorpus &parallelCorpus, const double &alpha)
{

    unsigned currentRound = 0;
    while (currentRound < iterations_)
    {
        cerr << "iteration " << currentRound << "..." << endl;
        if (currentRound == 0)
        {
            BatchInitializeTranslationTableEntries(parallelCorpus);
        }
        TrainWithBatches(parallelCorpus);
        if (currentRound == 0)
        {
            ttable.SetInitialized();
        }
        ttable.NormalizeVB(alpha);
        currentRound++;
    }
}