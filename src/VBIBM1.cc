#include "VBIBM1.h"

VBIBM1::VBIBM1():IBM1()
{
	
}

VBIBM1::~VBIBM1()
{
}

void VBIBM1::TrainModel(ParallelCorpus &parallelCorpus, const bool forward, const int iterations, const double alpha) {
	forward_ = forward;
	iterations_ = iterations; 
	int currentRound = 0;
	while (currentRound < iterations_) {  
		cerr<< "iteration "<< currentRound << "..."<<endl;
		if(currentRound==0 && ! tt_initialized){
			BatchInitializeTranslationTableEntries(parallelCorpus);
			tt_initialized = true; 
			ttable.Freeze();
		}
		TrainWithBatches(parallelCorpus);
		ttable.NormalizeVB(alpha); 
		currentRound++;  
	}
}