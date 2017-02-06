#ifndef VBIBM1_H
#define VBIBM1_H
#include "IBM1.h"
class VBIBM1 :public IBM1 { 
	public:
		VBIBM1();
		~VBIBM1();
		void TrainModel(ParallelCorpus &parallelCorpus, const bool forward, const int iterations, const double alpha); 
};

#endif // VBIBM1_H
