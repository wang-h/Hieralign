//　Copyright 2013 by Hao Wang  
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef IBM1_H
#define IBM1_H  
#include "Sentence.h" 
#include "ParallelCorpus.h"
#include "Alignment.h"
#include "TranslationTable.h"
  
class IBM1 {
	protected:   
		unsigned    iterations_;
		bool        forward_;    
		bool        with_null_;
		 
		size_t thread_buffer_size_ = 5000; 
		
		
		void BatchInitializeTranslationTableEntries(ParallelCorpus &parallelCorpus); 
		inline void AddTranslationOptionsInBatch(vector<vector<unsigned>>& insert_buffer);
		
		
		void TrainWithBatches(ParallelCorpus &parallelCorpus); 
		void UpdateParamsUsingBatch(const vector<SentencePair>& buffer);
		
		
	public:   
		IBM1(); 
		~IBM1();  
		
		TranslationTable ttable;  
		AlignmentList linksList;  
		void Config(const bool& forward, const unsigned& iterations, const bool& withNull, const size_t& size_source_vocab);
		void TrainModel(ParallelCorpus &parallelCorpus); 
		
		void ViterbiAlign(ParallelCorpus &corpus); 
		inline void ViterbiAlignSentencePair(const Sentence &source, const Sentence &target, Alignment *links);
		
		void EstimateViterbiProb(ParallelCorpus &corpus, const AlignmentList& link_list, const bool& symmetrized=false); 
};
#endif // IBM1_H
  