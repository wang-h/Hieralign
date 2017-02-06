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
		int    iterations_;
		bool   forward_;   
		bool   tt_initialized   = false;  
		double threshold_       = 1e-9;  
		size_t thread_buffer_size_ = 5000; 
		size_t num_target_vocab_ = 0;
		size_t num_source_vocab_ = 0;
		
		
		
		void TrainWithBatches(ParallelCorpus &parallelCorpus);
		
		void InitializeTranslationTableEntries(ParallelCorpus &parallelCorpus);
		
		void BatchInitializeTranslationTableEntries(ParallelCorpus &parallelCorpus); 
		inline void AddTranslationOptionsInBatch(vector<vector<unsigned>>& insert_buffer);
		
		void UpdateParamsUsingBatch(const vector<SentencePair>& buffer, const int lid);
		
		
	public:   
		IBM1(); 
		~IBM1();  
		
		TranslationTable ttable;  
		AlignmentList linksList;  
		
		void Config(const size_t& thread_buffer_size, const size_t& num_source_vocab, const size_t& num_target_vocab); 
		
		void TrainModel(ParallelCorpus &parallelCorpus, const bool forward=true, const int iteration=5); 
		
		void ViterbiAlign(ParallelCorpus &corpus); 
		inline void ViterbiAlignSentencePair(const Sentence &source, const Sentence &target, Alignment *links);
		
		void EstimateViterbiProb(ParallelCorpus &corpus, const AlignmentList& link_list, const bool& symmetrized=false); 
		void clear();
};
#endif // IBM1_H
  