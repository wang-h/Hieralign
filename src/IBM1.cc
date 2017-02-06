#include "IBM1.h" 
IBM1::IBM1() {  
} 
IBM1::~IBM1() {
}

void IBM1::Config(const size_t& thread_buffer_size, const size_t& num_source_vocab, const size_t& num_target_vocab) { 
	thread_buffer_size_ = thread_buffer_size;  
	num_target_vocab_   = num_target_vocab; 
	num_source_vocab_   = num_source_vocab; 
	if (!forward_)
		ttable.SetMaxF(static_cast<unsigned>(num_source_vocab));
	else
		ttable.SetMaxF(static_cast<unsigned>(num_target_vocab));	
}
void IBM1::TrainModel(ParallelCorpus &parallelCorpus, const bool forward, const int iterations) {
	forward_ = forward;
	iterations_ = iterations; 
	int currentRound = 0;
	while (currentRound < iterations_) { 
		cerr<< "iteration "<< currentRound << "..."<<endl;
		if(currentRound==0 && ! tt_initialized){
			//InitializeTranslationTableEntries(parallelCorpus);
			BatchInitializeTranslationTableEntries(parallelCorpus);
			tt_initialized = true; 
			ttable.Freeze();
		}
		TrainWithBatches(parallelCorpus);
		ttable.Normalize(); 
	 
		currentRound++;  
	}
}



void IBM1::TrainWithBatches(ParallelCorpus &parallelCorpus) {   
	vector<SentencePair> buffer; 
	for (size_t i = 0; i < parallelCorpus.size(); i++){ 
		SentencePair sentence_pair = parallelCorpus[i];  
		if (i % 5000 == 0 && i!=0){
			cerr << '.';  
		}
		if (i % 100000 == 0 && i!=0) { 
			cerr << "\t[" << i << "]\n" << flush;  
		}
		if (i == (parallelCorpus.size()-1)) { 
			cerr << "\nTotal:            \t[" << i << "]\n" << flush;  
		} 
		buffer.push_back(sentence_pair);
		if (buffer.size() >= thread_buffer_size_ || i == (parallelCorpus.size()-1)) { 
			UpdateParamsUsingBatch(buffer, i); 
			buffer.clear();
		} 
	} 
}

inline void IBM1::AddTranslationOptionsInBatch(vector<vector<unsigned> >& insert_buffer) {
	ttable.SetMaxF_count(insert_buffer.size()+1);
	#pragma omp parallel for schedule(dynamic)
	for (unsigned f = 0; f < insert_buffer.size(); f++) {
		for (unsigned e : insert_buffer[f]) {
			ttable.Insert(f, e);
		}
		insert_buffer[f].clear();
	}
}
void IBM1::BatchInitializeTranslationTableEntries(ParallelCorpus &parallelCorpus){
	vector<vector<unsigned>> insert_buffer; 
	double insert_buffer_items = 0.0;
	for (size_t i = 0; i < parallelCorpus.size(); i++){ 
		SentencePair sentence_pair = parallelCorpus[i];  
		Sentence *src=sentence_pair.first;
		Sentence *trg=sentence_pair.second; 
		if (!forward_)
			swap(src,trg);
		for (const unsigned f_i : *src) { 
			if (f_i >= insert_buffer.size()) {
				insert_buffer.resize(f_i+1);
			} 
			for (size_t j=1;j<trg->size();j++) { 
				const unsigned& e_j = (*trg)[j];  
				insert_buffer[f_i].push_back(e_j);
			}
			insert_buffer_items += trg->size();
			if (insert_buffer_items >= thread_buffer_size_*100 || i == (parallelCorpus.size()-1)) { 
				AddTranslationOptionsInBatch(insert_buffer);
				insert_buffer.clear();
				insert_buffer_items = 0.0;
			}
		}  
		  
	}
}  
void IBM1::InitializeTranslationTableEntries(ParallelCorpus &parallelCorpus){
	// single thread, very slow
	for (size_t i = 0; i < parallelCorpus.size(); i++){ 
		SentencePair sentence_pair = parallelCorpus[i];  
		Sentence *src=sentence_pair.first;
		Sentence *trg=sentence_pair.second; 
		if (!forward_){
			swap(src,trg); 
		} 
		for (const unsigned f_i : *src) { 
			for (size_t j=1;j<trg->size();j++) { 
				const unsigned& e_j = (*trg)[j];  
				ttable.Insert(f_i, e_j);
			} 
		}  
	}
}  
void IBM1::UpdateParamsUsingBatch(const vector<SentencePair>& buffer, const int lid){
	#pragma omp parallel for schedule(dynamic)
	for (int line_idx = 0; line_idx < static_cast<int>(buffer.size()); line_idx++) { 
		Sentence *src=buffer[line_idx].first;
		Sentence *trg=buffer[line_idx].second;
		if (!forward_)
			swap(src,trg);
		vector<double> probs(src->size());
		if (src->size() <= 1 || trg->size() <= 1) {
			cerr << "#Warning: empty line in " << lid+line_idx <<  endl; 
			continue;
		}
		for (size_t j = 1; j < trg->size(); j++) {
			const unsigned& e_j = (*trg)[j];
			double sum = 0;
			double prob_a_i = 1.0 / (src->size());  // uniform (model 1) 
			for (size_t i = 0; i <src->size(); i++) {  
				const unsigned& f_i = (*src)[i];
				probs[i] = ttable.prob(f_i, e_j) * prob_a_i;
				sum += probs[i];
			}
			for (unsigned i = 0; i <src->size(); i++) {
				const double p = probs[i] / sum;
				const unsigned& f_i = (*src)[i];
				ttable.Increment(f_i, e_j, p);   
			}
		}
		
	} 
}

void IBM1::clear() { 
	linksList.clear(); 
	ttable.clear();   
	tt_initialized = false;
} 
 
void IBM1::EstimateViterbiProb(ParallelCorpus &parallelCorpus, const AlignmentList& link_list, const bool& symmetrized){  

	ttable.ClearCounts();   
	CHECK(parallelCorpus.size()==link_list.size(),"#ERROR, not same size"); 
	#pragma omp parallel for schedule(dynamic) 
	for (size_t k = 0; k < parallelCorpus.size(); k++){   
		SentencePair sentence_pair = parallelCorpus[k];   
		Sentence *src=sentence_pair.first;
		Sentence *trg=sentence_pair.second; 
		const Alignment& links=link_list[k]; 
		if (!forward_){
			swap(src,trg); 
		}    
		for(auto const& point: links){
			// alignment is in e-f direction
			 
			unsigned i = point.second;
			unsigned j = point.first; 
			if (symmetrized)
				swap(i,j); 
			const unsigned f = (*src)[i];   
			const unsigned e = (*trg)[j]; 
			ttable.Increment(f, e, 1.0);
		}  
	} 
	ttable.Normalize();
}
 
void IBM1::ViterbiAlign(ParallelCorpus &parallelCorpus){ 
	linksList.resize(parallelCorpus.size());
	int count=0;  
	#pragma omp parallel for schedule(dynamic) reduction(+:count)
	for (size_t i = 0; i < parallelCorpus.size(); i++){  
		Alignment links; 
		SentencePair sentence_pair = parallelCorpus[i];   
		Sentence *src=sentence_pair.first;
		Sentence *trg=sentence_pair.second; 
		if (!forward_){
			swap(src,trg); 
		}  
		ViterbiAlignSentencePair(*src, *trg, &links);
		linksList[i]= links;
		count++; 
	}  
}
inline void IBM1::ViterbiAlignSentencePair(const Sentence &sourceSentence, const Sentence &targetSentence, Alignment *links){ 
	for (size_t j = 1; j < targetSentence.size(); j++){
		const unsigned e     = targetSentence[j]; 
		unsigned best_candidate = 0;
		double best_score =-DBL_MAX; 
		for (size_t i = 0; i < sourceSentence.size(); i++){ 
			const unsigned f = sourceSentence[i];     
			const double   score = ttable.prob(f,e);  
			if (score > best_score) {
				best_score = score;
				best_candidate = i;  
			}    
		} 
		// target_index-source-index
		// (j-1, best_candidate-1) is the true indices
		links->emplace_back(j, best_candidate); 
		
	} 
}