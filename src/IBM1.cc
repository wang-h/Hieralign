#include "IBM1.h" 
IBM1::IBM1() {  
} 
IBM1::~IBM1() {
}

void IBM1::Config(const bool& forward, const unsigned& iterations, const bool& withNull, const size_t& size_source_vocab) {
	forward_ = forward;
	iterations_ = iterations; 
	with_null_ = withNull;
	ttable.SetMaxF(size_source_vocab);
	//size_source_vocab {0:unk, 1:NULL, 2:vocab}
}
void IBM1::TrainModel(ParallelCorpus &parallelCorpus) { 
	unsigned currentRound = 0;
	while (currentRound < iterations_) { 
		cerr<< "iteration "<< currentRound << " ..."<<endl;
		if(currentRound==0){ 
			BatchInitializeTranslationTableEntries(parallelCorpus); 
		}
		TrainWithBatches(parallelCorpus);
		if(currentRound==0){
			ttable.SetInitialized();
		}
		ttable.Normalize();  
		currentRound++;  
	}
}

void IBM1::BatchInitializeTranslationTableEntries(ParallelCorpus &parallelCorpus){ 
	vector<vector<unsigned>> insert_buffer; 
	double insert_buffer_items = 0.0;
	size_t last= parallelCorpus.size()-1;
	for (size_t i = 0; i < parallelCorpus.size(); i++){ 
		SentencePair sentence_pair = parallelCorpus[i];  
		Sentence *src=sentence_pair.first;
		Sentence *trg=sentence_pair.second; 
		if (!forward_)
			swap(src,trg);
		// in training, f will never be 0. f["NULL"]=1 at least;  
		for (const unsigned& e : *trg) {  
			if(with_null_){  
				if (1 >= insert_buffer.size()) 
					insert_buffer.resize(2);
				insert_buffer[1].push_back(e);
			}
			insert_buffer_items += trg->size();
			for (unsigned& f : *src) { 
				if (f >= insert_buffer.size())
					insert_buffer.resize(f+1);
				insert_buffer[f].push_back(e); 
			}
		}
		if (insert_buffer_items >= thread_buffer_size_*100 || i== last) { 
			AddTranslationOptionsInBatch(insert_buffer);
			insert_buffer.clear();
			insert_buffer_items = 0.0; 
		} 
	}   
}  

inline void IBM1::AddTranslationOptionsInBatch(vector<vector<unsigned> >& insert_buffer) { 
	#pragma omp parallel for schedule(dynamic)
	for (unsigned f = 0; f < insert_buffer.size(); f++) {
		for (unsigned e : insert_buffer[f]) {
			ttable.Insert(f, e);
		}
		insert_buffer[f].clear();
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
		if (buffer.size() >= thread_buffer_size_ ) { 
			UpdateParamsUsingBatch(buffer); 
			buffer.clear();
		} 
	}
	UpdateParamsUsingBatch(buffer); 
	buffer.clear();
}


void IBM1::UpdateParamsUsingBatch(const vector<SentencePair>& buffer ){
	#pragma omp parallel for schedule(dynamic)
	for (int line_idx = 0; line_idx < static_cast<int>(buffer.size()); line_idx++) { 
		Sentence *src=buffer[line_idx].first;
		Sentence *trg=buffer[line_idx].second;
		if (!forward_)
			swap(src,trg);
		vector<double> probs(src->size()+1);
		if (src->size() < 1 || trg->size() < 1) { 
			continue;
		} 
		for (const unsigned& e: *trg) { 
			double sum = 0.0;
			const double m = (with_null_) ? src->size()+1:src->size();
			double prob_a_i = 1.0 /m;  // uniform (model 1) 
			if(with_null_){
				probs[0] = ttable.Prob(1, e) * prob_a_i;
				sum += probs[0];
			}
			for (size_t i = 0; i <src->size(); i++) {  
				const unsigned& f = (*src)[i];
				probs[i+1] = ttable.Prob(f, e) * prob_a_i;
				sum += probs[i+1]; 
			} 
				
			for (unsigned i = 0; i <src->size(); i++) {
				const double p = probs[i+1] / sum;
				const unsigned& f = (*src)[i];
				ttable.Increment(f, e, p);   
			}
			if(with_null_){
				const double p0 = probs[0] / sum;
				ttable.Increment(1, e, p0);   
			}
		}
		
	} 
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
		if(!forward_){
			swap(src,trg); 
		} 
		const double m = src->size(); 
		for(auto const& point: links){
			// alignment is in e-f direction
			unsigned j = point.first; 
			unsigned i = point.second;
			if(symmetrized)
				swap(i,j); 
			const unsigned e = (*trg)[j]; 
			if (i>=m){
				ttable.Increment(1, e, 1.0);
			}  
			else{
				const unsigned f = (*src)[i];    
				ttable.Increment(f, e, 1.0); 
			}
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
	const unsigned null_index = static_cast<unsigned>(sourceSentence.size());
	for (size_t j = 0; j < targetSentence.size(); j++){
		const unsigned e     = targetSentence[j]; 
		unsigned best_candidate = 0;
		double best_score =-DBL_MAX; 
		if(with_null_){  
			const double score = ttable.Prob(1, e);  
			if (score > best_score) {
				best_score = score;
				best_candidate = null_index;  
			} 
		}
		for (size_t i = 0; i < sourceSentence.size(); i++){ 
			const unsigned f = sourceSentence[i];     
			const double   score = ttable.Prob(f,e);  
			if (score > best_score) {
				best_score = score;
				best_candidate = i;  
			}    
		}  
		// target_index-source-index
		// (j, best_candidate) is the true indices
		links->emplace_back(j, best_candidate); 
		
	} 
}