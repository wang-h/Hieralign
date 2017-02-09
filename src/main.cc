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

#include <stdio.h> 
#include <string.h> 
#include "Utils.h"
#include "Alignment.h" 
#include "SimpleWordWrapper.h"
#include "ParallelCorpus.h" 
#include "TranslationTable.h"
#include "IBM1.h"
#include "VBIBM1.h"
#include "Heuristic.h"
#include "Hieralign.h"  
 
int main(int argc, char **argv)
{
	string input      = ""     ;
	string train      = ""     ; 
	string model_file = ""     ;
	string iterations = "5"    ;
	string prior      = "0.01" ;
	string beamsize   = "10"   ;
	string theta      = "1"    ; 
	string delta      = "0"    ; 
	string threshold  = "1e-7" ;  
	string p0         = "1e-4" ;  
	string model      = "ibm1_vbh"; 
	for (int i = 1; i < argc; i++) {
			if (      strcmp(argv[i], "--input"     ) == 0 && i + 1 < argc) {
				input = argv[++i]; 
			}else if (strcmp(argv[i], "--train"     ) == 0 && i + 1 < argc) {
				train = argv[++i]; 
			}else if (strcmp(argv[i], "--model"     ) == 0 && i + 1 < argc) {
				model = argv[++i]; 
			}else if (strcmp(argv[i], "--model_file") == 0 && i + 1 < argc) {
				model_file = argv[++i]; 
			}else if (strcmp(argv[i], "--theta"     ) == 0 && i + 1 < argc) {
				theta = argv[++i];  
			}else if (strcmp(argv[i], "--delta"     ) == 0 && i + 1 < argc) {
				delta = argv[++i];  
			}else if (strcmp(argv[i], "--prior"     ) == 0 && i + 1 < argc) {
				prior = stod(argv[++i]);  
			}else if (strcmp(argv[i], "--threshold" ) == 0 && i + 1 < argc) {
				threshold = argv[++i];
			}else if (strcmp(argv[i], "--p0"        ) == 0 && i + 1 < argc) {
				p0 = argv[++i];
			}else if (strcmp(argv[i], "--beamsize"  ) == 0 && i + 1 < argc) {
				beamsize = argv[++i];   
			}else if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
				iterations = atoi(argv[++i]); 
			} else {
				cerr << "usage: " << argv[0] << " [options]...\n"; 
				cerr << "version: 0.2\n"; 
				cerr << "1. as an aligner:\n\n";  
				cerr << "  --train      <string>  : input file, need to align.            \t["        << input     << "]\n";   
				cerr << "  --iterations <integer> : number of training iterations.        \t[default="<< iterations<< "]\n"; 
				cerr << "  --beamsize   <integer> : beam size during parsing.             \t[default="<< beamsize  << "]\n";  
				cerr << "  --threshold  <double>  : threshold to filter [f,e].            \t[default="<< threshold << "]\n"; 
				cerr << "  --p0         <double>  : probability for null, [f,e].          \t[default="<< p0        << "]\n"; 
				cerr << "  --model      <string>  : model used,                           \t["        << model     << "]\n"  
					 << "                           e.g., 1. ibm1,      2. vbibm1     \n"
					 << "                                 3. ibm1_vb,   4. vbibm1_vb  \n" 
					 << "                                 5. ibm1_vbh,  6. vbibm1_vbh \n"; 
				cerr << "  --prior      <double>  : hyper-parameter of Dirichlet prior.   \t[default="<< prior     << "]\n"
					 << "                           only works when using Variational Bayes.\n";
				cerr << "  --theta      <double>  : hyper parameter for translation model.\t[default="<< theta     << "]\n";
				cerr << "  --delta      <double>  : hyper parameter for distortion model. \t[default="<< delta     << "]\n\n";  
 
				cerr << "2. train the model in advance:\n\n"; 
				cerr << "  --train      <string>  : input training file.                 \t["         << train      << "]\n"; 
				cerr << "  --model_file <string>  : write model to file.                 \t["         << model_file << "]\n";  
				
				cerr << "  --iterations <integer> : number of training iterations.       \t[default=" << iterations << "]\n";  
				cerr << "  --prior      <double>  : hyper-parameter of Dirichlet prior.  \t[default=" << prior      << "]\n"
					 << "                           only works when using Variational Bayes.\n";
				cerr << "  --threshold  <double>  : threshold to filter tt[f,e].         \t[default=" << threshold  << "]\n"; 
				cerr << "  --model      <string>  : model used,                          \t[default=" << model      << "]\n" 
					 << "                           e.g., 1. ibm1,      2. vbibm1    \n"
					 << "                                 3. ibm1_vb,   4. vbibm1_vb \n" 
					 << "                                 5. ibm1_vbh,  6. vbibm_vbh \n\n"; 
				
	 
				cerr << "3. align bitext as an online aligner:\n\n"; 
				cerr << "  --input      <string>  : input file, need to align.           \t["         << input      << "]\n"; 
				cerr << "  --model_file <string>  : load model from file.                \t["         << model_file << "]\n"; 
				cerr << "  --beamsize   <integer> : beam size during parsing.             \t[default="<< beamsize  << "]\n"; 
				cerr << "  --theta      <double>  : hyper parameter for translation model.\t[default="<< theta     << "]\n";
				cerr << "  --delta      <double>  : hyper parameter for distortion model. \t[default="<< delta     << "]\n";
				cerr << "  --p0         <double>  : probability for null, [f,e].          \t[default="<< p0        << "]\n"; 
				
				return 1;
			}
	}
	CHECK(!train.empty()||!input.empty(), "Flag  --train or --input should not be empty");  
	SimpleWordWrapper sw2id, tw2id;
	ParallelCorpus trainCorpus,testCorpus; 
	TranslationTable tt;  
	Hieralign aligner(beamsize); 
	
	if (!train.empty()){ 
		//load the training corpus into memory. 
		ifstream input_train(train);
		trainCorpus.ReadParallelCorpus(&input_train, &sw2id, &tw2id);
		input_train.close(); 
		sw2id.Frozen();
		tw2id.Frozen();
		//need initialize the model 
		cerr <<  "Training:       \t"  << train             << endl;
		cerr <<  "lines:          \t[" << trainCorpus.size()<<  "]" << endl; 
		cerr <<  "source vocabs:  \t[" << sw2id.size()      <<  "]" << endl;
		cerr <<  "target vocabs:  \t[" << tw2id.size()      <<  "]" << endl;
		cerr <<  "threshold:      \t[" << threshold         <<  "]" << endl;
		
		if(!model.compare("ibm1")){
			// # case 1 
			IBM1 ibm1_forward, ibm1_backward;  
			cerr<<  "# initializing IBM model1 probabilities using EM algorithm."<<endl; 
			
			ibm1_forward.Config(true, 5, true, sw2id.size()); 
			//forward=true, iterations=5, withNull=true, num_source_words
			cerr << "# training IBM1 in forward direction ..." <<endl; 
			ibm1_forward.TrainModel( trainCorpus);  //training 5 forward loops. 
			 
			ibm1_backward.Config(false, 5, true, tw2id.size()); 
			//forward=false, iterations=5, withNull=true, num_target_words
			cerr << "# training IBM1 in backward direction ..." <<endl;
			ibm1_backward.TrainModel(trainCorpus); //training 5 backward loops. 
 
			ibm1_forward.ttable.Merge(ibm1_backward.ttable);  
			tt=ibm1_forward.ttable; 
			tt.Freeze();
		}
		else if(!model.compare("vbibm1")){
			// # case 2 
			VBIBM1 vbibm1_forward, vbibm1_backward;  
			cerr<<  "# initializing IBM model1 probabilities using Variational Bayes."<<endl; 
			
			vbibm1_forward.Config(true, 5, true, sw2id.size()); 
			cerr << "# training VBIBM1 in forward direction ..." <<endl; 
			vbibm1_forward.TrainModel(trainCorpus, stod(prior)); 
			
			vbibm1_backward.Config(false, 5, true, tw2id.size()); 
			cerr << "# training VBIBM1 in backward direction ..." <<endl; 
			vbibm1_backward.TrainModel(trainCorpus, stod(prior)); 
			
			vbibm1_forward.ttable.Merge(vbibm1_backward.ttable); 
			tt=vbibm1_forward.ttable;
			tt.Freeze();
		} 
		else if (!model.compare("ibm1_vb")){ 
			// # case 3 
			IBM1 ibm1_vb_forward, ibm1_vb_backward;  
			cerr<<  "# initializing IBM model1 Viterbi probabilities using EM algorithm."<<endl; 
			
			cerr << "# training IBM1 Viterbi in forward direction ..." <<endl; 
			ibm1_vb_forward.Config(true, 5, true, sw2id.size()); 
			ibm1_vb_forward.TrainModel(trainCorpus); 
			ibm1_vb_forward.ViterbiAlign(trainCorpus);
			ibm1_vb_forward.EstimateViterbiProb(trainCorpus,  ibm1_vb_forward.linksList);  
			
			
			cerr << "# training IBM1 Viterbi in backward direction ..." <<endl; 
			ibm1_vb_backward.Config(false, 5, true, sw2id.size()); 
			ibm1_vb_backward.TrainModel(trainCorpus); 
			ibm1_vb_backward.ViterbiAlign(trainCorpus);
			ibm1_vb_backward.EstimateViterbiProb(trainCorpus, ibm1_vb_backward.linksList);  
			
			ibm1_vb_forward.ttable.Merge(ibm1_vb_backward.ttable); 
			tt=ibm1_vb_forward.ttable;
			tt.Freeze();
		}
		else if (!model.compare("vbibm1_vb")){ 
			// # case 4 
			VBIBM1 vbibm1_vb_forward, vbibm1_vb_backward; 
			vbibm1_vb_forward.Config(true, 5, true, sw2id.size());  
			cerr<<  "# initializing IBM1 Viterbi probabilities using Variational Bayes."<<endl; 
			
			cerr << "# training VBIBM1 Viterbi in forward direction ..." <<endl; 
			vbibm1_vb_forward.TrainModel(trainCorpus, stod(prior)); 
			vbibm1_vb_forward.ViterbiAlign(trainCorpus);
			vbibm1_vb_forward.EstimateViterbiProb(trainCorpus, vbibm1_vb_forward.linksList);  
			
			vbibm1_vb_backward.Config(false, 5, true, sw2id.size()); 
			cerr << "# training VBIBM1 Viterbi in backward direction ..." <<endl; 
			vbibm1_vb_backward.TrainModel(trainCorpus, stod(prior)); 
			vbibm1_vb_backward.ViterbiAlign(trainCorpus);
			vbibm1_vb_backward.EstimateViterbiProb(trainCorpus, vbibm1_vb_backward.linksList);  
			
			vbibm1_vb_forward.ttable.Merge(vbibm1_vb_backward.ttable); 
			tt=vbibm1_vb_forward.ttable;
			tt.Freeze();
		}
		else if (!model.compare("ibm1_vbh")){ 
			// # case 5 
			IBM1 ibm1_vbh_forward, ibm1_vbh_backward;  
			cerr<<  "# initializing IBM1 Viterbi+heuristic probabilities using EM algorithm."<<endl; 
			ibm1_vbh_forward.Config(true, 5, true, sw2id.size());  
			cerr << "# training IBM1 Viterbi+heuristic in forward direction ..." <<endl; 
			ibm1_vbh_forward.TrainModel(trainCorpus); 
			ibm1_vbh_forward.ViterbiAlign(trainCorpus);  
			
			ibm1_vbh_backward.Config(false, 5, true, sw2id.size());
			cerr << "# training IBM1 Viterbi+heuristic in backward direction ..." <<endl; 
			ibm1_vbh_backward.TrainModel(trainCorpus); 
			ibm1_vbh_backward.ViterbiAlign(trainCorpus);   
			
			AlignmentList symmetric_list;  
			SymmetrizeAlignment(ibm1_vbh_forward.linksList, ibm1_vbh_backward.linksList, &symmetric_list); 
			ibm1_vbh_forward.EstimateViterbiProb( trainCorpus, symmetric_list); 
			ibm1_vbh_backward.EstimateViterbiProb(trainCorpus, symmetric_list, true); 
			ibm1_vbh_forward.ttable.Merge(ibm1_vbh_backward.ttable);   
			tt=ibm1_vbh_forward.ttable;  
			tt.Freeze();
		}
		else if (!model.compare("vbibm1_vbh")){ 
			// # case 6 
			VBIBM1 vbibm1_vbh_forward, vbibm1_vbh_backward; 
			
			cerr<<  "# initializing IBM1 Viterbi+heuristic probabilities using Variational Bayes."<<endl; 
			vbibm1_vbh_forward.Config(true, 5, true, sw2id.size());  
			cerr << "# training VBIBM1 Viterbi+heuristic in forward direction ..." <<endl; 
			vbibm1_vbh_forward.TrainModel( trainCorpus, stod(prior)); 
			vbibm1_vbh_forward.ViterbiAlign(trainCorpus);  
			
			vbibm1_vbh_backward.Config(false, 5, true, sw2id.size());
			cerr << "# training VBIBM1 Viterbi+heuristic in backward direction ..." <<endl; 
			vbibm1_vbh_backward.TrainModel(trainCorpus, stod(prior)); 
			vbibm1_vbh_backward.ViterbiAlign(trainCorpus);  
			
			AlignmentList symmetric_list;
			SymmetrizeAlignment(vbibm1_vbh_forward.linksList, vbibm1_vbh_backward.linksList, &symmetric_list);
			
			vbibm1_vbh_forward.EstimateViterbiProb( trainCorpus, symmetric_list);
			vbibm1_vbh_backward.EstimateViterbiProb(trainCorpus, symmetric_list, true);
			vbibm1_vbh_forward.ttable.Merge(vbibm1_vbh_backward.ttable); 
			tt=vbibm1_vbh_forward.ttable; 
			tt.Freeze();
		}
		else{
			cerr<<  "#ERROR! unknown model type."<<endl;
			
		}
		if (!model_file.empty()){  
				ofstream dump_model(model_file); 
				tt.WriteTranslationTable(&dump_model, sw2id, tw2id, stod(threshold), false);
		} 
	}  
	
	if (!input.empty()){
		//load the input corpus into memory. 
		cerr <<  "# reading input ..."<<endl; 
		ifstream input_test(input); 
		testCorpus.ReadParallelCorpus(&input_test, &sw2id, &tw2id); 
		input_test.close(); 
		cerr <<  "Input:           \t"  << input                    << endl;
		cerr <<  "lines:           \t[" << testCorpus.size()<<  "]" << endl; 
		cerr <<  "Beam size:       \t[" << beamsize         <<  "]" << endl;  
		if (!model_file.empty()){ 
			
			cerr <<"Model file:     \t" << model_file               << endl;  
			ifstream input_model(model_file);
			tt.LoadTranslationTable(&input_model, sw2id, tw2id); 
			tt.Freeze();
		}  
	}
	
	
	if (!input.empty()){
		aligner.LoadTranslationTable(tt); 
		cerr<<  "# start to align input bitext ..."<<endl;
		aligner.setHyperParameters(stod(theta), stod(delta), stod(p0)); 
		aligner.Align(testCorpus);
	}
	else if (model_file.empty()){
		aligner.LoadTranslationTable(tt); 
		cerr<<  "# start to align training bitext ..."<<endl;
		aligner.setHyperParameters(stod(theta), stod(delta), stod(p0)); 
		aligner.Align(trainCorpus); 
	}
	return 0;
}

