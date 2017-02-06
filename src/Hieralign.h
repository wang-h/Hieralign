//　Copyright 2017 by Hao Wang  
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
#ifndef HIERALIGN_H
#define HIERALIGN_H 

#include "Utils.h"
#include "Alignment.h"
#include "Matrix.h"
#include "ParallelCorpus.h" 
#include "TranslationTable.h"

class Hieralign 
{
private:
		int beam_width = 10;  
		
		struct ParserAction{
			int x;
			int y;
			bool option; 
			double score;
			ParserAction(int x_, int y_, bool option_, double score_=0.0) : x(x_), y(y_), option(option_), score(score_){};
			bool operator==(const ParserAction& rhs) const {
				return  (x == rhs.x  && y == rhs.y && option == rhs.option) ; 
			} 
			bool operator<(const ParserAction& rhs) const {
					return (score > rhs.score);
			}
		}; 
		struct TranslationEntry{
			unsigned x;
			unsigned y;
			double score;
			TranslationEntry(int x_, int y_, double score_) : x(x_), y(y_), score(score_){};
		};
		struct ParserBlock{
			int i1;
			int j1;
			int i2;
			int j2; 
			//type {0, 1}
			//0: terminal, 1: nonterminal
			ParserBlock(int i1_, int j1_, int i2_, int j2_): i1(i1_), j1(j1_), i2(i2_), j2(j2_){}
			bool operator==(const ParserBlock &rhs) const {
				return (rhs.i1 == i1 && rhs.j1 == j1 && rhs.i2 == i2 && rhs.j2 == j2);
			}
		};
		
		struct ParserState{
			public:
				double score;
				bool   terminal;
				vector<ParserBlock> stack;  // Stack of open blocks.
				vector<ParserAction> actions;  // History of actions. 
				explicit ParserState(int ls, int lt): score(0.0), terminal(false) {
					if (ls > 1 && lt> 1) {
						stack.emplace_back(0, 0, ls, lt);
					} 
				} 
//				ParserState(const ParserState &state, const ParserAction &action)
//					: score(state.score), stack(state.stack), actions(state.actions), terminal(false){
//					Advance(action,  );
//				} 
				ParserState(const ParserState &state)
					: score(state.score), stack(state.stack), actions(state.actions), terminal(state.terminal){ 
				} 
				void Advance(const ParserAction &action){
					//score *=  action.score;
					//score *= action.score;
					score += log(action.score); 
					const ParserBlock block = stack.back();
					stack.pop_back();
					
					const int i = action.x;
					const int j = action.y;
					const bool option = action.option;
					if (option){
						//inverted
						if (i - block.i1 >= 2 && block.j2-j >= 2) 
							stack.emplace_back(block.i1, j, i, block.j2);
//						else if (i - block.i1 == 1 || block.j2-j == 1){
//							double max_ = DBL_MIN;
//							for(int p=block.i1; p<i; p++)
//								for(int q=j; q<block.j2; q++){
//									const double score_ =softMatrix[p][q]; 
//									if (max_ < score_)
//										max_=score_;
//							}
//							score += log(max_); 
//						} 
						if (block.i2 - i >= 2 && j-block.j1 >= 2)
							stack.emplace_back(i, block.j1, block.i2, j);
//						else if (block.i2 - i == 1 || j-block.j1 == 1){
//							double max_ = DBL_MIN;
//							for(int p=i; p<block.i1; p++)
//								for(int q=block.j1; q<j; q++){
//									const double score_ = softMatrix[p][q]; 
//									if (max_ < score_)
//										max_=score_;
//							}
//							score += log(max_); 
//						} 
					}
					else{
						if (i-block.i1 >= 2 && j-block.j1 >= 2) 
							stack.emplace_back(block.i1, block.j1, i, j);
//						else if (i-block.i1 == 1 || j-block.j1 == 1){
//							double max_ = DBL_MIN;
//							for(int p=block.i1; p<i; p++)
//								for(int q=block.j1; q<j; q++){
//									const double score_ = softMatrix[p][q]; 
//									if (max_ < score_)
//										max_=score_;
//							}
//							score += log(max_); 
//						} 
						if (block.i2- i >= 2 && block.j2-j >= 2) 
							stack.emplace_back(i, j, block.i2, block.j2);
//						else if (block.i2- i == 1 || block.j2-j == 1){
//							double max_ = DBL_MIN;
//							for(int p=i; p<block.i2; p++)
//								for(int q=j; q<block.j2; q++){
//									const double score_ = softMatrix[p][q]; 
//									if (max_ < score_)
//										max_=score_;
//							}
//							score += log(max_); 
//						} 
					}
					if (stack.empty())
						terminal = true;
					actions.push_back(action);
				} 
				bool operator<(const ParserState &rhs) const {
					return (score > rhs.score);
				};
		};
		  
		typedef priority_queue<ParserState> Agenda; 
		typedef priority_queue<ParserAction> BestParserActions; 
		
		TranslationTable ttable; 
		
		
		inline void AddTranslationOptionsInBatch(vector<vector<pair<unsigned, double>>>& insert_buffer);
		
		void Parse(const Matrix &accumMatrix, const int ls, const int lt, 
								vector<ParserAction> *bestActions);
								
		void SearchBestPartition(const Matrix &accumMatrix, const int i1,  const int j1, 
												const int i2,  const int j2, BestParserActions *bestParserActions);
		void ActionToAlignment(const vector<ParserAction> &bestActions, const int ls, const int lt,  Alignment *links);									
		void BuildAccumMatrix(const Matrix &softMatrix, Matrix *accumMatrix);
		void Partitionize(const Matrix &accumMatrix, const int i1, const int j1, const int i2, const int j2, Alignment *links);
		void BuildSoftAlignmentMatrix(const SentencePair &sentencePair, Matrix *softMatrix);
		static void FmeasureXY(const Matrix &accumMatrix, const int i1, const int j1, const int i2, const int j2, const int i3, const int j3,
							double *score, double *score_);
		static void Ncut(const Matrix &accumMatrix, const int i1, const int j1, const int i2, const int j2, const int i3, const int j3,
							double *score, double *score_);  
		
		double thread_buffer_size=10000;
	public:
		Hieralign(const string& beamsize); 
		~Hieralign(); 
		
		double sigma_theta;
		double sigma_delta; 
		double sigma_threshold;
		void setHyperParameters(const double theta = 1, const double delta = 0, const double threshold = 0.0001); 
		void Align(ParallelCorpus &parallelCorpus); 
		void LoadTranslationTable(ifstream *file, const SimpleWordWrapper &sw2id, const SimpleWordWrapper &tw2id);
		void LoadTranslationTable(const TranslationTable &other);
		void LoadTranslationTableWithBatches(ifstream *file, const SimpleWordWrapper &sw2id, const SimpleWordWrapper &tw2id);
		void PrintAlignmentList(const AlignmentList &linksList);
};

#endif // CUTNALIGNER_H
