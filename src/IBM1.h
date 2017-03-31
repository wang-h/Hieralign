#ifndef IBM1_H
#define IBM1_H
#include "Sentence.h"
#include "ParallelCorpus.h"
#include "Alignment.h"
#include "TranslationTable.h"

class IBM1
{
  protected:
    unsigned iterations_;
    bool forward_;
    bool with_null_;

    size_t thread_buffer_size_ = 5000;

    void BatchInitializeTranslationTableEntries(ParallelCorpus &parallelCorpus);
    inline void AddTranslationOptionsInBatch(vector<vector<unsigned>> &insert_buffer);

    void TrainWithBatches(ParallelCorpus &parallelCorpus);
    void UpdateParamsUsingBatch(const vector<SentencePair> &buffer);

  public:
    IBM1();
    ~IBM1();

    TranslationTable ttable;
    AlignmentList linksList;
    void Config(const bool &forward, const unsigned &iterations, const bool &withNull, const size_t &size_source_vocab);
    void TrainModel(ParallelCorpus &parallelCorpus);

    void ViterbiAlign(ParallelCorpus &corpus);
    inline void ViterbiAlignSentencePair(const Sentence &source, const Sentence &target, Alignment *links);

    void EstimateViterbiProb(ParallelCorpus &corpus, const AlignmentList &link_list, const bool &symmetrized = false);
};
#endif // IBM1_H
