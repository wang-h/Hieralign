#ifndef PARALLELCORPUS_H
#define PARALLELCORPUS_H
#include "Utils.h"
#include "Sentence.h"
#include "SimpleWordWrapper.h"

class ParallelCorpus
{
  private:
    SentenceList sourceSentences;
    SentenceList targetSentences;

    inline const Sentence EncodeSentenceWithWrapper(const string &line, SimpleWordWrapper *w2id)
    {
        Sentence sentence;
        const vector<string> tokens = Split(line, " ");
        CHECK(tokens.size() > 0, "Invalid line: " + line);
        for (size_t j = 0; j < tokens.size(); j++)
            sentence.push_back(w2id->encode(tokens[j]));
        return sentence;
    };

  public:
    ParallelCorpus(){};
    ~ParallelCorpus(){};
    size_t size_ = 0;
    void ReadParallelCorpus(ifstream *file, SimpleWordWrapper *sw2id, SimpleWordWrapper *tw2id)
    {
        for (string line; !getline(*file, line).eof();)
        {
            const vector<string> sentence_pair = Split(line, "\t");
            CHECK(sentence_pair.size() == 2, "#ERROR! invalid line: " + line);
            sourceSentences.push_back(EncodeSentenceWithWrapper(sentence_pair[0], sw2id));
            targetSentences.push_back(EncodeSentenceWithWrapper(sentence_pair[1], tw2id));
        }
        CHECK(source_size() == target_size(), "#ERROR! sentences are not equal.");
        size_ = source_size();
        //cerr << "training corpus size:"  << " [" << size_ << "]" << endl;
    };
    void ReadCorpus(ifstream *file, SimpleWordWrapper *sw2id, SimpleWordWrapper *tw2id, bool isSource = true)
    {
        SentenceList *sentences = (isSource) ? &sourceSentences : &targetSentences;
        SimpleWordWrapper *w2id = (isSource) ? sw2id : tw2id;
        for (string line; !getline(*file, line).eof();)
        {
            sentences->push_back(EncodeSentenceWithWrapper(line, w2id));
        }
        if (!isSource)
        {
            CHECK(source_size() == target_size(), "#ERROR! sentences are not equal.");
            size_ = source_size();
            //cerr << "training corpus size:"  << " [" << size_ << "]" << "source vocab size:" << " [" << source_vocab_size() << "]" << "target vocab size:"  << " [" << target_vocab_size() << "]" << endl;
        }
    }
    size_t size() const
    {
        return size_;
    }
    size_t source_size() const
    {
        return sourceSentences.size();
    }
    size_t target_size() const
    {
        return targetSentences.size();
    }
    SentencePair operator[](int i)
    {
        return make_pair(&sourceSentences[i], &targetSentences[i]);
    }
};

#endif // PARALLELCORPUS_H
