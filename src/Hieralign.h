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
    double sigma_filter = 0.0;
    unsigned beam_width = 10;
    double sigma_p0 = 1e-6;

    struct ParserAction
    {
        unsigned x;
        unsigned y;
        bool option;
        double score;
        ParserAction(unsigned x_, unsigned y_, bool option_, double score_ = 0.0) : x(x_), y(y_), option(option_), score(score_){};
        bool operator==(const ParserAction &rhs) const
        {
            return (x == rhs.x && y == rhs.y && option == rhs.option);
        }
        bool operator<(const ParserAction &rhs) const
        {
            return (score > rhs.score);
        }
    };

    struct TranslationEntry
    {
        unsigned x;
        unsigned y;
        double score;
        TranslationEntry(unsigned x_, unsigned y_, double score_) : x(x_), y(y_), score(score_){};
    };

    struct ParserBlock
    {
        unsigned i1;
        unsigned j1;
        unsigned i2;
        unsigned j2;
        //type {0, 1}
        //0: terminal, 1: nonterminal
        ParserBlock(unsigned i1_, unsigned j1_, unsigned i2_, unsigned j2_) : i1(i1_), j1(j1_), i2(i2_), j2(j2_) {}
        bool operator==(const ParserBlock &rhs) const
        {
            return (rhs.i1 == i1 && rhs.j1 == j1 && rhs.i2 == i2 && rhs.j2 == j2);
        }
    };

    struct ParserState
    {
      public:
        double score;
        bool terminal;
        vector<ParserBlock> stack;    // Stack of open blocks.
        vector<ParserAction> actions; // History of actions.
        explicit ParserState(unsigned m, unsigned n) : score(0.0), terminal(false)
        {
            if (m > 1 && n > 1)
            {
                stack.emplace_back(0, 0, m, n);
            }
        }
        ParserState(const ParserState &state)
            : score(state.score), terminal(state.terminal), stack(state.stack), actions(state.actions)
        {
        }
        void Advance(const ParserAction &action)
        {
            //score *=  action.score;
            //score *= action.score;
            score += log(action.score);
            const ParserBlock block = stack.back();
            stack.pop_back();

            const unsigned i = action.x;
            const unsigned j = action.y;
            const bool option = action.option;
            if (option)
            {
                //inverted option=true
                if (i - block.i1 >= 2 && block.j2 - j >= 2)
                    stack.emplace_back(block.i1, j, i, block.j2);
                if (block.i2 - i >= 2 && j - block.j1 >= 2)
                    stack.emplace_back(i, block.j1, block.i2, j);
            }
            else
            {
                if (i - block.i1 >= 2 && j - block.j1 >= 2)
                    stack.emplace_back(block.i1, block.j1, i, j);
                if (block.i2 - i >= 2 && block.j2 - j >= 2)
                    stack.emplace_back(i, j, block.i2, block.j2);
            }
            if (stack.empty())
                terminal = true;
            actions.push_back(action);
        }
        bool operator<(const ParserState &rhs) const
        {
            return (score > rhs.score);
        };
    };

    typedef priority_queue<ParserState> Agenda;
    typedef priority_queue<ParserAction> BestParserActions;

    TranslationTable ttable;

    void BuildSoftAlignmentMatrix(const SentencePair &sentencePair,
                                  const unsigned &m, const unsigned &n, Matrix *softMatrix);

    void Parse(const Matrix &accumMatrix, const unsigned &m, const unsigned &n,
               vector<ParserAction> *bestActions);

    void SearchBestPartition(const Matrix &accumMatrix, BestParserActions *bestParserActions,
                             const unsigned &i1, const unsigned &j1, const unsigned &i2, const unsigned &j2);

    void ActionToAlignment(const vector<ParserAction> &bestActions, const unsigned &m,
                           const unsigned &n, Alignment *links, const Matrix &softMatrix);

    void BuildAccumMatrix(const Matrix &softMatrix, Matrix *accumMatrix);

    void Partitionize(const Matrix &accumMatrix, const unsigned i1, const unsigned j1,
                      const unsigned i2, const unsigned j2, Alignment *links);

    static void FmeasureXY(const Matrix &accumMatrix, const unsigned &i1, const unsigned &j1,
                           const unsigned &i2, const unsigned &j2, const unsigned &i3, const unsigned &j3,
                           double *score, double *score_);

    static void Ncut(const Matrix &accumMatrix, const unsigned &i1, const unsigned &j1,
                     const unsigned &i2, const unsigned &j2, const unsigned &i3, const unsigned &j3,
                     double *score, double *score_);

    double thread_buffer_size = 10000;

  public:
    Hieralign(const string &beamsize);
    ~Hieralign();

    double sigma_theta;
    double sigma_delta;
    double sigma_threshold;
    void PrintAlignmentList(const AlignmentList &linksList);
    void setFilter(const double &filter = 0.0);
    void setHyperParameters(const double &theta = 3, const double &delta = 5, const double &p0 = 1e-6);
    void Align(ParallelCorpus &parallelCorpus);
    void LoadTranslationTable(const TranslationTable &other);
    void unsignedAlignmentList(const AlignmentList &linksList);
};

#endif // HIERALIGN_H
