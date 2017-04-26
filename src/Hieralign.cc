#include "Hieralign.h"
Hieralign::Hieralign(const string &beamsize)
{
    beam_width = stoi(beamsize);
}

Hieralign::~Hieralign()
{
}
void Hieralign::setFilter(const double &filter)
{
    sigma_filter = filter;
}
void Hieralign::setHyperParameters(const double &theta, const double &delta, const double &p0)
{
    sigma_theta = theta;
    sigma_delta = delta;
    sigma_p0 = p0;
    cerr << "Parameters:     \t[theta: " << sigma_theta << ", delta: " << sigma_delta << ", p0: " << sigma_p0 << "]" << endl;
}
void Hieralign::LoadTranslationTable(const TranslationTable &other)
{
    ttable = other;
}

void Hieralign::Align(ParallelCorpus &parallelCorpus)
{
    AlignmentList link_list(parallelCorpus.size());

    int count = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : count)
    for (size_t i = 0; i < parallelCorpus.size(); i++)
    {
        if (i % 10000 == 0 && i != 0)
        {
            cerr << '.';
        }
        if (i % 100000 == 0 && i != 0)
        {
            cerr << " \t[" << i << "]\n"
                 << flush;
        }
        Alignment links;
        SentencePair sentence_pair = parallelCorpus[i];
        const unsigned m = static_cast<int>(sentence_pair.first->size());
        const unsigned n = static_cast<int>(sentence_pair.second->size());
        if (m == 1 || n == 1)
        {
            for (unsigned i = 0; i < m; i++)
                for (unsigned j = 0; j < n; j++)
                    links.emplace_back(i, j);
        }
        else
        {
            Matrix softMatrix(m, n);
            BuildSoftAlignmentMatrix(sentence_pair, m, n, &softMatrix);
            //softMatrix.PrintMatrix();
            Matrix accumMatrix(m + 1, n + 1);
            BuildAccumMatrix(softMatrix, &accumMatrix);
            vector<ParserAction> bestActions;
            Parse(accumMatrix, m, n, &bestActions);
            ActionToAlignment(bestActions, m, n, &links, softMatrix);
        }
        link_list[i] = links;
        count += 1;
    }
    cerr << "Total:          \t[" << count << "]" << flush;
    cerr << "\n# outputing result ..." << endl;
    PrintAlignmentList(link_list);
}

void Hieralign::ActionToAlignment(const vector<ParserAction> &bestActions,
                                  const unsigned &m, const unsigned &n, Alignment *links,
                                  const Matrix &softMatrix)
{

    vector<ParserBlock> stack;
    stack.emplace_back(0, 0, m, n);
    for (size_t k = 0; k < bestActions.size(); k++)
    {
        const ParserAction &action = bestActions[k];
        const ParserBlock block = stack.back();
        stack.pop_back();
        if (action.option)
        {
            //inverted
            if (action.x - block.i1 >= 2 && block.j2 - action.y >= 2)
                stack.emplace_back(block.i1, action.y, action.x, block.j2);
            else if (action.x - block.i1 == 1 || block.j2 - action.y == 1)
                for (unsigned i = block.i1; i < action.x; i++)
                    for (unsigned j = action.y; j < block.j2; j++)
                    {
                        if (softMatrix(i, j) > sigma_filter)
                            (*links).emplace_back(i, j);
                    }
            if (block.i2 - action.x >= 2 && action.y - block.j1 >= 2)
                stack.emplace_back(action.x, block.j1, block.i2, action.y);
            else if (block.i2 - action.x == 1 || action.y - block.j1 == 1)
                for (unsigned i = action.x; i < block.i2; i++)
                    for (unsigned j = block.j1; j < action.y; j++)
                    {
                        if (softMatrix(i, j) > sigma_filter)
                            (*links).emplace_back(i, j);
                    }
        }
        else
        {
            if (action.x - block.i1 >= 2 && action.y - block.j1 >= 2)
                stack.emplace_back(block.i1, block.j1, action.x, action.y);
            else if (action.x - block.i1 == 1 || action.y - block.j1 == 1)
                for (unsigned i = block.i1; i < action.x; i++)
                    for (unsigned j = block.j1; j < action.y; j++)
                    {
                        if (softMatrix(i, j) > sigma_filter)
                            (*links).emplace_back(i, j);
                    }
            if (block.i2 - action.x >= 2 && block.j2 - action.y >= 2)
                stack.emplace_back(action.x, action.y, block.i2, block.j2);
            else if (block.i2 - action.x == 1 || block.j2 - action.y == 1)
                for (unsigned i = action.x; i < block.i2; i++)
                    for (unsigned j = action.y; j < block.j2; j++)
                    {
                        if (softMatrix(i, j) > sigma_filter)
                            (*links).emplace_back(i, j);
                    }
        }
    }
}

void Hieralign::PrintAlignmentList(const AlignmentList &linksList)
{
    for (size_t k = 0; k < linksList.size(); k++)
    {
        const Alignment &links = linksList[k];
        if (links.empty())
        {
            cout << endl;
        }
        else
        {
            for (size_t i = 0; i < links.size(); i++)
            {
                if (i == links.size() - 1)
                    cout << links[i].first << "-" << links[i].second << endl;
                else
                    cout << links[i].first << "-" << links[i].second << " ";
            }
        }
    }
}

void Hieralign::BuildSoftAlignmentMatrix(const SentencePair &sentence_pair,
                                         const unsigned &m, const unsigned &n, Matrix *softMatrix)
{

    const Sentence *src = sentence_pair.first;
    const Sentence *trg = sentence_pair.second;
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < m; i++)
    {
        for (unsigned j = 0; j < n; j++)
        {
            const unsigned sw = (*src)[i];
            const unsigned tw = (*trg)[j];
            const double theta = log(max(ttable.safe_Prob(sw, tw), sigma_p0)) / sigma_theta;
            double score;
            if (sigma_delta < 0.001)
            {
                score = (1.0 - sigma_p0) * exp(theta);
            }
            else
            {
                const double offset = abs(static_cast<double>(i) / m - static_cast<double>(j) / n);
                const double delta = log(max(1.0 - offset, sigma_p0)) / sigma_delta;
                score = (1.0 - sqrt(sigma_p0)) * exp(theta + delta);
            }
            (*softMatrix)(i, j) = score;
        }
    }
}
void Hieralign::BuildAccumMatrix(const Matrix &softMatrix, Matrix *accumMatrix)
{
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < softMatrix.width(); i++)
    {
        for (size_t j = 0; j < softMatrix.height(); j++)
        {
            (*accumMatrix)(i + 1, j + 1) = (*accumMatrix)(i, j + 1) + (*accumMatrix)(i + 1, j) - (*accumMatrix)(i, j) + softMatrix(i, j);
        }
    }
}

void Hieralign::Parse(const Matrix &accumMatrix, const unsigned &m, const unsigned &n,
                      vector<ParserAction> *bestActions)
{
    //const int maxblocksize = 3;
    Agenda old_agenda;
    Agenda new_agenda;
    //top is the smallest
    // Add the initial parser state to the agenda.
    old_agenda.push(ParserState(m, n));
    double current_score = DBL_MAX;
    for (unsigned transition = 0; transition < min(m, n); transition++)
    {
        if (old_agenda.empty())
            continue;
        while (!old_agenda.empty())
        {
            const ParserState &state = old_agenda.top();
            if (state.terminal)
            {
                if (state.score < current_score)
                {
                    current_score = state.score;
                    *bestActions = state.actions;
                }

                old_agenda.pop();
                continue;
            }
            const ParserBlock &block = state.stack.back();
            BestParserActions bestParserActions;
            SearchBestPartition(accumMatrix, &bestParserActions, block.i1, block.j1, block.i2, block.j2);
            int top = 0;
            while (!bestParserActions.empty() && top++ < 3)
            {
                const ParserAction &action = bestParserActions.top();
                ParserState new_state = ParserState(state);
                new_state.Advance(action);
                const double threshold = new_agenda.empty() ? DBL_MAX : new_agenda.top().score;
                if (new_agenda.size() <= static_cast<size_t>(beam_width))
                {
                    new_agenda.push(new_state);
                }
                else
                {
                    if (new_state.score <= threshold)
                    {
                        new_agenda.pop();
                        new_agenda.push(new_state);
                    }
                }
                bestParserActions.pop();
            }
            old_agenda.pop();
        }
        if (!new_agenda.empty())
        {
            old_agenda.swap(new_agenda);
        }
        else
        {
            break;
        }
    }
    // Incremental top-down parsing with beam search.
}
void Hieralign::SearchBestPartition(const Matrix &accumMatrix, BestParserActions *bestParserActions,
                                    const unsigned &i1, const unsigned &j1, const unsigned &i2, const unsigned &j2)
{

    for (unsigned i = i1 + 1; i < i2; i++)
    {
        for (unsigned j = j1 + 1; j < j2; j++)
        {
            double score, score_;
            //            const int x = i-i1;
            //            const int y = i2-i;
            //            const int p = j-j1;
            //            const int q = j2-j;
            //if(y-x < 3 )
            //FmeasureXY(accumMatrix, i1, j1, i, j, i2, j2, &score, &score_);
            Ncut(accumMatrix, i1, j1, i, j, i2, j2, &score, &score_);
            //if (x <= 3 && abs((i2-i)-(j2-j))<=3)
            bestParserActions->emplace(i, j, false, score);
            //else if(|| (i1+5< i<i2-5 && j1+5< j<j2-5)))
            //if (abs((i-i1)-(j2-j))<=3 && abs((i2-i)-(j-j1))<=3 || (i1+5< i<i2-5 && j1+5< j<j2-5))
            bestParserActions->emplace(i, j, true, score_);
        }
    }
}
void Hieralign::FmeasureXY(const Matrix &accumMatrix, const unsigned &i1, const unsigned &j1,
                           const unsigned &i2, const unsigned &j2, const unsigned &i3, const unsigned &j3,
                           double *score, double *score_)
{
    const double Pi1j1 = accumMatrix(i1, j1);
    const double Pi1j2 = accumMatrix(i1, j2);
    const double Pi1j3 = accumMatrix(i1, j3);
    const double Pi2j1 = accumMatrix(i2, j1);
    const double Pi2j2 = accumMatrix(i2, j2);
    const double Pi2j3 = accumMatrix(i2, j3);
    const double Pi3j1 = accumMatrix(i3, j1);
    const double Pi3j2 = accumMatrix(i3, j2);
    const double Pi3j3 = accumMatrix(i3, j3);
    const double WAB = Pi2j2 + Pi1j1 - Pi2j1 - Pi1j2;
    const double W_A_B = Pi3j3 + Pi2j2 - Pi2j3 - Pi3j2;
    const double W_AB = Pi2j3 + Pi1j2 - Pi1j3 - Pi2j2;
    const double WA_B = Pi3j2 + Pi2j1 - Pi3j1 - Pi2j2;
    const double W_ABA_B = WA_B + W_AB;
    const double W_A_BAB = WAB + W_A_B;

    *score = 2 - (WAB / (W_ABA_B + WAB + WAB) + W_A_B / (W_ABA_B + W_A_B + W_A_B));
    *score_ = 2 - (W_AB / (W_A_BAB + W_AB + W_AB) + WA_B / (W_A_BAB + WA_B + WA_B));
    //    *score  = 2- WAB /(WAB + W_ABA_B) + W_A_B/(W_A_B + W_ABA_B);
    //    *score_ = 2- W_AB/(W_AB + W_A_BAB) + WA_B/(WA_B +  W_A_BAB);
}
void Hieralign::Ncut(const Matrix &accumMatrix, const unsigned &i1, const unsigned &j1,
                     const unsigned &i2, const unsigned &j2, const unsigned &i3, const unsigned &j3,
                     double *score, double *score_)
{
    const double Pi1j1 = accumMatrix(i1, j1);
    const double Pi1j2 = accumMatrix(i1, j2);
    const double Pi1j3 = accumMatrix(i1, j3);
    const double Pi2j1 = accumMatrix(i2, j1);
    const double Pi2j2 = accumMatrix(i2, j2);
    const double Pi2j3 = accumMatrix(i2, j3);
    const double Pi3j1 = accumMatrix(i3, j1);
    const double Pi3j2 = accumMatrix(i3, j2);
    const double Pi3j3 = accumMatrix(i3, j3);
    const double WXY = Pi2j2 + Pi1j1 - Pi2j1 - Pi1j2;
    const double W_X_Y = Pi3j3 + Pi2j2 - Pi2j3 - Pi3j2;
    const double W_XY = Pi2j3 + Pi1j2 - Pi1j3 - Pi2j2;
    const double WX_Y = Pi3j2 + Pi2j1 - Pi3j1 - Pi2j2;
    //case 1: A (X, Y) B (_X, _Y)
    const double cutAB = W_XY + WX_Y;
    const double assoAV = W_XY + WX_Y + WXY + WXY;
    const double assoBV = W_XY + WX_Y + W_X_Y + W_X_Y;
    *score = cutAB / assoAV + cutAB / assoBV;

    //case 2: A (X, _Y) B (_X, Y)
    const double cutAB_ = WXY + W_X_Y;
    const double assoAV_ = WXY + W_X_Y + W_XY + W_XY;
    const double assoBV_ = WXY + W_X_Y + WX_Y + WX_Y;
    *score_ = cutAB_ / assoAV_ + cutAB_ / assoBV_;
}
