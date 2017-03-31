<<<<<<< HEAD
//　Modified from the original file of Chris Dyer
=======
//　Copyright 2013 by Hao Wang, modified from the original code of Chris Dyer
>>>>>>> 696df8a42e20db9cf0fd617211ef71d9aace1264
// provided in https://github.com/clab/fast_align/blob/master/src/ttables.h

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

#ifndef TRANSLATION_TABLE_H
#define TRANSLATION_TABLE_H

typedef vector<W2Double>::const_iterator const_iterator;

static double digamma(double x)
{
    double result = 0, xx, xx2, xx4;
    for (; x < 7; ++x)
        result -= 1 / x;
    x -= 1.0 / 2.0;
    xx = 1.0 / x;
    xx2 = xx * xx;
    xx4 = xx2 * xx2;
    result += log(x) + (1. / 24.) * xx2 - (7.0 / 960.0) * xx4 + (31.0 / 8064.0) * xx4 * xx2 - (127.0 / 30720.0) * xx4 * xx4;
    return result;
}

class TranslationTable
{
  public:
    TranslationTable() : initialized_(false), frozen_(false){};
    TranslationTable(const TranslationTable &other) : initialized_(other.initialized_), frozen_(other.frozen_)
    {
        tt = other.tt;
        counts = other.counts;
    }
    ~TranslationTable(){};
    inline double Prob(const unsigned &f, const unsigned &e) const
    {
        return initialized_ ? tt[f].find(e)->second : 1e-9;
    }

    inline double safe_Prob(const unsigned &f, const unsigned &e) const
    {
        if (f < static_cast<unsigned>(tt.size()))
        {
            const W2Double &cpd = tt[f];
            const W2Double::const_iterator it = cpd.find(e);
            if (it == cpd.end())
                return 1e-7;
            return it->second;
        }
        else
            return 1e-7;
    }

    inline void Insert(const unsigned f, const unsigned e)
    {
        // NOT thread safe
        if (f >= counts.size())
            counts.resize(f + 1);
        counts[f][e] = 0;
    }
    inline void SetMaxF(const unsigned f)
    {
        // NOT thread safe
        if (f >= counts.size())
        {
            counts.resize(f + 1);
        }
    }
    inline void safe_SetValue(const unsigned f, const unsigned e, const double x)
    {
        if (f >= tt.size())
            tt.resize(f + 1);
        tt[f][e] = x;
    }
    inline void Increment(const unsigned f, const unsigned e, const double x)
    {
        // NOT thread safe
        counts[f].find(e)->second += x;
    }
    void Normalize()
    {
        CHECK(!frozen_, "#ERROR! tt has been initialized.");
        tt.swap(counts);
#pragma omp parallel for schedule(dynamic)
        for (unsigned f = 0; f < tt.size(); f++)
        {
            double total = 0;
            W2Double &cpd = tt[f];
            for (W2Double::iterator it = cpd.begin(); it != cpd.end(); ++it)
                total += it->second;
            if (!total)
                total = 1;
            for (W2Double::iterator it = cpd.begin(); it != cpd.end(); ++it)
                it->second /= total;
        }
        ClearCounts();
    }
    void NormalizeVB(const double alpha)
    {
        CHECK(!frozen_, "#ERROR! tt has been initialized.");
        tt.swap(counts);
#pragma omp parallel for schedule(dynamic)
        for (unsigned f = 0; f < tt.size(); f++)
        {
            double total = 0;
            W2Double &cpd = tt[f];
            for (W2Double::iterator it = cpd.begin(); it != cpd.end(); ++it)
                total += it->second + alpha;
            if (!total)
                total = 1;
            const double digamma_total = digamma(total);
            for (W2Double::iterator it = cpd.begin(); it != cpd.end(); ++it)
                it->second = exp(digamma(it->second + alpha) - digamma_total);
        }
        ClearCounts();
    }
    void SetInitialized()
    {
        CHECK(!initialized_, "#ERROR! tt has been initialized.");
        if (!initialized_)
        {
            tt.resize(counts.size());
            for (unsigned i = 0; i < counts.size(); i++)
            {
                tt[i] = counts[i];
            }
        }
        initialized_ = true;
    }
    void Freeze()
    {
        CHECK(!frozen_, "#ERROR! tt has frozen.");
        counts.clear();
        frozen_ = true;
    }
    const_iterator begin() const { return tt.begin(); }
    const_iterator end() const { return tt.end(); }
    W2Double operator[](int i) const
    {
        return tt[i];
    }
    size_t size() const
    {
        return tt.size();
    }
    void clear()
    {
        counts.clear();
        tt.clear();
        frozen_ = false;
        initialized_ = false;
    }

    void Merge(TranslationTable &other, const bool &reverse = true)
    {
        CHECK(!frozen_, "#ERROR! tt has frozen.");
        counts.clear();
//require this and other have same entries.
#pragma omp parallel for schedule(dynamic)
        for (unsigned f = 0; f < tt.size(); f++)
        {
            // remove NULL from F list
            W2Double &cpd = tt[f];
            for (W2Double::iterator it = cpd.begin(); it != cpd.end(); it++)
            {
                const unsigned e = it->first;
                const double score = (reverse) ? other.safe_Prob(e, f) : other.safe_Prob(f, e);
                it->second = sqrt(it->second * score);
            }
        }
    }

    void WriteTranslationTable(ofstream *file, const SimpleWordWrapper &sw2id, const SimpleWordWrapper &tw2id,
                               const double &threshold, const bool &keepNull)
    {
        CHECK(frozen_, "#ERROR! tt has not frozen.");
        if ((*file).is_open())
        {
            for (unsigned f = (keepNull) ? 0 : 1; f < tt.size(); f++)
            {
                for (auto const &e : tt[f])
                {
                    const double score = e.second;
                    if (score >= threshold)
                        (*file) << sw2id.decode(f) << " " << tw2id.decode(e.first) << " " << score << "\n";
                }
            }
            (*file).close();
        }
    }
    void LoadTranslationTable(ifstream *file, SimpleWordWrapper &sw2id, SimpleWordWrapper &tw2id)
    {
        CHECK(!frozen_, "#ERROR! tt has frozen.");
        int i;
        int effective = 0;
        string line;
        for (i = 0; !getline(*file, line).eof(); i++)
        {
            const vector<string> tokens = Split(line, " ");
            CHECK(tokens.size() == 3, "#ERROR! model format wrong: " + line);
            const string &f_str = tokens[0];
            const string &e_str = tokens[1];
            unsigned f = sw2id.encode(f_str);
            unsigned e = tw2id.encode(e_str);
            if (f && e)
            {
                safe_SetValue(f, e, stod(tokens[2]));
                effective++;
            }
        }
        initialized_ = true;
        (*file).close();
        cerr << "# of entries:          \t[" << i << "]" << endl;
        cerr << "# of effective entries:\t[" << effective << "]" << endl;
    }
    void ClearCounts()
    {
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < counts.size(); i++)
            for (auto &cnt : counts[i])
                cnt.second = 0.0;
    }

  private:
    W2WDouble tt;
    W2WDouble counts;
    double threshold_ = 1e-6;
    bool initialized_; // Disallow new e,f pairs to be added to counts
    bool frozen_;      // Disallow new e,f pairs to be added to counts
};

#endif // TRANSLATION_TABLE_H
