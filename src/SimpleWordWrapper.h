#ifndef SIMPLEWORDWRAPPER_H
#define SIMPLEWORDWRAPPER_H
#include "Utils.h"
class SimpleWordWrapper
{
  public:
    SimpleWordWrapper() : unk_("<UNK>")
    {
        words.reserve(20000);
        encode("NULL"); // dict["NULL "] = 1;
    };
    ~SimpleWordWrapper(){};
    typedef unordered_map<string, unsigned, hash<string>> HASH_MAP_TYPE;
    inline unsigned max() const { return words.size(); }
    inline unsigned encode(const string &word)
    {
        HASH_MAP_TYPE::iterator i = dict.find(word);
        if (i == dict.end())
        {
            if (frozen)
                return 0; //for unknown word
            words.push_back(word);
            dict[word] = words.size();
            return words.size();
        }
        else
            return i->second;
    }
    inline const string &decode(const unsigned id) const
    {
        if (id == 0)
            return unk_;
        return words[id - 1];
    }
    size_t size() const
    {
        return words.size();
    }
    void Frozen()
    {
        frozen = true;
    }

  private:
    string unk_;
    bool add_null_;
    bool frozen = false;
    vector<string> words;
    HASH_MAP_TYPE dict;
    HASH_MAP_TYPE dict_;
};
#endif // SIMPLEWORDWRAPPER_H