 #ifndef SIMPLEWORDWRAPPER_H
#define SIMPLEWORDWRAPPER_H 
#include "Utils.h"
class SimpleWordWrapper{
	public: 
		SimpleWordWrapper() : b0("<bad0>") {
			words.reserve(20000);
		};
		~SimpleWordWrapper(){};
		typedef unordered_map<string, unsigned, hash<string> > HASH_MAP_TYPE; 
		inline unsigned max() const { return words.size(); } 
		inline unsigned encode(const string& word, bool frozen = false) { 
			HASH_MAP_TYPE::iterator i = dict.find(word);
			if (i == dict.end()) {
				// dict["NULL "] = 1; 
				if (frozen)
					return 0;
				words.push_back(word);
				dict[word] = words.size();
				return words.size();
			}
			else 
				return i->second; 
		} 
		inline const string& decode(const unsigned id) const {
			if (id == 0) return b0;
			return words[id-1];
		}
		size_t size()const{
			return words.size();
		} 
		inline unsigned contain(const string& word) const{ 
			HASH_MAP_TYPE::const_iterator i = dict.find(word);
			if (i == dict.end()) 
				return 0;  
			return i->second; 
		} 
	private:
		string b0;
		vector<string> words;
		HASH_MAP_TYPE dict;
		HASH_MAP_TYPE dict_;
};
#endif // SIMPLEWORDWRAPPER_H