#ifndef HIERAGLIGN_UTILS_H
#define HIERAGLIGN_UTILS_H
#include <omp.h>
#include <vector>  
#include <string> 
#include <set>  
#include <unordered_map> 
#include <iostream> 
#include <fstream>  
#include <cmath>
#include <cassert>  
#include <map> 
#include <fstream> 
#include <algorithm>  
#include <limits> 
#include <ctime>
#include <cmath> 
#include <cfloat> 
#include <queue>
#include <stack> 
using namespace std;

#define CHECK(condition, message) do { \
	if (!(condition)) { \
		cerr << "ERROR(" << __FILE__ << ":" << __LINE__ << ") " << (message) << endl; \
		abort(); \
	} \
} while (0);
 

inline vector<string> Split(const string &str, const string &delim) {
	vector<string> result;
	string::size_type p, q;
	for (p = 0; (q = str.find(delim, p)) != string::npos;
				p = q + delim.size()) {
		result.emplace_back(str, p, q - p);
	}
	result.emplace_back(str, p);
	return result;
} 

typedef unordered_map<unsigned, double> W2Double;
typedef vector<W2Double> W2WDouble; 
#endif // HIERAGLIGN_UTILS_H
