#ifndef ALIGNMENT_SCORE_TYPE_FLAG
#define ALIGNMENT_SCORE_TYPE_FLAG
typedef double ALIGNMENT_SCORE_TYPE;
#endif

#ifndef MZ_DIFFERENCE_THRESHOLD_FLAG
#define MZ_DIFFERENCE_THRESHOLD_FLAG
#define MZ_DIFFERENCE_THRESHOLD (0.05)
#endif

#ifndef RT_DIFFERENCE_THRESHOLD_FLAG
#define RT_DIFFERENCE_THRESHOLD_FLAG
#define RT_DIFFERENCE_THRESHOLD (1)
#endif

#ifndef SCORING_FUNCTION_FLAG
#define SCORING_FUNCTION_FLAG
#include "feature_alignment.hpp"
extern ALIGNMENT_SCORE_TYPE SCORING_FUNCTION(FeatureAlignment &a, FeatureAlignment &b);
#endif

#ifndef RECURSION_THRESHOLD_FLAG
#define RECURSION_THRESHOLD_FLAG
#define RECURSION_THRESHOLD (10000)
#endif

#ifndef AFFINITY_BONUS_FLAG
#define AFFINITY_BONUS_FLAG
#define AFFINITY_BONUS (0.1)
#endif

#ifndef ALIGNMENT_HPP_FLAG
#define ALIGNMENT_HPP_FLAG

#include <cstdlib>
#include <string>
#include <vector>
#include "sample.hpp"
#include "feature_alignment.hpp"

using namespace std;

class Sample;

struct State{
	int x, y, z;
	State(): x(0), y(0), z(0){}
	State(int i, int j, int k): x(i), y(j), z(k){}
	bool operator==(const State &s) const{
		return x == s.x && y == s.y && z == s.z;
	}
	bool operator<(const State &s) const{
		return (x < s.x) || (x == s.x && y < s.y) || (x == s.x && y == s.y && z < s.z);
	}
};

class Alignment{
public:
	vector<FeatureAlignment> sample1, sample2;
	string sampleName;

	Alignment(Sample &sample1, Sample &sample2);
	ALIGNMENT_SCORE_TYPE fullAlignmentScore();
	int maxAlignmentSizePossible();
	ALIGNMENT_SCORE_TYPE localAlignment(State left, State right);
	vector<vector<ALIGNMENT_SCORE_TYPE> > forwardAlignment(State left, State right);
	vector<vector<ALIGNMENT_SCORE_TYPE> > reverseAlignment(State left, State right);
	void recursiveAlignment(State left, State right);
	ALIGNMENT_SCORE_TYPE fullAlignment();
	void write(FILE* csvFile);
	void write(const string csvFileName);
	void augmentingAlignment(int initSize1 = 200, int initSize2 = 200, int multiplier = 2, int windowSize = 10, int fullAlignmentRound = 1);
	ALIGNMENT_SCORE_TYPE alignmentScore();
	int alignmentCount();
	
private:
	bool isSortedByRt = false, isSortedByMz = false, isSortedByIntensity = false;
	
	void sortByRt();
	void sortByMz();
	void sortByIntensity();
};

#endif
