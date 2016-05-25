#ifndef FEATURE_SCORE_TYPE_FLAG
#define FEATURE_SCORE_TYPE_FLAG
#define FEATURE_SCORE_TYPE double
#endif

#ifndef FEATURE_SCORE_FUNCTION_FLAG
#define FEATURE_SCORE_FUNCTION_FLAG
#include "feature.hpp"
extern FEATURE_SCORE_TYPE FEATURE_SCORE_FUNCTION(Feature &a);
#endif

#ifndef SAMPLE_HPP_FLAG
#define SAMPLE_HPP_FLAG

#include <string>
#include <vector>
#include <algorithm>
#include "feature.hpp"
#include "alignment.hpp"

using namespace std;

class Sample{
public:
	string sampleName, runName;
	vector<Feature> features;
	
	Sample(){}
	Sample(FILE* csvFile, const string SampleName = "", const string RunName = "");
	Sample(const string csvFileName, const string SampleName = "", const string RunName = "");
	Sample(Sample &sample, int numRunsWithTopIntensity);
	Sample(Sample &sample);
	void write(FILE* csvFile);
	void write(const string csvFileName);
	void writeAlignment(Sample &sample, FILE* csvFile);
	void writeAlignment(Sample &sample, const string csvFileName);
	void sortByRt();
	void sortByMz();
	void sortByIntensity();
	void normalizeRt();
	ALIGNMENT_SCORE_TYPE alignmentScore(Sample &sample);
	int maxAlignmentSizePossible(Sample &sample);
	Sample unionAlign(Sample &sample);
	Sample intersectAlign(Sample &sample);
	vector<pair<Feature, Feature> > topAlignedFeatures(Sample &sample);
	vector<pair<Feature, Feature> > topAlignedFeatures(Sample &sample, int maxNumber);
	
private:
	bool isSortedByRt, isSortedByMz, isSortedByIntensity;
};

#endif

