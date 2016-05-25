#ifndef NO_ALIGNMENT_FLAG
#define NO_ALIGNMENT_FLAG
#define NO_ALIGNMENT (-1)
#endif

#ifndef FEATURE_ALIGNMENT_HPP_FLAG
#define FEATURE_ALIGNMENT_HPP_FLAG

#include <vector>
#include "feature.hpp"

struct FeatureAlignment{
	static bool rtCompare(const FeatureAlignment &feature1, const FeatureAlignment &feature2);
	static bool mzCompare(const FeatureAlignment &feature1, const FeatureAlignment &feature2);
	static bool intensityCompare(const FeatureAlignment &feature1, const FeatureAlignment &feature2);
	
	RT_TYPE rt;
	MZ_TYPE mz;
	INTENSITY_TYPE intensity;
	int charge;
	double quality;
	int alignedWith, mz_order, rt_order, intensity_order;
	
	FeatureAlignment(const Feature &feature); 
	Feature toFeature();
};

#endif
