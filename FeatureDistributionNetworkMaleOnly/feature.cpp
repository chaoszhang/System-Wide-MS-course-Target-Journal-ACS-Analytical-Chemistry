#include "feature.hpp"

bool Feature :: rtCompare(const Feature &feature1, const Feature &feature2){
	if (feature1.rt == feature2.rt) return feature1.mz < feature2.mz;
	return feature1.rt < feature2.rt;
}

bool Feature :: mzCompare(const Feature &feature1, const Feature &feature2){
	return feature1.mz < feature2.mz;
}

bool Feature :: intensityCompare(const Feature &feature1, const Feature &feature2){
	return feature1.intensity > feature2.intensity;
}
