#include "feature_alignment.hpp"

bool FeatureAlignment :: rtCompare(const FeatureAlignment &feature1, const FeatureAlignment &feature2){
	if (feature1.rt == feature2.rt) return feature1.mz < feature2.mz;
	return feature1.rt < feature2.rt;
}

bool FeatureAlignment :: mzCompare(const FeatureAlignment &feature1, const FeatureAlignment &feature2){
	return feature1.mz < feature2.mz;
}

bool FeatureAlignment :: intensityCompare(const FeatureAlignment &feature1, const FeatureAlignment &feature2){
	return feature1.intensity > feature2.intensity;
}

FeatureAlignment :: FeatureAlignment(const Feature &f): rt(f.rt), mz(f.mz), intensity(f.intensity), charge(f.charge),
	quality(f.quality), alignedWith(NO_ALIGNMENT){} 
	
Feature FeatureAlignment :: toFeature(){
	Feature f;
	f.rt = rt;
	f.mz = mz;
	f.intensity = intensity;
	f.charge = charge;
	f.quality = quality;
	return f;
}
