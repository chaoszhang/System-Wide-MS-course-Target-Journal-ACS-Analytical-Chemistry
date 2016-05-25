#ifndef RT_TYPE_FLAG
#define RT_TYPE_FLAG
typedef double RT_TYPE;
#endif

#ifndef MZ_TYPE_FLAG
#define MZ_TYPE_FLAG
typedef double MZ_TYPE;
#endif

#ifndef INTENSITY_TYPE_FLAG
#define INTENSITY_TYPE_FLAG
typedef double INTENSITY_TYPE;
#endif

#ifndef FEATURE_HPP_FLAG
#define FEATURE_HPP_FLAG

struct Feature{
	static bool rtCompare(const Feature &feature1, const Feature &feature2);
	static bool mzCompare(const Feature &feature1, const Feature &feature2);
	static bool intensityCompare(const Feature &feature1, const Feature &feature2);
	
	RT_TYPE rt;
	MZ_TYPE mz;
	INTENSITY_TYPE intensity;
	int charge;
	double quality; 
};

#endif

