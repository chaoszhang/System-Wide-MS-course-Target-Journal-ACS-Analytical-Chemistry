#include <iostream>
#include <fstream>
#include <map>
#include "sample.hpp"

using namespace std;

ALIGNMENT_SCORE_TYPE SCORING_FUNCTION(FeatureAlignment &a, FeatureAlignment &b){
	if (a.charge != b.charge) return -1e9;
	MZ_TYPE mzDifference = (a.mz > b.mz) ? a.mz - b.mz : b.mz - a.mz;
	return pow(a.intensity * b.intensity, 0.5 * a.quality * b.quality) * (0.051 / (0.001 + mzDifference) - 1);
}

FEATURE_SCORE_TYPE FEATURE_SCORE_FUNCTION(Feature &a){
	return pow(log(a.intensity * a.quality), 3);
}

int main(int argc, char** argv) {
	//spec-00003.mzXML.csv
	//Sample s1("spec-00003.mzXML.csv"), s2("spec-00004.mzXML.csv");
	ifstream fin("list.txt");
	ofstream fout("Matrix.txt");
	string s;
	vector<Sample*> samples;
	int c = 0;
	while(fin >> s) {
		cerr << "input: " << samples.size() << endl;
		samples.push_back(new Sample(s + ".csv"));
		cerr << s + ".csv" << endl;
		cerr << samples.back()->features.size() << endl;
	}
	for (Sample *s1: samples){
		for (Sample *s2: samples){
			c++;
			cerr << "align: " << c << "/" << samples.size() * samples.size() << endl;
			fout << s1->alignmentScore(*s2) << "\t";
		}
		fout << endl;
	}
	for (Sample *s1: samples) delete s1;
	//Sample s3(s1, 30);
	//Sample s4(s2, 30);
	//s3.write("input1.csv");
	//s4.write("input2.csv");
	//cerr << "sizes = " << s1.features.size() << ", " << s2.features.size() << " score = " << s1.alignmentScore(s2) << endl;
	//cerr << "sizes = " << s3.features.size() << ", " << s4.features.size() << " score = " << s3.alignmentScore(s4) << endl;
	//cerr << "sizes = " << s1.features.size() << ", " << s2.features.size() << " score = " << s1.maxAlignmentSizePossible(s2) << endl;
	//cerr << "sizes = " << s3.features.size() << ", " << s4.features.size() << " score = " << s3.maxAlignmentSizePossible(s4) << endl;
	//s1.writeAlignment(s2, "test_result1.csv");
	/*
	cerr << s1.alignmentScore(s2) << endl;
	cerr << s3.alignmentScore(s4) << endl;
	s3.sortByRt();
	for (Feature &f: s3.features)
		cerr << f.mz << "\t";
	cerr << endl;
	s4.sortByRt();
	for (Feature &f: s4.features)
		cerr << f.mz << "\t";
	cerr << endl << endl;
	*/
	return 0;
}
