#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <string>
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

map<string, string> unionfind;
map<string, int> counter;
string Find(string s){
	if (unionfind[s] == "") return s;
	return unionfind[s] = Find(unionfind[s]);
}
void Union(string s1, string s2){
	if (Find(s1) == Find(s2)) return;
	if (Find(s1) == "Control") {
		counter[Find(s1)] += counter[Find(s2)] + 1;
		unionfind[Find(s2)] = Find(s1);
	}
	else {
		counter[Find(s2)] += counter[Find(s1)] + 1;
		unionfind[Find(s1)] = Find(s2);
	}
}

int main(int argc, char** argv) {
	//spec-00003.mzXML.csv
	//Sample s1("spec-00003.mzXML.csv"), s2("spec-00004.mzXML.csv");
	ifstream fMatrix("MatrixWithLabel.tsv");
	ofstream fSif("CommonFeatureForFemaleOnly.sif");
	string name, path, headLine;
	vector<Sample*> samples;
	vector<vector<double> > matrix;
	map<string, set<string> > has;
	//map<string, set<string> > alignedWith;
	set<string> output; 
	getline(fMatrix, headLine);
	stringstream fList(headLine);
	while(fList >> path) samples.push_back(new Sample(path + ".csv", path, path));
	cerr << samples.size();
	for (Sample *s1: samples){
		fMatrix >> s1->sampleName;
		matrix.push_back({});
		for (Sample *s2: samples){
			double score;
			fMatrix >> score;
			matrix.back().push_back(score);
		}
	}
	for (int i = 0; i < samples.size(); i++){
		for (int j = i + samples.size() / 3; j < samples.size(); j++){
			Sample *s1 = samples[i], *s2 = samples[j];
			if (s1->sampleName.find("Person") == string::npos || s2->sampleName.find("Person") == string::npos) continue;
			if (matrix[i][j] > 0.15){
				int num;
				if ((s1->sampleName.find("_M_") != string::npos || s1->sampleName.find("Male") != string::npos)
					|| (s2->sampleName.find("_M_") != string::npos || s2->sampleName.find("Male") != string::npos)) num = 500;
				else num = 100;
				cerr << s1->sampleName << " " << s2->sampleName << endl; 
				auto res = s1->topAlignedFeatures(*s2, num);
				for (auto &e: res){
					string feature1 = s1->sampleName + "_RT=" + to_string(e.first.rt) + "_MZ=" + to_string(e.first.mz);
					string feature2 = s2->sampleName + "_RT=" + to_string(e.second.rt) + "_MZ=" + to_string(e.second.mz);
					has[s1->sampleName].insert(feature1);
					has[s2->sampleName].insert(feature2);
					//alignedWith[feature1].insert(feature2);
					//alignedWith[feature2].insert(feature1);
					Union(feature1, feature2);
					if (s1->sampleName.find("_M_") != string::npos || s1->sampleName.find("Male") != string::npos) Union(feature1, "Control");
					if (s2->sampleName.find("_M_") != string::npos || s2->sampleName.find("Male") != string::npos) Union(feature2, "Control");
				}
			}
		}
	}
	for (Sample *s1: samples) delete s1;
	for (auto &from: has){
		for (auto &to: from.second){
			if (Find(to) != "Control" && counter[Find(to)] > 3 && from.first.find("_M_") == string::npos
				&& from.first.find("Male") == string::npos) output.insert(from.first + " mf " + Find(to));
		}
	}
	for (auto &s: output)
		fSif << s << endl;
	/*
	for (auto &from: alignedWith){
		for (auto &to: from.second){
			if (Find(from.first) != "Control" && Find(to) != "Control") fSif << from.first << " AlignedWith " << to << endl;
		}
	}
	*/
	return 0;
}
