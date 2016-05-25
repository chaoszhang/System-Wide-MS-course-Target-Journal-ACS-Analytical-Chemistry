#include <cstring>
#include <map>
#include <iostream>
#include "sample.hpp"

using namespace std;

const int STRING_BUFFER_SIZE = 1024;

Sample :: Sample(FILE* file, const string sName, const string rName):
	sampleName(sName), runName(rName), isSortedByRt(false), isSortedByMz(false), isSortedByIntensity(false){
	char buffer[STRING_BUFFER_SIZE], dummy[STRING_BUFFER_SIZE];
	Feature temp;
	int charge;
	double rt, mz, intensity, quality, width;
	while(!feof(file) && fgets(buffer, STRING_BUFFER_SIZE, file)){
		if (!strcmp(buffer, "\n")) continue;
		if (!strcmp(buffer, "")) break;
		sscanf(buffer, "%s", dummy);
		if (strcmp(dummy, "FEATURE")) continue;
		sscanf(buffer, "%s\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\n", dummy, &rt, &mz, &intensity, &charge, &width, &quality);
		temp.rt = (RT_TYPE) rt;
		temp.mz = (MZ_TYPE) mz;
		temp.intensity = (INTENSITY_TYPE) intensity;
		temp.charge = charge;
		temp.quality = quality;
		features.push_back(temp);
	}
	sortByRt();
	//normalizeRt();
	fclose(file);
}

Sample :: Sample(const string fileName, const string sName, const string rName):
	Sample(fopen(fileName.c_str(), "r"), sName, rName){}

Sample :: Sample(Sample &s, const int num):
	sampleName(s.sampleName), runName(s.runName), isSortedByRt(false), isSortedByMz(true), isSortedByIntensity(false){
	s.sortByMz();
	FEATURE_SCORE_TYPE totalScore = 0, curScore = 0, targetScore;
	int cnt = 0, best = 0;
	for (Feature &f: s.features)
		totalScore += FEATURE_SCORE_FUNCTION(f);
	targetScore = totalScore / num;
	for (int i = 0; i < s.features.size(); i++){
		curScore += FEATURE_SCORE_FUNCTION(s.features[i]);
		if (FEATURE_SCORE_FUNCTION(s.features[i]) > FEATURE_SCORE_FUNCTION(s.features[best])) best = i;
		if (curScore > targetScore){
			features.push_back(s.features[best]);
			best = i + 1;
			cnt++; 
			if (cnt < num) targetScore = curScore + (totalScore - curScore) / (num - cnt);
			else break;
		}
	}
	if (cnt < num && best < s.features.size()) features.push_back(s.features[best]);
}

Sample :: Sample(Sample &s): Sample(s, s.features.size()){}

void Sample :: sortByRt(){
	if (isSortedByRt) return;
	sort(features.begin(), features.end(), Feature :: rtCompare);
	isSortedByRt = true;
	isSortedByMz = false;
	isSortedByIntensity = false;
}

void Sample :: sortByMz(){
	if (isSortedByMz) return;
	sort(features.begin(), features.end(), Feature :: mzCompare);
	isSortedByRt = false;
	isSortedByMz = true;
	isSortedByIntensity = false;
}

void Sample :: sortByIntensity(){
	if (isSortedByIntensity) return;
	sort(features.begin(), features.end(), Feature :: intensityCompare);
	isSortedByRt = false;
	isSortedByMz = false;
	isSortedByIntensity = true;
}

void Sample :: write(FILE* file){
	fprintf(file, "#rt,mz,intensity\n");
	for (Feature feature : features)
		fprintf(file, "%lf,%lf,%lf\n", (double) feature.rt, (double) feature.mz, (double) feature.intensity);
	fclose(file);
}

void Sample :: write(const string fileName){
	write(fopen(fileName.c_str(), "w"));
}

void Sample :: normalizeRt(){
	int size = features.size();
	for (int i = 0; i < size; i++){
		features[i].rt = (RT_TYPE) i;
	}
}

ALIGNMENT_SCORE_TYPE Sample :: alignmentScore(Sample &sample){
	Alignment alignment(*this, sample);
	alignment.augmentingAlignment();
	return alignment.alignmentScore();
}

int Sample :: maxAlignmentSizePossible(Sample &sample){
	Alignment alignment(*this, sample);
	return alignment.maxAlignmentSizePossible();
}

void Sample :: writeAlignment(Sample &sample, FILE* csvFile){
	Alignment alignment(*this, sample);
	alignment.augmentingAlignment();
	alignment.write(csvFile);
}

void Sample :: writeAlignment(Sample &sample, const string csvFileName){
	Alignment alignment(*this, sample);
	alignment.augmentingAlignment();
	//alignment.fullAlignment();
	alignment.write(csvFileName);
} 

Sample Sample :: unionAlign(Sample &sample){
	Sample result;
	Alignment alignment(*this, sample);
	alignment.augmentingAlignment();
	map<RT_TYPE, RT_TYPE> rt2to1;
	Feature f;
	for (int i = 0; i < alignment.sample1.size(); i++){
		if (alignment.sample1[i].alignedWith == NO_ALIGNMENT){
			f.rt = alignment.sample1[i].rt;
			f.mz = alignment.sample1[i].mz;
			f.intensity = alignment.sample1[i].intensity;
			f.charge = alignment.sample1[i].charge;
			f.quality = alignment.sample1[i].quality;
			result.features.push_back(f);
		}
		else {
			int k = alignment.sample1[i].alignedWith;
			f.rt = alignment.sample1[i].rt;
			f.mz = (alignment.sample1[i].mz * alignment.sample1[i].intensity + alignment.sample2[k].mz * alignment.sample2[k].intensity)
				/ (alignment.sample1[i].intensity + alignment.sample2[k].intensity);
			f.intensity = alignment.sample1[i].intensity + alignment.sample2[k].intensity;
			f.charge = alignment.sample1[i].charge;
			f.quality = (alignment.sample1[i].quality > alignment.sample2[k].quality) ? alignment.sample1[i].quality
				: alignment.sample2[k].quality;
			result.features.push_back(f);
			rt2to1[alignment.sample2[k].rt] = alignment.sample1[i].rt;
		}
	}
	
	int t = 0;
	for (int j = 0; j < alignment.sample2.size(); j++){
		if (alignment.sample2[j].alignedWith == NO_ALIGNMENT){
			auto low = rt2to1.lower_bound(alignment.sample2[j].rt), high = low--;
			if (low == rt2to1.end()) low = high;
			if (high == rt2to1.end()) high = low;
			for (int k = t; k < j; k++){
				f.rt = ((alignment.sample2[j].rt - (*low).first) * (*high).second + ((*high).first - alignment.sample2[j].rt) * (*low).second)
					/ ((*high).first - (*low).first);
				f.mz = alignment.sample2[k].mz;
				f.intensity = alignment.sample2[k].intensity;
				f.charge = alignment.sample2[k].charge;
				f.quality = alignment.sample2[k].quality;
				result.features.push_back(f);
			} 
			t = j + 1;
		}
	}
	return result;
}

Sample Sample :: intersectAlign(Sample &sample){
	Sample result;
	Alignment alignment(*this, sample);
	alignment.augmentingAlignment();
	Feature f;
	for (int i = 0; i < alignment.sample1.size(); i++){
		if (alignment.sample1[i].alignedWith != NO_ALIGNMENT){
			int k = alignment.sample1[i].alignedWith;
			f.rt = alignment.sample1[i].rt;
			f.mz = (alignment.sample1[i].mz * alignment.sample1[i].intensity + alignment.sample2[k].mz * alignment.sample2[k].intensity)
				/ (alignment.sample1[i].intensity + alignment.sample2[k].intensity);
			f.intensity = alignment.sample1[i].intensity + alignment.sample2[k].intensity;
			f.charge = alignment.sample1[i].charge;
			f.quality = (alignment.sample1[i].quality > alignment.sample2[k].quality) ? alignment.sample1[i].quality
				: alignment.sample2[k].quality;
			result.features.push_back(f);
		}
	}
	return result;
}

vector<pair<Feature, Feature> > Sample :: topAlignedFeatures(Sample &sample, int n){
	vector<pair<Feature, Feature> > result;
	Alignment alignment(*this, sample);
	alignment.augmentingAlignment();
	map<ALIGNMENT_SCORE_TYPE, int> dict;
	for (int i = 0; i < alignment.sample1.size(); i++){
		if (alignment.sample1[i].alignedWith != NO_ALIGNMENT){
			int j = alignment.sample1[i].alignedWith;
			dict[-SCORING_FUNCTION(alignment.sample1[i], alignment.sample2[j])] = i;
		}
	}
	int k = 0;
	for (auto e: dict){
		if (k++ == n) break;
		result.push_back({alignment.sample1[e.second].toFeature(), alignment.sample2[alignment.sample1[e.second].alignedWith].toFeature()});
	}
	cerr << dict.size() << " " << result.size() << endl;
	return result;
}

vector<pair<Feature, Feature> > Sample :: topAlignedFeatures(Sample &sample){
	return topAlignedFeatures(sample, features.size());
}
