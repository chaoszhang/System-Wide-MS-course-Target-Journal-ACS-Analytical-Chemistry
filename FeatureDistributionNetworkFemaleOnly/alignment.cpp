#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <queue>
#include <tuple>
#include <cmath> 
#include "alignment.hpp"

using namespace std;

Alignment :: Alignment(Sample &s1, Sample &s2){
	s1.sortByRt();
	s2.sortByRt();
	isSortedByRt = true;
	for (Feature &f: s1.features)
		sample1.push_back(FeatureAlignment(f));
	for (Feature &f: s2.features)
		sample2.push_back(FeatureAlignment(f));
	if (s1.sampleName == s2.sampleName) sampleName = s1.sampleName;
	sortByRt();
	for (int i = 0, n = sample1.size(); i < n; i++)
		sample1[i].rt_order = i;
	for (int i = 0, n = sample2.size(); i < n; i++)
		sample2[i].rt_order = i;
	sortByMz();
	for (int i = 0, n = sample1.size(); i < n; i++)
		sample1[i].mz_order = i;
	for (int i = 0, n = sample2.size(); i < n; i++)
		sample2[i].mz_order = i;
	sortByIntensity();
	for (int i = 0, n = sample1.size(); i < n; i++)
		sample1[i].intensity_order = i;
	for (int i = 0, n = sample2.size(); i < n; i++)
		sample2[i].intensity_order = i;
}

void Alignment :: sortByRt(){
	if (isSortedByRt) return;
	sort(sample1.begin(), sample1.end(), FeatureAlignment :: rtCompare);
	sort(sample2.begin(), sample2.end(), FeatureAlignment :: rtCompare);
	isSortedByRt = true;
	isSortedByMz = false;
	isSortedByIntensity = false;
}

void Alignment :: sortByMz(){
	if (isSortedByMz) return;
	sort(sample1.begin(), sample1.end(), FeatureAlignment :: mzCompare);
	sort(sample2.begin(), sample2.end(), FeatureAlignment :: mzCompare);
	isSortedByRt = false;
	isSortedByMz = true;
	isSortedByIntensity = false;
}

void Alignment :: sortByIntensity(){
	if (isSortedByIntensity) return;
	sort(sample1.begin(), sample1.end(), FeatureAlignment :: intensityCompare);
	sort(sample2.begin(), sample2.end(), FeatureAlignment :: intensityCompare);
	isSortedByRt = false;
	isSortedByMz = false;
	isSortedByIntensity = true;
}

const int MAXZ = 1 << RT_DIFFERENCE_THRESHOLD, MASK = MAXZ - 1;

class ScoreBinaryIndexedTree{
	inline static int lowBit(int i){
		return -i & i;
	}
	
public:
	const static ALIGNMENT_SCORE_TYPE oo;
	ScoreBinaryIndexedTree(int sz): size(sz + 1 + RT_DIFFERENCE_THRESHOLD){
		for (int k = 0; k < MAXZ; k++){
			scores[k].resize(size);
			states[k].resize(size);
			for (int j = 0; j < size; j++){
				scores[k][j] = -oo;
			}
		}
	}
	
	void set(State st, ALIGNMENT_SCORE_TYPE value){
		int j = st.y, k = MASK & st.z;
		for (; k != 0 && j < size && value > scores[k][j]; j++, k = MASK & (k << 1)){
			scores[k][j] = value;
			states[k][j] = st;
		}
		if (k != 0) return;
		for (; j < size && value > scores[0][j]; j += lowBit(j)){
			scores[0][j] = value;
			states[0][j] = st;
		}
	}
	
	pair<State, ALIGNMENT_SCORE_TYPE> get(State st){
		if (st.z != 0) return {states[st.z][st.y], scores[st.z][st.y]};
		int res = 0;
		for (int i = st.y; i != 0; i -= lowBit(i)){
			if (scores[0][i] > scores[0][res]) res = i;
		}
		return {states[0][res], scores[0][res]};
	}
private:
	const int size;
	vector<ALIGNMENT_SCORE_TYPE> scores[MAXZ];
	vector<State> states[MAXZ];
};
const ALIGNMENT_SCORE_TYPE ScoreBinaryIndexedTree :: oo = 1e20;

class ReverseScoreBinaryIndexedTree{
	const static ALIGNMENT_SCORE_TYPE oo;
	
	inline static int lowBit(int i){
		return -i & i;
	}
	
public:
	ReverseScoreBinaryIndexedTree(int sz): size(sz + 1 + RT_DIFFERENCE_THRESHOLD){
		for (int k = 0; k < MAXZ; k++){
			scores[k].resize(size);
			states[k].resize(size);
			for (int j = 0; j < size; j++){
				scores[k][j] = -oo;
			}
		}
	}
	
	void set(State st, ALIGNMENT_SCORE_TYPE value){
		int j = st.y, k = MASK & st.z;
		for (; k != 0 && j >= 0 && value > scores[k][j]; j--, k >>= 1){
			scores[k][j] = value;
			states[k][j] = st;
		}
		if (k != 0) return;
		for (; j >= 0 && value > scores[0][j]; j -= lowBit(j)){
			scores[0][j] = value;
			states[0][j] = st;
		}
	}
	
	pair<State, ALIGNMENT_SCORE_TYPE> get(State st){
		if (st.z != 0) return {states[st.z][st.y], scores[st.z][st.y]};
		if (st.y == 0) return {states[st.z][st.y], scores[st.z][st.y]};
		int res = st.y;
		for (int i = st.y; i < size; i += lowBit(i)){
			if (scores[0][i] > scores[0][res]) res = i;
		}
		return {states[0][res], scores[0][res]};
	}
private:
	const int size;
	vector<ALIGNMENT_SCORE_TYPE> scores[MAXZ];
	vector<State> states[MAXZ];
};
const ALIGNMENT_SCORE_TYPE ReverseScoreBinaryIndexedTree :: oo = 1e20;

inline int mzFind(MZ_TYPE mz, vector<int> &mzOrderS2, vector<FeatureAlignment> &sample2){
	int left = -1, right = mzOrderS2.size();
	while (left + 1 < right){
		int mid = (left + right) >> 1;
		if (sample2[mzOrderS2[mid]].mz < mz) left = mid;
		else right = mid;
	}
	return right;
}

struct MzOrderS2Sort{
	const vector<int> &mzOrderS2;
	const vector<FeatureAlignment> &sample2;
	MzOrderS2Sort(const vector<int> &mzs2, const vector<FeatureAlignment> &s2): mzOrderS2(mzs2), sample2(s2){}
	bool operator()(const int i, const int j){
		return sample2[i].mz_order < sample2[j].mz_order;
	}
};

ALIGNMENT_SCORE_TYPE Alignment :: fullAlignmentScore(){
	//return localAlignment(State(0, 0, MASK), State(sample1.size(), sample2.size(), MASK));
	vector<int> mzOrderS2(sample2.size());
	MzOrderS2Sort mzOrderS2Sort(mzOrderS2, sample2); 
	sortByRt();
	for (int i = 0, n = sample2.size(); i < n; i++){
		mzOrderS2[i] = i;
	}
	sort(mzOrderS2.begin(), mzOrderS2.end(), mzOrderS2Sort);
	ScoreBinaryIndexedTree score(sample2.size());
	score.set(State(0, 0, MASK), 0);
	for (int i = 0, n = sample1.size(), m = sample2.size(); i < n; i++){
		vector<pair<State, ALIGNMENT_SCORE_TYPE> > newStates;
		for (int mzIndex = mzFind(sample1[i].mz - MZ_DIFFERENCE_THRESHOLD, mzOrderS2, sample2),
			mzMax = mzFind(sample1[i].mz + MZ_DIFFERENCE_THRESHOLD, mzOrderS2, sample2); mzIndex < mzMax; mzIndex++){
			int j = mzOrderS2[mzIndex];
			ALIGNMENT_SCORE_TYPE curScore = SCORING_FUNCTION(sample1[i], sample2[j]);
			for (int k = 0; k < MAXZ; k++){
				for (int d = 0; d < RT_DIFFERENCE_THRESHOLD && j + d < m; d++){
					int bit = 1 << d;
					if (k & bit) continue;
					State from(i + 1, j + 1 + d, k), to(i + 1, j + 1 + d, k | bit);
	
					newStates.push_back({to, score.get(from).second + curScore});
				}
			}
		}
		for (pair<State, ALIGNMENT_SCORE_TYPE> s : newStates){
			score.set(s.first, s.second);
		}
	}
	int max_score = 0;
	State cur(sample1.size() + RT_DIFFERENCE_THRESHOLD, sample2.size() + RT_DIFFERENCE_THRESHOLD, 0);
	State st = score.get(cur).first;
	cerr << st.x << ", " << st.y << ", " << st.z << endl;
	return score.get(cur).second;
}

int Alignment :: maxAlignmentSizePossible(){
	sortByMz();
	int cnt = 0, i = 0, j = 0, n = sample1.size(), m = sample2.size();
	while (i < n && j < m){
		if (sample1[i].mz + MZ_DIFFERENCE_THRESHOLD < sample2[j].mz) i++;
		else if (sample2[j].mz + MZ_DIFFERENCE_THRESHOLD < sample1[i].mz) j++;
		else{ i++; j++; cnt++;}
	}
	return cnt;
}

ALIGNMENT_SCORE_TYPE Alignment :: localAlignment(State left, State right){
	//cerr << "l: (" << left.x << "," << left.y << "), (" << right.x << "," << right.y << ")" << endl;
	int minY = (left.y - RT_DIFFERENCE_THRESHOLD > 0) ? left.y - RT_DIFFERENCE_THRESHOLD : 0;
	vector<int> mzOrderS2(right.y - minY);
	MzOrderS2Sort mzOrderS2Sort(mzOrderS2, sample2); 
	sortByRt();
	for (int i = minY, n = right.y; i < n; i++){
		mzOrderS2[i - minY] = i;
	}
	sort(mzOrderS2.begin(), mzOrderS2.end(), mzOrderS2Sort);
	ScoreBinaryIndexedTree score(sample2.size());
	score.set(left, 0);
	map<State, tuple<State, ALIGNMENT_SCORE_TYPE, int> > backtrack;
	for (int i = left.x, n = right.x, m = right.y; i < n; i++){
		vector<pair<State, ALIGNMENT_SCORE_TYPE> > newStates;
		for (int mzIndex = mzFind(sample1[i].mz - MZ_DIFFERENCE_THRESHOLD, mzOrderS2, sample2),
			mzMax = mzFind(sample1[i].mz + MZ_DIFFERENCE_THRESHOLD, mzOrderS2, sample2); mzIndex < mzMax; mzIndex++){
			int j = mzOrderS2[mzIndex];
			ALIGNMENT_SCORE_TYPE curScore = SCORING_FUNCTION(sample1[i], sample2[j]);
			for (int k = 0; k < MAXZ; k++){
				for (int d = 0; d < RT_DIFFERENCE_THRESHOLD && j + d < m; d++){
					int bit = 1 << d;
					if (k & bit) continue;
					State from(i + 1, j + 1 + d, k), to(i + 1, j + 1 + d, k | bit);
					ALIGNMENT_SCORE_TYPE tempScore = score.get(from).second + curScore;
					newStates.push_back({to, tempScore});
					if (backtrack.count(to) == 0 || (backtrack.count(to) && get<1>(backtrack[to]) < tempScore))
						backtrack[to] = make_tuple(score.get(from).first, tempScore, j);
				}
			}
		}
		for (pair<State, ALIGNMENT_SCORE_TYPE> s : newStates){
			score.set(s.first, s.second);
		}
	}
	ALIGNMENT_SCORE_TYPE result = score.get(right).second;
	State cur = score.get(right).first;
	while (cur.x != left.x){
		sample1[cur.x - 1].alignedWith = get<2>(backtrack[cur]);
		sample2[get<2>(backtrack[cur])].alignedWith = cur.x - 1;
		cur = get<0>(backtrack[cur]);
	}
	return result;
}

vector<vector<ALIGNMENT_SCORE_TYPE> > Alignment :: forwardAlignment(State left, State right){
	//cerr << "f: (" << left.x << "," << left.y << "), (" << right.x << "," << right.y << ")" << endl;
	int minY = (left.y - RT_DIFFERENCE_THRESHOLD > 0) ? left.y - RT_DIFFERENCE_THRESHOLD : 0;
	vector<int> mzOrderS2(right.y - minY);
	MzOrderS2Sort mzOrderS2Sort(mzOrderS2, sample2); 
	sortByRt();
	for (int i = minY, n = right.y; i < n; i++){
		mzOrderS2[i - minY] = i;
	}
	sort(mzOrderS2.begin(), mzOrderS2.end(), mzOrderS2Sort);
	ScoreBinaryIndexedTree score(sample2.size());
	score.set(left, 0);
	for (int i = left.x, n = right.x, m = right.y; i < n; i++){
		vector<pair<State, ALIGNMENT_SCORE_TYPE> > newStates;
		for (int mzIndex = mzFind(sample1[i].mz - MZ_DIFFERENCE_THRESHOLD, mzOrderS2, sample2),
			mzMax = mzFind(sample1[i].mz + MZ_DIFFERENCE_THRESHOLD, mzOrderS2, sample2); mzIndex < mzMax; mzIndex++){
			int j = mzOrderS2[mzIndex];
			ALIGNMENT_SCORE_TYPE curScore = SCORING_FUNCTION(sample1[i], sample2[j]);
			for (int k = 0; k < MAXZ; k++){
				for (int d = 0; d < RT_DIFFERENCE_THRESHOLD && j + d < m; d++){
					int bit = 1 << d;
					if (k & bit) continue;
					State from(i + 1, j + 1 + d, k), to(i + 1, j + 1 + d, k | bit);
					ALIGNMENT_SCORE_TYPE tempScore = score.get(from).second + curScore;
					newStates.push_back({to, tempScore});
				}
			}
		}
		for (pair<State, ALIGNMENT_SCORE_TYPE> s : newStates){
			score.set(s.first, s.second);
		}
	}
	vector<vector<ALIGNMENT_SCORE_TYPE> > result;
	for (int j = left.y; j <= right.y; j++){
		result.push_back({});
		for (int k = 0; k < MAXZ; k++){
			State cur(right.x, j, k);
			result.back().push_back(score.get(cur).second);
		}
	}
	return result;
}

vector<vector<ALIGNMENT_SCORE_TYPE> > Alignment :: reverseAlignment(State left, State right){
	//cerr << "b: (" << left.x << "," << left.y << "), (" << right.x << "," << right.y << ")" << endl;
	int minY = (left.y - RT_DIFFERENCE_THRESHOLD > 0) ? left.y - RT_DIFFERENCE_THRESHOLD : 0;
	vector<int> mzOrderS2(right.y - minY);
	MzOrderS2Sort mzOrderS2Sort(mzOrderS2, sample2); 
	sortByRt();
	for (int i = minY, n = right.y; i < n; i++){
		mzOrderS2[i - minY] = i;
	}
	sort(mzOrderS2.begin(), mzOrderS2.end(), mzOrderS2Sort);
	ReverseScoreBinaryIndexedTree score(sample2.size());
	score.set(right, 0);
	for (int i = right.x - 2, n = left.x - 1, m = right.y; i >= n; i--){
		vector<pair<State, ALIGNMENT_SCORE_TYPE> > newStates;
		for (int mzIndex = mzFind(sample1[i].mz - MZ_DIFFERENCE_THRESHOLD, mzOrderS2, sample2),
			mzMax = mzFind(sample1[i].mz + MZ_DIFFERENCE_THRESHOLD, mzOrderS2, sample2); mzIndex < mzMax; mzIndex++){
			int j = mzOrderS2[mzIndex];
			ALIGNMENT_SCORE_TYPE curScore = SCORING_FUNCTION(sample1[i], sample2[j]);
			for (int k = 0; k < MAXZ; k++){
				for (int d = 0; d < RT_DIFFERENCE_THRESHOLD && j + d < m; d++){
					int bit = 1 << d;
					if (k & bit) continue;
					State from(i + 1, j + 1 + d, k), to(i + 1, j + 1 + d, k | bit);
					ALIGNMENT_SCORE_TYPE tempScore = score.get(from).second + curScore;
					newStates.push_back({to, tempScore});
				}
			}
		}
		for (pair<State, ALIGNMENT_SCORE_TYPE> s : newStates){
			score.set(s.first, s.second);
		}
	}
	vector<vector<ALIGNMENT_SCORE_TYPE> > result;
	for (int j = left.y; j <= right.y; j++){
		result.push_back({});
		for (int k = 0; k < MAXZ; k++){
			State cur(left.x, j, k);
			result.back().push_back(score.get(cur).second);
		}
	}
	return result;
}

void Alignment :: recursiveAlignment(State left, State right){
	//cerr << "r: (" << left.x << "," << left.y << "), (" << right.x << "," << right.y << ")" << endl;
	if (right.x - left.x < RECURSION_THRESHOLD) {
		localAlignment(left, State(right.x - 1, right.y, MASK & ~right.z));
		return;
	}
	int midX = (left.x + right.x) / 2;
	State midL(midX, right.y, 0), midR(midX + 1, left.y, 0);
	vector<vector<ALIGNMENT_SCORE_TYPE> > forward = forwardAlignment(left, midL);
	vector<vector<ALIGNMENT_SCORE_TYPE> > reverse = reverseAlignment(midR, right);
	int jValue = 0, k1Value = 0, k2Value = 0, fValue, rValue;
	ALIGNMENT_SCORE_TYPE t = -(ScoreBinaryIndexedTree :: oo);
	for (int j = 0; j < forward.size(); j++){
		vector<ALIGNMENT_SCORE_TYPE> &f = forward[j], &r = reverse[j];
		for (int k1 = 0; k1 < f.size(); k1++){
			for (int k2 = 0; k2 < r.size(); k2++){
				if (k1 & k2) continue;
				if (f[k1] + r[k2] > t){
					jValue = j + left.y;
					k1Value = k1;
					k2Value = k2;
					t = f[k1] + r[k2];
					fValue = f[k1];
					rValue = r[k2];
				}
			}
		}
	}
	recursiveAlignment(left, State(midX + 1, jValue, MASK & ~k1Value));
	recursiveAlignment(State(midX, jValue, MASK & ~k2Value), right);
}

ALIGNMENT_SCORE_TYPE Alignment :: fullAlignment(){
	vector<vector<ALIGNMENT_SCORE_TYPE> > forward = forwardAlignment(
		State(0, 0, MASK), State(sample1.size(), sample2.size(), 0));
	int jValue = 0, kValue = 0;
	ALIGNMENT_SCORE_TYPE t = -(ScoreBinaryIndexedTree :: oo), score = 0;
	for (int j = 0; j < forward.size(); j++){
		vector<ALIGNMENT_SCORE_TYPE> &temp = forward[j];
		for (int k = 0; k < temp.size(); k++){
			if (temp[k] > t){
				jValue = j;
				kValue = k;
				t = temp[k];
			}
		}
	}
	recursiveAlignment(State(0, 0, MASK), State(sample1.size() + 1, jValue, MASK & ~kValue));
	for (FeatureAlignment &s: sample1){
		if (s.alignedWith != NO_ALIGNMENT){
			score += SCORING_FUNCTION(s, sample2[s.alignedWith]);
		}
	}
	return score;
}

void Alignment :: write(FILE* file){
	fprintf(file, "#rt1\tmz1\tintensity1\trt2\tmz2\tintensity2\n");
	for (FeatureAlignment &feature : sample1)
		if (feature.alignedWith != NO_ALIGNMENT)
			fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", (double) feature.rt, (double) feature.mz, (double) feature.intensity,
				(double) sample2[feature.alignedWith].rt, (double) sample2[feature.alignedWith].mz, (double) sample2[feature.alignedWith].intensity);
	for (FeatureAlignment &feature : sample1)
		if (feature.alignedWith == NO_ALIGNMENT)
			fprintf(file, "%lf\t%lf\t%lf\t-1\t-1\t-1\n", (double) feature.rt, (double) feature.mz, (double) feature.intensity);
	for (FeatureAlignment &feature : sample2)
		if (feature.alignedWith == NO_ALIGNMENT)
			fprintf(file, "-1\t-1\t-1\t%lf\t%lf\t%lf\n", (double) feature.rt, (double) feature.mz, (double) feature.intensity);
	fclose(file);
}

void Alignment :: write(const string fileName){
	write(fopen(fileName.c_str(), "w"));
}

class ScoreRecursiveBinaryIndexedTree{
	inline static int lowBit(int i){
		return -i & i;
	}
	
public:
	const static ALIGNMENT_SCORE_TYPE oo;
	ScoreRecursiveBinaryIndexedTree(int sz): size(sz + 1), states(size, make_tuple(0, 0, -oo)){
		std::get<2>(states[0]) = 0;
	}
	
	void set(int x, int y, ALIGNMENT_SCORE_TYPE value){
		for (int j = y; j < size && value > std::get<2>(states[j]); j += lowBit(j)){
			states[j] = make_tuple(x, y, value);
		}
	}
	
	tuple<int, int, ALIGNMENT_SCORE_TYPE> get(int y){
		int res = 0;
		for (int j = y; j != 0; j -= lowBit(j)){
			if (std::get<2>(states[j]) > std::get<2>(states[res])) res = j;
		}
		return states[res];
	}
private:
	const int size;
	vector<tuple<int, int, ALIGNMENT_SCORE_TYPE> > states;
};
const ALIGNMENT_SCORE_TYPE ScoreRecursiveBinaryIndexedTree :: oo = 1e20;

class FeatureAlignmentPointer{
	FeatureAlignment* pointer;

public:
	FeatureAlignmentPointer(FeatureAlignment &feature): pointer(&feature){}
	FeatureAlignment &operator()(){
		return *pointer;
	}
	
	bool operator < (FeatureAlignmentPointer &other){
		if (pointer->rt < other.pointer->rt) return true;
		if (pointer->rt > other.pointer->rt) return false;
		return pointer < other.pointer;
	}
};

void augmentingPath(vector<FeatureAlignment> &s1, vector<FeatureAlignment> &s2, int sampleSize1, int sampleSize2, int windowSize){
	vector<FeatureAlignmentPointer> sample1, sample2;
	vector<int> s1rtOrderS2;
	int s1Front = 0, s1Back = 0;
	set<int> rtOrderS2;
	map<int, int> intensityS1toRtS1, rtS1toRtS2;
	set<pair<MZ_TYPE, int> > mzOrderS2;
	ScoreRecursiveBinaryIndexedTree scoreBIT(sampleSize2);
	vector<tuple<int, int, ALIGNMENT_SCORE_TYPE> > update;
	vector<tuple<int, int, ALIGNMENT_SCORE_TYPE, int, int> > backtrack;
	
	for (int i = 0; i < sampleSize1; i++)
		sample1.push_back(s1[i]);
	for (int j = 0; j < sampleSize2; j++)
		sample2.push_back(s2[j]);
	sort(sample1.begin(), sample1.end());
	sort(sample2.begin(), sample2.end());
	s1rtOrderS2.push_back(-1);
	for (int i = 0; i < sampleSize1; i++){
		if (sample1[i]().alignedWith == NO_ALIGNMENT) continue;
		intensityS1toRtS1[s2[sample1[i]().alignedWith].alignedWith] = i;
	}
	for (int j = 0; j < sampleSize2; j++){
		if (sample2[j]().alignedWith == NO_ALIGNMENT) continue;
		rtS1toRtS2[intensityS1toRtS1[sample2[j]().alignedWith]] = j;
	}
	for (int i = 0; i < sampleSize1; i++){
		if (sample1[i]().alignedWith == NO_ALIGNMENT) continue;
		s1rtOrderS2.push_back(rtS1toRtS2[i]);
	}
	s1rtOrderS2.push_back(sampleSize2);
	rtS1toRtS2[sampleSize1] = sampleSize2;
	s1Front = (windowSize < s1rtOrderS2.size()) ? windowSize - 1: s1rtOrderS2.size() - 1;
	for (int i = 0; i <= s1rtOrderS2[s1Front]; i++){
		rtOrderS2.insert(rtS1toRtS2[i]);
	}
	for (int j = 0; j < *(rtOrderS2.rbegin()); j++){
		if (sample2[j]().alignedWith == NO_ALIGNMENT) mzOrderS2.insert({sample2[j]().mz, j});
	}
		
	for (int i = 0; i < sampleSize1; i++){
		if (sample1[i]().alignedWith == NO_ALIGNMENT){
			MZ_TYPE mzMin = sample1[i]().mz - MZ_DIFFERENCE_THRESHOLD, mzMax = sample1[i]().mz + MZ_DIFFERENCE_THRESHOLD;
			for (auto it = mzOrderS2.lower_bound({mzMin, -1}); it != mzOrderS2.end() && (*it).first < mzMax; it++){
				int j = (*it).second;
				tuple<int, int, ALIGNMENT_SCORE_TYPE> state = scoreBIT.get(j);
				ALIGNMENT_SCORE_TYPE score = get<2>(state) + SCORING_FUNCTION(sample1[i](), sample2[j]());
				update.push_back(make_tuple(i + 1, j + 1, score));
				backtrack.push_back(make_tuple(i + 1, j + 1, score, get<0>(state), get<1>(state)));
			}
			while (update.size()){
				scoreBIT.set(get<0>(update.back()), get<1>(update.back()), get<2>(update.back()));
				update.pop_back();
			}
		}
		else{
			if (s1Front < s1rtOrderS2.size() - 1){
				s1Front++;
				if (s1rtOrderS2[s1Front] < *(rtOrderS2.begin())){
					int k = *(rtOrderS2.begin());
					for (int j = s1rtOrderS2[s1Front] + 1; j < k; j++){
						if (sample2[j]().alignedWith == NO_ALIGNMENT) mzOrderS2.insert({sample2[j]().mz, j});
					}
				}
				if (s1rtOrderS2[s1Front] > *(rtOrderS2.rbegin())){
					int k = s1rtOrderS2[s1Front];
					for (int j = *(rtOrderS2.rbegin()) + 1; j < k; j++){
						if (sample2[j]().alignedWith == NO_ALIGNMENT) mzOrderS2.insert({sample2[j]().mz, j});
					}
				}
				rtOrderS2.insert(s1rtOrderS2[s1Front]);
			}
			MZ_TYPE mzMin = sample1[i]().mz - MZ_DIFFERENCE_THRESHOLD, mzMax = sample1[i]().mz + MZ_DIFFERENCE_THRESHOLD;
			for (auto it = mzOrderS2.lower_bound({mzMin, -1}); it != mzOrderS2.end() && (*it).first < mzMax; it++){
				int j = (*it).second;
				tuple<int, int, ALIGNMENT_SCORE_TYPE> state = scoreBIT.get(j);
				ALIGNMENT_SCORE_TYPE score = get<2>(state) + SCORING_FUNCTION(sample1[i](), sample2[j]())
					- SCORING_FUNCTION(sample1[i](), s2[sample1[i]().alignedWith]);
				if (sample2[j]().intensity_order == sample1[i]().alignedWith)
					score = get<2>(state) + SCORING_FUNCTION(sample1[i](), sample2[j]())
						+ AFFINITY_BONUS * SCORING_FUNCTION(sample1[i](), s2[sample1[i]().alignedWith]);
				update.push_back(make_tuple(i + 1, j + 1, score));
				backtrack.push_back(make_tuple(i + 1, j + 1, score, get<0>(state), get<1>(state)));
			}
			while (update.size()){
				scoreBIT.set(get<0>(update.back()), get<1>(update.back()), get<2>(update.back()));
				update.pop_back();
			}
			if (s1Back <= s1Front - 2 * windowSize){
				rtOrderS2.erase(s1rtOrderS2[s1Back]);
				if (s1rtOrderS2[s1Back] < *(rtOrderS2.begin())){
					int k = *(rtOrderS2.begin());
					for (int j = s1rtOrderS2[s1Back] + 1; j < k; j++){
						if (sample2[j]().alignedWith == NO_ALIGNMENT) mzOrderS2.erase({sample2[j]().mz, j});
					}
				}
				if (s1rtOrderS2[s1Back] > *(rtOrderS2.rbegin())){
					int k = s1rtOrderS2[s1Back];
					for (int j = *(rtOrderS2.rbegin()) + 1; j < k; j++){
						if (sample2[j]().alignedWith == NO_ALIGNMENT) mzOrderS2.erase({sample2[j]().mz, j});
					}
				}
				s1Back++;
			}
		}
	}
	int x = get<0>(scoreBIT.get(sampleSize2)), y = get<1>(scoreBIT.get(sampleSize2));
	//cerr << x << y << endl;
	while(backtrack.size()){
		if (get<0>(backtrack.back()) == x && get<1>(backtrack.back()) == y){
			if (sample1[x - 1]().alignedWith != NO_ALIGNMENT){
				s2[sample1[x - 1]().alignedWith].alignedWith = NO_ALIGNMENT;
			}
			sample1[x - 1]().alignedWith = sample2[y - 1]().intensity_order;
			sample2[y - 1]().alignedWith = sample1[x - 1]().intensity_order;
			x = get<3>(backtrack.back());
			y = get<4>(backtrack.back());
		}
		backtrack.pop_back();
	}
}

void Alignment :: augmentingAlignment(int initSize1, int initSize2, int multiplier, int windowSize, int fullAlignmentRound){
	sortByIntensity();
	int currentSize1 = (initSize1 < sample1.size()) ? initSize1 : sample1.size();
	int currentSize2 = (initSize2 < sample2.size()) ? initSize2 : sample2.size();
	while (currentSize1 < sample1.size() || currentSize2 < sample2.size()){
		augmentingPath(sample1, sample2, currentSize1, currentSize2, windowSize);
		currentSize1 = (currentSize1 * multiplier < sample1.size()) ? currentSize1 * multiplier : sample1.size();
		currentSize2 = (currentSize2 * multiplier < sample2.size()) ? currentSize2 * multiplier : sample2.size();
	} 
	for (int i = 0; i < fullAlignmentRound; i++){
		augmentingPath(sample1, sample2, currentSize1, currentSize2, windowSize);
	}
}

ALIGNMENT_SCORE_TYPE Alignment :: alignmentScore(){
	//cerr << "in" << endl;
	ALIGNMENT_SCORE_TYPE score = 0, score1 = 0, score2 = 0;
	vector<ALIGNMENT_SCORE_TYPE> scores1, scores2;
	int num = (sample1.size() < sample2.size()) ? sample1.size() : sample2.size();
	for (FeatureAlignment &feature : sample1)
		if (feature.alignedWith != NO_ALIGNMENT)
			score += SCORING_FUNCTION(feature, sample2[feature.alignedWith]);
	for (FeatureAlignment &feature : sample1)
		scores1.push_back(SCORING_FUNCTION(feature, feature));
	for (FeatureAlignment &feature : sample2)
		scores2.push_back(SCORING_FUNCTION(feature, feature));
	sort(scores1.begin(), scores1.end());
	sort(scores2.begin(), scores2.end()); 
	for (int i = sample1.size() - 1; i >= (int)sample1.size() - num; i--) score1 += scores1[i];
	for (int i = sample2.size() - 1; i >= (int)sample2.size() - num; i--) score2 += scores2[i];
	//cerr << "out" << endl;
	return score / sqrt(score1 * score2);
}

int Alignment :: alignmentCount(){
	int cnt = 0;
	for (FeatureAlignment &feature : sample1)
		if (feature.alignedWith != NO_ALIGNMENT)
			cnt++;
	return cnt;
}
