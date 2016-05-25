#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

map<int, int> unionfind;
int Find(int s){
	if (unionfind.count(s) == 0) unionfind[s] = s;
	if (unionfind[s] == s) return s;
	return unionfind[s] = Find(unionfind[s]);
}
void Union(int s1, int s2){
	if (unionfind.count(s1) == 0) unionfind[s1] = s1;
	if (unionfind.count(s2) == 0) unionfind[s2] = s2;
	if (Find(s1) == Find(s2)) return;
	unionfind[Find(s1)] = Find(s2);
}

struct Similarity{
	int id1, id2;
	double score;
	Similarity(int i, int j, double k): id1(i), id2(j), score(k){}
	bool operator < (Similarity &other){
		if (score > other.score) return true;
		if (score < other.score) return false;
		if (id1 < other.id1) return true;
		if (id1 > other.id1) return false;
		if (id2 < other.id2) return true;
		return false;
	}
};

int main(int argc, char** argv) {
	ifstream fin("MatrixWithLabel.tsv");
	ofstream fout("EnvironmentMinimumSpanningTree.sif");
	
	string headline, s;
	vector<string> name, path;
	vector<vector<double> > matrix;
	vector<int> mostSimilarID;
	vector<double> bestSimilarityScore;
	vector<bool> validPerson, validEnvironment;
	getline(fin, headline);
	stringstream ss(headline);
	vector<Similarity> sim;
	
	while (ss >> s) path.push_back(s);
	for (int i = 0; i < path.size(); i++){
		fin >> s;
		name.push_back(s);
		matrix.push_back({});
		mostSimilarID.push_back(-1);
		bestSimilarityScore.push_back(-1);
		for (int j = 0; j < path.size(); j++){
			double t;
			fin >> t;
			matrix[i].push_back(t);
			if (i != j && bestSimilarityScore[i] < matrix[i][j]){
				mostSimilarID[i] = j;
				bestSimilarityScore[i] = matrix[i][j];
			}
		}
	}
	for (int i = 0; i < path.size(); i++){
		if (name[i].find("Person") != string::npos && (name[mostSimilarID[i]].find("Person") != string::npos
			|| name[mostSimilarID[i]].find("Environment") != string::npos))
			validPerson.push_back(true);
		else validPerson.push_back(false);
		if (name[i].find("Environment") != string::npos && (name[mostSimilarID[i]].find("Person") != string::npos
			|| name[mostSimilarID[i]].find("Environment") != string::npos))
			validEnvironment.push_back(true);
		else validEnvironment.push_back(false);
	}
	for (int i = 0; i < path.size(); i++){
		if (validEnvironment[i] == false) continue;
		for (int j = i + 1; j < path.size(); j++){
			if (validEnvironment[j] == false) continue;
			sim.push_back(Similarity(i, j, matrix[i][j]));
		}
	}
	sort(sim.begin(), sim.end());
	for (auto &e: sim){
		if (Find(e.id1) == Find(e.id2)) continue;
		fout << name[e.id1] << " ee " << name[e.id2] << endl;
		Union(e.id1, e.id2);
	}
	for (int i = 0; i < path.size(); i++){
		if (validPerson[i] == false) continue;
		int closestEnvironment = -1;
		double closestSimilarity = -1;
		for (int j = 0; j < path.size(); j++){
			if (validEnvironment[j] == false) continue;
			if (closestSimilarity < matrix[i][j]){
				closestEnvironment = j;
				closestSimilarity = matrix[i][j];
			}
		}
		fout << name[i] << " pe " << name[closestEnvironment] << endl;
	}
	return 0;
}
