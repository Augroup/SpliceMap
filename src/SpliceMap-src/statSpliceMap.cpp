

#include "statSpliceMap.h"
#include "SpliceMap_utils.h"

int main (int argc, char * const argv[]) {
    
	vector<string> filename_list;
	int lowbound_indep_read=15;  
	string temp = "";
	int num_read = 0;
	
	ifstream junction_file;
	ifstream new_junction_file;
	
	int junction_index = 0;
	int new_junction_index = 0;
	
	if (argc >= 2) {
		for (int i = 1; i<argc; i++) {
			if (argv[i][0] == 'N') {
				temp = argv[i];
				lowbound_indep_read = atoi(temp.substr(1).c_str());
			}else {
				filename_list.push_back(argv[i]);
				num_read++;
			}

		}
	}else {
		cout << "usage: ./statSpliceMap junction.bed [newjun.bed]" << endl;
		exit(1);
	}
	
	junction_file.open(filename_list[0].c_str(), ios::in);
	
	if (!junction_file.is_open()) {
		cerr << "ERROR: failed to open junction file -- " << filename_list[0] << endl;
		exit(1);
	}
	
	junction_index = file_index(filename_list[0]);
	if (junction_index < 0) {
		cerr << "ERROR: failed to open junction file extension wrong -- " << filename_list[0] << endl;
		cerr << "Only supports .bed and _ref.txt files" << endl;
		exit(2);
	}
	
	if(num_read == 2){
		new_junction_file.open(filename_list[1].c_str(), ios::in);
		
		if (!new_junction_file.is_open()) {
			cerr << "ERROR: failed to open new junction file -- " << filename_list[1] << endl;
			exit(1);
		}
		
		new_junction_index = file_index(filename_list[1]);
		
		if (new_junction_index < 0) {
			cerr << "ERROR: failed to open new junction file extension wrong -- " << filename_list[0] << endl;
			cerr << "Only supports .bed and _ref.txt files" << endl;
			exit(2);
		}
	}

	
	
	
	cout << "\n####Total####alljunction##########" << endl;
	count_read(junction_file, false, junction_index);
	junction_file.close();
	junction_file.clear();
	
	if(num_read == 2){
		cout << "\n####Total####newjunction##########" << endl;
		count_read(new_junction_file, false, new_junction_index);
		new_junction_file.close();
		new_junction_file.clear();
	}
	
	
	
	junction_file.open(filename_list[0].c_str(), ios::in);
	
	if(num_read == 2){
		new_junction_file.open(filename_list[1].c_str(), ios::in);
	}
	
	cout << "\n####nNR####alljunction#########" << endl;
	count_indep_read(junction_file,false,10,junction_index);
	junction_file.close();
	junction_file.clear();
	
	if (num_read == 2) {
		cout << "\n####nNR####newjunction#########" << endl;
		count_indep_read(new_junction_file,false,lowbound_indep_read,new_junction_index);
		new_junction_file.close();
		new_junction_file.clear();
	}
	
	if (num_read == 2) {
		new_junction_file.open(filename_list[1].c_str(), ios::in);
		print_idep_read(new_junction_file, lowbound_indep_read, new_junction_index);
		new_junction_file.close();
		new_junction_file.clear();
	}
	
    return 0;
}






void print_idep_read(ifstream &junction, int lowerbound, int index)
{

	string line = "";
	string temp = "";

	int left;
	int right;
	
	while (!junction.eof()) {
		getline(junction, line);
		trim2(line);
		
		if(line.length() != 0){
			if (line.substr(0, 5).compare("track") == 0) {
				continue;
			}
			
			vector<string> line_list = split(line, '\t');
			
			
			temp = line_list[index];
			
			
			
			
			left = (int)temp.find('_')+1;
			right = (int)temp.find(']',left);
			int num_indep = atoi(temp.substr(left, right-left).c_str());
			
			if (num_indep>lowerbound) {
				cout << line << endl;
			}
		}
	}
}

void count_indep_read(ifstream &junction,bool debug, int lowerbound, int index)
{
	int sum = 0;
	string line = "";
	string temp = "";
	map<int,int> dict;
	int left;
	int right;
	
	while (!junction.eof()) {
		getline(junction, line);
		trim2(line);
		
		if(line.length() != 0){
			if (line.substr(0, 5).compare("track") == 0) {
				continue;
			}
			
			vector<string> line_list = split(line, '\t');
			
			
			temp = line_list[index];
			
			
			
			
			left = (int)temp.find('_')+1;
			right = (int)temp.find(']',left);
			int num_indep = atoi(temp.substr(left, right-left).c_str());
			
			sum++;
			
			pair<map<int,int>::iterator,bool> out;
			out = dict.insert(pair<int,int>(num_indep,1)); 
			if (!out.second) {
				(out.first)->second++;
			}
			if (num_indep>lowerbound && debug) {
				cout << line << endl;
			}
		}
	}
	
	map<int,int>::iterator it = dict.begin();
	while (it != dict.end()) {
		
		cout << it->first << ":\t" << it->second << "\t\t" <<  
			100*((double)it->second)/sum << "%" << endl; 
		
		it++;
	}
	
}

void count_read(ifstream &junction, bool debug, int index)
{
	int c[7] = {0,0,0,0,0,0,0};
	int total = 0;
	int num_juntion = 0;
	int sum = 0;
	string temp = "";
	string line = "";
	
	while (!junction.eof()) {
		getline(junction, line);
		trim2(line);
		if(line.length() != 0){
			if (line.substr(0, 5).compare("track") == 0) {
				continue;
			}
			
			vector<string> line_list = split(line, '\t');
			
			
			temp = line_list[index];
			int left = (int)temp.find('(')+1;
			int right = (int)temp.find(')',left);
			int depth = atoi(temp.substr(left, right-left).c_str());
			
			num_juntion++;
			
			total = total + depth;
			
			if (depth == 1) {
				c[0]++;
			}else if(depth > 1 && depth <= 5){
				c[1]++;
			}else if(depth > 5 && depth <= 20){
				c[2]++;
				if(debug){
					cout << "5-20 :" << line << endl;
				}
			}else if(depth > 20 && depth <= 50){
				c[3]++;
				if(debug){
					cout << "20-50 :" << line << endl;
				}
			}else if(depth > 50 && depth <= 200){
				c[4]++;
				if(debug){
					cout << "50-200 :" << line << endl;
				}
			}else if(depth > 200 && depth <= 1000){
				c[5]++;
				if(debug){
					cout << "200-1000 :" << line << endl;
				}
			}else if(depth > 1000){
				c[6]++;
				if(debug){
					cout << "1000+ :" << line << endl;
				}
			}
			
		}
		

	}
	
	for (int i = 0; i<7; i++) {
		sum = sum + c[i];
	}
	
	double ave = ((double)total)/num_juntion;
	
	int col1 = 13;
	int col2 = 12;
	int col3 = 6;
	
	cout << setw(col1) << "total: " << total << endl;
	cout << setw(col1) << "ave: " << ave << endl;
	
	cout << setw(col1) << "1: " 
	<< setiosflags(ios::left) << setw(col2)<< c[0] 
	<< resetiosflags(ios::left)
	<< setw(col3) << setiosflags(ios::fixed) << setprecision(2) << 100.0*(((double)c[0]))/sum 
	<< "%"  << endl;
	
	cout << setw(col1) << "2-5: " 
	<< setiosflags(ios::left) << setw(col2)<< c[1] 
	<< resetiosflags(ios::left)
	<< setw(col3)  << 100.0*(((double)c[1]))/sum 
	<< "%"  << endl;
	
	cout << setw(col1) << "6-20: " 
	<< setiosflags(ios::left) << setw(col2)<< c[2] 
	<< resetiosflags(ios::left)
	<< setw(col3)  << 100.0*(((double)c[2]))/sum 
	<< "%"  << endl;
	
	cout << setw(col1) << "21-50: " 
	<< setiosflags(ios::left) << setw(col2)<< c[3] 
	<< resetiosflags(ios::left)
	<< setw(col3)  << 100.0*(((double)c[3]))/sum 
	<< "%"  << endl;
	
	cout << setw(col1) << "51-200: " 
	<< setiosflags(ios::left) << setw(col2)<< c[4] 
	<< resetiosflags(ios::left)
	<< setw(col3) << 100.0*(((double)c[4]))/sum 
	<< "%" << endl;
	
	cout << setw(col1) << "201-1000: " 
	<< setiosflags(ios::left) << setw(col2)<< c[5] 
	<< resetiosflags(ios::left)
	<< setw(col3)  << 100.0*(((double)c[5]))/sum 
	<< "%"  << endl;
	
	cout << setw(col1) << "1000+: " 
	<< setiosflags(ios::left) << setw(col2)<< c[6] 
	<< resetiosflags(ios::left)
	<< setw(col3)  << 100.0*(((double)c[6]))/sum
	<< "%"  << endl;
}





int file_index(string filename)
{
	int index = -1;
	string::size_type len = filename.length();
	if (filename.substr(len-4).compare(".bed") == 0) {
		index = 3;
	}
	else if(filename.substr(len-8).compare("_ref.txt") == 0){
		index = 0;
	}
	
	return index;
}




