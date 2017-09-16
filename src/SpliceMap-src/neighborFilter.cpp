

#include "neighborFilter.h"



int main (int argc, char * const argv[]) 
{
	int LIM = 80;
	
	string line;
	ifstream in_file;
	ofstream out_file;
	
	
    map<string,vector<int> > good_point_chr; 
	
	if (argc < 4 || argc > 5) {
		print_usage_and_exit();
	}
	
	string sam_filename = argv[1];
	string in_filename = argv[2];
	string out_filename = argv[3];
	
	
	if (argc == 5) {
		LIM = atoi(argv[4]);
	}
	
	cout << "Filtering using limit = " << LIM << endl;
	
	
	if (in_filename.substr(in_filename.length()-4).compare(".bed")!=0 ||
		out_filename.substr(out_filename.length()-4).compare(".bed")!=0) {
		cout << "ERROR: I'm sorry, this utility is only for bed files";
		exit(2);
	}

	in_file.open(sam_filename.c_str(), ios::in);
	
	if (!in_file.is_open()) {
		cout << "ERROR: I cannot open the SAM file, " << sam_filename << endl;
		exit(1);
	}
	

	
	
	cout << "Reading SAM file... " << endl;
	
	pair<map<string,vector<int> >::iterator,bool> good_pt_it;
	
	while (!in_file.eof()) {
		string curr_chr = "";
		string temp_curr_chr = "";
		int curr_left_pos;
		string curr_cigar;
		size_t tab_loc = 0;
		size_t tab_loc2 = 0;
		
		getline(in_file, line);
		trim2(line);
		
		if (line.length() == 0) {
			continue;
		}
		
		
		if (line[0] == '@') {
			continue;
		}
		
		tab_loc = line.find('\t'); 
		tab_loc = line.find('\t',tab_loc+1); 
		tab_loc2 = line.find('\t', tab_loc+1); 
		temp_curr_chr = line.substr(tab_loc+1, tab_loc2-tab_loc-1);
		
		if (curr_chr.length() == 0) {
			curr_chr = temp_curr_chr;
			good_pt_it = good_point_chr.insert(pair<string,vector<int> >(curr_chr,vector<int>()));
		}else if (temp_curr_chr.compare(curr_chr) != 0){
			
			sort(((good_pt_it.first)->second).begin(),((good_pt_it.first)->second).end());
			good_pt_it = good_point_chr.insert(pair<string,vector<int> >(curr_chr,vector<int>()));
		}

		
		tab_loc = tab_loc2;  
		tab_loc2 = line.find('\t', tab_loc+1); 
		curr_left_pos = atoi(line.substr(tab_loc+1, tab_loc2-tab_loc-1).c_str());
		
		tab_loc = tab_loc2;  
		tab_loc2 = line.find('\t', tab_loc+1); 
		tab_loc = tab_loc2;  
		tab_loc2 = line.find('\t', tab_loc+1); 
		
		curr_cigar = line.substr(tab_loc+1, tab_loc2-tab_loc-1);
		
		
		if (curr_cigar.find('N') == string::npos) {
			((good_pt_it.first)->second).push_back(curr_left_pos);
			((good_pt_it.first)->second).push_back(get_cigar_right_pt(curr_left_pos, curr_cigar));
		}
		
	}
	
	in_file.close();
	in_file.clear();
	
	cout << "Filtering... " << endl;
	
	in_file.open(in_filename.c_str(), ios::in);
	out_file.open(out_filename.c_str(), ios::out);
	
	while (!in_file.eof()) {
		getline(in_file, line);
		trim2(line);
		
		if(line.length() == 0){
			continue;
		}
		
		
		if (line.substr(0, 5).compare("track") == 0) {
			out_file << line << "\n";
			continue;
		}
		
		
		vector<string> line_list = split(line, '\t');
		string chr_name = line_list[0];
		vector<string> thickness_str = split(line_list[10],',');
		pair<int,int> thickness = pair<int,int>(atoi(thickness_str[0].c_str()),atoi(thickness_str[1].c_str()));
		int left_pt = atoi(line_list[1].c_str()) + thickness.first;
		int right_pt = atoi(line_list[2].c_str()) - thickness.second;
		
		int nNR_begin = (int)line_list[3].find('_') + 1;
		int nNR_end = (int)line_list[3].find(']');
		int nNR = atoi(line_list[3].substr(nNR_begin,nNR_end-nNR_begin).c_str());
		
		
		
		
		if (nNR > 1) {
			out_file << line << "\n";
			continue;
		}
		
		
		
		map<string,vector<int> >::iterator good_pt_list = good_point_chr.find(chr_name);
		
		
		
		vector<int>::iterator close_left_pt_it = (upper_bound((good_pt_list->second).begin(), (good_pt_list->second).end(), left_pt) - 1);
		vector<int>::iterator close_right_pt_it = (lower_bound((good_pt_list->second).begin(), (good_pt_list->second).end(), right_pt) + 1);
		
		
		
		bool good = false;
		
		if (close_left_pt_it != (good_pt_list->second).begin()) {
			close_left_pt_it--;
			
			

			
			
			if (left_pt - *close_left_pt_it < LIM) {
				good = true;
			}
		}
		
		if (close_right_pt_it != (good_pt_list->second).end()) {
			
			
			
			
			if (*close_right_pt_it - right_pt < LIM) {
				good = true;
			}
		}
		
		
		if (good) {
			out_file << line << "\n";
		}
		
	}
	
	
	in_file.close();
	in_file.clear();
	out_file.close();
	out_file.clear();
	
	return 0;
}

int get_cigar_right_pt(int left_pos, string &cigar)
{
	int result = left_pos-1;
	string temp_str = "";
	unsigned int index = 0;
	
	while (index < cigar.length()) {
		char curr_chr = cigar[index];
		if (isdigit(curr_chr)) {
			temp_str += curr_chr;
		}else {
			switch (curr_chr) {
				case 'S':
					
					break;
				case 'N':
					result += atoi(temp_str.c_str());
					break;
				case 'M':
					result += atoi(temp_str.c_str());
					break;
				default:
					cout << "ERROR: Malformed CIGAR string: " << cigar << endl;
					exit(2);
					break;
			}
			
			temp_str = "";
		}

		
		index++;
	}
	
	return result;
}

void print_usage_and_exit()
{
	cout << "Usage: ./neighborFilter good_hits.sam infile.bed outfile.bed [limit]" << endl;
	cout << "SAM file must be sorted by chromosome position" << endl;
	cout << "Filters out all junctions with non-redundant reads = 1 and no reads close to it, within the limit. " << endl;
	cout << "The default limit is 80nt" << endl;
	cout << "See also: " << WEBSITE << endl;
	exit(1);
}




