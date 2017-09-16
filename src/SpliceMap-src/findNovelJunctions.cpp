

#include "findNovelJunctions.h"


int main (int argc, char * const argv[]) {
	
	string primary_filename = "refFlat.txt";
	string bed_filename = "junction.bed";
	
	string temp = "";
	string line = "";
	
	chrdict_t primary_jundict;
	
	ifstream primary_file;
	int primary_file_type = 0; 
	ifstream bed_file;
	
	ofstream newjun_file;
	
	struct timeval tv;
	struct timeval start_tv;
	
	gettimeofday(&start_tv, NULL);
	
	if (argc >= 3 && argc <=3) {
		primary_filename = argv[1];
		bed_filename = argv[2];
		
		
		
	} else {
		cout << "usage: ./findNovelJunction refFlat.txt junction.bed" << endl;
		cout << " or    ./findNovelJunction refFlat.bed junction.bed" << endl;
		cout << "Finds the novel juctions that exist in junction.bed, but not in the refFlat file" << endl;
		
		
		exit(1);
	}
	
	if(primary_filename.substr(primary_filename.length()-4).compare(".bed") == 0){
		primary_file_type = 1;
	}else if(primary_filename.substr(primary_filename.length()-4).compare(".txt") == 0){
		primary_file_type = 0;
	}else {
		cout << "ERROR: I'm sorry, only ref (.txt) or bed (.bed) files are accepted." << endl;
		exit(2);
	}
	
	primary_file.open(primary_filename.c_str(),ios::in);
	cout << "Loading primary file: " << primary_filename << endl;
	
	if(!primary_file.is_open()){
		cerr << "ERROR: I'm sorry, could not open the primary file, " << primary_filename << endl;
		exit(1);
	}
	
	load_bed_ref_file(primary_file, primary_jundict, primary_file_type);
	
	primary_file.close();
	primary_file.clear();

	
	
	int total_primary_junctions = 0;
	int total_chromosomes = 0;
	
	chrdict_t::iterator chrdict_it = primary_jundict.begin();
	
	while (chrdict_it != primary_jundict.end()) {
		
		
		
		beginpos_t::iterator begin_it = chrdict_it->second.begin();
		
		while (begin_it != chrdict_it->second.end()) {
			endpos_t::iterator  end_it = begin_it->second.begin();
			
			while (end_it != begin_it->second.end()) {
				
				total_primary_junctions++;
				end_it++;
			}
			
			begin_it++;
		}
		
		chrdict_it++;
		total_chromosomes++;
	}
	
	cout << "Total chromosomes loaded: " << total_chromosomes << endl;
	cout << "Total junctions loaded: " << total_primary_junctions << endl;
	
	
	
	bed_file.open(bed_filename.c_str(), ios::in);
	
	cout << "Loading bed file: " << bed_filename << endl;
	
	if(!bed_file.is_open()){
		cerr << "ERROR: I'm sorry, could not open the bed file, " << bed_filename << endl;
		exit(1);
	}
	
	temp = bed_filename+".new.bed";
	newjun_file.open(temp.c_str(),ios::out);
	
	
	
	
	while (!bed_file.eof()) {
		getline(bed_file, line);
		trim2(line);
		
		int indicator = 0;
		
		if (line.length() == 0){
			continue;
		}
		
		
		if (line.substr(0, 5).compare("track") == 0) {
			newjun_file << line << endl;
			
		}else {
			vector<string> bedline_list = split(line);
			
			string line_chr = bedline_list[0];
			
			chrdict_t::iterator primarychr_it = primary_jundict.find(line_chr);
			
			if (primarychr_it != primary_jundict.end()) {
				vector<string> thickness = split(bedline_list[10], ',');
				int leftpos = atoi(bedline_list[1].c_str()) + atoi(thickness[0].c_str());
				int rightpos = atoi(bedline_list[2].c_str()) - atoi(thickness[1].c_str());
				
				
				
				beginpos_t::iterator begin_it = (primarychr_it->second).find(leftpos);
				
				if (begin_it != (primarychr_it->second).end()) {
					endpos_t::iterator finish_it = (begin_it->second).find(rightpos);
					 
					if (finish_it != (begin_it->second).end()) {
						indicator = 1;
					}else {
						newjun_file << line << endl;
					}

					
				}else {
					newjun_file << line << endl;
				}
				
			}else {
				newjun_file << line << endl;
			}
			
			
			
		}
		
	}
	
	
	newjun_file.close();
	newjun_file.clear();
	bed_file.close();
	bed_file.clear();
	
	
	
	gettimeofday(&tv, NULL);
	cout << "______________________" << endl;
	cout<<"Total Novel Junction Execution Time: "<<diffclock(start_tv,tv)<<" s."<<endl;
	
    return 0;
}






