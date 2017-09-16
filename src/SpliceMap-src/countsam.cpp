

#include "countsam.h"

using namespace std;

int main (int argc, char * const argv[]) 
{
	string line;
	string temp;
	
	

	string infile_name;
	string reffile_name;
	
	ifstream infile;
	ifstream reffile;
	
	
	
	
	
	
	
	
	
	
	
	map<int,uint8_t> read_store; 
	
	
	
	
	
	
	
	map<string, vector<pair<int,int > > > read_dict;
	
	if (argc != 3) {
		print_usage_and_exit();
	}
	
	
	infile_name = argv[1];
	reffile_name = argv[2];
	
	
	
	
	cout << "Loading SAM file: " << infile_name << endl; 
	
	infile.open(infile_name.c_str(), ios::in);
	
	if (!infile.is_open()) {
		cout << "ERROR: I'm sorry, could not open the infile, " << infile_name << endl;
		exit(1);
	}
	

	
	while (!infile.eof()) {
		getline(infile, line);
		trim2(line);
		if (line.length() == 0) {
			continue;
		}
		
		
		if (line[0] == '@') {
			continue;
		}
		
		vector<string> line_list = split(line,'\t');
		
		size_t start_loc = line_list[0].rfind('[');
		size_t end_loc = line_list[0].rfind(']');
		int read_idx = (int)atoi(line_list[0].substr(start_loc+1,end_loc-start_loc-1).c_str());
		string chr_name = line_list[2];
		int left_pos = (int)atoi(line_list[3].c_str());
		string cigar = line_list[5];
		
		size_t n_loc = line.find("XS:");
		if (n_loc != string::npos) {  

			read_store[read_idx] |= 0x04; 
			
		}else {
			
			
			pair<map<string, vector<pair<int,int> > >::iterator, bool> out = read_dict.insert(pair<string, vector<pair<int,int> > >(chr_name,vector<pair<int,int> >()));
			
			(out.first)->second.push_back(pair<int,int>(left_pos,read_idx));
		}

	}
	
	infile.close();
	infile.clear();
	
	
	
	map<string, vector<pair<int,int> > >::iterator read_dict_it = read_dict.begin();
	while (read_dict_it != read_dict.end()) {
		
		sort(read_dict_it->second.begin(), read_dict_it->second.end());
		
		read_dict_it++;
	}
	
	cout << "Number of chromosomes loaded: " << read_dict.size() << endl;
	
	
	
	cout << "Reading the annotations file: " << reffile_name << endl;
	
	reffile.open(reffile_name.c_str(), ios::in);
	
	while (!reffile.eof()) {
		getline(reffile, line);
		trim2(line);
		
		if (line.length() == 0){
			continue;
		}
		
		
		
		
		
		
		
		vector<string> refline_list = split(line);
	
		
		
		
		vector<string> temp_list = split(refline_list[9], ',');
		list<int> exon_start_list;
		for (vector<string>::iterator it = temp_list.begin(); it!=temp_list.end(); it++) {
			exon_start_list.push_back(atoi(it->c_str()));
		}
		
		
		
		temp_list = split(refline_list[10], ',');
		list<int> exon_end_list;
		for (vector<string>::iterator it = temp_list.begin(); it!=temp_list.end(); it++) {
			exon_end_list.push_back(atoi(it->c_str()));
		}
		
		
		
		
		
		string line_chr = refline_list[2]; 
		

		list<int>::iterator start_it = exon_start_list.begin();
		list<int>::iterator end_it = exon_end_list.begin(); 
		
		map<string,vector<pair<int,int> > >::iterator chr_read_it = read_dict.find(line_chr);
		
		
		
		while (start_it != exon_start_list.end()) {
			
			vector<pair<int,int> >::iterator read_start_it 
			= lower_bound(chr_read_it->second.begin(), chr_read_it->second.end(), pair<int,int>(*start_it,INT32_MIN));
			read_start_it++;
			
			vector<pair<int,int> >::iterator read_end_it 
			= lower_bound(chr_read_it->second.begin(), chr_read_it->second.end(), pair<int,int>(*end_it,INT32_MAX));
			read_end_it++;
			
			while (read_start_it != read_end_it) {
				uint_fast32_t curr_read_idx = (int)abs(read_start_it->second);
				read_start_it->second = -curr_read_idx; 

				
				read_store[curr_read_idx] |= 0x02; 
				
				read_start_it++;
			}
			
			
			start_it++;
			
			if (start_it == exon_start_list.end()) {
				continue;
			}
			
			read_start_it 
			= lower_bound(chr_read_it->second.begin(), chr_read_it->second.end(), pair<int,int>((*end_it)+1,INT32_MIN));
			read_start_it++;
			
			read_end_it 
			= lower_bound(chr_read_it->second.begin(), chr_read_it->second.end(), pair<int,int>((*start_it)-1,INT32_MAX));
			read_end_it++;
			
			while (read_start_it != read_end_it) {
				uint_fast32_t curr_read_idx = (int)abs(read_start_it->second);
				read_start_it->second = -curr_read_idx; 
				
				read_store[curr_read_idx] |= 0x08; 
				
				read_start_it++;
			}
			
			end_it++;
			
		}
		
		
		
	}
	
	reffile.close();

	
	
	map<string, vector<pair<int,int > > >::iterator full_read_dict_it = read_dict.begin();
	
	while (full_read_dict_it != read_dict.end()) {
		
		vector<pair<int,int > >::iterator full_read_it = full_read_dict_it->second.begin();
		
		while (full_read_it != full_read_dict_it->second.end()) {
			
			int_fast32_t curr_read_idx = full_read_it->second;
			
			if (curr_read_idx > 0){ 

				read_store[curr_read_idx] |= 0x01; 
			}
			
			full_read_it++;
		}
		
		full_read_dict_it++;
	}
	
	
	
	
	uint_fast32_t total_reads_mapped = (uint_fast32_t)read_store.size();
	uint_fast32_t exon_mapped = 0;
	uint_fast32_t junction_mapped = 0;
	uint_fast32_t intron_mapped = 0;
	uint_fast32_t inter_mapped = 0;
	uint_fast32_t exon_intron_mapped = 0;
	uint_fast32_t exon_junction_mapped = 0;
	uint_fast32_t intron_junction_mapped = 0;
	
	map<int,uint8_t>::iterator read_store_it = read_store.begin();
	
	while (read_store_it != read_store.end()) {
		
		if ((read_store_it->second) & 0x02) {
			exon_mapped++;
		}
		
		if ((read_store_it->second) & 0x04) {
			junction_mapped++;
		}
		
		if ((read_store_it->second) & 0x08) {
			intron_mapped++;
		}
		
		if ((read_store_it->second) & 0x01) {
			inter_mapped++;
		}
		
		if ((read_store_it->second) & 0x02 && (read_store_it->second) & 0x08) {
			exon_intron_mapped++;
		}
		
		if ((read_store_it->second) & 0x02 && (read_store_it->second) & 0x04) {
			exon_junction_mapped++;
		}
		
		if ((read_store_it->second) & 0x08 && (read_store_it->second) & 0x04) {
			intron_junction_mapped++;
		}
		
		read_store_it++;
	}
	
	
	cout << "Total reads mapped:        " << total_reads_mapped << endl;
	
	cout << "Total mapped to exons:     " << exon_mapped;
	output_percentage(cout,total_reads_mapped,exon_mapped);
	cout << endl;
	
	cout << "Total mapped to junctions: " << junction_mapped;
	output_percentage(cout,total_reads_mapped,junction_mapped);
	cout << endl;
	
	cout << "Total mapped to introns:   " << intron_mapped;
	output_percentage(cout,total_reads_mapped,intron_mapped);
	cout << endl;
	
	cout << "Total mapped to exons and junctions:   " << exon_junction_mapped;
	output_percentage(cout,total_reads_mapped,exon_junction_mapped);
	cout << endl;
	
	cout << "Total mapped to exons and introns:     " << exon_intron_mapped;
	output_percentage(cout,total_reads_mapped,exon_intron_mapped);
	cout << endl;
	
	cout << "Total mapped to introns and junctions: " << intron_junction_mapped;
	output_percentage(cout,total_reads_mapped,intron_junction_mapped);
	cout << endl;

	cout << "Total mapped to intergenic/unannotated: " << inter_mapped;
	output_percentage(cout,total_reads_mapped,inter_mapped);
	cout << endl;
	
	return 0;
}

void output_percentage(ostream &out, uint_fast32_t total, uint_fast32_t num)
{
	out << " (" << setprecision(2) << (num*100.0/(double)total) << resetiosflags (ios_base::fixed) << "%)";
}

void print_usage_and_exit()
{
	cout << "usage: countsam infile.sam refFlat.txt" << endl;
	cout << "Output is the number of reads mapped to exons/junctions" << endl;
	exit(2);
}



