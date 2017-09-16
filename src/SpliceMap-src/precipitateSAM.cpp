

#include "precipitateSAM.h"


using namespace std;

int main (int argc, char * const argv[]) 
{
	string temp;
	string line;
	
	ifstream SAM_file;
	string SAM_file_name;
	
	string outfile_base;
	
	map<uint_fast32_t,pair<float,float> > curr_coverage; 
	
	if (argc != 3) {
		print_usage_and_exit();
	}
	
	SAM_file_name = argv[1];
	outfile_base = argv[2];
	
	
	cout << "Reading SAM file and outputting coverage... " << endl;
	
	string coverage_all_name = outfile_base+"_all.wig";
	string coverage_up_name = outfile_base+"_up.wig";
	string coverage_down_name = outfile_base+"_down.wig";
	
	ofstream coverage_all;
	ofstream coverage_up;
	ofstream coverage_down;
	
	coverage_all.open(coverage_all_name.c_str(), ios::out);
	if (!coverage_all.is_open()) {
		cout << "ERROR: could not write to coverage file, " << coverage_all_name << endl;
		exit(1);
	}
	
	coverage_up.open(coverage_up_name.c_str(), ios::out);
	if (!coverage_up.is_open()) {
		cout << "ERROR: could not write to coverage file, " << coverage_up_name << endl;
		exit(1);
	}
	
	coverage_down.open(coverage_down_name.c_str(), ios::out);
	if (!coverage_down.is_open()) {
		cout << "ERROR: could not write to coverage file, " << coverage_down_name << endl;
		exit(1);
	}
	
	
	coverage_all << "track\ttype=bedGraph\tname=\"SpliceMap all coverage\"" << '\n';
	coverage_up << "track\ttype=bedGraph\tname=\"SpliceMap unique coverage\"" << '\n';
	coverage_down << "track\ttype=bedGraph\tname=\"SpliceMap non-unique coverage\"" << '\n';
		
	
	out_group_t all_g = {0,0,0.0f};   
	out_group_t up_g = {0,0,0.0f};
	out_group_t down_g = {0,0,0.0f};
	
	
	
	SAM_file.open(SAM_file_name.c_str(), ios::in);
	
	if (!SAM_file.is_open()) {
		cout << "ERROR: could not open SAM file, " <<  SAM_file_name << endl;
		exit(1);
	}
	
	string curr_chr = "";
	
	while (!SAM_file.eof()) {
		getline(SAM_file, line);
		
		trim2(line);
		
		if (line.length() == 0 || line[0] == '@') { 
			continue;
		}
		
		vector<string> line_list = split(line, '\t');  
		
		if (curr_chr == "") {
			curr_chr = line_list[2];
		}else if(line_list[2] != curr_chr){
			
			
			
			write_coverage(coverage_up,up_g,0,0,curr_chr);
			
			
			write_coverage(coverage_down,down_g,0,0,curr_chr);
			
			
			write_coverage(coverage_all,all_g,0,0,curr_chr);
			
			
			
			all_g.start = 0;
			all_g.end = 0;
			all_g.val = 0.0f;
			up_g.start = 0;
			up_g.end = 0;
			up_g.val = 0.0f;
			down_g.start = 0;
			down_g.end = 0;
			down_g.val = 0.0f;
			
			
			curr_chr = line_list[2];
		}

		
		int num_mapped = 1; 
		
		for (int i = 11; i<(int)line_list.size(); i++) {
			
			if (line_list[i].length() > 5 && line_list[i][0] == 'N' && line_list[i][1] == 'H' && line_list[i][2] == ':' && line_list[i][3] == 'i' && line_list[i][4] == ':') {

				num_mapped = atoi(line_list[i].substr(5).c_str());
				
				i = (int)line_list.size();
				
			}
		}
		
		float up_cover   = 0;
		float down_cover = 0;
		
		if (num_mapped == 1) {
			up_cover = 1.0f;
		}else if(num_mapped == 0){
			cout << "ERROR in SAM file at line:" << endl;
			cout << line << endl;
			exit(2);
		}else {
			down_cover = 1.0f/num_mapped;
		}
		
		
		
		
		string cigar = line_list[5];
		uint_fast32_t curr_start = atoi(line_list[3].c_str());  
		uint_fast32_t start_loc = curr_start;
		size_t read_ptr = 0;
		
		while (read_ptr != string::npos) {
			size_t temp_read_ptr = cigar.find_first_of("SMN", read_ptr);

			if(temp_read_ptr != string::npos){
				if(cigar[temp_read_ptr] == 'M'){
					uint_fast32_t range = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
					
					
					
					
					for(uint_fast32_t i = start_loc;i<start_loc+range;i++){
						pair<map<uint_fast32_t,pair<float,float> >::iterator, bool> out 
						= curr_coverage.insert(pair<uint_fast32_t,pair<float,float> >(i,pair<float,float>(up_cover,down_cover)));
						
						if (!out.second) {
							(out.first)->second.first += up_cover;
							(out.first)->second.second += down_cover;
						}
					}
					
					start_loc = start_loc + range;
				}else if(cigar[temp_read_ptr] == 'N'){ 
					uint_fast32_t range = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
					

					start_loc = start_loc + range;
				} else if (read_ptr == 0) { 
					
				}
				
				read_ptr = temp_read_ptr + 1;
				
				if (read_ptr == cigar.length()) {
					break;
				}
			}else {
				read_ptr = temp_read_ptr;
			}
			
			
		}
		
		
		
		if (curr_coverage.begin()->first < curr_start) {
			
			
			map<uint_fast32_t,pair<float,float> >::iterator curr_coverage_it = curr_coverage.begin();
			while (curr_coverage_it->first < curr_start) {
				uint_fast32_t curr_loc = curr_coverage_it->first;
				float curr_up = curr_coverage_it->second.first;
				float curr_down = -(curr_coverage_it->second.second);
				float curr_all = curr_up-curr_down; 
				
				
				
				if(curr_up != 0){
					write_coverage(coverage_up,up_g,curr_up,curr_loc,curr_chr);
				}
				
				
				if(curr_down != 0){
					write_coverage(coverage_down,down_g,curr_down,curr_loc,curr_chr);
				}
				
				
				write_coverage(coverage_all,all_g,curr_all,curr_loc,curr_chr);
				
				
				curr_coverage_it++;
			}
			
			curr_coverage.erase(curr_coverage.begin(),curr_coverage_it);
		}
		
		
		
	}
	
	
	write_coverage(coverage_up,up_g,0,0,curr_chr);
	write_coverage(coverage_down,down_g,0,0,curr_chr);
	write_coverage(coverage_all,all_g,0,0,curr_chr);
	
	
	SAM_file.close();
	SAM_file.clear();
	
	coverage_all.close();
	coverage_all.clear();
	coverage_up.close();
	coverage_up.clear();
	coverage_down.close();
	coverage_down.clear();
	
	return 0;
}


void write_coverage(ofstream &coverage, out_group_t &g, float curr_val,uint_fast32_t curr_loc, string chr)
{
	if (curr_val != g.val) {
		if(g.start > 0){
			coverage << chr << '\t' << g.start << '\t' << g.end << '\t' << g.val << '\n';
		}
			
		
		g.start = curr_loc;
		g.end   = curr_loc;
		g.val   = curr_val;
	}else {
		if (g.start == 0) {
			g.start = curr_loc;
			g.end   = curr_loc;
			g.val   = curr_val;
		}else {
			if (curr_loc == (g.end + 1)) {
				g.end = curr_loc;
			}else { 
				
				if(g.start > 0){
					coverage << chr << '\t' << g.start << '\t' << g.end << '\t' << g.val << '\n';
				}
				
				
				g.start = curr_loc;
				g.end   = curr_loc;
				g.val   = curr_val;
			}
			
		}
		
	}
}

void print_usage_and_exit()
{
	cout << "Usage: ./precipitateSAM good_hits.sam junction" << endl;
	cout << "For internal SpliceMap use only!" << endl;
	cout << "SAM file must be sorted by coordinate" << endl;
	cout << "it will output junction.bed, junction_all.wig, junction_up.wig, junction_down.wig" << endl;
	exit(1);
}







