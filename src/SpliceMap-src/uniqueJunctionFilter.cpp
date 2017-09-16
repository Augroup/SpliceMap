

#include "uniqueJunctionFilter.h"

int main (int argc, char * const argv[]) 
{
	string line;
	ifstream in_file;
	ofstream out_file;
	
	if (argc != 3) {
		print_usage_and_exit();
	}
	
	
	string in_filename = argv[1];
	string out_filename = argv[2];
	
	
	if (in_filename.substr(in_filename.length()-4).compare(".bed")!=0 ||
		out_filename.substr(out_filename.length()-4).compare(".bed")!=0) {
		cout << "ERROR: I'm sorry, this utility is only for bed files";
		exit(2);
	}
	 
	
	
	
	in_file.open(in_filename.c_str(), ios::in);
	
	
	if (!in_file.is_open()) {
		cout << "ERROR: I cannot open the input file, " << in_filename << endl;
		exit(1);
	}
	
	out_file.open(out_filename.c_str(), ios::out);
	
	
	getline(in_file, line);
	out_file << line << "\n";
	
	
	while (!in_file.eof()) {
		getline(in_file, line);
		int nNR_start = (int)line.find('_') + 1;
		int nNR_end = (int)line.find(']',nNR_start);
		
		int nNR = atoi(line.substr(nNR_start, nNR_end - nNR_start).c_str());
		
		
		if (nNR == 1){
			int UR_start = nNR_end + 2;
			int UR_end = (int)line.find('/',UR_start);
			
			int UR = atoi(line.substr(UR_start, UR_end - UR_start).c_str());
			
			if(UR > 0){
				out_file << line << '\n';
			}
			
		}else{ 
			out_file << line << '\n';
		}

	}
	
	
	out_file << flush;
	
	in_file.close();
	out_file.close();
	
	
	return 0;
}



void print_usage_and_exit()
{
	cout << "Usage: ./uniqueJunctionFilter infile.bed outfile.bed" << endl;
	cout << "Filters out all junctions with no unique reads and non-redundant reads = 1" << endl;
	cout << "See also: " << WEBSITE << endl;
	exit(1);
}



