

#include "randomJunctionFilter.h"



int main (int argc, char * const argv[]) 
{
	srand((unsigned)time(0)); 
	
	string line;
	ifstream in_file;
	ofstream out_file;
	
	int percentage;
	
	if (argc != 4) {
		print_usage_and_exit();
	}
	
	
	string in_filename = argv[1];
	string out_filename = argv[2];
	
	
	if (in_filename.substr(in_filename.length()-4).compare(".bed")!=0 ||
		out_filename.substr(out_filename.length()-4).compare(".bed")!=0) {
		cout << "ERROR: I'm sorry, this utility is only for bed files";
		exit(2);
	}
	
	percentage = atoi(argv[3]);
	
	if (percentage < 0 || percentage > 99) {
		cout << "ERROR: Please enter a percentage between 0 and 99 inclusive";
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
			int lucky_number = 1 + int(99*(rand()/(RAND_MAX + 1.0)));
			
			if(lucky_number <= percentage){
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
	cout << "Usage: ./randomJunctionFilter infile.bed outfile.bed percentage" << endl;
	cout << "Keeps a random percentage of junctions with non-redundant reads = 1" << endl;
	cout << "See also: " << WEBSITE << endl;
	exit(1);
}





