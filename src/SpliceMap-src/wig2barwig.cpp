

#include "wig2barwig.h"


using namespace std;


int main (int argc, char * const argv[]) 
{
	
	
	string line;
	string temp;
	
	string infile_name;
	string outfile_name;
	
	ifstream infile;
	ofstream outfile;
	
	if (argc != 3) {
		cout << "usage: ./wig2barwig infile.wig outfile.barwig" << endl;
		cout << "Coverts the wig file to a barwig file for use in cisgenome browser" << endl;
		exit(2);
	}
	
	infile_name = argv[1];
	outfile_name = argv[2];
	
	
	infile.open(infile_name.c_str(), ios::in);
	
	if (!infile.is_open()) {
		cout << "ERROR: I'm sorry, cannot open the infile, " << infile_name << endl;
		exit(1);
	}
	
	bool track_line = false;
	getline(infile, line); 
	
	if (line.find("track") == string::npos) {
		cout << "Waring: The infile is possibly malformed " << infile_name << endl;
	}else {
		track_line = true;
	}
	
	cout << "Loading and processing file: " << infile_name << endl;

	outfile.open(outfile_name.c_str(), ios::out);
	
	
	while (!infile.eof()) {
		
		
		if (track_line) {
			getline(infile, line);
		}else {
			track_line = true;
		}

		trim2(line);
		if (line.length() == 0) {
			continue;
		}
		
		vector<string> line_list = split(line, '\t');
		int start = atoi(line_list[1].c_str());
		int end = atoi(line_list[2].c_str());
		float val = (float)atof(line_list[3].c_str());
		
		
		for(int i = start; i<=end;i++){
			
			outfile << line_list[0] << "\t" << i << "\t" << val << "\n";
		}

	}
	
	infile.close();
	infile.clear();
	
	
	
	
	
	
	
	
	outfile.close();
	outfile.clear();
	
	return 0;
}



