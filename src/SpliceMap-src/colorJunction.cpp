

#include "colorJunction.h"

using namespace std;

int main (int argc, char * const argv[]) {

	const int new_depth = 50; 
	
	ifstream infile;
	ifstream newfile;
	ofstream outfile;
	
	string infile_name;
	string newfile_name = "";
	string outfile_name;
	
	string temp;
	string line;
	
	chrdict_t new_jundict;
	
	
	
	if (argc >= 2 && argc <= 3) {
		infile_name = argv[1];
		for (int i = 2; i<argc; i++) {
			newfile_name = argv[i];
		}
	}else {
		print_usage_and_exit();
	}
	
	
	if (infile_name.substr(infile_name.rfind('.')).compare(".bed") != 0) {
		cout << "ERROR: I'm sorry, only .bed files are supported for input file, " << infile_name << endl;
		print_usage_and_exit();
	}
	if (newfile_name.length() > 0 && newfile_name.substr(newfile_name.rfind('.')).compare(".bed") != 0) {
		cout << "ERROR: I'm sorry, only .bed files are supported for newjun file, " << newfile_name << endl;
		print_usage_and_exit();
	}
	
	outfile_name = infile_name.substr(0, infile_name.rfind('.')) + "_color.bed";
	cout << "Output filename: " << outfile_name << endl;
	
	
	if(newfile_name.length() > 0){
		newfile.open(newfile_name.c_str(), ios::in);
		
		if(!newfile.is_open()){
			cout << "ERROR: I'm sorry, cannot open the newjun file, " << newfile_name << endl;
			exit(1);
		}
		
		load_bed_ref_file(newfile, new_jundict, 1); 
		
		newfile.close();
		newfile.clear();
	}
	
	infile.open(infile_name.c_str(), ios::in);
	if(!infile.is_open()){
		cout << "ERROR: I'm sorry, cannot open the input file, " << infile_name << endl;
		exit(1);
	}
	
	outfile.open(outfile_name.c_str(), ios::out);
	
	while (!infile.eof()) {
		getline(infile, line);
		trim2(line);
		
		if(line.length() == 0){
			continue;
		}
		
		
		if (line.substr(0, 5).compare("track") == 0) {
			outfile << line << " itemRgb=\"On\"\n";
			continue;
		}
		
		string color;
		
		vector<string> line_list = split(line, '\t');
		
		
		int A = atoi(line_list[1].c_str()) + atoi(line_list[4].c_str()) - new_depth;
        int B = atoi(line_list[2].c_str()) - atoi(line_list[4].c_str()) + new_depth;
		
		string str_show;
		
		int start = (int)line_list[3].find('_') + 1; 
		int len = (int)line_list[3].find(']') - start;
		int indep_read = atoi(line_list[3].substr(start,len).c_str());
		
		str_show = line_list[3]; 
		
		
		
		
		
		color = color_read(indep_read);
		
		if(newfile_name.length() > 0){
			
			string line_chr = line_list[0];
			chrdict_t::iterator primarychr_it = new_jundict.find(line_chr);
			
			
			if (primarychr_it != new_jundict.end()) {
				vector<string> thickness = split(line_list[10], ',');
				int leftpos = atoi(line_list[1].c_str()) + atoi(thickness[0].c_str());
				int rightpos = atoi(line_list[2].c_str()) - atoi(thickness[1].c_str());
	
				beginpos_t::iterator begin_it = (primarychr_it->second).find(leftpos);
				
				if (begin_it != (primarychr_it->second).end()) {
					endpos_t::iterator finish_it = (begin_it->second).find(rightpos);
					
					if (finish_it != (begin_it->second).end()) {
						
						
						
						if (indep_read == 1) {
							color = "255,192,192";
						}else {
							color = "255,0,0";
						}
					}
				}				
			}
		}
		
		
		outfile << line_list[0] << '\t' << A << '\t' << B << '\t' << str_show << '\t' << new_depth << '\t'
		<< line_list[5] << '\t' << A << '\t' << B << '\t' << color << "\t2\t" << new_depth << ','
		<< new_depth << "\t0," << (B-A-new_depth) << "\n";
		
	}
	
	outfile << flush;
	
	infile.close();
	outfile.close();
	
	
	return 0;
}

void print_usage_and_exit()
{
	cout << "Usage: ./colorJunction infile.bed [newjun.bed]" << endl;
	cout << "  infile.bed -- This file will be coloured and outputed as infile_color.bed" << endl;
	cout << "  newjun.bed -- This is file of new junctions from \'findNovelJunctions\', they will be highlighted with color" << endl;
	cout << "The junctions are also drawn so that they are more visible" << endl;
	cout << "See also: " << WEBSITE << endl;
	exit(1);
}

string color_read(int indep_read)
{
	string color = "0,255,0";
	
    if (indep_read == 1){
        color = "192,192,192";
	}else if (indep_read > 1 && indep_read <= 5){
        color = "96,96,96";
	}else if (indep_read > 5){
        color = "0,0,0";
	}
	
	return color;
}

