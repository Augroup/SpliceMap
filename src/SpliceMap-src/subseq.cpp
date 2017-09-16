

#include "subseq.h"

using namespace std;

int main (int argc, char * const argv[]) 
{
	
	
	
	string genome_filename;
	ifstream genome_file;
	good_t good_var;
	
	string chr_line = "";
	string chr_name;
	
	string line;
	string temp;
	
	if (argc != 4) {
		cout << "usage: ./subseq chr_file start end" << endl;
		cout << "Displays the chromosone nucleotides between locations 'start' and 'end' inclusive." << endl;
		exit(2);
	}
	
	genome_filename = argv[1];
	good_var.a = atoi(argv[2]);
	good_var.b = atoi(argv[3]);
	good_var.c = 0;
	good_var.d = 0;

	
	
	genome_file.open(genome_filename.c_str(), ios::in);
	
	if (!genome_file.is_open()) {
		cout << "ERROR: I'm sorry, the genome file cannot be opened. " << endl;
		exit(1);
	}
	
	getline(genome_file, chr_name);
	chr_name.erase(0,1);
	cout << "Chromosone Name: " << chr_name << endl;
	
	while (!genome_file.eof()) {
		getline(genome_file, line);
		
		trim2(line);
		
		if(line.length() == 0)
			continue;
		
		chr_line.append(line);
	}
	
	cout << "Chromosone Length: " << chr_line.length() << endl;

	cout << "Parameters: [" << good_var.a << "," << good_var.b << "]" << endl;
	
	genome_file.close();
	
	
	
	
	
	temp = chr_line.substr(good_var.a-1, good_var.b-good_var.a+1);
	
	if(!make_DNA_upper(temp)){
		cout << "ERROR: Unexpected character in genome:" << endl;
		cout << temp << endl;
		exit(2);
	}
	cout << temp << endl;
	compleseq(temp);
	cout << temp << endl;
	
	
	return 0;
	
}




