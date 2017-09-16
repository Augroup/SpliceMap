#include "bar.h"
#include "chr_region.h"

int main(int argc, char *argv[]){
	bar mybar;

	if (argc == 2) {
		if (!convert_from_text_to_bar(argv[1], string(argv[1])+".bar")) {
			return 1;
		}
	} else if (argc == 3) {
		if (!mybar.read_from_file(argv[1])) {
			cout << "error read file" << endl;
			return 1;
		}

		if (!mybar.write_to_text_file(argv[2])) {
			cout << "error write file" << endl;
			return 1;
		}
	} else if (argc == 4) {
		string chr;
		int start, end;
		if (!parse_region(argv[3], chr, start, end)) {
			cout << "bad region" << endl;
			return 1;
		}
		vector<vector<bar_column> > data;
		if (!mybar.read_from_file_region(argv[1], chr, start, end, data)) {
			cout << "error read file" << endl;
			return 1;
		}
		FILE *file = fopen(argv[2], "wt");
		if (!file) {
			cout << "error open file" << argv[2] << endl;
			return 1;
		}
		if (!mybar.dump_data(file, chr, data)) {
			cout << "error write file" << endl;
			return 1;
		}
		fclose(file);
	} else {
		cout << "barloader <txt file>" << endl;
		cout << "barloader <bar file> <txt file>" << endl;
		cout << "barloader <bar file> <txt file> <chr region>" << endl;
		return 1;
	}

	return 0;
}
