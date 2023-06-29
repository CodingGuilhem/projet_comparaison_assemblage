#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <gkampi.h>

using namespace std;

struct anchors {
    string sequence;
    int position;
};

struct anchors extract_anchors(string fasta_name)
{
    ifstream fasta_file(fasta_name);
    if (fasta_file.is_open())
    {
        string line;
        while (getline(fasta_file, line))
        {
            if (line[0] == '>')
            {
                string position = line;
                
            }
        }
        getline(fasta_file, line);
        
        fasta_file.close();
        
    }
    else
    {
        cerr << "Unable to open file" << endl;
        exit(1);
    }
}
int 
main()
{
    
    string arg = "gkampi -k 100 -P -o fatsq.fasta --dump=ascii -f -i  test-1000_reads.fastq.gz";
    int result = system(&arg[0]);

    if (result == 0){

    }
    else { cerr<< "Error during the file generation"<< endl;}
    
}

