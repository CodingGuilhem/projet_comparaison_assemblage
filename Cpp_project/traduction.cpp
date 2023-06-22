#include <iostream>
#include <cstdlib>

using namespace std;


int 
main()
{

    string arg = "gkampi -k 100 -P -o fatsq.fasta --dump=ascii -f -i  test-1000_reads.fastq.gz";
    int result = system(&arg[0]);

    if (result == 0){

    }
    else { cerr<< "Error during the file generation"<< endl;}
    
}

