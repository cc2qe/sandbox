#include "bam_select.hpp"

using namespace std;

int main(int argc, char **argv) {
  bam1_t *bamrec = bam_init1();
  samfile_t *samfile;

  samfile = samopen(argv[1], "rb", NULL);

  while (samread(samfile, bamrec) > 0) {
    if (testFlag(bamrec->core.flag, BAM_FPROPER_PAIR)) {
	cout << bam1_qname(bamrec) << endl;
      }
  }
 
}
