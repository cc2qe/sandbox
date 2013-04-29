#include <iostream>
#include <fstream>

#include "kstring.h"
#include "sam.h"

using namespace std;

char compBase (char base) {
	switch ( base ) {
		case 'A':
			return 'T';
			break;
		case 'C':
			return 'G';
			break;
		case 'G':
			return 'C';
			break;
		case 'T':
			return 'A';
			break;
		case 'N':
			return 'N';
			break;
                case 'a':
                        return 't';
                        break;
                case 'c':
                        return 'g';
                        break;
                case 'g':
                        return 'c';
                        break;
                case 't':
                        return 'a';
                        break;
		default:
			return 'N';
			break;
	}
}

int min (int a, int b) {
	if (a < b) { return a; }
	else { return b; }
}

int main( int argc, char *argv[] )
{
	if ( argc != 3 ) /* argc should be 2 for correct execution */
	{
		cout << "\nProgram: unmappedToFastq\n\
Version: 0.0.1\n\
Contact: Colby Chiang <cchiang3@gmail.com>\n\
Info: Finds unaligned reads in a paired BAM file,\n\
  clips the unaligned end(s) to a given length,\n\
  and outputs the pair in fastq. (BAM file must be name-sorted)" << endl;
		cout << "\n  usage: " << argv[0] << " [clipTo] [bamFile]\n" << endl;
	}

	else
	{
	
		bam1_t *bamrec = bam_init1();
		samfile_t *samfile;
		if ( (samfile = samopen(argv[2], "rb", NULL)) == 0) {
			cout << "error: could not open bam file" << endl;
			return 1;
		}
		
		// open output stream
		samfile_t *fpBothMapped;
		fpBothMapped = samopen("bothMapped.bam", "wb", samfile->header);
		
		ofstream fq1;
		fq1.open("1_ex.fq");
		
		ofstream fq2;
		fq2.open("2_ex.fq");
		
		int flag;
		int maxBases = atoi(argv[1]);
		
		while (samread(samfile, bamrec) > 0) {
	
			flag = bamrec->core.flag;
//			cout << flag << endl;
			
			bool secondRead;
			if (flag % 256 < 128) {
				secondRead = false;
			}
			else {
				secondRead = true;
			}
			
			bool revStrand;
			if (flag % 32 < 16) {
				revStrand = false;
			}
			else {
				revStrand = true;
			}
			
			bool isAligned;
			if (flag % 8 < 4) {
				isAligned = true;
			}
			else {
				isAligned = false;
			}
			
			bool mateAligned;
			if (flag % 16 < 8) {
				mateAligned = true;
			}
			else {
				mateAligned = false;
			}
			
			bool isProper;
			if (flag % 4 < 2) {
				isProper = false;
			}
			else {
				isProper = true;
			}
			
			if (isAligned && mateAligned) {
				samwrite(fpBothMapped, bamrec);
			}
			
			else {
				if (!secondRead)
				{
					fq1 << "@" << bam1_qname(bamrec) << endl;
					
					int clipTo;
					if ( ! isAligned ) {
						clipTo = min(bamrec->core.l_qseq, maxBases);
					}
					else { clipTo = bamrec->core.l_qseq; }
					
					if (!revStrand) {
											
						int i;
						for (i = 0; i < clipTo; ++i) {
							fq1 << bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)] ;
						}
						fq1 << "\n+" << endl;
				
						for (i = 0; i < clipTo; ++i) {
							fq1 << (char) (bam1_qual(bamrec)[i] + 33) ;
						}
						fq1 << endl;
					}
					
					else {
						int i;
						for (i = clipTo - 1; i >= 0 ; --i) {
							fq1 << compBase(bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)]) ;
						}
						fq1 << "\n+" << endl;
						
						for (i = clipTo - 1; i >= 0; --i) {
							fq1 << (char) (bam1_qual(bamrec)[i] + 33) ;
						}
						fq1 << endl;
					}
				}
				
				else if (secondRead)
				{
					fq2 << "@" << bam1_qname(bamrec) << endl;
	
					int clipTo;
					if ( ! isAligned ) {
						clipTo = min(bamrec->core.l_qseq, maxBases);
					}
					else { clipTo = bamrec->core.l_qseq; }
					
					if (!revStrand) {
												
						int i;
						for (i = 0; i < clipTo; ++i) {
							fq2 << bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)] ;
						}
						fq2 << "\n+" << endl;
				
						for (i = 0; i < clipTo; ++i) {
							fq2 << (char) (bam1_qual(bamrec)[i] + 33) ;
						}
						fq2 << endl;
					}
					
					else {
						int i;
						for (i = clipTo - 1; i >= 0 ; --i) {
							fq2 << compBase(bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)]) ;
						}
						fq2 << "\n+" << endl;
						
						for (i = clipTo - 1; i >= 0; --i) {
							fq2 << (char) (bam1_qual(bamrec)[i] + 33) ;
						}
						fq2 << endl;
					}
				}
			}
		}
		samclose(fpBothMapped);
		
		fq1.close();
		fq2.close();
	}
	return 0;
}

