#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double get_X(double *observed,
             double *expected)
{
  int i;
  double sum = 0;
  double numer;
  for (i = 0; i < 27; ++i) {
    if (expected[i] != 0) {
      numer = (observed[i] - expected[i]);
      sum += (numer*numer)/expected[i];
    }
  }
  return sum;
}
void get_expected(double *rates_1,
                  double *rates_2,
                  double *rates_3,
                  int num_samples,
                  double *d_expected)
{
  int i,j,k;
  for (i= 0; i < 3; ++i)
    for (j= 0; j < 3; ++j)
      for (k= 0; k < 3; ++k) {
	d_expected[9*i + 3*j + 1*k] = rates_1[i] * 
	    rates_2[j] *
	    rates_3[k] * 
	  num_samples;
      }
}

void *get_observed(int *loci_1,
                   int *loci_2,
                   int *loci_3,
                   int num_samples,
                   double *d_observed)
{
  int i;
  for (i = 0; i < 27; ++i) 
    d_observed[i]=0;
  
  int genotype;
  for (i = 0; i < num_samples; ++i) {
    genotype = 9 * loci_1[i] +
      3 * loci_2[i] +
      1 * loci_3[i];
    d_observed[genotype] += 1;
  }
}

void get_rates(int *loci, 
               int num_samples,
               double *rates)
{
  int num_hom_ref = 0,
    num_het = 0,
    num_hom_alt = 0;
  
  int i,loci_i;
  for (i = 0; i < num_samples; ++i) {
    loci_i = loci[i];
    if (loci_i == 0) 
      ++num_hom_ref;
    else if (loci_i == 1)
      ++num_het;
    else
      ++num_hom_alt;
  }

  int total = num_hom_ref + num_het + num_hom_alt;
  
  rates[0] = ((double)num_hom_ref) / ((double)total);
  rates[1] = ((double)num_het) / ((double)total);
  rates[2] = ((double)num_hom_alt) / ((double)total);
}
 
int main (int argc, char **argv)
{
  if (argc != 6) {
    printf("usage %s: <genotypes file> <samples file> <num samples> <num loci> <min distance>\n", argv[0]);
    return 1;
  }
  
  char *file_name = argv[1];
  char *samples_file_name = argv[2];
  int num_samples = atoi(argv[3]);
  int num_loci = atoi(argv[4]);
  int min_distance = atoi(argv[5]);
  int max_line = 5000;
  char *sep = "\t";
  
  FILE *f = fopen(file_name, "rt");
  
  char line[max_line];
  
  char *chrArr[num_loci];
  int locArr[num_loci];
  char *geneArr[num_loci];


  // array of arrays containing genotype info for each sample
  // at each locus
  int *M[num_loci];
  
  // chr1  69510 OR4F5 0.65  0.64  0.32  0.87  0.69  2
  int j = 0;
  while (fgets(line, max_line, f) != NULL) {
    char *chr = strtok(line, sep);
    int loc = atoi(strtok(NULL, sep));
    char *gene = strtok(NULL, sep);
    
    char *rate_1 = strtok(NULL, sep);
    char *rate_2 = strtok(NULL, sep);
    char *rate_3 = strtok(NULL, sep);
    char *rate_4 = strtok(NULL, sep);
    char *rate_5 = strtok(NULL, sep);

    chrArr[j] = strdup(chr);
    locArr[j] = loc;
    geneArr[j] = strdup(gene);

    int *loci = (int *) malloc(num_samples * sizeof(int));
    
    int i = 0;

    char *tok = strtok(NULL,sep);
    while ((tok != NULL) && (i < num_samples)) {
      loci[i] = atoi(tok);
      tok = strtok(NULL,sep);
      ++i;
    }
    
    M[j] = loci;
    ++j;
  }
  fclose(f);
  
  int i,k;
  double rates_1[3], rates_2[3], rates_3[3];
  double expected[27], observed[27];
  double x;
  
  double *rates[num_loci];
  for (i = 0; i < num_loci; ++i) {
    double *rate = (double *) malloc(3 * sizeof(double));
    get_rates(M[i], num_samples, rate);
    rates[i] = rate;
  }
  
  for (i = 0; i < num_loci; ++i) {
    for (j = i + 1; j < num_loci; ++j) {
      for (k = j + 1; k < num_loci; ++k) {
	// only calc chi-square if loci are each separated by
	// minimum distance
	if (strcmp(chrArr[i],chrArr[j]) == 0 && abs(locArr[i] - locArr[j]) < min_distance ||
	    strcmp(chrArr[i],chrArr[k]) == 0 && abs(locArr[i] - locArr[k]) < min_distance ||
	    strcmp(chrArr[j],chrArr[k]) == 0 && abs(locArr[j] - locArr[k]) < min_distance) {
	  continue;
	}
	
	get_expected(rates[i],
		     rates[j],
		     rates[k],
		     num_samples,
		     expected);
	
	get_observed(M[i],
		     M[j],
		     M[k],
		     num_samples,
		     observed);
	
	x = get_X(observed,expected);

	int l = 0;
	//for (l = 0; l < 27; ++l) {
	printf("%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%.0f|%.1f|%.1f\t%f\n",
	       chrArr[i],locArr[i],geneArr[i],rates[i][0],rates[i][1],rates[i][2],
	       chrArr[j],locArr[j],geneArr[j],rates[j][0],rates[j][1],rates[j][2],
	       chrArr[k],locArr[k],geneArr[k],rates[k][0],rates[k][1],rates[k][2],
	       observed[l],expected[l],observed[0]-expected[l],
	       x);
	//}
      }
    }
  }
  return 0;
}
