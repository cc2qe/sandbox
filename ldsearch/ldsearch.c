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
                  int num_genotype,
                  double *d_expected)
{
  int i,j,k;
  for (i= 0; i < 3; ++i)
    for (j= 0; j < 3; ++j)
      for (k= 0; k < 3; ++k) {
	d_expected[9*i + 3*j + 1*k] = rates_1[i] * 
	  rates_2[j] *
	  rates_3[k] * 
	  num_genotype;
      }
}

void *get_observed(int *loci_1,
                   int *loci_2,
                   int *loci_3,
                   int num_genotype,
                   double *d_observed)
{
  int i;
  for (i = 0; i < 27; ++i) 
    d_observed[i]=0;
  
  int geno_type;
  for (i = 0; i < num_genotype; ++i) {
    geno_type = 9 * loci_1[i] +
      3 * loci_2[i] +
      1 * loci_3[i];
    d_observed[geno_type] += 1;
  }
}

void get_rates(int *loci, 
               int num_genotype,
               double *rates)
{
  int num_homo_ref = 0,
    num_hetro_ref = 0,
    num_homo_alt = 0;
  
  int i,loci_i;
  for (i = 0; i < num_genotype; ++i) {
    loci_i = loci[i];
    if (loci_i == 0) 
      ++num_homo_ref;
    else if (loci_i == 1)
      ++num_hetro_ref;
    else
      ++num_homo_alt;
  }

  int total = num_homo_ref + num_hetro_ref + num_homo_alt;
  
  rates[0] = ((double)num_homo_ref) / ((double)total);
  rates[1] = ((double)num_hetro_ref) / ((double)total);
  rates[2] = ((double)num_homo_alt) / ((double)total);
}
 
int main (int argc, char **argv)
{
  if (argc != 4) {
    printf("usage %s: <file> <num genotypes> <num loci>\n", argv[0]);
    return 1;
  }
  
  char *file_name = argv[1];
  int num_genotype = atoi(argv[2]);
  int num_loci = atoi(argv[3]);
  int max_line = 5000;
  char *sep = "\t";
  
  FILE *f = fopen(file_name, "rt");
  
  char line[max_line];
  
  int *M[num_loci];
  
  //chr1  69510   OR4F5   0.65    0.64    0.32    0.87    0.69    2
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
    
    int *loci = (int *) malloc(num_genotype * sizeof(int));
    
    int i = 0;
    
    char *tok = strtok(NULL,sep);
    while ((tok != NULL) && (i < num_genotype)) {
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
  int total = 0;
  
  double *rates[num_loci];
  for (i = 0; i < num_loci; ++i) {
    double *rate = (double *) malloc(3 * sizeof(double));
    get_rates(M[i], num_genotype, rate);
    rates[i] = rate;
  }
  
  for (i = 0; i < num_loci; ++i) {
    for (j = i + 1; j < num_loci; ++j) {
      for (k = j + 1; k < num_loci; ++k) {
	get_expected(rates[i],
		     rates[j],
		     rates[k],
		     num_genotype,
		     expected);
	
	get_observed(M[i],
		     M[j],
		     M[k],
		     num_genotype,
		     observed);
	
	x = get_X(observed,expected);
	printf("%f\n",x);
	++total;
	if (total==1000000)
	  return 0;
      }
    }
  }
  return 0;
}
