#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

double get_X(double observed,
             double expected)
{
  double numer;
  if (expected != 0) {
    numer = (observed - expected);
    return (numer*numer)/expected;
  }
  else return 0;
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

void parseArgs ()
{
  return;
}

int usage()
{
  fprintf(stderr,
	  "usage: ldsearch [options] <file> <samples>\n\n"
	  "authors: Ryan Layer and Colby Chiang\n"
	  "description: finds loci in linkage disequilibrium based on\n"
	  "  genotype frequencies in a set of individuals\n"
	  "\n"
	  "positional arguments:\n"
	  "  file                   tab-delimited input file of genotypes\n"
	  "  samples                tab-delimited file of sample names, superpopulations,\n"
	  "                           and subpopulations\n"
	  "\n"
	  "optional arguments:\n"
	  "  -h, --help             show this help and exit\n"
	  "  -s, --num_samples      number of samples in file\n"
	  "  -l, --num_loci         number of loci in file\n"
	  "  -d, --min_distance     minimum distance between loci\n"
	  "  -x, --min_chi_sum      minimum chi-squared sum to print\n"
	  "\n"
	  );
  return 1;
}

int main (int argc, char **argv)
{
  int min_distance = 0;
  int num_samples;
  int num_loci;
  double min_chi_sum = 0;
  char *file_name;
  char *samples_file_name;

  int i;
  int c;
  opterr = 0;

  while ((c = getopt(argc, argv, "hd:s:l:x:")) != -1) {
    switch (c) {
    case 'h':
      return usage();
    case 'd':
      min_distance = atoi(optarg);
      break;
    case 's':
      num_samples = atoi(optarg);
      break;
    case 'l':
      num_loci = atoi(optarg);
      break;
    case 'x':
      min_chi_sum = atoi(optarg);
      break;
    case '?':
      if (optopt == 'c')
	fprintf(stderr, "Option -%c requires an argument\n", optopt);
      else if (isprint(optopt))
	fprintf(stderr, "Unknown option '-%c'\n", optopt);
      else
	fprintf(stderr, "Unknown option character '\\x%x'\n", optopt);
      return 1;
    default:
      abort();
    }
  }

  // parse the positional arguments
  samples_file_name = argv[argc-1];
  file_name = argv[argc-2];
  
  if (argc < 2) {
    return usage();
  }
  
  int max_line = 5000; // the maximum length of a line to read in
  char *sep = "\t";
  
  FILE *f = fopen(file_name, "rt");
  FILE *s = fopen(samples_file_name, "rt");

  // make an array from the samples file
  char line[max_line];
  char **sample = (char **) malloc(3 * sizeof(char*));
  int j = 0;
  while (fgets(line, max_line, s) != NULL) {
    char *sample_name = strtok(line, sep);
    char *sample_subpop = strtok(NULL, sep);
    char *sample_superpop = strtok(NULL, sep);

    // copy the string into an array that will be retained
    // through the end of the run
    sample[0] = strdup(sample_name);
    sample[1] = strdup(sample_subpop);
    sample[2] = strdup(sample_superpop);

    ++j;
  }
  fclose(s);

  // array of arrays containing genotype info for each sample
  // at each locus
  char *chrArr[num_loci];
  int posArr[num_loci];
  char *geneArr[num_loci];
  int *M[num_loci];
  
  // chr1  69510 OR4F5 0.65  0.64  0.32  0.87  0.69  2
  j = 0;
  while (fgets(line, max_line, f) != NULL) {
    char *chr = strtok(line, sep);
    int pos = atoi(strtok(NULL, sep));
    char *gene = strtok(NULL, sep);
    
    char *rate_1 = strtok(NULL, sep);
    char *rate_2 = strtok(NULL, sep);
    char *rate_3 = strtok(NULL, sep);
    char *rate_4 = strtok(NULL, sep);
    char *rate_5 = strtok(NULL, sep);

    chrArr[j] = strdup(chr);
    posArr[j] = pos;
    geneArr[j] = strdup(gene);

    int *locus_gts = (int *) malloc(num_samples * sizeof(int));
    
    i = 0;
    char *tok = strtok(NULL,sep);
    while ((tok != NULL) && (i < num_samples)) {
      locus_gts[i] = atoi(tok);
      tok = strtok(NULL,sep);
      ++i;
    }
    
    M[j] = locus_gts;
    ++j;
  }
  fclose(f);

  int k,l;
  double rates_1[3], rates_2[3], rates_3[3];
  double expected[27], observed[27], chi[27];
  double chi_sum;
  
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
	if (strcmp(chrArr[i],chrArr[j]) == 0 && abs(posArr[i] - posArr[j]) < min_distance ||
	    strcmp(chrArr[i],chrArr[k]) == 0 && abs(posArr[i] - posArr[k]) < min_distance ||
	    strcmp(chrArr[j],chrArr[k]) == 0 && abs(posArr[j] - posArr[k]) < min_distance) {
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

	// calculate chi values for each cell and the chi_sum value for the trio
	chi_sum = 0;
	for (l = 0; l < 27; ++l) {
	  chi[l] = get_X(observed[l],expected[l]);
	  chi_sum += chi[l];
	}

	if (chi_sum >= min_chi_sum) {
	  for (l = 0; l < 27; ++l) {
	    printf("%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%.0f|%.1f|%.1f\t%f\t%f\n",
	  	   chrArr[i],posArr[i],geneArr[i],rates[i][0],rates[i][1],rates[i][2],
		   chrArr[j],posArr[j],geneArr[j],rates[j][0],rates[j][1],rates[j][2],
		   chrArr[k],posArr[k],geneArr[k],rates[k][0],rates[k][1],rates[k][2],
		   observed[l],expected[l],observed[l]-expected[l],
		   chi[l],chi_sum);
	  }
	}
      }
    }
  }
  return 0;
}
