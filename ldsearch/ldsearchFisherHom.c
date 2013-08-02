#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

double log_choose(double n, double k)
{
  int d;
  double r = 0;

  // swap for efficiency if k is more than half of n
  if (k * 2 > n) {
    k = n - k;
  }

  for (d = 1; d <= k; ++d) {
    r += log10(n--);
    r -= log10(d);
  }
  return r;
}

double log_fisher(double *obs)
{
  double log_p_val;
  double p_val;
  log_p_val =
    log_choose(obs[0] + obs[1], obs[0]) +
    log_choose(obs[2] + obs[3], obs[2]) -
    log_choose(obs[0] + obs[1] + obs[2] + obs[3], obs[0] + obs[2]);

  return log_p_val;
}

void get_expected(double *rates_1,
                  double *rates_2,
                  int num_inf,
                  double *d_expected)
{
  int i,j;
  for (i= 0; i < 3; ++i)
    for (j= 0; j < 3; ++j) {
      int genotype;
      genotype = 3 * i + j;
      if (genotype == 0) {
	d_expected[0] = rates_1[i] *
	  rates_2[j] *
	  num_inf;
      }
      else if (genotype == 2) {
        d_expected[1] = rates_1[i] *
          rates_2[j] *
          num_inf;
      }
      else if (genotype == 6) {
        d_expected[2] = rates_1[i] *
          rates_2[j] *
          num_inf;
      }
      else if (genotype == 8) {
        d_expected[3] = rates_1[i] *
          rates_2[j] *
          num_inf;
      }
    }
}

void get_observed(int *locus_1,
                   int *locus_2,
                   int num_samples,
                   double *d_observed,
		   int *multi_informative)
{
  int i;
  for (i = 0; i < 4; ++i) 
    d_observed[i] = 0;

  int genotype;
  for (i = 0; i < num_samples; ++i) {
    // only assess samples that are informative at all k loci
    if (locus_1[i] >= 0 && locus_2[i] >= 0) {
      // multi_informative is the number of samples that are informative
      // at ALL k loci
      *multi_informative += 1;

      // bitwise representation of genotype
      genotype = 3 * locus_1[i] +
	1 * locus_2[i];

      if (genotype == 0) {
	d_observed[0] += 1;
      }
      else if (genotype == 2) {
	d_observed[1] += 1;
      }
      else if (genotype == 6) {
	d_observed[2] += 1;
      }
      else if (genotype == 8) {
	d_observed[3] += 1;
      }
    }
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
    else if (loci_i == 2)
      ++num_hom_alt;
    // else it's -1 (uninformative)
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
	  "  -h                     show this help and exit\n"
	  "  -s NUM_SAMPLES         number of samples in file\n"
	  "  -l NUM_LOCI            number of loci in file\n"
	  "  -k SET_SIZE            number of loci in each set\n"
	  "                           (currently 2 no matter what you put, sucka)\n"
	  "  -d MIN_DISTANCE        minimum distance between loci\n"
	  "  -p MAX_LOG_P           maximum log p-value to print (default: 0)\n"
	  "  -i START_I             for parallelizing: the starting i value to iterate on\n"
	  "  -z BLOCK_SIZE          for parallelizing: the number of loci to iterate through\n"
	  "  -b                     brief but faster output\n"
	  "\n"
	  );
  return 1;
}

int main (int argc, char **argv)
{
  int min_distance = 0;
  int num_samples;
  int num_loci;
  int set_size;
  double max_log_p = 0;
  int brief = 0;
  char *file_name;
  char *samples_file_name;
  int start_i = 0;
  int block_size = -1;

  int i;
  int c;
  opterr = 0;

  while ((c = getopt(argc, argv, "hd:s:l:k:p:e:i:z:b")) != -1) {
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
    case 'k':
      set_size = 2;
      break;
    case 'p':
      max_log_p = atoi(optarg);
      break;
    case 'i':
      start_i = atoi(optarg);
      break;
    case 'z':
      block_size = atoi(optarg);
      break;
    case 'b':
      brief = 1;
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
  
  // if the user didn't specify a block size then just make it the entire file
  if (block_size == -1)
    block_size = num_loci;

  // the maximum chars per line of file to read in. (should be more than double the number of samples)
  int max_line = 10000;
  char *sep = "\t";
  
  FILE *f = fopen(file_name, "rt");
  FILE *s = fopen(samples_file_name, "rt");

  // make an array from the samples file
  char line[max_line];
  char **sample = (char **) malloc(3 * sizeof(char*));
  int j = 0;
  while (fgets(line, max_line, s) != NULL) {
    char *sample_name = strtok(line, sep);
    char *sample_ethn = strtok(NULL, sep);
    char *sample_subpop = strtok(NULL, sep);

    // copy the string into an array that will be retained
    // through the end of the run
    //sample[0] = strdup(sample_name);
    //sample[1] = strdup(sample_ethn);
    //sample[2] = strdup(sample_subpop);

    ++j;
  }
  fclose(s);

  // for (j = 0; j < num_samples; ++j) {
  //   printf("%d\n", j);   
  // }
  
  // array of arrays containing genotype info for each sample
  // at each locus
  char **chrArr = (char **) malloc(num_loci * sizeof(char*));
  int *posArr = (int *) malloc(num_loci * sizeof(int));
  char **geneArr = (char **) malloc(num_loci * sizeof(char*));
  char ** rsIdArr = (char **) malloc(num_loci * sizeof(char*));
  int *num_informative = (int *) malloc(num_loci * sizeof(int));
  int **M = (int **) malloc(num_loci * sizeof(int*));
  
  // chr1  69510 OR4F5 0.65  0.64  0.32  0.87  0.69  2
  j = 0;
  while (fgets(line, max_line, f) != NULL) {
    char *chr = strtok(line, sep);
    int pos = atoi(strtok(NULL, sep));
    char *rsId = strtok(NULL, sep);
    char *gene = strtok(NULL, sep);
    int inf = 0;
    
    char *rate_1 = strtok(NULL, sep);
    char *rate_2 = strtok(NULL, sep);
    char *rate_3 = strtok(NULL, sep);
    char *rate_4 = strtok(NULL, sep);
    char *rate_5 = strtok(NULL, sep);

    chrArr[j] = strdup(chr);
    posArr[j] = pos;
    geneArr[j] = strdup(gene);
    rsIdArr[j] = strdup(rsId);

    int *locus_gts = (int *) malloc(num_samples * sizeof(int));
    
    i = 0;
    char *tok = strtok(NULL,sep);
    while ((tok != NULL) && (i < num_samples)) {
      locus_gts[i] = atoi(tok);
      tok = strtok(NULL,sep);

      if (locus_gts[i] != -1) {
	++inf;
      }
      ++i;
    }
    
    // store the genotypes into the locus x gt matrix
    M[j] = locus_gts;
    // store the number of informative samples at locus
    num_informative[j] = inf;
    ++j;
  }
  fclose(f);

  // generate array of genotypes (eg 00, 01, 20) so we don't have to
  // call the dec to base function for every locus
  int m_gts[4];
  if (! brief) {
    m_gts[0] = 0;
    m_gts[1] = 2;
    m_gts[2] = 20;
    m_gts[3] = 22;
  }
  
  int k,l;
  double observed[4];

  // make sure that it won't blow through all the loci available
  int end_i = start_i + block_size;
  if (end_i > num_loci)
    end_i = num_loci;
  
  for (i = start_i; i < end_i; ++i) {
    for (j = i + 1; j < num_loci; ++j) {
      // only calc chi-square if loci are each separated by
      // minimum distance
      if (strcmp(chrArr[i],chrArr[j]) == 0 && abs(posArr[i] - posArr[j]) < min_distance) {
	continue;
      }	
      
      // number of samples at are informative at all loci in k
      int num_multi_informative = 0;
      
      get_observed(M[i],
		   M[j],
		   num_samples,
		   observed,
		   &num_multi_informative);

      double log_p = log_fisher(observed);
      if (log_p <= max_log_p) {
	if (brief) {
	  printf("%d\t%d\t%f\n",
		 i,j,
		 log_p);
	}
	
	else {
	  double *rate_1 = (double *) malloc(3 * sizeof(double));
	  double *rate_2 = (double *) malloc(3 * sizeof(double));
	  get_rates(M[i], num_samples, rate_1);
	  get_rates(M[j], num_samples, rate_2);

	  double expected[4];	  
	  get_expected(rate_1,
		       rate_2,
		       num_multi_informative,
		       expected);

	  for (l = 0; l < 4; ++l) {
	    printf("%s\t%d\t%s\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%s\t%.3f\t%.3f\t%.3f\t%02d\t%.0f|%.1f|%.1f\t%f\t%f\n",
		   chrArr[i],posArr[i],rsIdArr[i],geneArr[i],rate_1[0],rate_1[1],rate_1[2],
		   chrArr[j],posArr[j],rsIdArr[j],geneArr[j],rate_2[0],rate_2[1],rate_2[2],
		   m_gts[l],
		   observed[l],expected[l],observed[l]-expected[l],
		   log_p, pow(10, log_p));
	  }
	  
	}
      }
    }
  }
  
  for (j = 0; j < num_loci; ++j) {
    free(chrArr[j]);
    free(geneArr[j]);
    free(rsIdArr[j]);
    free(M[j]);
  }
  free(posArr);
  free(sample);

  return 0;
}
