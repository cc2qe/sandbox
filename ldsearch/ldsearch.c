#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

double get_X(double observed,
             double expected,
	     double min_exp)
{
  double numer;
  if (expected > min_exp && expected != 0) {
    numer = (observed - expected);
    return (numer*numer)/expected;
  }
  else return 0;
}

void get_expected(double *rates_1,
                  double *rates_2,
                  double *rates_3,
                  int num_inf,
                  double *d_expected)
{
  int i,j,k;
  for (i= 0; i < 3; ++i)
    for (j= 0; j < 3; ++j)
      for (k= 0; k < 3; ++k) {
	d_expected[9*i + 3*j + 1*k] = rates_1[i] * 
	  rates_2[j] *
	  rates_3[k] * 
	  num_inf;
      }
}

void get_observed(int *locus_1,
                   int *locus_2,
                   int *locus_3,
                   int num_samples,
                   double *d_observed,
		   int *multi_informative)
{
  int i;
  for (i = 0; i < 27; ++i) 
    d_observed[i] = 0;

  int genotype;
  for (i = 0; i < num_samples; ++i) {
    // only assess samples that are informative at all k loci
    if (locus_1[i] >= 0 && locus_2[i] >= 0 && locus_3[i] >= 0) {
      // multi_informative is the number of samples that are informative
      // at ALL k loci
      *multi_informative += 1;

      // bitwise (tripwise?) representation of genotype
      genotype = 9 * locus_1[i] +
	3 * locus_2[i] +
	1 * locus_3[i];
      d_observed[genotype] += 1;
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

int decToBase(int x,
	      int base)
{
  int i = 0;
  int digit;
  int x_b = 0;

  while(x != 0) {
    digit = x % base;
    x = x / base;

    int tenPower = 1;
    int j;
    for (j = 0; j < i; ++j) {
      tenPower = tenPower * 10;
    }
    x_b += (digit * tenPower);

    ++i;
  }

  return x_b;
}

void parseArgs ()
{
  return;
}

int usage()
{
  fprintf(stderr,
	  "usage: ldsearch [options] <file> <samples>\n\n"
	  "ldsearch\n"
	  "authors: Ryan Layer and Colby Chiang\n"
	  "description: reports loci in linkage disequilibrium based on\n"
	  "  genotype frequencies in a set of individuals\n"
	  "\n"
	  "positional arguments:\n"
	  "  file                   tab-delimited input file of genotypes\n"
	  "  samples                tab-delimited file of sample names, superpopulations,\n"
	  "                           and subpopulations\n"
	  "\n"
	  "optional arguments:\n"
	  "  -h                     show this help and exit\n"
<<<<<<< Updated upstream
	  "  -s NUM_SAMPLES         number of samples in file\n"
	  "  -l NUM_LOCI            number of loci in file\n"
	  "  -k SET_SIZE            number of loci in each set\n"
	  "                           (currently 3 no matter what you put, sucka)\n"
	  "  -d MIN_DISTANCE        minimum distance between loci\n"
	  "  -x MIN_CHI             minimum chi-squared sum to print\n"
	  "  -e MIN_EXP             conservative chi-square, cells only contribute\n"
	  "                           if the expected freq is greater than MIN_EXP\n"
	  "                           (default: 5)\n"
	  "  -b                     brief but faster output\n"
=======
	  "  -s NUM_SAMP            number of samples in file\n"
	  "  -l NUM_LOCI            number of loci in file\n"
	  "  -k SET_SIZE            number of loci in set of interactions\n"
	  "                           (currently 3 no matter what you put, sucka)\n"
	  "  -d MIN_DIST            minimum distance between loci\n"
	  "  -x MIN_CHI             minimum chi-squared sum to print\n"
>>>>>>> Stashed changes
	  "\n"
	  );
  return 1;
}

int main (int argc, char **argv)
{
  int set_size = 3;
  int min_distance = 0;
  int num_samples;
  int num_loci;
  int set_size;
  double min_chi_sum = 0;
  int yates = 0;
  double min_exp = 5.0;
  int brief = 0;
  char *file_name;
  char *samples_file_name;

  int i;
  int c;
  opterr = 0;

<<<<<<< Updated upstream
  while ((c = getopt(argc, argv, "hd:s:l:k:x:e:b")) != -1) {
=======
  while ((c = getopt(argc, argv, "hk:d:s:l:x:")) != -1) {
>>>>>>> Stashed changes
    switch (c) {
    case 'h':
      return usage();
    case 'k':
      // set_size = atoi(optarg);
      break;
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
      set_size = 3;
      break;
    case 'x':
      min_chi_sum = atoi(optarg);
      break;
    case 'y':
      yates = 1;
      break;
    case 'e':
      min_exp = atof(optarg);
      break;
    case 'b':
      brief = 1;
      break;
    case '?':
      if (optopt == 'd' || optopt == 's' || optopt == 'l' || optopt == 'x')
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
  
  // the maximum chars per line of file to read in. (should be more than double the number of samples)
  int max_line = 10000;
  char *sep = "\t";
  
  FILE *f = fopen(file_name, "rt");
  FILE *s = fopen(samples_file_name, "rt");

  // make an array from the samples file
  char line[max_line];
  char *sample_subpops[num_loci];
  int j = 0;
  while (fgets(line, max_line, s) != NULL) {
    char *sample_name = strtok(line, sep);
<<<<<<< Updated upstream
    char *sample_ethn = strtok(NULL, sep);
=======
    char *sample_eth = strtok(NULL, sep);
>>>>>>> Stashed changes
    char *sample_subpop = strtok(NULL, sep);

    // copy the string into an array that will be retained
    // through the end of the run
<<<<<<< Updated upstream
    //sample[0] = strdup(sample_name);
    //sample[1] = strdup(sample_ethn);
    //sample[2] = strdup(sample_subpop);
=======
    sample_subpops[j] = strndup(sample_subpop,strlen(sample_subpop)-1);
>>>>>>> Stashed changes

    ++j;
  }
  fclose(s);

<<<<<<< Updated upstream
  // for (j = 0; j < num_samples; ++j) {
  //   printf("%d\n", j);   
  // }
  
=======
  // test some things out about the samples files
  // get metrics on the subpopulations

  // make list of the available populations
  char **pops = (char **) malloc(5 * sizeof(char*));
  pops[0] = "AFR";
  pops[1] = "AMR";
  pops[2] = "ASN";
  pops[3] = "EUR";
  pops[4] = "ALL";
  

  //  int *afr_indices;
  //  int *amr_indices;
  //  int *asn_indices;
  //  int *eur_indices;
  //  printf("%d\n", countSamples(samples,"AFR",num_samples));

  //  for (j = 0; j < num_samples; ++j) {
  // printf("%s\n", sample_subpops[j]);
  // }

>>>>>>> Stashed changes
  // array of arrays containing genotype info for each sample
  // at each locus
  char *chrArr[num_loci];
  int posArr[num_loci];
  char *geneArr[num_loci];
  int num_informative[num_loci];
  int *M[num_loci];
  
  // chr1  69510 OR4F5 0.65  0.64  0.32  0.87  0.69  2
  j = 0;
  while (fgets(line, max_line, f) != NULL) {
    char *chr = strtok(line, sep);
    int pos = atoi(strtok(NULL, sep));
    char *gene = strtok(NULL, sep);
    int inf = 0;
    
    char *rate_1 = strtok(NULL, sep);
    char *rate_2 = strtok(NULL, sep);
    char *rate_3 = strtok(NULL, sep);
    char *rate_4 = strtok(NULL, sep);
    char *rate_5 = strtok(NULL, sep);

    // copy strings to array that will be retained
    // through the end of the run
    chrArr[j] = strdup(chr);
    posArr[j] = pos;
    geneArr[j] = strdup(gene);

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

  // generate array of genotypes (eg 000, 012, 202) so we don't have to
  // call the dec to base function for every locus
  int m_gts[27];
  if (! brief) {
    for (j = 0; j < 27; ++j) {
      m_gts[j] = decToBase(j, 3);
    }
  }

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
<<<<<<< Updated upstream
	
	// number of samples at are informative at all loci in k
	int num_multi_informative = 0;
	
=======

	// store the expected number in expected
	get_expected(rates[i],
		     rates[j],
		     rates[k],
		     num_samples,
		     expected);

	// store the observed number in observed
>>>>>>> Stashed changes
	get_observed(M[i],
		     M[j],
		     M[k],
		     num_samples,
		     observed,
		     &num_multi_informative);
      
        get_expected(rates[i],
                     rates[j],
                     rates[k],
		     num_multi_informative,
                     expected);

	// calculate chi values for each cell and the chi_sum value for the trio
	chi_sum = 0;
	for (l = 0; l < 27; ++l) {
	  chi[l] = get_X(observed[l],expected[l], min_exp);
	  chi_sum += chi[l];
	}

	if (chi_sum >= min_chi_sum) {
	  if (brief) {
	    printf("%d\t%d\t%d\t%f\n",
		   i,j,k,
		   chi_sum);
	  }

	  else {
	    for (l = 0; l < 27; ++l) {
	      printf("%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%03d\t%.0f|%.1f|%.1f\t%f\t%f\n",
		     chrArr[i],posArr[i],geneArr[i],rates[i][0],rates[i][1],rates[i][2],
		     chrArr[j],posArr[j],geneArr[j],rates[j][0],rates[j][1],rates[j][2],
		     chrArr[k],posArr[k],geneArr[k],rates[k][0],rates[k][1],rates[k][2],
		     m_gts[l],
		     observed[l],expected[l],observed[l]-expected[l],
		     chi[l],chi_sum);
	    }
	  }
	}
      }
    }
  }
  
  for (j = 0; j < num_loci; ++j) {
    free(chrArr[j]);
    free(geneArr[j]);
  }

  return 0;
}
