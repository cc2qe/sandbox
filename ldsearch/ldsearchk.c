#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int powi(int x, int y)
{
  int product = 1;
  int i;
  for (i = 0; i < y; ++i) {
    product *= x;
  }
  return product;
}

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

void get_expected(double **rates,
		  int *locus_indices,
		  int set_size,
		  int num_gt_combos,
                  int num_inf,
                  double *d_expected)
{
  int i, j;
  for (i = 0; i < num_gt_combos; ++i) {
    int itr = i;
    double exp_rate = 1;
    
    for (j = 0; j < set_size; ++j) {
      int p = powi(3,(set_size - j - 1));
      exp_rate *= rates[locus_indices[j]][itr / p];
      itr = itr % p;
    }
    d_expected[i] = exp_rate * num_inf;
  }
}

void get_observed(int **M,
		  int *locus_indices,
		  int set_size,
		  int num_gt_combos,
		  int num_samples,
		  int *multi_informative,
		  double *d_observed)
{
  int i, j;
  for (i = 0; i < num_gt_combos; ++i) 
    d_observed[i] = 0;

  int genotype;
  for (i = 0; i < num_samples; ++i) {
    // only assess samples that are informative at all k loci
    int checkInf = 1;
    for (j = 0; j < set_size; ++j) {
      checkInf *= (M[locus_indices[j]][i] + 1);
    }
    if (checkInf != 0) {
      *multi_informative += 1;

      // bitwise (tripwise?) representation of genotype
      genotype = 0;
      for (j = 0; j < set_size; ++j) {
	genotype += powi(3,j) * M[locus_indices[set_size-j-1]][i];
      }
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
      tenPower *= 10;
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
	  "                           (currently 3 no matter what you put, sucka)\n"
	  "  -d MIN_DISTANCE        minimum distance between loci\n"
	  "  -x MIN_CHI             minimum chi-squared sum to print\n"
	  "  -e MIN_EXP             conservative chi-square, cells only contribute\n"
	  "                           if the expected freq is greater than MIN_EXP\n"
	  "                           (default: 5)\n"
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
  double min_chi_sum = 0;
  int yates = 0;
  double min_exp = 5.0;
  int brief = 0;
  char *file_name;
  char *samples_file_name;

  int i;
  int c;
  opterr = 0;

  while ((c = getopt(argc, argv, "hd:s:l:k:x:e:b")) != -1) {
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
      set_size = atoi(optarg);
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
  char *chrArr[num_loci];
  int posArr[num_loci];
  char *geneArr[num_loci];
  int num_informative[num_loci];
  int *M[num_loci];
  int num_gt_combos = powi(3,set_size);
  
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
  char m_gts[num_gt_combos][set_size + 1];
  if (! brief) {
    for (j = 0; j < num_gt_combos; ++j) {
      char myGt[9];
      sprintf(myGt, "%08d", decToBase(j, 3));
      sprintf(m_gts[j], "%s", &myGt[8-set_size]);
    }
  }

  int k,l;
  double rates_1[3], rates_2[3], rates_3[3];
  double expected[num_gt_combos], observed[num_gt_combos], chi[num_gt_combos];
  double chi_sum;
  
  double *rates[num_loci];
  for (i = 0; i < num_loci; ++i) {
    double *rate = (double *) malloc(set_size * sizeof(double));
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

	// create an array of the [set_size] locus indices that will be interrogated
	int locus_indices[set_size];
	locus_indices[0] = i;
	locus_indices[1] = j;
	locus_indices[2] = k;
	
	
	// number of samples at are informative at all loci in k
	int num_multi_informative = 0;
	
	get_observed(M,
		     locus_indices,
		     set_size,
		     num_gt_combos,
		     num_samples,
		     &num_multi_informative,
		     observed);
      
        get_expected(rates,
		     locus_indices,
		     set_size,
		     num_gt_combos,
		     num_multi_informative,
                     expected);

	// calculate chi values for each cell and the chi_sum value for the trio
	chi_sum = 0;
	for (l = 0; l < num_gt_combos; ++l) {
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
	    for (l = 0; l < num_gt_combos; ++l) {
	      printf("%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%.0f|%.1f|%.1f\t%f\t%f\n",
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
