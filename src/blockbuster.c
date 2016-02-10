#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>


/* PROTOTYPES */
struct  read;
void    read_bed_file(char*);
void    assignReadsToBlocks(struct read*);
void    writeSuperGaussian(struct read*, double*, int);
void    writeBlocks(struct read*, double bThreshold);
int     assignReads(struct read*, int, int, int);
double  gaussian(double, double, double);
int     getHighestPeak(double*, int);
int     getRest(struct read*);
double  stddev(double*, double*, int);
void    sortArray(double*, int, int);
int     double_cmp(const void *a, const void *b);
void    writeHeader(char*);
void    printUsage(void);
void    freeList(struct read*);
char    *append(char*, char);

/* GLOBAL VARIABLES */
double pi = 3.14159265;
int    clusterStart = -1; 
int    clusterEnd = -1;
double clusterHeight = 0;
int    readCount = 0;
int    tagCount = 0;
char   *clusterChrom = "x";
char   *clusterStrand = "x";
int    clusterCounter = 0;

/* USER VARIABLES */
double minblockheight = 1;
double sizescale = 0.5;
int    merge = 0;
int    distance = 30;
int    minClusterHeight = 50;
int    print = 1;
double tagFilter = 0;
char   *sep = "\t";
int    type = 1;
char   *chr = NULL;
char   *blockheight = "abs";






/* MAIN FUNCTION */
int main(int argc, char *argv[])
{  
	// CHECK OPTIONAL ARGUMENTS
	if((argc) == 1){printUsage(); exit(0);}
	
	int i; int check = -1;
	for (i = 1; i < argc; i++)  
	{	
        if (strcmp(argv[i], "-scale") == 0)  
        {sizescale = atof(argv[i+1]); check = i+1;}
		
		if (strcmp(argv[i], "-minBlockHeight") == 0)
		{minblockheight = atof(argv[i+1]); check = i+1;}
		
		if (strcmp(argv[i], "-merge") == 0)
		{merge = atoi(argv[i+1]); check = i+1;}
		
		if (strcmp(argv[i], "-distance") == 0)
		{distance = atoi(argv[i+1]); check = i+1;}
		
		if (strcmp(argv[i], "-minClusterHeight") == 0)
		{minClusterHeight = atoi(argv[i+1]); check = i+1;}
		
		if (strcmp(argv[i], "-tagFilter") == 0)
		{tagFilter = atoi(argv[i+1]); check = i+1;}
		
		if (strcmp(argv[i], "-print") == 0)
		{print = atoi(argv[i+1]); check = i+1;}
		
		if (strcmp(argv[i], "-format") == 0)
		{type = atoi(argv[i+1]); check = i+1;}

		if (strcmp(argv[i], "-chr") == 0)
		{chr = argv[i+1]; check = i+1;}
		
		if (strcmp(argv[i], "-blockHeight") == 0)
		{blockheight = argv[i+1]; check = i+1;}

		if (strcmp(argv[i], "-h") == 0)
		{printUsage(); exit(0);}
	}
	if((argc-1) == check){printUsage(); exit(0);}
	
	read_bed_file(argv[argc-1]);
	return(0);
}





/*****************
 * STRUCTURES
 ****************/

struct read
{
	char *chrom;
	int start;
	int end;
	int block;
	double height;
	char *id;
	char *strand;
	struct read *next;
};


/***************** 
 *   FUNCTIONS
 ****************/

void read_bed_file(char *file)
/* READ BED FILE */
{		
	struct read *thisCluster = NULL;
	
	int  lastEnd = -1;
	char *lastChrom = "x";
	char *lastStrand = "x";
	double threshold;

	FILE *f;
	f = fopen(file, "r");
	if(!f)
	{
		printf("cannot open %s\n", file);
	}
	else
	{
		char c; int header = 0; int e=0;
		while((c=getc(f)) != EOF) {

			// if "#" at the beginning of line -> header
			if(c == '#') header = 1;
			ungetc(c,f);
			
			// if the end of line is hit -> parse line
			if(header != 1){
				char chrom[50], id[100], strand[5], info[100];
				int  start, end;
				double height, freq = -1;

				if(type == 1) {
					int tmp1, tmp2, tmp3;
					char str1[30], str2[50], str3[50];
					e = fscanf(f, "%s %d %d %s %lf %s %d %d %s %d %s %s", chrom, &start, &end, id, &height, strand, &tmp1, &tmp2, str1, &tmp3, str2, str3);
					//e = fscanf(f, "%s %d %d %s %lf %s", chrom, &start, &end, id, &height, strand);
					//printf("%s\t%d\t%d\t%s\t%lf\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n", chrom, start, end, id, height, strand, tmp1, tmp2, str1, tmp3, str2, str3);
					//if(isdigit(start) == 1 || isdigit(end) == 1 || strand == NULL){printf("wrong file format\n"); printUsage(); exit(0);}
					if(start < 0 || end < 0 || end < start){printf("wrong file format at %d %d\n", start, end); printUsage(); exit(0);}
				}
				else if(type == 2){
					int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
					char str1[10]; char str2[200];
					//e = fscanf(f, "%s %d %d %lf %d %d %d %d %d %d %s %d %d %s %s %lf", info, &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &tmp7, &tmp8, &tmp9, strand, &start, &end, chrom, str2, &freq); 
					//printf("%s\t%d\t%d\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%lf\n", info, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, strand, start, end, chrom, str2, freq);
					e = fscanf(f, "%s %s %d %d %d %d %d %d %d %s %d %d %s %s %lf", str1, info, &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &tmp7, strand, &start, &end, chrom, str2, &freq); 

					// split id|freq
					char *q = strtok(info, "|"); int j = 0;
					while(q != 0)
					{
						if(j == 0){
							strcpy(id, q);
							strcat(id, "_");
							char sFreq[10];
							sprintf(sFreq, "%d", (int) freq);
							strcat(id, sFreq);
							//printf("%lf\t%s\n", freq, sFreq);
						}
						if(j == 1) height = atof(q);
						j++;
						q = strtok(NULL, "|");
					}
					if(height != -1){height = height / freq;}
					else{height = 1 / freq;}

					//if(isdigit(start) == 1 || isdigit(end) == 1 || strand == NULL){printf("wrong file format\n"); printUsage(); exit(0);}
					if(start < 0 || end < 0 || end < start){printf("wrong file format at %d %d\n", start, end); printUsage(); exit(0);}
				}

				if(height >= tagFilter && (chr==NULL || strcmp(chrom, chr)==0)){
					if(strcmp(chrom, lastChrom) != 0 || strcmp(strand, lastStrand) != 0 || (start - lastEnd) > distance)
					{					
						if(clusterHeight > minClusterHeight && (clusterEnd-clusterStart)<=1000000){
							// ANALYZE CLUSTER
							assignReadsToBlocks(thisCluster);
							if(strcmp(blockheight, "rel")==0)
								threshold = (clusterHeight * minblockheight)/100;
							else
								threshold = minblockheight;
					
							writeBlocks(thisCluster, threshold);
						}
						
						//  linked list for next cluster
						freeList(thisCluster);
						thisCluster = NULL;
						
						// reset cluster dimensions
						clusterStart = start; clusterEnd = end;					
						clusterHeight = height;
						tagCount = 1;
					}
					else{
						// update cluster dimensions
						if(clusterStart > start) clusterStart = start;
						if(clusterEnd < end) clusterEnd = end;
						clusterHeight += height;
						tagCount++;
					}	
					
					// add read in linked list
					struct read *thisRead = NULL;
					thisRead = (struct read*)malloc(sizeof(struct read));
					thisRead->chrom = malloc(strlen(chrom)+1);
					strcpy(thisRead->chrom, chrom); 
					thisRead->id = malloc(strlen(id)+1);
					strcpy(thisRead->id, id);
					thisRead->strand = malloc(strlen(strand)+1);
					strcpy(thisRead->strand, strand);	
					thisRead->start  = start;
					thisRead->end    = end;
					thisRead->height = height;
					thisRead->block  = -1;
					thisRead->next   = thisCluster; 
					thisCluster      = thisRead; 
					
					clusterChrom = realloc(NULL, strlen(chrom)+1); 
					strcpy(clusterChrom, chrom);
					clusterStrand = realloc(NULL, strlen(strand)+1);
					strcpy(clusterStrand, strand);
					lastChrom = realloc(NULL, strlen(chrom)+1); 
					strcpy(lastChrom, chrom);
					lastStrand = realloc(NULL, strlen(strand)+1);
					strcpy(lastStrand, strand);
					lastEnd    = end;	
				}
			}
		} 
		if(clusterHeight > minClusterHeight &&  (clusterEnd-clusterStart)<=1000000){
			assignReadsToBlocks(thisCluster);
			if(strcmp(blockheight, "rel")==0)
				threshold = (clusterHeight * minblockheight)/100;
			else
				threshold = minblockheight;
			
			writeBlocks(thisCluster, threshold);
		}
		freeList(thisCluster);
	}
	fclose(f);
}

char *append(char *oldstring, char c)
{
    int result;
    char *newstring;
    result = asprintf(&newstring, "%s%c", oldstring, c);
    if (result == -1) newstring = NULL;
    return newstring;
}

void writeHeader(char *file)
{	
	time_t t;
	t = time(NULL);
	printf("# blockbuster result file generated %s# query file: %s\n# scale: %.1f, minblockheight: %lf, mergeDistance: %i\n# block_number\tchromosome\tstart_of_block\tend_of_block\tstrand\treadIDs\n", ctime(&t), file, sizescale, minblockheight, merge);
}

void printUsage(void)
{ 
	printf("\nusage: blockbuster.x [OPTIONS] <file>\n");
	printf("Detect blocks of overlapping reads using a gaussian distribution approach\n\n");
	printf("[OPTIONS]\n");
	printf("-format <int>              file format of input file (default: 1)\n");
	printf("                            (1) bed\n");
	printf("                            (2) segemehl-output\n");
	printf("-distance <int>            minimum distance between two clusters (default: 30)\n");
	printf("-minClusterHeight <double> minimum height of a cluster (default: 50)\n");
	printf("-minBlockHeight <double>   minimum height of a block (default: 2)\n");
	printf("-scale <int>               scale stddev for a single read (default: 0.2)\n");
	printf("-merge <int>               merge reads with almost similar means (default: 0)\n");
	printf("-tagFilter <int>           skip tags with expression smaller than this value (default: 0)\n");
	printf("-print <int>               print out: (1) blocks (2) reads (default: 1)\n");
	printf("-chr <string>              chromosome (default: all)\n");
	printf("-blockHeight <string>      minBlockHeight threshold (default: abs)\n");
	printf("                            abs: absolute read count\n");
	printf("                            rel: relative read count with respect to expression of cluster\n");
	printf("[COMPLAINT DEPARTMENT]\n");
	printf("Please be nice when complaining to david@bioinf.uni-leipzig.de\n");
	printf("or steve@bioinf.uni-leipzig.de\n\n");
}

void assignReadsToBlocks(struct read *anchor)
/* ASSIGN READS TO BLOCKS */
{
	int blockCount = 1;
	
	// create a double array with clusterSize entries for the superGaussian distribution
	int clusterSize = (clusterEnd - clusterStart);
	double distrib[clusterSize];
	
	int old = 1;
	int new = 0;
			
	// run through sorted peaks
	while(old != new)
	{		
		old = getRest(anchor);

		// clean distribution array
		int p = 0;
		for(p = 0; p < clusterSize; p++){distrib[p] = 0;}
		
		// write distribution
		writeSuperGaussian(anchor, distrib, clusterSize);
		//	for(p = 0; p < clusterSize; p++){printf("x: %i\ty: %f\n", p, distrib[p]);}
		int highestPeak = getHighestPeak(distrib, clusterSize);

		// assign reads to the highest peak
		int sum = assignReads(anchor, highestPeak, readCount, blockCount);
		if(sum != 0) blockCount++;
		
		new = getRest(anchor);
	}
}

double gaussian(double x, double mean, double variance)
/* CALCULATE THE GAUSSIAN */
{
	return (1/(variance * sqrt(2 * pi))) * exp(-1 * (pow((x - mean),2))/pow((2 * variance),2));
}
 

int getHighestPeak(double *distrib, int size)
/* GET THE HIGHEST PEAK IN THE DISTRIBUTION ARRAY */
{
	int i;
	double max = 0;
	int result = 0;
	for(i = 0; i < size; i++)
	{
		if(distrib[i] > max){result = i; max = distrib[i];}
	}
	distrib[result] = 0;
	return result;
}

int getRest(struct read *anchor)
/* GET THE AMOUNT OF READS NOT ASSIGNED TO BLOCKS */
{
	int sum = 0;	
	while (anchor != NULL)
	{
		if(anchor->block == -1){sum++;}
		anchor = anchor->next;
	}
	return (sum);
}

void writeSuperGaussian(struct read *anchor, double *distrib, int clusterHeight)
/* CALCULATE THE GAUSSIAN DISTRIBUTIONS OF ALL READS AND SUM THEM UP */
{
	while (anchor != NULL)
	{
		if(anchor->block == -1)
		{
			double mean = ((anchor->start + anchor->end) / 2) - (double) clusterStart;
			double variance = sizescale * (abs(anchor->end - anchor->start)/2);
		
			int i = 0;
			for(i = 0; i <= 2 * variance; i++)
			{
				double x = mean + i;
				double y = anchor->height * gaussian(x, mean, variance);
				if((int) x < clusterHeight) distrib[(int) x] += y;
				if((int) (mean -i) > 0) distrib[(int)(mean - i)] += y;
			}						
		}
		anchor = anchor->next;
	}
}

int assignReads(struct read *anchor, int highestPeak, int clusterSize, int blockCount)
/* ASSIGN READS TO A BLOCK */
{	
	double* readMeans;
	double* readHeights;
	readMeans = malloc(sizeof(double) * tagCount);
	readHeights = malloc(sizeof(double) * tagCount);
	
	int meanCounter = 0;
	int p;
	for(p = 0; p < tagCount; p++){readMeans[p] = -1; readHeights[p] = -1;}
	
	int counterNew = 0; int counterOld = -1;

	while(counterOld != counterNew) 
	{	
		double dev = stddev(readMeans, readHeights, tagCount);
		
		counterOld = counterNew;
		struct read *start = anchor;
		while (start != NULL)
		{
			if(start->block == -1){
				double mean = ((start->start + start->end) / 2) - clusterStart;
				double variance = sizescale * (abs(start->end - start->start)/2);
				
				if(((mean-variance-dev) <= highestPeak && (mean+variance+dev) >= highestPeak) || (mean >= (highestPeak - merge) && mean <= (highestPeak + merge)))
				{
//					printf("id: %s\tstart: %i\tend: %i\tpeak: %i\tmean: %f\tvariance: %f\tdev: %f\tblock: %i\n", start->id, start->start, start->end, highestPeak, mean, variance, dev, blockCount);

					readMeans[meanCounter] = mean;
					readHeights[meanCounter] = start->height;
					meanCounter++;
					start->block = blockCount;
					counterNew++;
				}
//				printf("id: %s\tstart: %i\tend: %i\tpeak: %i\tmean: %f\tvariance: %f\tdev: %f\tblock: %i\n", start->id, start->start, start->end, highestPeak, mean, variance, dev, start->block);
			}
			start = start->next;
		}
	}
	free(readMeans);
	free(readHeights);
	return counterNew;
}

double stddev(double *readMeans, double *readHeights, int size)
/* CALCULATE THE STANDARD DEVIATION */
{
	double sum = 0;
	int counter = 0;
	
	int i;
	for(i = 0; i < size; i++)
	{
		if(readMeans[i] != -1) 
		{
			int j;
			for(j = 0; j < readHeights[i]; j++) 
			{
				sum += readMeans[i];
				counter++;
			}
		}
	}
	
	if(counter == 0) return 0;
	
	double mean = sum / (double) counter;
	
	sum = 0;
	for(i = 0; i < size; i++)
	{
		if(readMeans[i] != -1) 
		{
			int j;
			for(j = 0; j < readHeights[i]; j++) 
			{
				sum += pow(((readMeans[i]) - mean), 2);
			}
//			printf("value: %f\tmean: %f\tsum: %f\n", means[i], mean, sum);
		}
	}
	return sqrt(sum / counter);
}

void writeBlocks(struct read *anchor, double threshold)
/* WRITE THE READS THAT ARE ASSIGNED TO A BLOCK TO STDOUT */
{
	
	int    thisBlock         = 0; 
	double blockHeight       = 0;	
	int    blockNb           = 0;
	int    thisTagCount      = 0;
	double thisClusterHeight = 0;
	int    absTagCount       = 0;
	double absClusterHeight  = 0;
	int    size              = 1;
	int    thisBlockStart    = -1;
	int    thisBlockEnd      = -1;
	int    thisClusterStart  = -1;
	int    thisClusterEnd    = -1;

	// get cluster information
	while(size > 0) 
	{
		thisBlock++;
		
		// reset variables
		size = 0;
		blockHeight = 0;
		thisClusterHeight = 0;
		thisTagCount = 0;
		thisBlockStart = -1;
		thisBlockEnd = -1;
		
		// run through linked list of reads
		struct read *start = anchor;
		while (start != NULL)
		{
			// if current read is in thisBlock
			if(start->block == thisBlock)
			{
				if(thisBlockStart == -1){thisBlockStart = start->start; thisBlockEnd = start->end;}
				if(thisBlockStart > start->start){thisBlockStart = start->start;}
				if(thisBlockEnd < start->end){thisBlockEnd = start->end;}
				size++; 
				blockHeight += start->height;	
				thisClusterHeight += start->height;
				thisTagCount++;
			}
			start = start->next;
		}		

		// check if block is high enough
		if(blockHeight >= threshold && size > 0)
		{
			if(thisClusterStart == -1){thisClusterStart = thisBlockStart; thisClusterEnd = thisBlockEnd;}
			if(thisClusterStart > thisBlockStart){thisClusterStart = thisBlockStart;}
			if(thisClusterEnd < thisBlockEnd){thisClusterEnd = thisBlockEnd;}
			blockNb++;
			absClusterHeight += thisClusterHeight;
			absTagCount += thisTagCount;
		} 
	}

	if(blockNb > 0 && absClusterHeight > minClusterHeight)
	{
		clusterCounter++;
		
		// print header
		printf(">cluster_%i\t%s\t%i\t%i\t%s\t%.2f\t%i\t%i\n", clusterCounter, clusterChrom, thisClusterStart, thisClusterEnd, clusterStrand, absClusterHeight, absTagCount, blockNb);
		
		// print blocks
		if(print == 1)
		{
			thisBlock = 0;
			size = 1;
			int writeBlock = 0;
			while(size > 0) 
			{
				thisBlock++;
				size = 0;		
				double thisBlockHeight = 0;
				int    thisBlockTags   = 0;
				int    thisBlockStart  = -1;
				int    thisBlockEnd    = -1;
				struct read *start = anchor;
				while (start != NULL)
				{
					if(start->block == thisBlock)
					{
						if(thisBlockStart == -1){thisBlockStart = start->start; thisBlockEnd = start->end;}
						if(start->start < thisBlockStart) thisBlockStart = start->start;
						if(start->end > thisBlockEnd) thisBlockEnd = start->end;
						thisBlockHeight += start->height;
						thisBlockTags++;
						size++;
					}
					start = start->next;
				}
				if(thisBlockHeight >= threshold && size > 0)
				{
					writeBlock++;
					printf("%i\t%s\t%i\t%i\t%s\t%.2f\t%i\n", writeBlock, clusterChrom, thisBlockStart, thisBlockEnd, clusterStrand, thisBlockHeight, thisBlockTags);
				}
			}
		}
		
		// print tags
		if(print == 2)
		{
			thisBlock = 0;
			size = 1;
			int    writeBlock = 0;
			while(size > 0) 
			{
				double thisBlockHeight = 0;
				thisBlock++;
				size = 0;				
				struct read *start = anchor;
				while (start != NULL)
				{
					if(start->block == thisBlock)
					{
						thisBlockHeight += start->height;
						size++;
					}
					start = start->next;
				}	
				if(thisBlockHeight >= threshold && size > 0)
				{
					writeBlock++;
					struct read *start = anchor;
					while (start != NULL)
					{
						if(start->block == thisBlock)
						{
							printf("%s\t%d\t%d\t%s\t%lf\t%s\t%i\n",start->chrom,start->start,start->end,start->id,start->height,start->strand,writeBlock);	
						}
						start = start->next;
					}
				}
			}
		}
	}
}

			   

void freeList(struct read *anchor)
{
	struct read *tmp;
	
	while (anchor != NULL) 
	{
		tmp = anchor->next;
//		free(anchor->id);
//		free(anchor->chrom);
//		free(anchor->strand);
//		free(anchor->next);
		free(anchor);
		anchor = tmp;
	}
}

			   
