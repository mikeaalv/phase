#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define THREDHOLD 0.001
#define UNI_READ_WIN_1 10
#define UNI_READ_WIN_2 10
#define UNI_RAO_21_1 0.5
#define UNI_RAO_21_2 0.25
#define CER_LOC 3
#define SIZE 200000
#define TEST_NUM  2.0E4// to reduce false positive, bonferroni correction {estimation}, which is somewhat too strict
#define AMOUNT 2000
#define Cyc_len_1 21
#define Cyc_len_2 21
//#define LIMITS_NUM 141

char data[SIZE][100] = {0};
long place[SIZE][3] = {0};
int uni_oth[AMOUNT] = {0}, uni_oth_nolim[AMOUNT][50] = {0}, sig[AMOUNT] = {0}, uni_oth_min[SIZE] = {0}, uni_oth_nolim_min[AMOUNT][50] = {0};

//the border are all included
//regard repeat sequences in a windows as distinct ones
//fit in situation as a lot of short reference sequences such as gene sequences rather than genome sequcnes.
// ./pandphase input cycle_num readamout
//format of input should be the ouput of preanalysis.c
//NS = not significant, NE = not enough, SN = significant but not strict, SL = significant but only loosely strict, SS = significant and strict
//all unique reads binding, n unique reads binding 21bp, k unique reads in regiester 21bp, amount unique reads in register 21bp add together, reads_loci reads amount figure.

double CC(long k, long n)//the combination, and optimization on read accuracy and limits
{
  long temp = n - k;
  double result = 1;
  while(k > 0)
    {
      result *= ((double)n/(double)k);
      n--;
      k--;
    }
  while(temp > 0)
    {
      result *= ((double)n/(double)temp);
      temp--;
      n--;
    }

  return result;
}



int main(int argc, char **argv)
{
  FILE * input, * output;
  int cycle, number = 0;
  input = fopen(argv[1], "r");
  output = fopen("Phasescore", "w");//file
  cycle = atoi(argv[2]);

  printf("the cycle number: %d\n",cycle);//cycle number
  while(!feof(input))
    {
      int st = 0;
      while(st < AMOUNT)//for sense or antisense
	{
	  sig[st] = 0;
	  st++;
	}
      number++;
      printf("%d:\n", number);//for the CDS or chromosome or groups number
      int i = -1;
      do
	{
	  i++;
	  fgets(data[i], 100, input);//data of a group

	 }while(data[i][0] != '/');//the end of a group

      int t = 0;
      char name[30];//the unique read number in a group
      while(t < i)
	{

	  int c = 0;
	  while(data[t][c] != '\t')//omit read name
	    c++;
	  c++;

	  int l = 0;
	  while(data[t][c] != '\t')
	    {
	      name[l] = data[t][c];//get group name
	      l++;
	      c++;
	    }
	  name[l] = '\0';

	  c++;
	  place[t][0] = atoi(&data[t][c]);//get read start
	  while(data[t][c] != '\t')
	    c++;
	  c++;
	  place[t][1] = atoi(&data[t][c]);//get read end
	  while(data[t][c] != '\t')
	    c++;
	  c++;
	  if(data[t][c] == '-')//for sense and antisense sequence
	    { place[t][0] += 2;
	      place[t][1] += 2;
	      sig[t] = 1;
	    }
	  c++;
	  place[t][2] = atoi(&data[t][c]);//amount of the read
	  t++;
	}

         fprintf(output, "\n**************************\nlocus: %s\t number of sequence: %d\n the score number is of the middle of the window (locus : score probablity )", name, i);

       long max = place[0][1], min = place[0][0];
       int max_loci = 0, min_loci = 0;
	  int c = 0;
	while(c < i)	    //find the mim and max of all the read
		{
	      if(place[c][0] < min)
		{
		  min = place[c][0];
		  min_loci = c;
		}
	      if(place[c][1] > max)
		{
		  max = place[c][1];
		  max_loci = c;
		}
		c++;
	}

	  long run_min = min, run_max = max;//use the max and min as the border of searching

	  if((run_max - run_min) > Cyc_len_2 * cycle)//can use to screen at least one window
    {

           while(run_min < run_max)//every window
	    {
	      int othe = 0;
	      long num_win = 0, num_map = 0, num_win_no21 = 0, amount_read = 0, num_win_min = 0, num_map_min = 0, num_win_no21_min = 0, amount_loci = 0;
	      int uni_reg[100] = {0} , uni_reg_num = 0, j = 0, uni_reg_min[100] = {0}, uni_reg_num_min = 0;
	      double amount_window = 0.0, amount_num = 0.0;
	      double pro = 0.0;
	      int size_i = 0;
	      while(size_i < AMOUNT)
		{
		  uni_oth[size_i] = 0;
		  uni_oth_min[size_i] = 0;
		  int i_50 = 0;
		  while(i_50 < 50)
		    {
		  uni_oth_nolim[size_i][i_50] = 0;
		  uni_oth_nolim_min[size_i][i_50] = 0;
		  i_50++;
		    }
		  size_i++;
		}

	      while(othe < i) // n, k :::every read
		{
		  if(place[othe][0] >= run_min && place[othe][0] < run_min + Cyc_len_2 * cycle)//in the window
		    {
		      if(place[othe][1] - place[othe][0] == Cyc_len_2)//21bp read
			{
			  if((place[othe][0] - run_min) % Cyc_len_2 == 0)//in phase
			    {
			      int sitp = (place[othe][0] - run_min) / Cyc_len_2;
			      if((uni_reg[sitp] == 0 && sig[othe] == 0) || (uni_reg_min[sitp] == 0 && sig[othe] == 1))
				{
				  amount_read += place[othe][2];//21 in phase amount unique
				}
			      if(sig[othe] == 0)//this site has been calculated
				uni_reg[sitp] = 1;
			      else
				uni_reg_min[sitp] = 1;
			    }
			      if(sig[othe] == 0 && uni_oth[place[othe][0] - run_min] == 0)//for all 21bp this site has been marked
                {
                    uni_oth[place[othe][0] - run_min] = 1;
                    amount_num += place[othe][2];
                }
			      else if(sig[othe] == 1 && uni_oth_min[place[othe][0] - run_min] == 0)
                {
				uni_oth_min[place[othe][0] - run_min] = 1;
                amount_num += place[othe][2];
                }
			}

			  int nolim_i =0, nolim_min_i = 0;
		      if(sig[othe] == 0)//for all reads and their amount
			{
			  while(uni_oth_nolim[place[othe][0] - run_min][nolim_i] != place[othe][1] - place[othe][0] && uni_oth_nolim[place[othe][0] - run_min][nolim_i] != 0)
			    nolim_i++;
			  if(uni_oth_nolim[place[othe][0] - run_min][nolim_i] == 0)
			    {
			    amount_window += place[othe][2];
			    num_win_no21++;

			    uni_oth_nolim[place[othe][0] - run_min][nolim_i] = place[othe][1] - place[othe][0];
			    }
			   }
		      else
			{
			  while(uni_oth_nolim_min[place[othe][0] - run_min][nolim_min_i] != place[othe][1] - place[othe][0] && uni_oth_nolim_min[place[othe][0] - run_min][nolim_min_i] != 0)
			    nolim_min_i++;
			   if(uni_oth_nolim_min[place[othe][0] - run_min][nolim_min_i] == 0)
			     {
			       amount_window += place[othe][2];
			       num_win_no21_min++;

			       uni_oth_nolim_min[place[othe][0] - run_min][nolim_min_i] = place[othe][1] - place[othe][0];
			     }
			    }
		    }
			if(place[othe][0] == run_min + Cyc_len_2 * cycle / 2)
			amount_loci += place[othe][2];
		  othe++;
		}

		int sip = 0;
	      while(sip < cycle)
		{
		  if(uni_reg[sip] == 1)
		    uni_reg_num++;
		sip ++;
		}
	      sip = 0;
	      while(sip < cycle * Cyc_len_2)
		{
		  if(uni_oth[sip] == 1)
		    num_win++;
		  sip ++;
		}
	       sip = 0;
	        while(sip < cycle)
		{
		  if(uni_reg_min[sip] == 1)
		    uni_reg_num_min++;
		sip ++;
		}
	      sip = 0;
	      while(sip < cycle * Cyc_len_2)
		{
		  if(uni_oth_min[sip] == 1)
		    num_win_min++;
		sip ++;
		}
	      sip = 0;

	       double score = 0;
	       score = (double)(uni_reg_num + uni_reg_num_min - 2) * log(1.0 + 10.0 * ((double)amount_read / (double)(1.0  + (double)amount_num - (double)amount_read)));
	       fprintf(output, "\nthe window at %ld -- %ld\t with all = %ld n = %ld, k = %d, amount\t%ld\treads_loci\t%ld\t", run_min, run_min + Cyc_len_2 * cycle, num_win_no21 + num_win_no21_min, num_win + num_win_min, uni_reg_num + uni_reg_num_min, amount_read, amount_loci);//all: all unique no limits on length; n: all unique and 21bp; k unique 21bp in phase; amount_read: 21bp in phase amount; read_loci: read amount at this site
	       if(uni_reg_num + uni_reg_num_min >= 3)//can calculate Phasiscore
	       {
	       if(score > 0)
		{
		  fprintf(output, "%ld: %lf\t", run_min + Cyc_len_2 * cycle / 2, score);//position and score
		}
		else
		  fprintf(output, "%ld: %lf\t", run_min + Cyc_len_2 * cycle / 2, 0.0);
	       }
	       else
		 {
		   fprintf(output, "\t\t");
		 }

	       num_map = uni_reg_num + uni_reg_num_min;
	       if(uni_reg_num + uni_reg_num_min >= 3)
		 {
		   while(num_map <= cycle * 2)
		     {
		       pro += CC(num_win + num_win_min - num_map, (Cyc_len_2 - 1) * cycle * 2) * CC(num_map, cycle * 2) / CC(num_win + num_win_min, Cyc_len_2 * cycle * 2);
		       num_map++;

		     }
		   fprintf(output,"%lf\t", pro);
		   if(pro < (0.001 / TEST_NUM))
		  {

		    if((num_win_no21 + num_win_no21_min >= UNI_READ_WIN_1) && ((double)(num_win + num_win_min) / (double)(num_win_no21+ num_win_no21_min)) >= UNI_RAO_21_1 && (uni_reg_num + uni_reg_num_min) >= CER_LOC)
		      {
			fprintf(output,"SS\n");
		      }
		    else if(((num_win_no21 + num_win_no21_min) >= UNI_READ_WIN_2) && ((double)(num_win + num_win_min) / (double)(num_win_no21 + num_win_no21_min)) >= UNI_RAO_21_2 && (uni_reg_num + uni_reg_num_min) >= CER_LOC)
		      {
			fprintf(output,"SL\n");
		      }
		    else
		      {
			fprintf(output,"SN\n");
		      }
		  }
		   else if(pro < 0.001)
		     {
		       if((num_win_no21 + num_win_no21_min >= UNI_READ_WIN_1) && ((double)(num_win + num_win_min) / (double)(num_win_no21+ num_win_no21_min)) >= UNI_RAO_21_1 && (uni_reg_num + uni_reg_num_min) >= CER_LOC)
		      {
			fprintf(output,"sS\n");
		      }
		    else if(((num_win_no21 + num_win_no21_min) >= UNI_READ_WIN_2) && ((double)(num_win + num_win_min) / (double)(num_win_no21 + num_win_no21_min)) >= UNI_RAO_21_2 && (uni_reg_num + uni_reg_num_min) >= CER_LOC)
		      {
			fprintf(output,"sL\n");
		      }
		    else
		      {
			fprintf(output,"sN\n");
		      }
		     }
		  else if((uni_reg_num + uni_reg_num_min) >= 3)
		  {
		    fprintf(output, "NS\n");
		  }
		 }
	       else
		 {
		   fprintf(output,"\t");
		   fprintf(output, "NE\n");//not enough data to predict

		 }
		run_min++;
	    }
    }
	  else
	    {
	      fprintf(output, "too little mapping in this area/sequence for a prediction\n");
	    }
	  fprintf(output,"\n//\n");
    //    if(number == 141)
      //  {exit(0);}
   //     printf("**%d**\n", number);
	    }
}
