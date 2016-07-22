#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// the preanalysis of the file;
// the file should be analyzed by perl before, to put amount in the end!
// please delete the last '\n'


char  locus[100000000][300] = {0};
int main(int argc, char **argv)
{
  FILE * input, *output, *media;
  input = fopen(argv[1], "r");
  media = fopen("median", "w+");
  output = fopen(argv[2], "w");

  long record = 0, length = 0, situ = 0;
  long amount_all = 0;

  while(!feof(input))
    {
      char line[1000] = {0};
      fgets(line, 1000, input);

      long i = 0;
      char temp[1000] = {0};


      while(line[i] != '\t')
	i++;
      long t = 0;
      i++;

      while(line[i] != '\t')
	{
	  temp[t] = line[i];
	  i++;
	  t++;
	}
      temp[t] = '\0';

      i++;
      long pla_i = 0;
      while(pla_i < 3)
	{
	   while(line[i] != '\t')
	     i++;
	   i++;
	  pla_i++;
	}
      amount_all += atoi(&line[i]);

      long j = 0;
      int flag = 1;
      if(situ != 0)
	{
	  while(locus[j][0] != 0)
	    {
	      if(strcmp(locus[j], temp) == 0)
		{
		  flag = 0;
		  break;
		}
	  j++;
	    }
	}

      length = strlen(line);

      if(line[length - 1] == '\n')
	{
	length = length - 1;
	}

      char group[1000] ;
      sprintf(group,"%ld\n", j);

      if(flag == 0)
	{
	  line[length] = '\t';
	  strcat(&line[length + 1], group);
	}
      else
	{
	  strcpy(locus[j], temp);
	  line[length] = '\t';
	  strcat(&line[length + 1], group);
	  record = j;

	  printf("%ld:\n", record);
	}
      fputs(line, media);

      situ++;
    }

  printf("{2}\n");

  fseek(media, 0, SEEK_SET);
  int l = 0;

  while(l <= record)
    {
      printf("%d\n", l);
      char line[1000] = {0};
      while(fgets(line, 1000, media) != NULL)
	{
	  length = strlen(line);
	  int fc = 0, spa = 0;
	  while(fc < 6)
	    {
	      while(line[spa] != '\t')
		spa++;
	      spa++;
	      fc++;
	    }
	  if(atoi(&line[spa]) == l)
	    {
	    fputs(line, output);
	    }
	}

      fseek(media, 0, SEEK_SET);
      l++;
      fputs("//\n", output);
    }
  printf("\nthe total amount : %ld\n", amount_all);
  fclose(input);
  fclose(media);
  fclose(output);
}
