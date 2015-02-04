/*
  io.c: I/O routines for RNAlila
  Last changed Time-stamp: <2015-02-04 18:15:05 mtw>
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <lila.h>

/**/
void
lila_parse_seq_struc(FILE *fp)
{
  char *line=NULL;
  
  line = get_line(fp);
  /* skip comment lines */
  while ((*line == '*')||(*line == '\0')||(*line == '>')) {
    free(line);
    line = get_line(fp);
  }
  lilass.sequence  = (char *) calloc (strlen(line)+1, sizeof(char));
  lilass.structure = (char *) calloc (strlen(line)+1, sizeof(char));
  assert(lilass.sequence != NULL); assert(lilass.structure != NULL);
  sscanf(line, "%s", lilass.sequence);
  free (line);
  line = get_line(fp);
  sscanf(line, "%s", lilass.structure);
  free (line);
  lilass.length = strlen(lilass.sequence);
}

/**/
void
lila_parse_sequence(FILE *fp)
{
  char *line=NULL;
  
  line = get_line(fp);
  /* skip comment lines */
  while ((*line == '*')||(*line == '\0')||(*line == '>')) {
    free(line);
    line = get_line(fp);
  }
  lilass.sequence  = (char *) calloc (strlen(line)+1, sizeof(char));
  assert(lilass.sequence != NULL);
  sscanf(line, "%s", lilass.sequence);
  free (line);
  lilass.length = strlen(lilass.sequence);
}
