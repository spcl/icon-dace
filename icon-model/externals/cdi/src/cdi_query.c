#include <string.h>
#include <stdlib.h>
#include "cdi.h"

static int
sum_found(int listSize, bool *listFound)
{
  int numFound = 0;
  for (int i = 0; i < listSize; ++i) numFound += listFound[i];
  return numFound;
}

static int
sum_not_found(int listSize, bool *listFound)
{
  return listSize - sum_found(listSize, listFound);
}

void
cdiQueryInit(CdiQuery *query)
{
  query->numEntries = 0;

  query->numNames = 0;
  query->names = NULL;
  query->namesFound = NULL;

  query->numCellidx = 0;
  query->cellidx = NULL;
  query->cellidxFound = NULL;

  query->numLevidx = 0;
  query->levidx = NULL;
  query->levidxFound = NULL;

  query->numStepidx = 0;
  query->stepidx = NULL;
  query->stepidxFound = NULL;
}

CdiQuery *
cdiQueryCreate(void)
{
  CdiQuery *query = (CdiQuery *) malloc(sizeof(CdiQuery));
  cdiQueryInit(query);
  return query;
}

void
cdiQueryDelete(CdiQuery *query)
{
  if (query)
    {
      if (query->numNames)
        {
          for (int i = 0; i < query->numNames; ++i) free(query->names[i]);
          free(query->names);
          free(query->namesFound);
        }

      if (query->numCellidx)
        {
          free(query->cellidx);
          free(query->cellidxFound);
        }

      if (query->numLevidx)
        {
          free(query->levidx);
          free(query->levidxFound);
        }

      if (query->numStepidx)
        {
          free(query->stepidx);
          free(query->stepidxFound);
        }

      cdiQueryInit(query);
      free(query);
    }
}

int
cdiQueryNumNames(const CdiQuery *query)
{
  return query ? query->numNames : 0;
}

int
cdiQueryNumCellidx(const CdiQuery *query)
{
  return query ? query->numCellidx : 0;
}

int
cdiQueryNumStepidx(const CdiQuery *query)
{
  return query ? query->numStepidx : 0;
}

int
cdiQueryNumEntries(const CdiQuery *query)
{
  return query ? query->numEntries : 0;
}

void
cdiQuerySetNames(CdiQuery *query, int numEntries, char **names)
{
  if (numEntries)
    {
      query->numEntries += numEntries;
      query->numNames = numEntries;
      query->namesFound = (bool *) calloc(numEntries, sizeof(bool));
      query->names = (char **) malloc(numEntries * sizeof(char *));
      for (int i = 0; i < numEntries; ++i) query->names[i] = strdup(names[i]);
    }
}

void
cdiQuerySetCellidx(CdiQuery *query, int numEntries, size_t *cellidx)
{
  if (numEntries)
    {
      query->numEntries += numEntries;
      query->numCellidx = numEntries;
      query->cellidxFound = (bool *) calloc(numEntries, sizeof(bool));
      query->cellidx = (size_t *) malloc(numEntries * sizeof(size_t));
      for (int i = 0; i < numEntries; ++i) query->cellidx[i] = cellidx[i];
    }
}

void
cdiQuerySetLevidx(CdiQuery *query, int numEntries, int *levidx)
{
  if (numEntries)
    {
      query->numEntries += numEntries;
      query->numLevidx = numEntries;
      query->levidxFound = (bool *) calloc(numEntries, sizeof(bool));
      query->levidx = (int *) malloc(numEntries * sizeof(int));
      for (int i = 0; i < numEntries; ++i) query->levidx[i] = levidx[i];
    }
}

void
cdiQuerySetStepidx(CdiQuery *query, int numEntries, int *stepidx)
{
  if (numEntries)
    {
      query->numEntries += numEntries;
      query->numStepidx = numEntries;
      query->stepidxFound = (bool *) calloc(numEntries, sizeof(bool));
      query->stepidx = (int *) malloc(numEntries * sizeof(int));
      for (int i = 0; i < numEntries; ++i) query->stepidx[i] = stepidx[i];
    }
}

size_t
cdiQueryGetCellidx(const CdiQuery *query, int index)
{
  return (index >= 0 && index < query->numCellidx) ? query->cellidx[index] : (size_t) -1;
}

CdiQuery *
cdiQueryClone(const CdiQuery *query)
{
  CdiQuery *queryOut = cdiQueryCreate();

  if (query)
    {
      cdiQuerySetNames(queryOut, query->numNames, query->names);
      cdiQuerySetCellidx(queryOut, query->numCellidx, query->cellidx);
      cdiQuerySetLevidx(queryOut, query->numLevidx, query->levidx);
      cdiQuerySetStepidx(queryOut, query->numStepidx, query->stepidx);
    }

  return queryOut;
}

static void
print_list_compact_int(int n, const int *list)
{
  for (int i = 0; i < n; ++i)
    {
      int value = list[i];
      printf(" %d", value);
      if ((i + 2) < n && (value + 1) == list[i + 1] && (value + 2) == list[i + 2])
        {
          printf("/to/");
          int last = list[++i];
          while ((i + 1) < n && (last + 1) == list[i + 1]) last = list[++i];
          printf("%d", last);
        }
    }
  printf("\n");
}

void
cdiQueryPrint(const CdiQuery *query)
{
  if (query)
    {
      if (query->numNames)
        {
          printf("Names:");
          for (int i = 0; i < query->numNames; ++i) printf(" %s", query->names[i]);
          printf("\n");
        }

      if (query->numCellidx)
        {
          printf("Cellidx:");
          for (int i = 0; i < query->numCellidx; ++i) printf(" %zu", query->cellidx[i]);
          printf("\n");
        }

      if (query->numLevidx)
        {
          printf("Levidx:");
          print_list_compact_int(query->numLevidx, query->levidx);
        }

      if (query->numStepidx)
        {
          printf("Stepidx:");
          print_list_compact_int(query->numStepidx, query->stepidx);
        }
    }
}

int
cdiQueryNumEntriesFound(const CdiQuery *query)
{
  int numEntriesFound = 0;

  if (query)
    {
      if (query->numNames) numEntriesFound += sum_found(query->numNames, query->namesFound);
      if (query->numCellidx) numEntriesFound += sum_found(query->numCellidx, query->cellidxFound);
      if (query->numLevidx) numEntriesFound += sum_found(query->numLevidx, query->levidxFound);
      if (query->numStepidx) numEntriesFound += sum_found(query->numStepidx, query->stepidxFound);
    }

  return numEntriesFound;
}

void
cdiQueryPrintEntriesNotFound(const CdiQuery *query)
{
  if (query)
    {
      int numEntriesNotFound = cdiQueryNumEntries(query) - cdiQueryNumEntriesFound(query);
      if (numEntriesNotFound > 0)
        {
          if (query->numNames)
            {
              if (sum_not_found(query->numNames, query->namesFound) > 0)
                {
                  printf("Name not found:");
                  for (int i = 0; i < query->numNames; ++i)
                    if (!query->namesFound[i]) printf(" %s", query->names[i]);
                  printf("\n");
                }
            }

          if (query->numCellidx)
            {
              if (sum_not_found(query->numCellidx, query->cellidxFound) > 0)
                {
                  printf("Grid cell index not found:");
                  for (int i = 0; i < query->numCellidx; ++i)
                    if (!query->cellidxFound[i]) printf(" %zu", query->cellidx[i]);
                  printf("\n");
                }
            }

          if (query->numLevidx)
            {
              if (sum_not_found(query->numLevidx, query->levidxFound) > 0)
                {
                  printf("Level index not found:");
                  for (int i = 0; i < query->numLevidx; ++i)
                    if (!query->levidxFound[i]) printf(" %d", query->levidx[i]);
                  printf("\n");
                }
            }

          if (query->numStepidx)
            {
              if (sum_not_found(query->numStepidx, query->stepidxFound) > 0)
                {
                  printf("Step index not found:");
                  for (int i = 0; i < query->numStepidx; ++i)
                    if (!query->stepidxFound[i]) printf(" %d", query->stepidx[i]);
                  printf("\n");
                }
            }
        }
    }
}

int
cdiQueryName(CdiQuery *query, const char *name)
{
  if (query && query->numNames && name && *name)
    {
      for (int i = 0; i < query->numNames; ++i)
        if (strcmp(name, query->names[i]) == 0)
          {
            query->namesFound[i] = true;
            return 0;
          }
    }

  return -1;
}

int
cdiQueryCellidx(CdiQuery *query, size_t cellidx)
{
  if (query && query->numCellidx)
    {
      for (int i = 0; i < query->numCellidx; ++i)
        if (query->cellidx[i] == cellidx)
          {
            query->cellidxFound[i] = true;
            return 0;
          }
    }

  return -1;
}

int
cdiQueryLevidx(CdiQuery *query, int levidx)
{
  if (query && query->numLevidx)
    {
      for (int i = 0; i < query->numLevidx; ++i)
        if (query->levidx[i] == levidx)
          {
            query->levidxFound[i] = true;
            return 0;
          }
    }

  return -1;
}

int
cdiQueryStepidx(CdiQuery *query, int stepidx)
{
  if (query && query->numStepidx)
    {
      for (int i = 0; i < query->numStepidx; ++i)
        if (query->stepidx[i] == stepidx)
          {
            query->stepidxFound[i] = true;
            return 0;
          }
    }

  return -1;
}
