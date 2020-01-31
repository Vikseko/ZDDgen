#include <stdio.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <regex.h>
#include <ctype.h>
#include <iterator>
#include <time.h>
#include <chrono>
#include <queue>
#include <unordered_map>
#include <utility>


/* Files included from CUDD package */
#include "/home/vsk/Diagrams/cudd-3.0.0/util/util.h"
#include "/home/vsk/Diagrams/cudd-3.0.0/cudd/cudd.h"

static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}

template <typename T1, typename T2>
struct less_second {
    typedef std::pair<T1, T2> type;
    bool operator ()(type const& a, type const& b) const {
        return a.second > b.second;
    }
};

void print_dd (DdManager *gbm, DdNode *dd, int n, int pr )
{
    printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(gbm)); /*Reports the number of live nodes in BDDs and ADDs*/
    printf("DdManager vars: %d | ", Cudd_ReadSize(gbm) ); /*Returns the number of BDD variables in existence*/
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(gbm) ); /*Returns the number of times reordering has occurred*/
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(gbm) ); /*Returns the memory in use by the manager measured in bytes*/
    Cudd_zddPrintDebug(gbm, dd, n, pr);  // Prints to the standard output a DD and its statistics: number of nodes, number of leaves, number of minterms.
}

void write_dd (DdManager *gbm, DdNode *dd, char* filename)
{
    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(filename,"w");
    DdNode **ddnodearray = (DdNode**)malloc(sizeof(DdNode*)); // initialize the function array
    ddnodearray[0] = dd;
    Cudd_zddDumpDot(gbm, 1, ddnodearray, NULL, NULL, outfile); // dump the function to .dot file
    free(ddnodearray);
    fclose (outfile); // close the file */
}

void SplitStringIntoVector (const std::string& line, std::vector<int>& learnt){
  std::stringstream ss(line);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> vstrings(begin, end);
  std::transform(vstrings.begin(), vstrings.end(), std::back_inserter(learnt),
               [](const std::string& str) { return std::stoi(str); });
}



int main (int argc, char *argv[])
{
    clock_t tStart = clock();
    clock_t tzddstart;
    clock_t tzddend;
    clock_t tzddtotal = 0;
    clock_t clausestart;
    clock_t clauseend;
    clock_t clausetotal = 0;
    clock_t unionstart;
    clock_t unionend;
    clock_t uniontotal = 0;
    clock_t crzddvarsstart;
    clock_t crnegzddvarsstart;
    clock_t crzddvarsend;
    std::string inputfilename;
    int maxvar = 0;
    int max_lit;
    std::vector<int> learnt;
    std::vector<int> learntemp;
    int litcount = 0;
    std::vector<std::vector<int>> learnts;
    std::vector<std::vector<int>> sets;
    std::vector<clock_t> uniontime;
    char filename[30];
    std::vector<int> variables;
    std::unordered_map<int, int> clause_hasher;

    if (argc > 1){
      inputfilename = argv[1];
    }else{
      inputfilename = "test.cnf";
    }
    std::cout << " CNF " << inputfilename <<std::endl;

    std::ifstream file(inputfilename.c_str(), std::ios::in);
    if(!file.is_open())
    {
      throw std::runtime_error(std::string("Can't open file ") + inputfilename);
      std::cout << " Can't open file " << inputfilename  <<std::endl;
    }
    learnts.clear();
    while(file.good() && !file.eof())
    {
      learnt.clear();
      std::string line;
      std::getline(file, line);
      if (line[0] != 'c' && line[0] != 'p' )
      {
        SplitStringIntoVector(line,learnt);
        if(!learnt.empty())
        {
          for (int i = 0; i < learnt.size(); ++i)
          {
            learnt[i] = learnt[i] * (-1);
          }
          learnt.pop_back();
          std::sort(learnt.begin(), learnt.end(), [](int o, int p) { return abs(o) < abs(p); });
          litcount += learnt.size();
          learnts.push_back(learnt);
          max_lit = *std::max_element(learnt.begin(), learnt.end(), abs_compare);
          maxvar = abs(std::max(abs(max_lit),abs(maxvar)));
        }
      }
    }
    file.close();

    std::cout << " Count of literals: " << litcount <<std::endl;
    std::cout << " Maximum number of x: " << maxvar << std::endl;
    std::cout << " Count of disjuncts: " << learnts.size() <<std::endl;
    std::cout << " ZDD start " <<std::endl;


    std::vector<clock_t> clausetime (learnts.size());
    
    for (int j = 0; j < learnts.size(); ++j) //making sets
    {
      learntemp.clear();
      for (int i = 0; i < learnts[j].size(); ++i)
      {
        if (learnts[j][i]>0)
        {
          learntemp.push_back((2*learnts[j][i])-1);
        }else{
          learntemp.push_back(((-2)*learnts[j][i]));
        }
      }

      std::sort(learntemp.begin(), learntemp.end(), [](int o, int p) { return abs(o) < abs(p); });
      sets.push_back(learntemp);
    }

    for (int i = 0; i < sets.size(); ++i) //hash_table
    {
      for (int j = 0; j < sets[i].size(); ++j)
      {
        auto it = clause_hasher.find(sets[i][j]);
        if(it != clause_hasher.end())
        {
          it->second++;
        }else{
          clause_hasher[sets[i][j]] = 1;
        }
      }
    }

    std::vector<std::pair<int, int> > mapcopy(clause_hasher.begin(), clause_hasher.end());
    std::sort(mapcopy.begin(), mapcopy.end(), less_second<int, int>()); // sorting vars by frequency

    std::vector<std::pair<int,int>> mapvars ; //map for reordering vars
    for (int i = 0; i < mapcopy.size(); ++i)
    {
      std::pair<int,int> tmpmap;
      tmpmap.first = i+1;
      tmpmap.second = mapcopy[i].first;
      mapvars.push_back(tmpmap);
    }

    tzddstart = clock();
    std::cout << " Start renaming vars" << std::endl;
    for (int i = 0; i < sets.size(); ++i)
    {
      for (int j = 0; j < sets[i].size(); ++j)
      {
        for (int k = 0; k < mapvars.size(); ++k)
        {
          if (abs(sets[i][j]) == mapvars[k].second)
          {
            if (sets[i][j]>0)
            {
              sets[i][j] = mapvars[k].first;
            }else{
              sets[i][j] = (-1)*mapvars[k].first;
            }
            break;
          }
        }
      }
    }
    for (int i = 0; i < sets.size(); ++i)
    {
      for (int j = 1; j <= mapvars.size(); j++){  //literal addition
        if(!(std::find(sets[i].begin(), sets[i].end(), j) != sets[i].end()) && !(std::find(sets[i].begin(), sets[i].end(), (-1)*j) != sets[i].end())) {
          sets[i].push_back(j*(-1));
        }
      }
    }

    DdManager *gbm; // Global BDD manager.
    gbm = Cudd_Init(0,mapvars.size(),CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0); // Initialize a new BDD manager. 
    DdNode *bdd, *disj_zdd, *var, *var1, *var2, *var3, *zdd, *f, *tmp_zdd, *tmp, *tmp_neg, *clause;

    std::cout << " ZDD variables count: " << Cudd_ReadZddSize(gbm) << std::endl;
    crzddvarsstart = clock();
    std::vector<DdNode *> zddvars (mapvars.size());
    for (int i = 0; i < mapvars.size(); ++i)
    {
      zddvars[i] = Cudd_zddIthVar(gbm, i);
      Cudd_Ref(zddvars[i]);
    }
    crnegzddvarsstart = clock();
    std::vector<DdNode *> zddnegvars (mapvars.size());
    for (int i = 0; i < mapvars.size(); ++i)
    {
      zddnegvars[i] = Cudd_zddChange(gbm, zddvars[i], i);
      Cudd_Ref(zddnegvars[i]);
    }
    std::cout << " ZDD variables created" << std::endl;
    crzddvarsend = clock();
    
    for (int j = 0; j < sets.size(); ++j) //zdd for all disjuncts, one by one
    {
      clausestart = clock();

      std::queue<DdNode *> clauseq;
      for (int i = 0; i < sets[j].size(); i++) 
      {
          if (sets[j][i]<0){
            clauseq.push(zddnegvars[abs(sets[j][i])-1]); 
          }else{
            clauseq.push(zddvars[abs(sets[j][i])-1]); 
          }
      }
      while (clauseq.size()>1)
      {
        var1 = clauseq.front();
        clauseq.pop();
        var2 = clauseq.front();
        clause = Cudd_zddIntersect(gbm, var1, var2);
        Cudd_Ref(clause);
        clauseq.pop();
        clauseq.push(clause);
      }
      clausetime[j] = clock() - clausestart;
      clausetotal += clausetime[j];
      unionstart = clock();
      if (j == 0)
      {
        zdd = clauseq.front();
        clauseq.pop();
        Cudd_Ref(zdd);
      }else{
        zdd = Cudd_zddUnion(gbm, clauseq.front(), zdd);
        clauseq.pop();
        Cudd_Ref(zdd);
      }
      unionend = clock();
      uniontotal += unionend - unionstart;
    }
    tzddend = clock();
    std::string outfilename = "./zdd/out_" + inputfilename;
    for (int i = 0; i < 3; ++i)
    {
      outfilename.erase( outfilename.end() - 1 );
    }
    outfilename = outfilename + "dot";

    outfilename.copy(filename,outfilename.size()+1);
    filename[outfilename.size()] = '\0';

    printf("Time taken for creating zdd variables: %.2fs\n", (double)(crnegzddvarsstart - crzddvarsstart)/CLOCKS_PER_SEC);
    printf("Time taken for cumputing complements of zdd variables: %.2fs\n", (double)(crzddvarsend - crnegzddvarsstart)/CLOCKS_PER_SEC);
    printf("Time taken for conjunction of all subsets: %.2fs\n", (double)(clausetotal)/CLOCKS_PER_SEC);
    printf("Time taken for union of all sets: %.2fs\n", (double)(uniontotal)/CLOCKS_PER_SEC);
    printf("Time taken for zdd: %.2fs\n", (double)(tzddend - tzddstart)/CLOCKS_PER_SEC);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    //pr = 0 : prints nothing
    //pr = 1 : prints counts of nodes and minterms
    //pr = 2 : prints counts + disjoint sum of products
    //pr = 3 : prints counts + list of nodes
    //pr > 3 : prints counts + disjoint sum of products + list of nodes 
    print_dd (gbm, zdd, 2, 1); /*Print the dd to standard output*/
    //sprintf(filename, "./zdd/graph.dot"); /*Write .dot filename to a string*/
    //write_dd(gbm, zdd, filename);  /*Write the resulting cascade dd to a file*/
    Cudd_Quit(gbm);
    return 0; 
}
