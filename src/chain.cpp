#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <random>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include "chaincmdline.h"
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <chrono>

//#define NDEBUG

#define MAXDEGREE 50 //max number of precincts incident with any precinct.
using namespace std;


//GLOBAL VARIABLES
int g_NUMDISTRICTS;  //number of districts in the state
char g_debuglevel;   //debug output level
////////////////////////



struct edge{
  int u;
  int j;
  edge(int a, int b){
    u=a;
    j=b;
  }
  edge(){}
  operator int () {   //assigns an unique integer to an edge
    return MAXDEGREE*this->u+j;
  }
};


void Tokenize( const char* szString, vector<std::string>&
	       vecstrTokens, const char* szSeparators, bool fNoEmpties){
  const char* pc;
  string strCur;
  bool fPush;

  if( !( pc = szString ) )
    return;

  fPush = false;
  while( true ) {
    strCur.clear( );
    if( fNoEmpties )
      for( ; *pc && strchr( szSeparators, *pc ); ++pc );
    if( !*pc ) {
      if( !fNoEmpties && fPush )
	vecstrTokens.push_back( strCur );
      return; }
    for( ; *pc && !strchr( szSeparators, *pc ); ++pc )
      strCur += *pc;
    if( (fPush = !!*pc) )
      pc++;
    vecstrTokens.push_back( strCur ); 
  } 
}


class precinct{
  bool initialized;
public:
  bool selfinitialized;
  int degree;
  int * neighbors; //list of INDICES (in pr[]) of neighbors of the precinct
  int * self; //my place in each of my neighbor's neighbor lists
  double * shared_perimeters; //list of shared perimeters with the neighbors
  double area;
  int population;
  int voteA;
  int voteB;
  int current_district;
  int original_district;
  string county;
  bool frozen;
  bool giant; //Am I in the intial component of the original district.
  friend void computeselves();

  precinct(){
    initialized=false;
    selfinitialized=false;
    giant=false;
    frozen=false;
  }

  precinct(int d){
    degree=d;
    neighbors=new int[degree];
    self=new int[degree];
    shared_perimeters= new double[degree];
    initialized=selfinitialized=true;
    frozen=false;
    giant=false;
  }

  precinct(vector<string> dataline, bool flip, bool use_counties){
    initialized=true;
    selfinitialized=false;
    vector<string> neighbors_string;
    vector<string> shared_perimeters_string;
    Tokenize(dataline[1].c_str(),neighbors_string,",",true);
    Tokenize(dataline[2].c_str(),shared_perimeters_string,",",true);
    degree=neighbors_string.size();
    if (shared_perimeters_string.size()!=degree){
      cerr << "ERROR: perimeter list differs in length from neighbor list?" << endl;
      cerr << "neighbor list was: "<< dataline[1] << endl; 
      exit(-1);
    }
    neighbors=new int[degree];
    shared_perimeters= new double[degree];
    
    for (int i=0; i<degree; i++){
      neighbors[i]=atoi(neighbors_string[i].c_str()); //subtract one because of indexing from 0
      shared_perimeters[i]=atof(shared_perimeters_string[i].c_str());
    }
    area=atof(dataline[4].c_str());
    population=atoi(dataline[5].c_str());
    if (flip){
      voteA=atoi(dataline[7].c_str());
      voteB=atoi(dataline[6].c_str());
    }
    else{
      voteA=atoi(dataline[6].c_str());
      voteB=atoi(dataline[7].c_str());
    }
    
    original_district=atoi(dataline[8].c_str())-1; //subtract one because of indexing from 0
    current_district=original_district;
    if (use_counties)
      county=dataline[9];
    frozen=false;
    giant=false;
    if (current_district>=g_NUMDISTRICTS){
      cerr << "ERROR: district number out of range"<<endl;
      exit(-1);
    }
  }

  precinct(const precinct& source){
    initialized=source.initialized;
    selfinitialized=source.selfinitialized;
    degree=source.degree;
    area=source.area;
    population=source.population;
    voteA=source.voteA;
    voteB=source.voteB;
    current_district=source.current_district;
    original_district=source.original_district;
    county=source.county;
    frozen=source.frozen;
    giant=source.giant;
    if (initialized){
      neighbors=new int[degree];
      shared_perimeters=new double[degree];
      for (int i=0; i<degree; i++){
	neighbors[i]=source.neighbors[i];
	shared_perimeters[i]=source.shared_perimeters[i];
      }
    }
    if (selfinitialized){
      self=new int[degree];
      for (int i=0; i<degree; i++)
	self[i]=source.self[i];
    }
  }
  precinct& operator= (const precinct& source){
    if (&source!=this){
      initialized=source.initialized;
      selfinitialized=source.selfinitialized;
      degree=source.degree;
      area=source.area;
      population=source.population;
      voteA=source.voteA;
      voteB=source.voteB;
      current_district=source.current_district;
      original_district=source.original_district;
      county=source.county;
      frozen=source.frozen;
      giant=source.giant;
      if (initialized){
	neighbors=new int[degree];
	shared_perimeters=new double[degree];
	for (int i=0; i<degree; i++){
	  neighbors[i]=source.neighbors[i];
	  shared_perimeters[i]=source.shared_perimeters[i];
	}
      }
    if (selfinitialized){
      self=new int[degree];
      for (int i=0; i<degree; i++)
	self[i]=source.self[i];
    }
    }
    return *this;
  }
  ~precinct(){
    if (initialized){
      delete[] neighbors;
      delete[] shared_perimeters;
    }
    if (selfinitialized)
      delete[] self;
  }
};


void computeselves(precinct * pr, int N){   //find edges back to me 
  //from my neighboring precints.
  //where multiple edges are present between a pair,
  //we don't need or ensure geometric correspondence.  Only 1-1 correspondence.

  int * boundary;//indices of precincts bordering boundary, in clockwise order
  int * bself; //boundary back pointers
  boundary=new int[N]; 
  bself=new int[N];

  for (int i=0; i<N; i++){
    boundary[i]=-1;  
    bself[i]=-1;
  }


  int boundarydegree=0;

  for (int i=0; i<N; i++){
    pr[i].selfinitialized=true;
    pr[i].self=new int[pr[i].degree];
    for (int j=0; j<pr[i].degree; j++){
      if (pr[i].neighbors[j]>-1){
	int setcount=0;
	for (int l=0; l<pr[pr[i].neighbors[j]].degree; l++){
	  if (pr[pr[i].neighbors[j]].neighbors[l]==i && pr[pr[i].neighbors[j]].shared_perimeters[l]==pr[i].shared_perimeters[j]){
	    setcount++;
	    pr[i].self[j]=l;
	  }
	}
	assert(setcount<2);
	if (setcount!=1){
	  cerr << "setcount is " << setcount<<endl;
	  cerr << "i is "<< i<<endl;
	  cerr << "j is "<< j<<endl;
	  cerr << "pr[i].neighbors[j] is " << pr[i].neighbors[j] << endl;
	}
	assert(setcount==1);
      }
      else{
	boundarydegree++;
	boundary[-pr[i].neighbors[j]-1]=i;  
	bself[-pr[i].neighbors[j]-1]=j;

	pr[i].self[j]= (-pr[i].neighbors[j]-1); //outside edges are labeled by order
	pr[i].neighbors[j]=N; //reset pointer to boundary precinct
      }
    }
  }
  
  assert(boundary[boundarydegree]<0);
  assert(boundary[boundarydegree-1]>=0);

  assert(bself[boundarydegree]<0);
  assert(bself[boundarydegree-1]>=0);


  pr[N]=precinct(boundarydegree);
  for (int j=0; j<pr[N].degree; j++){
    pr[N].neighbors[j]=boundary[j];
    pr[N].self[j]=bself[j];
    pr[N].shared_perimeters[j]=pr[pr[N].neighbors[j]].shared_perimeters[ pr[N].self[j] ];
  }
  pr[N].original_district=pr[N].current_district=-1; //dummy district
}

     //CODE to check for counties which belong to only one district and set
     //their precincts to frozen:
void freezedistrictsbycounty(precinct * pr, int N){   
  unordered_map<string,int> countyfreeze_map; //unique dists for frozen counties
  for (int i=0; i<N; i++){
    if (countyfreeze_map.count(pr[i].county)==0)
      countyfreeze_map[pr[i].county]=pr[i].current_district;
    else if (countyfreeze_map[pr[i].county]!=pr[i].current_district) 
      countyfreeze_map[pr[i].county]=-1; //not preserved
  }
    
  for (auto keyvaluepair : countyfreeze_map){
    //    cout << "first is " << keyvaluepair.first<< " and second is "<<keyvaluepair.second<<endl;
    if (keyvaluepair.second!=-1)
      cout << keyvaluepair.first << " County is preserved by initial districting"<<endl;
  }
    
  for (int i=0; i<N; i++){
    if (countyfreeze_map[pr[i].county]!=-1){
      assert(countyfreeze_map[pr[i].county]==pr[i].current_district);
      pr[i].frozen=true;
    }
    else
      assert(pr[i].frozen==false);
  }
}





template<class TYPE>  //type must be uniquely castable to (int) 
class rdpile{ //fast deletion of specific elements and fast sampling of random elements
private:
  TYPE * array;
  int max; //max number of TYPES to be stored
  int maxindex; //max size of (int) TYPE to be seen
  int * index; //index of locations
public:
  int count; //number of elements

  TYPE access(int i){
    assert (i>=0 && i<count);
    return array[i];
  }
  void insert(TYPE e){
    index[(int) e]=count;
    array[count++]=e;
  }
  void removeindex(int i){
    assert (i>=0 && i<count);
    array[i]=array[--count];
    index[(int) array[i]]=i;
  }
  void remove(TYPE e){
    int i=index[(int) e];
    assert (i>=0 && i<count);
    array[i]=array[--count];
    index[(int) array[i]]=i;
  }
  rdpile(int m, int M){
    max=m;
    maxindex=M;
    array=new TYPE[max];
    index=new int[maxindex];
    count=0;
  }
  rdpile(const rdpile& source){
    max=source.max;
    maxindex=source.maxindex;
    count=source.count;
    array=new TYPE[source.max];
    index=new int[source.maxindex];
    for (int i=0; i<count; i++){
      array[i]=source.array[i];
      index[(int) array[i]]=source.index[(int) array[i]];
    }
  }
  rdpile& operator= (const rdpile& source){
    if (&source!=this){
      max=source.max;
      maxindex=source.maxindex;
      count=source.count;
      array=new TYPE[source.max];
      index=new int[source.maxindex];
      for (int i=0; i<count; i++){
	array[i]=source.array[i];
	index[(int) array[i]]=source.index[(int) array[i]];
      }
    }
    return *this;
  }
  ~rdpile(){
    delete[] array;
    delete[] index;
  }
};



bool validpop(int pop, double avgpop, double popthresh){
  return pop<avgpop*(1+popthresh) && pop>avgpop*(1-popthresh);
}


double dcompact(double Area, double Perim){  //inverse of Polsby-Popper: bigger is worse
  return (pow(Perim,2)/(Area*4*M_PI));
}

double compactsum(double * Area, double * Perim, int arraylength){
  double sum=0;
  for (int i=0; i<arraylength; i++)
    sum+=dcompact(Area[i],Perim[i]);
  return sum;
}

double compactsum(int * Area, double * Perim, int arraylength){
  double sum=0;
  for (int i=0; i<arraylength; i++)
    sum+=dcompact(Area[i],Perim[i]);
  return sum;
}

double compactL2(double * Area, double * Perim, int arraylength){
  double sum=0;
  for (int i=0; i<arraylength; i++)
    sum+=pow(dcompact(Area[i],Perim[i]),2);
  return pow(sum,.5);
}

double compactL2(int * Area, double * Perim, int arraylength){
  double sum=0;
  for (int i=0; i<arraylength; i++)
    sum+=pow(dcompact(Area[i],Perim[i]),2);
  return pow(sum,.5);
}

double arraysum(double * array, int arraylength){
  double sum=0;
  for (int i=0; i<arraylength; i++)
    sum+=array[i];
  return sum;
}


double arrayL2(double * array, int arraylength){
  double sum=0;
  for (int i=0; i<arraylength; i++)
    sum+=pow(array[i],2);
  return pow(sum,.5);
}


int reps(double * Ashare){
  int reps=0;
  for (int k=0; k<g_NUMDISTRICTS; k++){
    if (Ashare[k]>.5)
      reps++;
  }
  return reps;
}

double piecewise_reps(double * Ashare){
  double reps=0;
  for (int k=0; k<g_NUMDISTRICTS; k++){
    if (Ashare[k]>=0.8) {
      reps += Ashare[k];
    }
    else if (Ashare[k]>=0.65) {
      reps += 8.0/3 * Ashare[k] - 4.0/3;
    }
    else if (Ashare[k]>=0.60) {
      reps += 7 * Ashare[k] - 4.15;
    }
    else if (Ashare[k]>=0.50) {
      reps += 0.5 * Ashare[k] - 0.25;
    }
    // implict 0 if x <0.5
  }
  return reps;
}

double smoothed_piecewise_reps(double * Ashare){
  double reps=0;
  for (int k=0; k<g_NUMDISTRICTS; k++){
    if (Ashare[k]>=0.7) {
      reps += 1;
        }
    else if (Ashare[k]>=0.55) {
      reps += (5.0 / 3) * Ashare[k] - 1.0 / 6;
    }
    else if (Ashare[k]>=0.45) {
      reps += 5 * Ashare[k] - 2;
    }
    else if (Ashare[k]>=0.3) {
      reps += 5.0 / 3 * Ashare[k] - 0.5;
        }
    // implict 0 if x <0.3
  }
  return reps;
}

int get_bin_width_and_count(int steps, double minval, double maxval, double * bin_width) {
  // https://stats.stackexchange.com/a/862
  //  Naive IQR (interquartile range) as just 25% - 75% of total range
  double rng = maxval - minval;
  *bin_width = (2 * (0.75*rng - 0.25*rng) / pow(steps, 0.33333));
  return rng / *bin_width;
}

double clip(double n, double lower, double upper) {
  return std::max(lower, std::min(n, upper));
}

int attr_bin(double attr_value, double bin_width, int bins) {
  return clip(floor(attr_value / bin_width), 0, bins);
}

double variance(double * Ashare)
{
  double meansq=0;
  double mean=0;
  for (int k=0; k<g_NUMDISTRICTS; k++){
    meansq+=pow(Ashare[k],2);
    mean+=Ashare[k];
  }
  meansq=meansq/g_NUMDISTRICTS;
  mean=mean/g_NUMDISTRICTS;
  return meansq-pow(mean,2);
}

double efficiency_gap(int * DvotesA, int * DvotesB)
{
  int total=0;
  for (int k=0; k<g_NUMDISTRICTS; k++){
    total+=DvotesA[k];
    total+=DvotesB[k];
  }
  int wastedA=0;
  int wastedB=0;
  for (int k=0; k<g_NUMDISTRICTS; k++){
    if (DvotesA[k]>DvotesB[k]){
      wastedA+=(DvotesA[k]-DvotesB[k])/2;
      wastedB+=DvotesB[k];
    }
    else{
      wastedB+=(DvotesB[k]-DvotesA[k])/2;
      wastedA+=DvotesA[k];
    }
  }
  return ((double) wastedA-wastedB)/total;
}

double median_mean(double * Ashare)
{
  double * AshareCpy;
  AshareCpy = new double[g_NUMDISTRICTS];
  copy(Ashare, Ashare+g_NUMDISTRICTS,AshareCpy);
  double mean,median,sum;
  int middle=(g_NUMDISTRICTS+1)/2; //top median
  nth_element(AshareCpy, AshareCpy+middle, AshareCpy+g_NUMDISTRICTS);
  nth_element(AshareCpy, AshareCpy+middle-1, AshareCpy+middle);
  if (!(g_NUMDISTRICTS%2))//for median with even g_NUMDISTRICTS
    median=(AshareCpy[middle]+AshareCpy[middle-1])/2;
  else  //median for odd g_NUMDISTRICTS
    median=AshareCpy[middle-1];
  sum=0;
  for (int k=0; k<g_NUMDISTRICTS; k++)
    sum+=AshareCpy[k];
  mean=sum/g_NUMDISTRICTS;
  delete [] AshareCpy;
  return mean-median;
}

void tosvg(char* svgfilename, char* inputsvgfilename, precinct* pr, int firstline, int N)
{
    ofstream districtingsvg(svgfilename);
    ifstream inputsvg(inputsvgfilename);
    if (!districtingsvg.good() || !inputsvg.good()){
    cerr << "ERROR with one of the svg files"<<endl;
    exit(-1);
    }
    string correct_header ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");


    mt19937_64 gen(314159265358979); //reinitialize so colors are consistent on each run
    std::uniform_int_distribution<int> randpercent(0,100);
    string * color;
    color=new string[g_NUMDISTRICTS];
    for (int i=0; i<g_NUMDISTRICTS; i++){
      int red=randpercent(gen);
      int green=randpercent(gen);
      int blue=randpercent(gen);
      char tmpcolor[25];
      sprintf(tmpcolor,"rgb(%d%%,%d%%,%d%%);",red,blue,green);
      color[i]=string(tmpcolor);
    }
    string line;
    getline (inputsvg,line); //file header
    districtingsvg << line << endl;
    if (correct_header.compare(line)){
      cerr << "ERROR: incorrect file header at pos " << correct_header.compare(line) << endl;
      exit(-1);
    }

    for (int i=0; i<firstline-2; i++){
      getline (inputsvg,line);
      districtingsvg << line << endl;
    }

    for (int i=0; i<N; i++){
      assert(pr[i].current_district>=0);
      assert(pr[i].current_district<g_NUMDISTRICTS);
      getline (inputsvg,line); //color line ignored
      districtingsvg << color[pr[i].current_district]<<endl;
      getline (inputsvg,line); 
      districtingsvg << line << endl;    
    }
    while (getline (inputsvg,line))
      districtingsvg << line;    
    districtingsvg.close();
    inputsvg.close();
}


void tofile(char* filename, precinct* pr, int N, bool use_counties)
{
  ofstream output(filename);
  if (!output.good()){
    cerr << "ERROR with precinct output file"<<endl;
    exit(-1);
  }
  if (use_counties)
    output << "precinctlistv02" <<endl;
  else
    output << "precinctlistv01" <<endl;
  output << N << endl;
  output << "\tnb\tsp\tunshared\tarea\tpop\tvoteA\tvoteB\tcongD"<<endl;
  for (int i=0; i<N; i++){
    output << i;

    output << '\t';

    for (int j=0; j<pr[i].degree-1; j++){
      if (pr[i].neighbors[j]<N)
	output << pr[i].neighbors[j] << ',';
      else
	output << (-pr[i].self[j]-1) << ',';
    }
    if (pr[i].neighbors[pr[i].degree-1]<N)
      output << pr[i].neighbors[pr[i].degree-1];
    else
      output << (-pr[i].self[pr[i].degree-1]-1);

    output << '\t';

    for (int j=0; j<pr[i].degree-1; j++)
      output << pr[i].shared_perimeters[j] << ',';
    output << pr[i].shared_perimeters[pr[i].degree-1];
    output << '\t';
    output << 0;
    output << '\t';
    output << pr[i].area;
    output << '\t';
    output << pr[i].population;
    output << '\t';
    output << pr[i].voteA;
    output << '\t';
    output << pr[i].voteB;
    output << '\t';
    output << pr[i].current_district+1;
    if (use_counties){
      output << '\t';
      output << pr[i].county;
    }
    output << endl;
  }
  output.close();
}
 



void chainstep(precinct * pr, rdpile<edge> & edgeset, edge e, int * DvotesA, int * DvotesB, double * Ashare, int * Dpop, double * Dareas, double * Dperims){  //add e.u(e.j) to district of e.u
  int Du=pr[e.u].current_district;
  int v=pr[e.u].neighbors[e.j];
  int Dv=pr[v].current_district;
  if (Du==Dv){
    cout <<"WTF! Du is "<<Du<<" and Dv is "<<Dv<<endl;
    exit(-1);
  }
  assert (Dv>=0); //not adding outside
  assert (Du!=Dv); 

  DvotesA[Du]+= pr[v].voteA;
  DvotesA[Dv]-= pr[v].voteA;
  DvotesB[Du]+= pr[v].voteB;
  DvotesB[Dv]-= pr[v].voteB;
  Ashare[Du]=((double) DvotesA[Du])/(DvotesA[Du]+DvotesB[Du]);
  Ashare[Dv]=((double) DvotesA[Dv])/(DvotesA[Dv]+DvotesB[Dv]);

  Dpop[Du]+= pr[v].population;
  Dpop[Dv]-= pr[v].population;

  Dareas[Du]+= pr[v].area;
  Dareas[Dv]-= pr[v].area;

  for (int l=0; l<pr[v].degree; l++){

    if (pr[v].neighbors[l]<0){  //outside
      Dperims[Du]+=pr[v].shared_perimeters[l]; //because it avoids Du
      Dperims[Dv]-=pr[v].shared_perimeters[l]; //because it avoids Dv


    }
    else if (pr[pr[v].neighbors[l]].current_district==Du){
      if (!pr[pr[v].neighbors[l]].frozen){
	edgeset.remove(edge(v,l));
	edgeset.remove(edge(pr[v].neighbors[l],pr[v].self[l]));
      }
      Dperims[Du]-=pr[v].shared_perimeters[l]; //because it intersects Du
      Dperims[Dv]-=pr[v].shared_perimeters[l]; //because it avoids Dv
    }
    else if (pr[pr[v].neighbors[l]].current_district==Dv){
      if (!pr[pr[v].neighbors[l]].frozen){
	edgeset.insert(edge(v,l));
	edgeset.insert(edge(pr[v].neighbors[l],pr[v].self[l]));
      }
      Dperims[Du]+=pr[v].shared_perimeters[l]; //because it avoids Du
      Dperims[Dv]+=pr[v].shared_perimeters[l]; //because it intersects Dv
    }
    else{ //avoids Du AND Dv
      Dperims[Du]+=pr[v].shared_perimeters[l]; //because it avoids Du
      Dperims[Dv]-=pr[v].shared_perimeters[l]; //because it avoids Dv
    }
  }

  pr[v].current_district=Du;

}


class HistogramPopulator{
private:
  double m_bin_width;
  double m_minval;
  double m_maxval;
  int m_bins;
  std::vector<int64_t> m_histogram;
public:

  HistogramPopulator(int64_t steps, double m_minval, double m_maxval) {
    m_bins = get_bin_width_and_count(steps, m_minval, m_maxval, &m_bin_width);
    m_histogram = std::vector<int64_t>(m_bins);
  }

  ~HistogramPopulator(){}

  int bins() {
    return m_bins;
  }
  double bin_width() {
    return m_bin_width;
  }
  void add_value(double value, int64_t increment){
    m_histogram[attr_bin(value, m_bin_width, m_bins)] += increment;
  }
  void save_file(std::string fname) {
    ofstream myfile;
    myfile.open(fname);
    myfile << "# nbins, bin_start, bin_width\n";
    myfile << m_bins << ", " << m_minval << ", "<< m_bin_width << "\n";
    myfile << "# counts in each bin\n";

    for (int j=0; j<=m_bins; j++){
      myfile << m_histogram[j] << ", ";
    }
    myfile << "\n";
    myfile.close();
  }
};



int main(int argc, char* argv[])
{
  gengetopt_args_info lineArgs;
  if (cmdline_parser(argc, argv, &lineArgs)) {
    cmdline_parser_print_help();
    return 1;
  }

  g_NUMDISTRICTS=lineArgs.numdists_arg;
  
  int numberfrozen=lineArgs.freeze_given;
  bool * frozen_districts;  //which districts are frozen
  frozen_districts=new bool[g_NUMDISTRICTS];

  for (int k=0; k<g_NUMDISTRICTS; k++)
    frozen_districts[k]=0;                       //zero them

  for (int i=0; i<numberfrozen; i++)
    frozen_districts[lineArgs.freeze_arg[i]-1]=1;//set to one the given frozen 
                                                 //(indexing from 0)
  
  double perimthresh=9999999999;
  double popperthresh=9999999999;
  double L1thresh=9999999999;
  double L2thresh=9999999999;
  bool PerimTest=false;
  bool L1Test=false;
  bool L2Test=false;
  bool PopperTest=false;

  std::chrono::time_point<std::chrono::system_clock> start_time = std::chrono::system_clock::now();
  float target_time_min=999999999;

  if (lineArgs.perimeter_given){
    PerimTest=true;
    perimthresh=lineArgs.perimeter_arg;
  }
  if (lineArgs.polsby_popper_given){
    PopperTest=true;     //district-by-district compactness
    popperthresh=lineArgs.polsby_popper_arg; //actually inverse of Polsby-Popper value
  }
  if (lineArgs.L1_compactness_given){
    L1Test=true;     
    L1thresh=lineArgs.L1_compactness_arg;
  }
  if (lineArgs.L2_compactness_given){
    L2Test=true;     
    L2thresh=lineArgs.L2_compactness_arg;
  }
  if (lineArgs.target_time_given){
    target_time_min = lineArgs.target_time_arg;
  }
  

  g_debuglevel=0;

  int N; //num precincts
  precinct * pr; //array of precincts


  ifstream myfile (lineArgs.filename_arg);
  string line;
  if (!myfile.good()){
    cerr << "ERROR with file"<<endl;
    exit(-1);
  }
  bool use_counties=false;
  string old_header ("precinctlistv01");
  string county_header ("precinctlistv02");
  getline (myfile,line); //file header
  if (old_header.compare(line)){
    if (county_header.compare(line)){
      cerr << "ERROR: incorrect file header at pos " << county_header.compare(line) << endl;
    exit(-1);
    }
    use_counties=true;
  }
  if (lineArgs.counties_flag && ! use_counties){
    cerr << "ERROR: preserving preserved counties required a file v02 file with county data!"<<endl;
    exit(-1);
  }
  getline(myfile,line); //#precints
  N=atoi(line.c_str());
  pr=new precinct[N+1];  

  
  cout << "We have "<<N<<" precincts."<<endl;
  //upper bound on crossing edges.
  int MAXEDGES=2*(2*N+(g_NUMDISTRICTS-numberfrozen)-6);
  
  getline(myfile,line);  //skip data header
 
  int I=0;
  int sum=0;
  while (getline(myfile,line) && I<N){
    vector<string> dataline;
    Tokenize(line.c_str(), dataline, "\t", true);
    if (dataline.size()!=9+use_counties){
      cerr << "ERROR: shouldn't there be "<<9+use_counties<<" chunks per line?"<<endl;
      cerr << "there are "<<dataline.size()<<" at I="<<I<<endl;
      exit(-1);
    }
    if (atoi(dataline[0].c_str())!=I){
      cerr << "ERROR: line numbers don't match"<<endl;
      exit(-1);
    }
    pr[I]=precinct(dataline,lineArgs.flip_flag,use_counties);

    sum+=pr[I].degree;
    I++;
  }

  // TODO(bojanserafimov): tidy up later
  // TODO(bojanserafimov): test this.
  // Replaces only the election results with information from another file.
  // The election results file is a CSV with N rows and 3 columns.
  // First column is the index, second and third are election results.
  // The flip flag applies.
  if (lineArgs.filename_election_results_given) {
    ifstream file_er(lineArgs.filename_election_results_arg);
    if (!file_er.good()) {
      cerr << "ERROR with election results file"<<endl;
      exit(-1);
    }
    int I=0;
    while (getline(file_er,line)){
      if (I >= N) {
        cerr << "ERROR: election results file is longer than it should be";
        exit(-1);
      }
      vector<string> dataline;
      Tokenize(line.c_str(), dataline, ", ", true);
      if (dataline.size()!=3) {
        cerr << "ERROR: shouldn't there be "<<3<<" chunks per line in election results file?"<<endl;
        cerr << "there are "<<dataline.size()<<" at I="<<I<<endl;
        exit(-1);
      }
      if (atoi(dataline[0].c_str())!=I){
        cerr << "ERROR: line numbers don't match"<<endl;
        cerr << dataline[0].c_str() << " " << I << endl;
        exit(-1);
      }
      if (lineArgs.flip_flag) {
        pr[I].voteA=atoi(dataline[2].c_str());
        pr[I].voteB=atoi(dataline[1].c_str());
      } else {
        pr[I].voteA=atoi(dataline[1].c_str());
        pr[I].voteB=atoi(dataline[2].c_str());
      }
      I++;
    }
    if (I != N) {
      cerr << "ERROR: election results file is shorter than it should be";
      exit(-1);
    }
  }

  // read congD from filename_wes_units
  if (lineArgs.filename_wes_units_given){
    ifstream file_wu(lineArgs.filename_wes_units_arg);
    if (!file_wu.good()) {
      cerr << "ERROR with wes units file"<<endl;
      exit(-1);
    }
    getline(file_wu,line);  //skip data header
    int I = 0;
    while (getline(file_wu, line)) {
      vector<string> dataline;
      Tokenize(line.c_str(), dataline, ",", true);
      if (dataline.size()!=2) {
        cerr << "ERROR: shouldn't there be "<<2<<" chunks per line in b2wid?"<<endl;
        cerr << "there are "<<dataline.size()<<" at I="<<I<<endl;
        exit(-1);
      }

      for (int i = 0; i < dataline.size(); i++) {
        if (dataline[i][0] == '\"') {
          dataline[i] = dataline[i].substr(1, dataline[i].size() - 2);
        }
      }
      int wid = atoi(dataline[0].c_str());
      int district = atoi(dataline[1].c_str()) - 1;
      pr[wid].original_district = district;
      pr[wid].current_district = district;
      I++;
    }
  }
 
  if (I<N || !myfile.eof() ){
    cerr << "ERROR: Line numbers don't match precinct count"<<endl;
    exit(-1);
  }

  if (lineArgs.counties_flag && use_counties)
    freezedistrictsbycounty(pr,N);

  //create backneighbor lists:
  computeselves(pr,N);


  int64_t * reps_histogram;
  reps_histogram = new int64_t[g_NUMDISTRICTS+1];

  for (int i=0; i<=g_NUMDISTRICTS; i++) {
    reps_histogram[i]=0;
  }

  int ** adjM;                                       //precinct adjacency matrix
  rdpile<edge> edgeset(MAXEDGES,MAXEDGES*MAXDEGREE);                    //set of edges on district boundaries

  adjM=new int*[N];
  for (int i=0; i<N; i++){
    adjM[i]=new int[N];
    for (int j=0; j<N; j++)
      adjM[i][j]=0;
  }

  for (int i=0; i<N; i++){
    for (int j=0; j<pr[i].degree; j++){
      int idistrict=pr[i].current_district;
      int jdistrict=pr[pr[i].neighbors[j]].current_district;
      if (pr[i].neighbors[j]<N){  //not state boundary connection
	adjM[i][pr[i].neighbors[j]]++;
	if (idistrict!=jdistrict &&
	    frozen_districts[idistrict]==0 &&
	    frozen_districts[jdistrict]==0 &&
	    !pr[i].frozen &&
	    !pr[pr[i].neighbors[j]].frozen
	    ){
	  edgeset.insert(edge(i,j)); 
	}
      }
    }
  }
  cout << "There are "<<edgeset.count<<"/2 boundary edges"<<endl;
  assert(edgeset.count%2==0);



  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (adjM[i][j]!=adjM[j][i]){
	cerr << "ERROR: graph is not undirected!"<< endl;
	cerr << "bad pair is "<<i<<","<<j<<"."<<endl;
	exit(-1);
      }
    }
  }

  /*     OBSOLETE FROM CONNECTIVITY CHECK BELOW
  for (int i=0; i<N; i++){
    bool isolated=true;
    for (int j=0; j<pr[i].degree; j++){
      if (pr[i].neighbors[j]<N && pr[i].current_district==pr[pr[i].neighbors[j]].current_district){
	isolated=false;
	break;
      }
    }
    if (isolated)
    cerr << i <<" is isolated"<<endl;
    assert(!isolated);
    }*/

  int count=0;

  if (g_debuglevel>1)
  cout << "connectivity check..."<<endl;
  for (int k=0; k<g_NUMDISTRICTS; k++){
    int localcount=0;
    queue <int> BFSqueue;
    for (int i=0; i<N; i++){
      if (pr[i].current_district==k){
	BFSqueue.push(i);
	pr[i].giant=true;
	count++; localcount++;
	break;
      }
    }
    while (!BFSqueue.empty()){
      int myindex=BFSqueue.front();
      BFSqueue.pop();
      assert(pr[myindex].giant==true);
      for (int j=0; j<pr[myindex].degree; j++){
	if (pr[pr[myindex].neighbors[j]].current_district==pr[myindex].current_district && pr[pr[myindex].neighbors[j]].giant==false){
	  pr[pr[myindex].neighbors[j]].giant=true; count++; localcount++;
	  BFSqueue.push(pr[myindex].neighbors[j]);
	}
      }
    }
    if (g_debuglevel>1)
      cerr << "district "<<k<< " giant has "<<localcount<<" precincts."<<endl;
  }
  
  bool passed=true;
  for (int i=0; i<N; i++){
    if (pr[i].giant==false){
      passed=false;
      cerr << "ERROR: Precinct "<<i<< " is not in its giant!"<<endl;
    }
  }

  if (passed){
    assert(count==N);
    if (g_debuglevel>1)
      cout << "passed connectivity check" <<endl;
  }
  else
    exit(-1);


  int64_t steps=pow(2,lineArgs.steps_arg);
  int64_t period=pow(2,lineArgs.period_arg);
  double popthresh=lineArgs.poperror_arg;

  int * DvotesA;
  int * DvotesB;
  int * Dpop;
  double * Ashare;

  DvotesA=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
  DvotesB=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
  Ashare= new double[g_NUMDISTRICTS];
  Dpop=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
  double * Dareas;
  double * Dperims;
  Dareas=new double[g_NUMDISTRICTS];   //district areas
  Dperims=new double[g_NUMDISTRICTS];  //district perimeters
  for (int k=0; k<g_NUMDISTRICTS; k++){
    DvotesA[k]=DvotesB[k]=Dpop[k]=0;
    Dareas[k]=Dperims[k];
  }


  mt19937_64 gen(314159265358979);

  HistogramPopulator median_mean_hist = HistogramPopulator(steps, 0.0, 0.1);
  HistogramPopulator eg_hist = HistogramPopulator(steps, 0, 0.8);
  HistogramPopulator pp_hist = HistogramPopulator(steps, 0, 0.8);
  HistogramPopulator seats_hist = HistogramPopulator(steps, 0, g_NUMDISTRICTS);
  HistogramPopulator smoothed_seats_hist = HistogramPopulator(steps, 0, g_NUMDISTRICTS);

  //filling in first step...
  int revisitations;
  {
    double p=((double) edgeset.count)/MAXEDGES;
    std::geometric_distribution<int> revisits(p);
    revisitations=1+revisits(gen);
    cout << "first edgeset.count is "<<edgeset.count<<endl;
    cout << "first p is "<<p<<endl;
    cout << "first revisitations is "<<revisitations<<endl;
  }
  int totalpop=0;
  for (int i=0; i<N; i++){
    Dpop[pr[i].current_district]+=pr[i].population;
    totalpop+=pr[i].population;
    DvotesA[pr[i].current_district]+= pr[i].voteA;
    DvotesB[pr[i].current_district]+= pr[i].voteB;
    Dareas[pr[i].current_district]+=  pr[i].area;
    for (int j=0; j<pr[i].degree; j++){
      if (pr[i].neighbors[j]>=0 && pr[i].current_district!=pr[pr[i].neighbors[j]].current_district){
	Dperims[pr[i].current_district]+=pr[i].shared_perimeters[j];
      }
      else if (pr[i].neighbors[j]<0){
	Dperims[pr[i].current_district]+=pr[i].shared_perimeters[j];
      }
    }
  }
  double avgpop = (double) totalpop/g_NUMDISTRICTS;
  int Avotes=0;
  int Bvotes=0;
  for (int k=0; k<g_NUMDISTRICTS; k++){
    Ashare[k]=((double) DvotesA[k])/(DvotesA[k]+DvotesB[k]);
    Avotes+=DvotesA[k];
    Bvotes+=DvotesB[k];
    if (g_debuglevel>0)
      cout << "district "<<k<<" has " <<DvotesA[k]<<" voteA and "<<DvotesB[k]<<" voteB.  share is "<<Ashare[k]<<endl;
  }
  cout << "A has "<<Avotes <<" votes"<<endl;
  cout << "B has "<<Bvotes <<" votes"<<endl;
  cout << "total population is "<<totalpop<<endl;
  bool OutVariance=lineArgs.variance_flag;
  bool OutMedianMean=lineArgs.median_mean_flag;
  bool OutEfficiencyGap=lineArgs.efficiency_gap_flag;
  bool OutSeats=lineArgs.seats_flag;
  bool OutHistogram=lineArgs.histogram_flag;

  double initial_variance=variance(Ashare);
  double initial_median_mean=median_mean(Ashare);
  double initial_efficiency_gap=efficiency_gap(DvotesA,DvotesB);
  int initial_seat_count=reps(Ashare);
  double initial_piecewise_seat_count=piecewise_reps(Ashare);
  double initial_smoothed_seat_count=smoothed_piecewise_reps(Ashare);
  


  if (OutVariance)
    cout << "initial variance is "<<initial_variance<<endl;
  if (OutMedianMean)
    cout << "initial median_mean is "<<initial_median_mean<<endl;
  if (OutEfficiencyGap)
    cout << "initial efficiency_gap is "<<initial_efficiency_gap<<endl;
  if (OutSeats) {
    cout << "initial seat count is "<<initial_seat_count<<endl;
    cout << "initial piecewise seat count is "<<initial_piecewise_seat_count<<endl;
    cout << "initial smoothed seat count is "<<initial_smoothed_seat_count<<endl;
  }
  
  
  //(filling in first step)
  int64_t variance_lessunusual=0;
  int64_t median_mean_lessunusual=0;
  int64_t efficiency_gap_lessunusual=0;
  int64_t seat_count_lessunusual=0;
  
  int64_t variance_moreunusual=revisitations; //itself
  int64_t median_mean_moreunusual=revisitations; //itself
  int64_t efficiency_gap_moreunusual=revisitations; //itself
  int64_t seat_count_moreunusual=revisitations; //itself
  reps_histogram[reps(Ashare)]+=revisitations;
  int64_t totalsteps=revisitations;



  if (PerimTest){
    for (int i=0; i<N; i++)
      if (arraysum(Dperims,g_NUMDISTRICTS)>perimthresh){
	cerr << "ERROR: initial districting violates perimeter criterion.  Current perimeter is "<< arraysum(Dperims,g_NUMDISTRICTS)  <<endl;
	exit(-1);
      }
  }
  if (PopperTest){
    for (int k=0; k<g_NUMDISTRICTS; k++){
      if (dcompact(Dareas[k],Dperims[k])>popperthresh){
	cerr << "ERROR: given district "<<k<<" does not satisfy compactness requirement!\n It has compactness value "<< dcompact(Dareas[k],Dperims[k]) <<"."<<endl;
	exit(-1);
      }
    }
  }
  if (L1Test){
    for (int i=0; i<N; i++)
      if (compactsum(Dareas,Dperims,g_NUMDISTRICTS)>L1thresh){
	cerr << "ERROR: initial districting violates L1 compactness criterion.  Current L1 norm is "<< compactsum(Dareas,Dperims,g_NUMDISTRICTS)  <<endl;
	exit(-1);
      }
  }
  if (L2Test){
    for (int i=0; i<N; i++)
      if (compactL2(Dareas,Dperims,g_NUMDISTRICTS)>L2thresh){
	cerr << "ERROR: initial districting violates L2 criterion.  Current L2 norm is "<< compactL2(Dareas,Dperims,g_NUMDISTRICTS)  <<endl;
	exit(-1);
      }
  }

 
  int outputcount=0;

  //MAIN LOOP
  for (int64_t i=0; i<steps; i++){
    float time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                               std::chrono::system_clock::now()- start_time).count() / 60.0;
    if (time_elapsed > target_time_min) break;

    uniform_int_distribution<> intdist(0,edgeset.count-1);
    int rindex=intdist(gen);
    edge e=edgeset.access(rindex);    //we'll try adding vertex u(j) to u's district
    int Du=pr[e.u].current_district;
    int v=pr[e.u].neighbors[e.j];
    int Dv=pr[v].current_district;

    
    int Dusegmentcount=0;
    int Dvsegmentcount=0;
    int precinct=pr[v].neighbors[0];

    int pointback=(pr[v].self[0]+1) % pr[precinct].degree;

    int lastincycle=pr[precinct].neighbors[pointback];
    int lastprecinct=lastincycle;
    int edgeindex=pr[precinct].self[pointback] ;
    
      
    int newprecinct=-1;
    int testcount=0;

      
    while ( ! ( newprecinct==pr[v].neighbors[0] && pr[lastprecinct].self[edgeindex]==pointback )  ){
      


      assert(precinct==pr[lastprecinct].neighbors[edgeindex]);
      testcount++;
      assert(testcount<100);
      int lastD=pr[lastprecinct].current_district;
      int newD=pr[precinct].current_district;



      
      if (lastD==Du && newD!=Du)
	Dusegmentcount++;
      else if (lastD==Dv && newD!=Dv)
	Dvsegmentcount++;

      edgeindex=(pr[lastprecinct].self[edgeindex]-1+pr[precinct].degree) % pr[precinct].degree;
      newprecinct=pr[precinct].neighbors[edgeindex];
      if (newprecinct==v){ //if I'm a neighbor of v, skip v for next step of cycle
	edgeindex=(edgeindex-1+pr[precinct].degree) % pr[precinct].degree;
	newprecinct=pr[precinct].neighbors[edgeindex];
      }


      lastprecinct=precinct;
      precinct=newprecinct;
    }


      
    assert(Dusegmentcount>0);
    if (Dvsegmentcount==0){
      cerr << "u is " <<e.u<<" and v is " <<v << endl;
    }
    assert(Dvsegmentcount>0);

    if (Dusegmentcount>1 || Dvsegmentcount>1){
      i--;
      continue;
    }

    
    //END CONNECTIVITY CHECK
    chainstep(pr, edgeset, e, DvotesA, DvotesB, Ashare, Dpop, Dareas, Dperims);


    bool LOOP=false;
    if (!validpop(Dpop[Du],avgpop,popthresh) || !validpop(Dpop[Dv],avgpop,popthresh)){
      LOOP=true;
      if (g_debuglevel>1)
	cout << "LOOP! for population violation"<<endl;
    }
    else if (PerimTest && arraysum(Dperims,g_NUMDISTRICTS)>perimthresh){
      LOOP=true;
      if (g_debuglevel>1)
	cout << "LOOP! for Perimeter violation"<<endl;
    }
    else if (PopperTest && (dcompact(Dareas[Du],Dperims[Du])>popperthresh || dcompact(Dareas[Dv],Dperims[Dv])>popperthresh)){
      LOOP=true;
      if (g_debuglevel>1)
	cout << "LOOP! for Polsby-Popper violation"<<endl;
    }
    else if (L1Test && compactsum(Dareas,Dperims,g_NUMDISTRICTS)>L1thresh){
      LOOP=true;
      if (g_debuglevel>1)
	cout << "LOOP! for compactness L1 violation"<<endl;
    }
    else if (L2Test && compactL2(Dareas,Dperims,g_NUMDISTRICTS)>L2thresh){
      LOOP=true;
      if (g_debuglevel>1)
	cout << "LOOP! for compactness L2 violation"<<endl;
    }

    if (LOOP){
      edge rev_e; //edge to reverse the move
      bool foundone=false;
      for (int l=0; l<pr[v].degree; l++){
	if (pr[pr[v].neighbors[l]].current_district==Dv){
	  rev_e=edge(pr[v].neighbors[l],pr[v].self[l]);
	  foundone=true;
	  break;
	}
      }
      assert(foundone);
      chainstep(pr, edgeset, rev_e, DvotesA, DvotesB, Ashare, Dpop, Dareas, Dperims);
    }


    double p=((double) edgeset.count)/MAXEDGES;
    std::geometric_distribution<int> revisits(p);
    revisitations=1+revisits(gen);
    totalsteps+=revisitations;
    if (OutVariance){
      if (variance(Ashare)>=initial_variance)
	variance_moreunusual+=revisitations;
      else
	variance_lessunusual+=revisitations;
    }
    double median_mean_val;
    double eg_val;
    int seats;
    double piecewise_seats;
    double smoothed_seats;

    if (OutMedianMean){
      median_mean_val = median_mean(Ashare);
      if (median_mean_val>=initial_median_mean)
	median_mean_moreunusual+=revisitations;
      else
	median_mean_lessunusual+=revisitations;
    }

    if (OutEfficiencyGap){
      eg_val = efficiency_gap(DvotesA,DvotesB);
      if (eg_val>=initial_efficiency_gap)
	efficiency_gap_moreunusual+=revisitations;
      else
	efficiency_gap_lessunusual+=revisitations;
    }

    if (OutSeats){
      seats = reps(Ashare);
      piecewise_seats = piecewise_reps(Ashare);
      smoothed_seats = smoothed_piecewise_reps(Ashare);
      if (seats<=initial_seat_count)
	seat_count_moreunusual+=revisitations;
      else
	seat_count_lessunusual+=revisitations;
    }

    if (OutHistogram){
      reps_histogram[reps(Ashare)]+=revisitations;

      if (OutMedianMean)
        median_mean_hist.add_value(median_mean_val, revisitations);

      if (OutEfficiencyGap)
        eg_hist.add_value(eg_val, revisitations);

      if (OutSeats){
        seats_hist.add_value(piecewise_seats, revisitations);
        smoothed_seats_hist.add_value(smoothed_seats, revisitations);
      }
    }

    
    if (!(i%period) || i==steps-1){
      
      outputcount++;
      if (lineArgs.stages_flag){
	if (lineArgs.svg_filename_given){
	  char svgfilename[100];
	  sprintf(svgfilename,"%s_%d.svg",lineArgs.svg_filename_arg,outputcount);
	  char inputsvgfilename[100];
	  sprintf(inputsvgfilename,"%s",lineArgs.inputsvg_filename_arg);
	  tosvg(svgfilename,inputsvgfilename,pr,lineArgs.svg_firstline_arg,N);
	}
      
	if (lineArgs.precinct_filename_given){
	  char precinctfilename[100];
	  sprintf(precinctfilename,"%s_%d",lineArgs.precinct_filename_arg,outputcount);
	  tofile(precinctfilename,pr,N,use_counties);
	}
      }

     
      cout << endl << "for i="<<i<<","<<endl;
      if (OutHistogram){
	for (int j=0; j<=g_NUMDISTRICTS; j++){
	  if (reps_histogram[j]>0){
	    cout << setw(2) << j << ": ";
	    for (int k=0; k<(50*reps_histogram[j])/(totalsteps); k++)
	      cout << "*";
	    for (int k=ceil(50*reps_histogram[j])/(totalsteps); k<50; k++)
	      cout << " ";
	    cout <<100*(double) reps_histogram[j]/(totalsteps)<<"%"<<endl;
	  }
	}

  // TODO: here's where histograms get saved
  if (OutMedianMean)
    median_mean_hist.save_file("median_mean_hist.txt");

  if (OutEfficiencyGap)
    eg_hist.save_file("eg_hist.txt");

  if (OutSeats){
    seats_hist.save_file("seats_hist.txt");
    smoothed_seats_hist.save_file("smoothed_seats_hist.txt");
  }

	cout << endl;
      }
      if (g_debuglevel>0){
	cout << "Dshare is ";
	for (int j=0; j<g_NUMDISTRICTS; j++)
	  cout << setprecision(2) << 100*Ashare[j] <<"%, ";
	cout << endl;
	cout << "edge count is "<<edgeset.count<<endl;
      }    
      if (OutVariance){
	cout << "--FOR VARIANCE--"<<endl;
	cout << "variance is "<<setprecision(5)<< variance(Ashare)<<endl;
	double V_ep=((double)variance_moreunusual) / totalsteps;
	double V_pv=sqrt(2*V_ep);
	cout << "moreunusual is " << variance_moreunusual<<endl;
	cout << "lessunusual is " << variance_lessunusual<<endl;
	cout << "ep="<<V_ep<<endl;
	cout << "p="<<V_pv<<endl << endl;
      }
      if (OutMedianMean){
	cout << "--FOR MEDIAN/MEAN--"<<endl;
	cout << "median_mean is "<< setprecision(5)<< median_mean(Ashare) << endl;
	double M_ep=((double)median_mean_moreunusual) / totalsteps;
	double M_pv=sqrt(2*M_ep);
	cout << "moreunusual is " << median_mean_moreunusual<<endl;
	cout << "lessunusual is " << median_mean_lessunusual<<endl;
	cout << "ep="<<M_ep<<endl;
	cout << "p="<<M_pv<<endl << endl;
      }
      if (OutEfficiencyGap){
	cout << "--FOR Efficiency Gap--"<<endl;
	cout << "efficiency_gap is "<< setprecision(5)<< efficiency_gap(DvotesA,DvotesB) << endl;
	double M_ep=((double)efficiency_gap_moreunusual) / totalsteps;
	double M_pv=sqrt(2*M_ep);
	cout << "moreunusual is " << efficiency_gap_moreunusual<<endl;
	cout << "lessunusual is " << efficiency_gap_lessunusual<<endl;
	cout << "ep="<<M_ep<<endl;
	cout << "p="<<M_pv<<endl << endl;
      }
      if (OutSeats){
	cout << "--FOR Seat Count--"<<endl;
	cout << "seat_count is "<< setw(2)<< reps(Ashare) << endl;
	double M_ep=((double)seat_count_moreunusual) / totalsteps;
	double M_pv=sqrt(2*M_ep);
	cout << "moreunusual is " << seat_count_moreunusual<<endl;
	cout << "lessunusual is " << seat_count_lessunusual<<endl;
	cout << "ep="<<M_ep<<endl;
	cout << "p="<<M_pv<<endl << endl;
      }

    }
  }


  if (lineArgs.svg_filename_given){
    char svgfilename[100];
    sprintf(svgfilename,"%s.svg",lineArgs.svg_filename_arg);
    char inputsvgfilename[100];
    sprintf(inputsvgfilename,"%s",lineArgs.inputsvg_filename_arg);
    tosvg(svgfilename,inputsvgfilename,pr,lineArgs.svg_firstline_arg,N);
  }

  if (lineArgs.precinct_filename_given){
    char precinctfilename[100];
    sprintf(precinctfilename,"%s",lineArgs.precinct_filename_arg);
    tofile(precinctfilename,pr,N,use_counties);
  }
}
