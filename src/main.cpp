#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <list>
#include <cstddef>
#include <cstdlib>
#include <algorithm>
#include <vector>
using namespace std;

const double kel = 273.15;
const double g = 9.81;
const double Rd = 287.04;
const double Rv = 461.5;
const double epsilon = Rd/Rv;
const double cp = 1005;
const double xlv = 2501000;
const double xls = 2834000;
const double cpv = 1870;
const double cpl = 4190;
const double cpi = 2106;
const double ttrip = 273.15;
const double eref = 611.2;
const double Gamma_d = g/cp;
const double pref = 611.65;
const double EC_T1 = 273.15;
const double EC_T2 = 253.15;
const double EC_L = 250;
const double sigma = 1.1;
const double Pr = 1.0 / 3.0;
const double alpha = 0.8;
const double ksq = 0.18;
  
class IndicesCollector;
class InfoCollector;
class Sounding;
class Cache;
double alcl(double ts, double tds, double ps, double tmp);

template <typename Typ>

Typ Get(list<Typ> *lista, int index){
  typename std::list<Typ>::iterator it = lista->begin();
  if(index<lista->size())std::advance(it, index);
  return *it;
}

double sign(double x){
  if (x > 0) return 1;
  else if (x < 0) return -1;
  else if (x == 0) return 0;
}

double computeMean(const vector<double>& arr){ 
    double sum = 0.0, mean, standardDeviation = 0.0;   
    int size = arr.size(); 
    for (int i = 0; i < size; ++i) { 
        sum += arr[i]; 
    } 
    mean = sum / size; 
    return mean; 
} 

double computeStandardDeviation(const vector<double>& arr){ 
    double sum = 0.0, mean, standardDeviation = 0.0;   
    int size = arr.size(); 
    for (int i = 0; i < size; ++i) { 
        sum += arr[i]; 
    } 
    mean = sum / size; 
    for (int i = 0; i < size; ++i) { 
        standardDeviation += pow(arr[i] - mean, 2); 
    } 
    return sqrt(standardDeviation / size); 
} 

double ESAT(double t)
{
  double temp = t + kel;
  return (pow(10, 23.832241 - 5.02808 * log10(temp) - 0.00000013816 * pow(10, 11.344 - 0.0303998 * temp) + 0.0081328 * pow(10, 3.49149 - 1302.8844 / temp) - 2949.076 / temp));
}

double W(double t, double p)
{
  return 622 * ESAT(t) / (p - ESAT(t));
}

double heaviside(double x){   
  return (sign(x) + 1)/2;
}

double compute_rsat(double T, double p, double iceflag){
  double qsat = 0;
  double omeg = ((T-EC_T1)/(EC_T2-EC_T1))*heaviside((T-EC_T1)/(EC_T2-EC_T1))*heaviside((1-(T-EC_T1)/(EC_T2-EC_T1))) + heaviside(-(1 - (T-EC_T1)/(EC_T2-EC_T1)));
  if(iceflag==0){
    double term1=(cpv-cpl)/Rv;
    double term2=(xlv-ttrip*(cpv-cpl))/Rv;
    double esl=exp((T-ttrip)*term2/(T*ttrip))*eref*pow((T/ttrip),(term1));
    qsat=epsilon*esl/(p-esl);         
  }
  if(iceflag==1){
    double term1=(cpv-cpl)/Rv;
    double term2=(xlv-ttrip*(cpv-cpl))/Rv;
    double esl_l=exp((T-ttrip)*term2/(T*ttrip))*eref*pow((T/ttrip),(term1));
    double qsat_l=epsilon*esl_l/(p-esl_l);
    term1=(cpv-cpi)/Rv;
    term2=(xls-ttrip*(cpv-cpi))/Rv;
    double esl_i=exp((T-ttrip)*term2/(T*ttrip))*eref*pow((T/ttrip),(term1));
    double qsat_i=epsilon*esl_i/(p-esl_i);
    qsat=(1-omeg)*qsat_l + (omeg)*qsat_i;
  }
  if(iceflag==2){
    double term1=(cpv-cpi)/Rv;
    double term2=(xls-ttrip*(cpv-cpi))/Rv;
    double esl=exp((T-ttrip)*term2/(T*ttrip))*eref*pow((T/ttrip),(term1));
    qsat=epsilon*esl/(p-esl);
  }
  return qsat;
}

double O(double t, double p)
{
  return (t + kel) * pow(1000.0 / p, 0.28541);
}

double wobf(double temp)
{
  double x = temp - 20;
  double pol = 0;
  double wbts = 0;
  if(x <= 0){
    pol = 1. + x * (-0.0088416604999999992 + x * (0.00014714143000000001 + x * (-9.6719890000000006e-07 + x * (-3.2607217000000002e-08 + x * (-3.8598072999999999e-10)))));
    wbts = 15.130000000000001/(pow(pol,4));
  } else {
    pol = 1. + x * (0.0036182989000000001 + x * (-1.3603273e-05 + x * (4.9618921999999997e-07 + x * (-6.1059364999999998e-09 + x * (3.9401550999999998e-11 + x * (-1.2588129e-13 + x * (1.668828e-16)))))));
    wbts = 29.93/(pow(pol,4)) + 0.95999999999999996 * x - 14.800000000000001;
  }
  return wbts;
}

double OS(double t, double p)
{
  double akap = 0.28541;
  double pt = (t+kel)*pow((1000/p),akap)-kel;
  return pt-wobf(pt)+wobf(t);
}

double TMR(double W, double p)
{
  if (W != -622)
  {
    double zmienna = log10(W * p / (622 + W));
    return pow(10, 0.0498646455 * zmienna + 2.4082965) - 280.23475 + 38.9114 * (pow(pow(10, 0.0915 * zmienna) - 1.2035, 2));
  }
  else return -kel ;
}

double TDA(double O, double p)
{
  return O * pow(p / 1000.0, 0.28541) - kel;
}

double TSA(double OS, double p)
{
  double thw = OS;
  double akap = 0.28541;
  double pwrp = pow(p/1000,akap);
  double tone = (thw + kel) * pwrp - kel;
  double eone = wobf(tone) - wobf(thw);
  double rate = 1;
  double dlt = 1;
  double ttwo = 0;
  double pt = 0;
  double etwo = 0;
  while(abs(dlt) > 0.10000000000000001) {
    ttwo = tone - eone * rate;
    pt = (ttwo + kel)/pwrp - kel;
    etwo = pt + wobf(ttwo) - wobf(pt) - thw;
    dlt = etwo * rate;
    rate = (ttwo - tone)/(etwo - eone);
    tone = ttwo;
    eone = etwo;
  }
  double result = (ttwo - dlt) + kel;
  return result-273.15;
}

double TW(double t, double d, double p,  double *OW)
{
  double w = W(d, p);
  double theta = O(t, p);
  double pii = p;
  double x = 0;
  
  for (int i = 1; i <= 10; i++)
  {
    x = 0.02 * (TMR(w, pii) - TDA(theta, pii));
    if (abs(x) <= 0.01) break;
    pii *= pow(2.0, x);
  }
  double ti = TDA(theta, pii);
  double lump = OS(ti, pii);
  double Wet = TSA(lump, p);
  (*OW) = TSA(lump, 1000);
  return Wet;
}

double OE(double ts, double tds, double ps)
{
  double syf = 0;
  double alift = alcl(ts, tds, ps, syf);
  double olift = O(ts, ps);
  double tlift = TDA(olift, alift) + kel;
  return olift * exp(2.6518986 * W(tds, ps) / tlift) - kel;
}

double tv(double t, double w)
{
  double tv = ((t+kel) * (((w/1000) + 0.622) / (0.622 + 0.622 * (w/1000))));
  return tv-kel;
}

double alcl(double ts, double tds, double ps, double tmp)
{
  double aw = W(tds, ps);
  double ao = O(ts, ps);
  double pi = ps;
  double x = 0;
  for (int i = 1; i < 20; i++)
  {
    x = 0.01 * (TMR(aw, pi) - TDA(ao, pi));
    if (abs(x) < 0.01) break;
    pi *= pow(2.0, x);
  }
  tmp = TMR(aw, pi);
  return pi;
}

class Vector{
private:
  double x;
  double y;
  double z;
  
public:

  double X(){
    return this->x;
  }
  double Y(){
    return this->y;
  }
  double Z(){
    return this->z;
  }
  Vector(double X, double Y, double Z);
  Vector(double a, double v);
  Vector();
  Vector operator=(Vector const& r);
  
  Vector operator+(Vector const& r);
  Vector operator-(Vector const& r);
  Vector operator*(double const& r);
  Vector operator/(double const& r);
  
  Vector operator+=(Vector const& r);
  Vector operator-=(Vector const& r);
  Vector operator*=(double const& r);
  Vector operator/=(double const& r);
  
  double scalar(Vector const& r);
  static Vector vec(Vector const& l, Vector const& r);
  
  double* toAV();
  
  double abs();
};

Vector::Vector(double X, double Y, double Z){
  this->x=X;
  this->y=Y;
  this->z=Z;
}
Vector::Vector(double a,double v){
  this->x = v * cos(M_PI * a / 180.0);
  this->y = v * sin(M_PI * a / 180.0);
  this->z=0;
}
Vector Vector::operator*(double const& r){
  return Vector(this->x*r,this->y*r,this->z*r);
}
Vector Vector::operator/(double const& r){
  return Vector(this->x/r,this->y/r,this->z/r);
}
Vector Vector::operator+(Vector const& r){
  return Vector(this->x+r.x,this->y+r.y,this->z+r.z);
}
Vector Vector::operator-(Vector const& r){
  return Vector(this->x-r.x,this->y-r.y,this->z-r.z);
}
Vector Vector::operator=(Vector const& r){
  x=r.x;
  y=r.y;
  z=r.z;
  return *this;
}
Vector Vector::operator*=(double const& r){
  Vector tmp= Vector(this->x*r,this->y*r,this->z*r);
  this->x=tmp.x;
  this->y=tmp.y;
  this->z=tmp.z;
  return *this;
}
Vector Vector::operator/=(double const& r){
  Vector tmp= Vector(this->x/r,this->y/r,this->z/r);
  this->x=tmp.x;
  this->y=tmp.y;
  this->z=tmp.z;
  return *this;
}
Vector Vector::operator+=(Vector const& r){
  Vector tmp= Vector(this->x+r.x,this->y+r.y,this->z+r.z);
  this->x=tmp.x;
  this->y=tmp.y;
  this->z=tmp.z;
  return *this;
}
Vector Vector::operator-=(Vector const& r){
  Vector tmp= Vector(this->x-r.x,this->y-r.y,this->z-r.z);
  this->x=tmp.x;
  this->y=tmp.y;
  this->z=tmp.z;
  return *this;
}
double Vector::scalar(Vector const& r){
  return this->x*r.x + this->y*r.y + this->z*r.z;
}
Vector Vector::vec(Vector const& l, Vector const& r){
  double a = l.y * r.z - l.z * r.y;
  double b = l.z * r.x - r.z * l.x;
  double c = l.x * r.y - l.y * r.x;
  return Vector(a,b,c);
}
Vector::Vector(){
  this->x=0;
  this->y=0;
  this->z=0;
}

double* Vector::toAV(){
  double *result = new double[2];
  result[1]=sqrt(this->x*this->x + this->y*this->y + this->z*this->z);
  double x = this->x;
  double y= this->y;
  double angle = atan2(y,x)*180/M_PI;
  if(angle<0)angle =360+angle;
  
  result[0]=angle;
  return result;
}
double Vector::abs(){
  double v = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
  return v;
}

double distance(Vector x0, Vector x1, Vector x2){
  Vector l = x2-x1;
  Vector r = x1-x0;
  Vector up = Vector::vec(l,r);
  double u = up.abs();
  Vector down=x2-x1;
  double d = down.abs();
  double position = sign((x1.X() - x2.X()) * (x0.Y() - x2.Y()) - (x1.Y() - x2.Y()) * (x0.X() - x2.X()));
  return (u/d) * position;
}

class InfoCollector
{
protected:
  double h0;
  double p0;
  bool iLine;
  
  Cache* cache;
  
  void setInitialLine(double p, double h)
  {
    if (!iLine)
    {
      this->p0 = p;
      this->h0 = h;
      this->iLine = true;
    }
  }
  virtual void putSpecificLine(int i, double p, double h, double t, double d, double a, double v)=0;
public: 
  void setSoundingCache(Cache *c){
    this->cache = c;
  }
  
  InfoCollector()
  {
    this->iLine = false;
  }
  
  void putLine(int i, double p, double h, double t, double d, double a, double v)
  {
    this->setInitialLine(p, h);
    this->putSpecificLine(i, p, h, t, d, a, v);
  }
};

class Cache{
  friend class IndicesCollector;
private:
  int *pindex;
  double *p;	
  
  int *hindex;
  double *h;
  
  void initp();
  void inith();
  
  double h0;
  
public:
  Cache();
  virtual ~Cache();
  int getPressureIndex(double pressure);
  int getHeightIndex(double height);
  double * getArray(int i, int* len);
  double getH0();
  void setH0(double val);
  friend void setPressureIndex(int parrayI, int index, Cache *C);
  friend void setHeightIndex(int harrayI, int index, Cache *C);
  
};

Cache::Cache(){
  this->pindex = new int[10];
  this->p = new double[10];
  
  this->h = new double[15];
  this->hindex = new int[15];
  
  this->initp();
  this->inith();
}
Cache::~Cache(){
  delete[](this->pindex);
  delete[](this->p);
  
  delete[](this->hindex);
  delete[](this->h);
}
void Cache::initp(){
  for(int i=0;i<10;i++){
    this->pindex[i]=-1;
  }
  this->p[0]=1000;
  this->p[1]=850;
  this->p[2]=800;
  this->p[3]=700;
  this->p[4]=600;
  this->p[5]=500;
  this->p[6]=400;
  this->p[7]=300;
  this->p[8]=200;
  this->p[9]=100;
}
void Cache::inith(){
  for(int i=0;i<15;i++){
    this->hindex[i]=-1;
  }
  this->h[0]=100;
  this->h[1]=250;
  this->h[2]=500;
  this->h[3]=750;
  this->h[4]=1000;
  this->h[5]=1500;
  this->h[6]=2000;
  this->h[7]=3000;
  this->h[8]=4000;
  this->h[9]=5000;
  this->h[10]=6000;
  this->h[11]=7000;
  this->h[12]=8000;
  this->h[13]=9000;
  this->h[14]=10000;
}
int Cache::getHeightIndex(double height){
  for(int i=0; i<15; i++){
    if(this->h[i]==height)return this->hindex[i];
  }
  return -1;
}
int Cache::getPressureIndex(double pressure){
  for(int i=0; i<10; i++){
    if(this->p[i]==pressure)return this->pindex[i];
  }
  return -1;
}
double * Cache::getArray(int i, int* len){
  if(i==0){
    (*len)=10;
    return this->p;
  }else{
    (*len)=15;
    return this->h;
  }
}
double Cache::getH0(){
  return this->h0;
}
void Cache::setH0(double val){
  this->h0=val;
}
void setHeightIndex(int harrayI, int index, Cache *C){
  C->hindex[harrayI]=index;
}
void setPressureIndex(int parrayI, int index, Cache *C){
  C->pindex[parrayI]=index;
  
}
class Kinematics:public InfoCollector{
  friend class IndicesCollector;
private:
  list<Vector> *vw;
  Vector mean0;
  Vector mean6;
  double nsurf;
  double nsix;
  Vector v0;
  Vector v1;
  Vector v2;
  Vector v3;
  Vector v4;
  Vector v5;
  Vector v6;
  Vector v7;
  Vector v8;
  Vector v9;
  Vector v10;
  Vector llj;
  Vector mean06;
  Vector mean02;
  Vector mean03;
  Vector mean13;
  Vector mean26;
  Vector mean020;

  Vector mean36;
  Vector mean69;
  Vector mean912;
  Vector mean16;
  
  double SR_500_RM;
  double SR_1000_RM;
  double SR_3000_RM;
  double SR_16_RM;
  double SR_36_RM;

  double SR_500_LM;
  double SR_1000_LM;
  double SR_3000_LM;
  double SR_16_LM;
  double SR_36_LM;

  double n16sr;
  double n36sr;

  double n500;
  double n1000;
  double n3000;
  
  double n26;
  double n020;
  double n13;

  double n36;
  double n69;
  double n912;
  double n16;

  double n2;
  double n3;
  double n6;
  Vector mean01eff;
  Vector mean01;
  double n1;
  double n1eff;
  Vector CorfidiA;
  Vector Corfidi_upwind;
  Vector Corfidi_downwind;
  double lasth;
  
  double srh100rm;
  double srh100lm;

  double srh100rm2;
  double srh100lm2;

  double srh250rm;
  double srh250lm;
  
  double srh500rm;
  double srh500lm;

  double srh500rm2;
  double srh500lm2;

  double srh01lm;
  double srh01rm;
  double srh01sm;
  double srh01lm_eff;
  double srh01rm_eff;
  double srh01sm_eff;
  
  double srh01lmf;
  double srh01rmf;
  double srh01smf;
  double srh01lmf_eff;
  double srh01rmf_eff;
  double srh01smf_eff;
  
  double srh13lm;
  double srh13rm;
  double srh13sm;
  
  double srh13lmf;
  double srh13rmf;
  double srh13smf;

  double shear_l2;
  double sw13rm2;
  double sw13lm2;
  double srh13rm2;
  double srh13lm2;
  double srh13sm2;

  double srh03lm;
  double srh03rm;
  double srh03sm;
  double srh03lm_eff;
  double srh03rm_eff;
  double srh03sm_eff;
  
  double srh03lmf;
  double srh03rmf;
  double srh03smf;
  double srh03lmf_eff;
  double srh03rmf_eff;
  double srh03smf_eff;
  
  double srh36lm;
  double srh36rm;
  double srh36sm;

  double sw100rm;
  double sw100lm;

  double sw100rm2;
  double sw100lm2;

  double sw500rm;
  double sw500lm;

  double sw500rm2;
  double sw500lm2;

  double sw01rm;
  double sw01lm;
  
  double sw03rm;
  double sw03lm;

  double shear100m;
  double shear500m;
  double shear1000m; 
  double shear3000m;

  double shear100m2;
  double shear500m2;

  double sw13rm;
  double sw13lm;
  
  double sw13rmf;
  double sw13lmf;
  double shear_l;
  
  
  void putMandatoryVectors(int i, double p, double h, double t, double d, double a, double v, Vector v_);
  void putMeanVectors(int i, double p, double h, double t, double d, double a, double v, Vector v_);
  void putLLJ(int i, double p, double h, double t,double d,double a,double v, Vector v_);
  void putSpecificLine(int i, double p, double h, double t, double d, double a, double v);
  void prepareSupercellVectors();
  void prepareCorfidiVectors();
  void doSRH(int i, double p, double h, double t, double d, double a,double v);
  void doSRH2(int i, double p, double h, double t, double d, double a,double v);
  void finishMeanVectors();
  Vector shear06();
  

public:
  double muheight;
  double mumlheight;
  Sounding *S;
  Kinematics();
  Vector rm;
  Vector lm;
  virtual ~Kinematics();
  void putSecondPhaseLine(int i, double p, double h, double t, double d, double a, double v)
  {
    doSRH(i, p, h, t, d, a, v);
    doSRH2(i, p, h, t, d, a, v);
    lasth = h;
  }
  void effvec(){
      if (n1eff != 0) mean01eff *= 1.0 / n1eff;
      else mean01eff = Vector(0, 0, 0);
  }

  void finishPhase1()
  {
    finishMeanVectors();
    prepareSupercellVectors();
    prepareCorfidiVectors();
    lasth = h0;
    
    srh100rm=0;
    srh100lm=0;

    srh100rm2=0;
    srh100lm2=0;

    srh250rm=0;
    srh250lm=0;
    
    srh500rm=0;
    srh500lm=0;

    srh500rm2=0;
    srh500lm2=0;

    srh01lm=0;
    srh01rm=0;
    srh01sm=0;
    srh01lm_eff=0;
    srh01rm_eff=0;
    srh01sm_eff=0;

    srh01lmf=0;
    srh01rmf=0;
    srh01smf=0;
    
    srh13lm=0;
    srh13rm=0;
    srh13sm=0;
    
    srh13lmf=0;
    srh13rmf=0;
    srh13smf=0;

    shear_l2=0;
    sw13rm2=0;
    sw13lm2=0;
    srh13rm2=0;
    srh13lm2=0;
    srh13sm2=0;
    
    srh03lm = 0;
    srh03rm = 0;
    srh03sm = 0;
    srh03lm_eff = 0;
    srh03rm_eff = 0;
    srh03sm_eff = 0;
    
    srh03lmf = 0;
    srh03rmf = 0;
    srh03smf = 0;
    srh03lmf_eff = 0;
    srh03rmf_eff = 0;
    srh03smf_eff = 0;
    
    srh36lm = 0;
    srh36rm = 0;
    srh36sm = 0;

    sw100rm = 0;
    sw100lm = 0;

    sw100rm2 = 0;
    sw100lm2 = 0;

    sw500rm = 0;
    sw500lm = 0;

    sw500rm2 = 0;
    sw500lm2 = 0;

    sw01rm = 0;
    sw01lm = 0;
    
    sw03rm = 0;
    sw03lm = 0;

    shear100m = 0;
    shear500m = 0;

    shear100m2 = 0;
    shear500m2 = 0;

    shear1000m = 0;
    shear3000m = 0;
    
    sw13rm = 0;
    sw13lm = 0;
    
    sw13lmf=0;
    sw13rmf=0;  
    shear_l=0;
    
    SR_500_RM=0;
    SR_1000_RM=0;
    SR_3000_RM=0;
    SR_16_RM=0;
    SR_36_RM=0;
    
    SR_500_LM=0;
    SR_1000_LM=0;
    SR_3000_LM=0;
    SR_16_LM=0;
    SR_36_LM=0;
    
  }
  
  void finish()
  {
    
  }
};

Kinematics::Kinematics(){
  this->mean01eff = Vector(0, 0, 0);
  this->mean01 = Vector(0, 0, 0);
  this->mean06 = Vector(0, 0, 0);
  this->mean02 = Vector(0,0,0);
  this->mean03 = Vector(0,0,0);
  this->mean13 = Vector(0,0,0);	

  this->mean36 = Vector(0,0,0);	
  this->mean69 = Vector(0,0,0);	
  this->mean912 = Vector(0,0,0);	
  this->mean16 = Vector(0,0,0);	
  
  this->mean0=Vector(0,0,0);
  this->mean6=Vector(0,0,0);  
  this->mean26 = Vector(0,0,0);
  this->mean020 = Vector(0,0,0);
  
  this->lasth = h0;
  this->vw = new list<Vector>();
  this->llj = Vector(0, 0, 0);
  n2=0;
  n3=0;
  nsurf=0;
  nsix=0;
  n6=0;
  n1=0;
  n1eff=0;
  n13=0;
  n500=0;
  n1000=0;
  n3000=0;
  
  srh01lm = 0;
  srh01rm = 0;
  srh01sm = 0;
  srh01lm_eff = 0;
  srh01rm_eff = 0;
  srh01sm_eff = 0;
  
  srh01lmf = 0;
  srh01rmf = 0;
  srh01smf = 0;
  srh01lmf_eff = 0;
  srh01rmf_eff = 0;
  srh01smf_eff = 0;
  
  srh13lm = 0;
  srh13rm = 0;
  srh13sm = 0;
  
  srh13lmf = 0;
  srh13rmf = 0;
  srh13smf = 0;
  
  sw13lmf = 0;
  sw13rmf = 0;

  shear_l = 0;
  
  srh03lm = 0;
  srh03rm = 0;
  srh03sm = 0;
  srh03lm_eff = 0;
  srh03rm_eff = 0;
  srh03sm_eff = 0;
  
  srh03lmf = 0;
  srh03rmf = 0;
  srh03smf = 0;
  srh03lmf_eff = 0;
  srh03rmf_eff = 0;
  srh03smf_eff = 0;
  
  srh36lm = 0;
  srh36rm = 0;
  srh36sm = 0;
  
  SR_500_RM=0;
  SR_1000_RM=0;
  SR_3000_RM=0;
  SR_16_RM=0;
  SR_36_RM=0;
  
  SR_500_LM=0;
  SR_1000_LM=0;
  SR_3000_LM=0;
  SR_16_LM=0;
  SR_36_LM=0;
  
  n26=0;
  n36=0;
  n69=0;
  n912=0;
  n16=0;  
  n020=0;
  n16sr=0;
  n36sr=0;
  
}
Kinematics::~Kinematics(){
  delete(this->vw);
}
void Kinematics::putMandatoryVectors(int i, double p, double h, double t, double d, double a, double v, Vector v_)
{
  Cache *c = this->cache;
  if (i == 0)
  {
    this->v0 = v_;
  }
  else
  {
    int i1000 = c->getHeightIndex(1000);
    int i2000 = c->getHeightIndex(2000);
    int i3000 = c->getHeightIndex(3000);
    int i4000 = c->getHeightIndex(4000);
    int i5000 = c->getHeightIndex(5000);
    int i6000 = c->getHeightIndex(6000);
    int i7000 = c->getHeightIndex(7000);
    int i8000 = c->getHeightIndex(8000);
    int i9000 = c->getHeightIndex(9000);
    int i10000 = c->getHeightIndex(10000);
    
    if (i == i1000)
    {
      this->v1 = v_;
    }
    else if (i == i2000)
    {
      this->v2 = v_;
    }
    else if (i == i3000)
    {
      this->v3 = v_;
    }
    else if (i == i4000)
    {
      this->v4 = v_;
    }
    else if (i == i5000)
    {
      this->v5 = v_;
    }
    else if (i == i6000)
    {
      this->v6 = v_;
    }
    else if (i == i7000)
    {
      this->v7 = v_;
    }
    else if (i == i8000)
    {
      this->v8 = v_;
    }
    else if (i == i9000)
    {
      this->v9 = v_;
    }
    else if (i == i10000)
    {
      this->v10 = v_;
    }
  }
  
}

void Kinematics::putMeanVectors(int i, double p, double h, double t, double d, double a, double v, Vector v_)
{
  if((fmod(abs(h-h0),100.0)==0.0)||(h==h0)){
    
    if(h-h0<=500){
      mean0+=v_;
      nsurf+=1;
    }
          
    if ((h - h0 <= 1000))
    {
      mean01 += v_;
      n1 += 1;
      
    }
    if(h-h0<=6000 &&h-h0>=2000){
      mean26+=v_;
      n26+=1;
    }
    if(t>=-40 && t<=-10){
      mean020+=v_;
      n020+=1;
    }
    
  }
  
  if((fmod(abs(h-h0),200.0)==0.0)||(h==h0)){
    
    if (h-h0<=2000){
      mean02+=v_;
      n2+=1;
      
    }
    
    if (h-h0<=3000){
      mean03+=v_;
      n3+=1;
      
    }
    
    if (h-h0<=3000&&h-h0>=1000){
      mean13+=v_;
      n13+=1;
    }

     if (h-h0<=6000&&h-h0>=3000){
      mean36+=v_;
      n36+=1;
    }

     if (h-h0<=9000&&h-h0>=6000){
      mean69+=v_;
      n69+=1;
    }

     if (h-h0<=12000&&h-h0>=9000){
      mean912+=v_;
      n912+=1;
    }

    if (h-h0<=6000&&h-h0>=1000){
      mean16+=v_;
      n16+=1;
    }

    if (h - h0 <= 6000)
    {
      mean06 += v_;
      n6 += 1;
      
      if(h-h0>=5400){
        mean6+=v_;
        nsix+=1;
      }
    }
  }
}

void Kinematics::putLLJ(int i, double p, double h, double t,double d,double a,double v, Vector v_)
{
  double cabs = llj.abs();
  double vabs = v_.abs();
  if (h-h0 <= 1500 && vabs >= cabs) llj = v_;
}

void Kinematics::putSpecificLine(int i, double p, double h, double t, double d, double a, double v)
{
  Vector v_ = Vector(a,v * 0.514444);
  this->vw->push_back(v_);
  this->putMandatoryVectors(i, p, h, t, d, a, v, v_);
  this->putMeanVectors(i, p, h, t, d, a, v, v_);
  this->putLLJ(i, p, h, t, d, a, v, v_);
  lasth = h;
}

Vector Kinematics::shear06(){
  return this->v6 - this->v0;
}

void Kinematics::prepareSupercellVectors()
{
  Vector meanwind = this->mean06;
  Vector tv = Vector(0, 0, 1);
  Vector dev = Vector(0, 0, 0);
  Vector tshear = this->mean6-this->mean0;
  dev = Vector::vec(tshear,tv);
  dev *= 7.5;
  dev *= 1.0 / tshear.abs();
  this->rm = meanwind - dev;
  this->lm = meanwind + dev;
  
}

void Kinematics::prepareCorfidiVectors()
{
  int i850 = this->cache->getPressureIndex(850);
  int i700 = this->cache->getPressureIndex(700);
  int i500 = this->cache->getPressureIndex(500);
  int i300 = this->cache->getPressureIndex(300);
  
  Vector v850 = Get(this->vw, i850);
  Vector v700 = Get(this->vw, i700);
  Vector v500 = Get(this->vw, i500);
  Vector v300 = Get(this->vw, i300);
  
  
  CorfidiA = (v850 + v700 + v500 + v300) * 0.25;
  Corfidi_upwind = CorfidiA - llj;
  Corfidi_downwind = CorfidiA + Corfidi_upwind;
}

void Kinematics::doSRH(int i, double p, double h, double t, double d, double a,double v)
{	
  Vector meanwind = this->mean06;
  
  if((fmod(abs(h-h0),100.0)==0.0)||(h==h0)){
    
    double h_MU = this->mumlheight;
    Vector v_ = Get(this->vw, i);
    
    if (h-h_MU >= 0 && h-h_MU <= 500)
    {
      mean01eff += v_;
      n1eff += 1;
    }
  }
  
 if ((size_t)i < vw->size()-1 && h-h0<=6000){
      
    std::list<Vector>::iterator it = vw->begin();
    std::advance(it, i);
    
    std::list<Vector>::iterator it2 = vw->begin();
    std::advance(it2, i+1);
    
    Vector v1 = *it;
    Vector v2 = *it2;
    
    double tmps1 = (v1.X() - rm.X()) * (v2.Y() - v1.Y()) - (v1.Y() - rm.Y()) * (v2.X() - v1.X());
    double tmps2 = (v1.X() - lm.X()) * (v2.Y() - v1.Y()) - (v1.Y() - lm.Y()) * (v2.X() - v1.X());
    double tmps3 = (v1.X() - meanwind.X()) * (v2.Y() - v1.Y()) - (v1.Y() - meanwind.Y()) * (v2.X() - v1.X());
    
    double SR_U_rm = ( (v1.Y()+v2.Y() ) / 2 ) - rm.Y();
    double SR_V_rm = ( (v1.X()+v2.X() ) / 2 ) - rm.X();
    
    double SR_U_lm = ( (v1.Y()+v2.Y() ) / 2 ) - lm.Y();
    double SR_V_lm = ( (v1.X()+v2.X() ) / 2 ) - lm.X();

    double SR_U_sm = ( (v1.Y()+v2.Y() ) / 2 ) - meanwind.Y();
    double SR_V_sm = ( (v1.X()+v2.X() ) / 2 ) - meanwind.X();

    Vector SR_vec_lm = Vector(SR_U_lm, SR_V_lm,0);
    Vector SR_vec_rm = Vector(SR_U_rm, SR_V_rm,0);
    Vector SR_vec_sm = Vector(SR_U_sm, SR_V_sm,0);
    
    double SR_M_lm = SR_vec_lm.abs(); 
    double SR_M_rm = SR_vec_rm.abs();
    double SR_M_sm = SR_vec_sm.abs();
    
    double VORT_U = -(v2.X()-v1.X()); // (h-lasth);
    double VORT_V = (v2.Y()-v1.Y()); // (h-lasth);	  
    
    double OMEGA_rm = (SR_U_rm*VORT_U+SR_V_rm*VORT_V) / (sqrt( (SR_U_rm*SR_U_rm) + (SR_V_rm*SR_V_rm) ) );
    double OMEGA_lm = (SR_U_lm*VORT_U+SR_V_lm*VORT_V) / (sqrt( (SR_U_lm*SR_U_lm) + (SR_V_lm*SR_V_lm) ) );
    
    double shear_layer = sqrt(((v2.X() - v1.X()) * (v2.X() - v1.X())) + ((v2.Y() - v1.Y()) * (v2.Y() - v1.Y())));
    
    shear_l += shear_layer;
    
    sw13rm += OMEGA_rm;
    sw13lm += OMEGA_lm;
    
    srh13rm += tmps1;
    srh13lm += tmps2;
    srh13sm += tmps3;
    
    // if((fmod(abs(h-h0),100.0)==0.0)||(h==h0)){
    
    double h_MU = this->mumlheight;
    Vector v_ = Get(this->vw, i);
    
    if (h-h_MU >= 0 && h-h_MU <= 1000)
    {
      srh01rm_eff += tmps1;
      srh01lm_eff += tmps2;
      srh01sm_eff += tmps3;
    }

    if (h-h_MU >= 0 && h-h_MU <= 3000)
    {
      srh03rm_eff += tmps1;
      srh03lm_eff += tmps2;
      srh03sm_eff += tmps3;
    }
    // }
    
    if((fmod(abs(h-h0),100.0)==0.0)||(h==h0)){
      if(h-h0<=500){
        SR_500_RM += SR_M_rm;
        SR_500_LM += SR_M_lm;
        n500+=1;
      }
      if(h-h0<=1000){
        SR_1000_RM += SR_M_rm;
        SR_1000_LM += SR_M_lm;
        n1000+=1;
      }
      
    }
    
    if((fmod(abs(h-h0),200.0)==0.0)||(h==h0)){
      if(h-h0<=3000){
        SR_3000_RM += SR_M_rm;
        SR_3000_LM += SR_M_lm;
        n3000+=1;
      }
      if(h-h0 >= 3000 && h-h0 <= 6000){
        SR_36_RM += SR_M_rm;
        SR_36_LM += SR_M_lm;
        n36sr+=1;
      }
      
     if(h-h0 >= 1000 && h-h0 <= 6000){
        SR_16_RM += SR_M_rm;
        SR_16_LM += SR_M_lm;
        n16sr+=1;
      }
    }
    
    if(h-h0 >= 3000 && h-h0 <= 6000){
      srh36rm = srh13rm;
      srh36lm = srh13lm;
      srh36sm = srh13sm;
    }
    
    if(h-h0 <= 3000){
      srh03rm = srh13rm;
      srh03lm = srh13lm;
      srh03sm = srh13sm;
      sw03rm = sw13rm;
      sw03lm = sw13lm;
      shear3000m = shear_l;
    }
    
    if (h-h0<=1000){
      srh01rm = srh13rm;
      srh01lm = srh13lm;
      srh01sm = srh13sm;
      sw01rm = sw13rm;
      sw01lm = sw13lm;
      shear1000m = shear_l;
    }	  
    
    if(h-h0<=500){
      srh500rm = srh13rm;
      srh500lm = srh13lm;
      sw500rm = sw13rm;
      sw500lm = sw13lm;
      shear500m = shear_l;
    }
    
    if(h-h0<=250){
      srh250rm = srh13rm;
      srh250lm = srh13lm;
    }
    
    if(h-h0<=100){
      srh100rm = srh13rm;
      srh100lm = srh13lm;
      sw100rm = sw13rm;
      sw100lm = sw13lm;
      shear100m = shear_l;

    }
  }
}

void Kinematics::doSRH2(int i, double p, double h, double t, double d, double a, double v)
{	
  Vector meanwind = this->mean06;
  
 if ((size_t)i < vw->size()-1 && h-h0<=6000){

    std::list<Vector>::iterator it = vw->begin();
    std::advance(it, i);
    
    std::list<Vector>::iterator it2 = vw->begin();
    std::advance(it2, i+1);
    
    Vector v1 = *it;
    Vector v2 = *it2;
   
    if(i==0){
      v1 = Vector(0,0,0);
    }

    double tmps1 = (v1.X() - rm.X()) * (v2.Y() - v1.Y()) - (v1.Y() - rm.Y()) * (v2.X() - v1.X());
    double tmps2 = (v1.X() - lm.X()) * (v2.Y() - v1.Y()) - (v1.Y() - lm.Y()) * (v2.X() - v1.X());
    double tmps3 = (v1.X() - meanwind.X()) * (v2.Y() - v1.Y()) - (v1.Y() - meanwind.Y()) * (v2.X() - v1.X());
    
    double SR_U_rm = ( (v1.Y()+v2.Y() ) / 2 ) - rm.Y();
    double SR_V_rm = ( (v1.X()+v2.X() ) / 2 ) - rm.X();
    
    double SR_U_lm = ( (v1.Y()+v2.Y() ) / 2 ) - lm.Y();
    double SR_V_lm = ( (v1.X()+v2.X() ) / 2 ) - lm.X();

    double SR_U_sm = ( (v1.Y()+v2.Y() ) / 2 ) - meanwind.Y();
    double SR_V_sm = ( (v1.X()+v2.X() ) / 2 ) - meanwind.X();

    Vector SR_vec_lm = Vector(SR_U_lm, SR_V_lm,0);
    Vector SR_vec_rm = Vector(SR_U_rm, SR_V_rm,0);
    Vector SR_vec_sm = Vector(SR_U_sm, SR_V_sm,0);
    
    double SR_M_lm = SR_vec_lm.abs(); 
    double SR_M_rm = SR_vec_rm.abs();
    double SR_M_sm = SR_vec_sm.abs();
    
    double VORT_U = -(v2.X()-v1.X()); // (h-lasth);
    double VORT_V = (v2.Y()-v1.Y()); // (h-lasth);	  
    
    double OMEGA_rm = (SR_U_rm*VORT_U+SR_V_rm*VORT_V) / (sqrt( (SR_U_rm*SR_U_rm) + (SR_V_rm*SR_V_rm) ) );
    double OMEGA_lm = (SR_U_lm*VORT_U+SR_V_lm*VORT_V) / (sqrt( (SR_U_lm*SR_U_lm) + (SR_V_lm*SR_V_lm) ) );
    
    double shear_layer = sqrt(((v2.X() - v1.X()) * (v2.X() - v1.X())) + ((v2.Y() - v1.Y()) * (v2.Y() - v1.Y())));
    
    shear_l2 += shear_layer;
    
    sw13rm2 += OMEGA_rm;
    sw13lm2 += OMEGA_lm;
  
    srh13rm2 += tmps1;
    srh13lm2 += tmps2;
    srh13sm2 += tmps3;
           
    if(h-h0<=500){
      srh500rm2 = srh13rm2;
      srh500lm2 = srh13lm2;
      sw500rm2 = sw13rm2;
      sw500lm2 = sw13lm2;
      shear500m2 = shear_l2;
    }
    
    if(h-h0<=100){
      srh100rm2 = srh13rm2;
      srh100lm2 = srh13lm2;
      sw100rm2 = sw13rm2;
      sw100lm2 = sw13lm2;
      shear100m2 = shear_l2;
    }
  }
}


void Kinematics::finishMeanVectors()
{
  if (n6 != 0) mean06 *= 1.0 / n6;
  else mean06 = Vector(0, 0, 0);
  if (n1 != 0) mean01 *= 1.0 / n1;
  else mean01 = Vector(0, 0, 0);

  if (n36 != 0) mean36 *= 1.0 / n36;
  else mean36 = Vector(0, 0, 0);

  if (n69 != 0) mean69 *= 1.0 / n69;
  else mean69 = Vector(0, 0, 0);

  if (n912 != 0) mean912 *= 1.0 / n912;
  else mean912 = Vector(0, 0, 0);

  if (n16 != 0) mean16 *= 1.0 / n16;
  else mean16 = Vector(0, 0, 0);

  if (n13 != 0) mean13 *= 1.0 / n13;
  else mean13 = Vector(0, 0, 0);
  if (n2 != 0) mean02 *= 1.0 / n2;
  else mean02 = Vector(0, 0, 0);
  if (n3 != 0) mean03 *= 1.0 / n3;
  else mean03 = Vector(0, 0, 0);
  if(nsurf!=0)mean0/=nsurf;
  else mean0=Vector(0,0,0);
  if(nsix!=0)mean6/=nsix;
  else mean6=Vector(0,0,0);
  
  mean26/=n26;
  mean020/=n020;
}

class LapseRate{
  friend class IndicesCollector;
private:
  
  int lclIndex;
  int vLclIndex;
  
  int lfcIndex;
  int vLfcIndex;
  
  int elIndex;
  int vElIndex;
  
  double cape;
  double cin;
  double to3cape;
  double to2cape;
  
  double vcape;
  double vcin;
  double vto3cape;
  double vto2cape;

  double vcin500;
  
  double os;
  double o;
  double w;
  double tcin;
  double tvcin;
  double vos;
  double vo;
  double vw;
  double tch;
  
  double middlecape;
  double coldcape; 
  double coldcapeTV; 
  double peakB;
  double peakB_M10;
  double peakB_3km;

  double starth;

  bool isSet;
  bool dcape_;
  
  double h0;
  
  double dcape;
  double dvcape;
  
  int i700index;
  
  void allocate();
  
  void free();
  
  void testSpecificLCL(int i,double p,double t,double tmr, double tda, int* lclInd_, int* lfcInd_, list<double>* curve_, double*os_);
  void doRest(int i, double p,double h, double t, double TSA, int* lfcInd_, double* cape_, double* to3, double* to2, double* cin_, int* elInd_, list<double>* curve_);
  void putClassicLine(int i, double p, double h, double t,double d, double a, double v);
  void putVirtualLine(int i, double p, double h, double t, double d, double a, double v);
  
  
public:
  double NCAPE;
  double NCAPE_N;
  void putNCAPE(int i, double value){
    if(i >= vLfcIndex && i <= vElIndex){
      NCAPE += value;
      NCAPE_N += 1;
    }
  } 
  list<double> *values;
  list<double> *virtualValues;
  double lasth;
  list<double> *getVirtualValues(){
    return this->virtualValues;
  }
  int startIndex;
  LapseRate();
  ~LapseRate();
  void setInitialConditions(int i, double p, double h, double t, double d, double a, double v, double h0);
  void setInitialW(double w_, double o_);
  void prepareForDCAPE();
  void putLine(int i, double p, double h, double t, double d, double a, double v);
  void finish();
  void seti700index(int value){
    i700index=value;
  }
};

void LapseRate::allocate(){
  NCAPE=0;
  NCAPE_N=0;
  values = new list<double>();
  virtualValues = new list<double>();
  cape = cin = to3cape = to2cape = vcape = vcin = vto3cape = vto2cape = os = o = w = vos=vo=vw=dcape=dvcape=0;
  vcin500 = 0;
  middlecape=0;
  coldcape=0;
  coldcapeTV=0;
  peakB=0;
  peakB_M10=0;
  peakB_3km=0;
  lclIndex = vLclIndex = lfcIndex = vLfcIndex = elIndex = vElIndex = -1;
  startIndex=-1;
  isSet = false;
  tcin = 0;
  tvcin = 0;
  i700index=-1;
  this->dcape_=false;
}
void LapseRate::setInitialW(double w_, double o_){
  this->w = w_;
  this->vw = w_;
  this->o=o_;
  this->vo=o_;
}
void LapseRate::free(){
  delete(this->values);
  delete(this->virtualValues);
}
LapseRate::LapseRate(){
  this->allocate();
}
LapseRate::~LapseRate(){
  this->free();
}
void LapseRate::setInitialConditions(int i, double p, double h, double t, double d, double a, double v, double h0){
  this->free();
  this->allocate();
  this->h0 = h0;
  
  this->isSet = true;
  this->o = O(t, p);
  this->w = W(d, p);
  
  this->vo = O(t, p);
  this->vw = W(d, p);

  this->starth=h;
  this->vcin500=0;
  this->startIndex=i;
  this->lasth = h;
  this->tch= W(d,p);
}
void LapseRate::testSpecificLCL(int i,double p,double t,double tmr, double tda, int* lclInd_, int* lfcInd_, list<double>* curve_, double*os_){
  if ((*lclInd_) == -1)
  {
    if (tmr < tda)
    {
      curve_->push_back(tda);
      (*lfcInd_)= -1;
    }
    else
    {
      (*lfcInd_)= -1;
      (*os_)= OS(tmr, p);
      (*lclInd_)= i;
      
    }
  }
}
void LapseRate::doRest(int i, double p,double h, double t, double TSA, int* lfcInd_, double* cape_, double* to3, double* to2, double* cin_, int* elInd_, list<double>* curve_){
  curve_->push_back(TSA);
  double dz = abs(h - lasth);
  double tcap = g * dz * (TSA - t) / (t + kel);
  if (TSA >= t)
  {
    if ((*lfcInd_) == -1)
    {
      (*lfcInd_) = i;
      
    }
    if ((*elInd_)!= -1)
    {
      (*elInd_ )= -1;
      cin += tcin;
      tcin = 0;
    }
    (*cape_ )+= tcap;
    if (h - h0 <3000) (*to3) = (*cape_);
    
  }
  else
  {
    if (i <= i700index&&this->dcape_)
    {
      dcape += tcap;
    }
    
    
    if ((*lfcInd_) == -1)
    {
      (*cin_ )+= tcap;
    }
    else
    {
      tcin += tcap;
      if ((*elInd_ )== -1) (*elInd_ )= i;
    }
  }
}
void LapseRate::putClassicLine(int i, double p, double h, double t,double d, double a, double v){
  double tmr = TMR(w, p);
  double tda = TDA(o, p);
  
  
  
  
  this->testSpecificLCL(i, p,  t,tmr, tda, &lclIndex, &lfcIndex, values, &os);
  if (lclIndex != -1)
  {
    double TSA_ = TSA(os, p);
    this->doRest(i, p, h, t, TSA_, &lfcIndex, &cape, &to3cape, &to2cape, &cin, &elIndex, values);
  }    else
  {
    
    double dz = abs(h - lasth);
    double tcap = g * dz * (tda - t) / (t + kel);
    double t_parcel = tda;
    double tch_ = W(t_parcel, p);
    double tw = W(d, p);
    double vt_parcel = tv(t_parcel, tch_);
    double t_ = tv(t, tw);
    double tvcap = g * dz * (vt_parcel - t_) / (t_ + kel);
    if (vLclIndex == -1)
    {
      if (t_parcel < t)
      {
        cin += tcap;
      }
      if (vt_parcel < t_)
      {
        vcin += tvcap;
      }
    }
    
    
  }
  
}
void LapseRate::prepareForDCAPE(){
  this->free();
  values = new list<double>();
  virtualValues = new list<double>();
  cape = cin = to3cape = to2cape = vcape = vcin = vto3cape = vto2cape = 0;
  dcape = 0;
  dvcape = 0;
  startIndex = 0;
  dcape_ = true;
  lasth = h0;
}
void LapseRate::putVirtualLine(int i, double p, double h, double t, double d, double a, double v)
{
  this->vLclIndex = lclIndex;
  double t_parcel = *(--(this->values->end()));
  if(this->vLclIndex!=-1)tch=W(t_parcel,p);
  
  double tw = W(d, p);
  double vt_parcel = tv(t_parcel, tch);
  this->virtualValues->push_back(vt_parcel);
  double t_ = tv(t, tw);
  double dz = abs(h - lasth);

  if (vLclIndex != -1) {  
  double ttt = (vt_parcel - t_) / (t_+273.15) * 9.81;
  if(ttt>this->peakB) {
    peakB = ttt;
  }

  if ( t<=-10&&t>=-40 ) {
  double ttt2 = (vt_parcel - t_) / (t_+273.15) * 9.81;
  if(ttt2>this->peakB_M10) {
    peakB_M10 = ttt2;
     }
   }

  if ( (h-h0) <= 3000) {
  double ttt3 = (vt_parcel - t_) / (t_+273.15) * 9.81;
  if(ttt3>this->peakB_3km) {
    peakB_3km = ttt3;
     }
   }
  }
  
  double tcap = g * dz * (vt_parcel - t_) / (t_ + kel);
 	
  if (vLclIndex != -1) { 
    if (vt_parcel >= t_)
    {
      {
        if (vLfcIndex == -1) vLfcIndex = i;
        if (vElIndex != -1)
        {
          vElIndex = -1;
          vcin += tvcin;
          tvcin = 0;
        }
        vcape += tcap;
        if (h - h0 < 3000) vto3cape = vcape;
        if (h - h0 < 2000) vto2cape = vcape;
        if(t<=-10&&t>=-40)middlecape+=tcap;
        if(t<= -10)coldcape+=tcap;
        if(vt_parcel<= -10)coldcapeTV+=tcap;
      }
    }
    else
    {
      if (vLfcIndex == -1)
      {
        vcin += tcap;
      }
      else
      {
        tvcin += tcap;
        if (vElIndex == -1)
          
          vElIndex = i;
        
        
      }
      if (i <= i700index && this->dcape_)
      {
        dvcape += tcap;
      }
    }
  }

  if( (h <= starth+4000) && (tcap < 0) ) {
       vcin500 += tcap;
   cout << vcin500 << " " << vcin << " " << (h-h0) << " " << vt_parcel << " " << t_ << endl;
   } 

}
void LapseRate::putLine(int i, double p, double h, double t, double d, double a, double v){
  if (i >= startIndex) {
    putClassicLine(i, p, h, t, d, a, v);
    putVirtualLine(i, p, h, t, d, a, v); 
    lasth = h;}
}
void LapseRate::finish(){
  if (cape ==0) cin = 0;
  if (vcape ==0) vcin = 0;
  if (dcape > 0) dcape = 0;
  if (dvcape > 0) dvcape = 0;
  dvcape = abs(dvcape);
  dcape = abs(dcape);
}

class Thermodynamics:public InfoCollector{
  friend class IndicesCollector;
public:

  double last_MSE0;
  double aggregated_MSE0;

  double meanLayerZHeight;
  double meanLayerBottom;
  double meanLayerTop;

  double n;
  double mp;
  double mh;
  double mt;
  double md;
  double mmr;

  double meanmostUnstableUP;

  double mumli1;
  double oe1;
  double mp1;
  double mh1;
  double mt1;
  double md1;
  double mmr1;
  double mo1;

  double mumli2;
  double oe2;
  double mp2;
  double mh2;
  double mt2;
  double md2;
  double mmr2;
  double mo2;

  double mumli3;
  double oe3;
  double mp3;
  double mh3;
  double mt3;
  double md3;
  double mmr3;
  double mo3;

  double mumli4;
  double oe4;
  double mp4;
  double mh4;
  double mt4;
  double md4;
  double mmr4;
  double mo4;

  double mumli5;
  double oe5;
  double mp5;
  double mh5;
  double mt5;
  double md5;
  double mmr5;
  double mo5;

  double mumli6;
  double oe6;
  double mp6;
  double mh6;
  double mt6;
  double md6;
  double mmr6;
  double mo6;

  double mumliLAST;
  double oeLAST;
  double mpLAST;
  double mhLAST;
  double mtLAST;
  double mdLAST;
  double mmrLAST;
  double moLAST;

  double mumliMAX;
  double oeMAX;
  double mpMAX;
  double mhMAX;
  double mtMAX;
  double mdMAX;
  double mmrMAX;
  double moMAX;

  double t0;
  double pwater;
  double pwater_eff;
  double lastp;
  double minTHTE;
  int minTHTEpos;
  double maxOE;
  double maxOE500;
  double mo;
  int zeropos;
  int wbzeropos;
  double zero;
  double wb0;
  
  int mintenpos;
  double minten;
  
  int min40pos;
  double min40;

  int min25pos;
  double min25;

  double mr1000;
  
  double lr01;
  double h01;
  
  double lr24;
  double h24;
  
  int _700;
  double lastt;
  double lasth;
  
  double t10;
  int p10;
  
  double downmr;
  double downmrn;
  double downo;
  double downon;
  
  double thetd;
  double thetn;
  
  double thet01d;
  double thet01n;
  
  double thet02d;
  double thet02n;
  
  double mthet;
  double mthetn;
  
  list<double>* wbt;
  list<double>* oe;
  list<double>* mixing;
  list<double>* virt;
  
  double meanhum1;
  double meand1b;
  
  double meanhum2;
  double meand2b;
  
  double meanhum25;
  double meand25;
  
  double meanhum36;
  double meand36;
  
  double meanhum14;
  double meand14;

  double meanhum850500;
  double meand850500;
  double meanhumMIDDLE;
  double meandMIDDLE;
  
  double meanmxr2;
  double meand2;  

  list<double>* MSE0;
  list<double>* MSE0_star;
  list<double>* MSE0_bar;
  list<double>* int_arg_MSE0;

  LapseRate* mostUnstable;
  LapseRate* mostU500;
  LapseRate* surfaceBased;
  LapseRate* meanLayer;
  LapseRate* downdraft;
  LapseRate* showalter;
  LapseRate* meanmostUnstable;
  
  void putMinTHTE(int i, double p, double h, double oe);
  void startConditions(int i, double p, double h, double t, double d, double a, double v, double oe);
  void putMaxTHTE(int i, double p, double h, double t, double d, double a, double v, double oe, double mr);
  void putMeanLayerParameters(int i, double p, double h, double t, double d, double a, double v,double mr);
  void determineDowndraft700(int i, double p, double h, double t, double d, double a, double v);
  void determineDowndraftByMinTHTE(int i, double p, double h, double t, double d, double a, double v);
  void putShowalter(int i, double p, double h, double t, double d, double a, double v);
  void putPWATER( int i, double p, double h, double t, double d, double a, double v);
  void putLowLapseRates(int i, double p, double h, double t, double d, double a,double v);
  void finishLowLapseRates();
  void ZeroPosStartingConditions(int i, double p, double h, double t, double d, double a, double v, double wbt);
  void putZeroPos(int i, double p, double h, double t, double d, double a, double v, double wbt);
  inline void putSpecificLine(int i, double p, double h, double t, double d, double a, double v);
  
  
  void setMlIndex(int i, double p, double h, double t, double d, double a, double v);
  void putMlLine(int i, double p, double h, double t, double d, double a, double v);
  Thermodynamics();
  virtual ~Thermodynamics();
  void prepareMeanLayer();
  void putMeanLine(int i, double p, double h, double t, double d, double a, double v);
  void finish();
  
  
  
};

void Thermodynamics::putMinTHTE(int i, double p, double h, double oe){
  
  if (h-h0 <= 4000 && oe < minTHTE){
    minTHTE = oe;
    minTHTEpos = i;
  }
}

// Sekcja ML nie wykorzystywana
void Thermodynamics::setMlIndex(int i, double p, double h, double t, double d, double a, double v){
  this->meanLayer->setInitialConditions(i, p, h, t,d, a, v,h0);
  this->meanLayer->setInitialW(mmr, mo);
}

void Thermodynamics::putMlLine(int i, double p, double h, double t, double d, double a, double v){
  this->meanLayer->putLine(i, p, h, t, d, a, v);
}
Thermodynamics::Thermodynamics(){

  last_MSE0 = 0;
  aggregated_MSE0 = 0;
  this->MSE0 = new list<double>();
  this->MSE0_star = new list<double>();
  this->MSE0_bar = new list<double>();
  this->int_arg_MSE0 = new list<double>();
  
  mr1000 = 0;
  meanLayerZHeight=500.0;
  n=0;
  mp=0;
  mh=0;
  mt=0;
  md=0;
  mo=0;
  mmr=0;

  meanmostUnstableUP=500;

  mumli1=0;
  oe1=0;
  mp1=0;
  mh1=0;
  mt1=0;
  md1=0;
  mmr1=0;
  mo1=0;

  mumli2=0;
  oe2=0;
  mp2=0;
  mh2=0;
  mt2=0;
  md2=0;
  mmr2=0;
  mo2=0;

  mumli3=0;
  oe3=0;
  mp3=0;
  mh3=0;
  mt3=0;
  md3=0;
  mmr3=0;
  mo3=0;
  
  mumli4=0;
  oe4=0;
  mp4=0;
  mh4=0;
  mt4=0;
  md4=0;
  mmr4=0;
  mo4=0;
  
  mumli5=0;
  oe5=0;
  mp5=0;
  mh5=0;
  mt5=0;
  md5=0;
  mmr5=0;
  mo5=0;
  
  mumli6=0;
  oe6=0;
  mp6=0;
  mh6=0;
  mt6=0;
  md6=0;
  mmr6=0;
  mo6=0;

  mumliLAST=0;
  oeLAST=-273.15;
  mpLAST=0;
  mhLAST=0;
  mtLAST=0;
  mdLAST=0;
  mmrLAST=0;
  moLAST=0;
  
  mumliMAX=0;
  oeMAX=-273.15;
  mpMAX=0;
  mhMAX=0;
  mtMAX=0;
  mdMAX=0;
  mmrMAX=0;
  moMAX=0;
  
  t0=0;
  t10=0;
  p10=0;
  pwater=0;
  pwater_eff=0;
  lastp=1000;
  minTHTE=0;
  minTHTEpos=0;
  maxOE=0;	
  zeropos=0;
  wbzeropos=0;
  zero=0;
  wb0=0;
  mr1000=0;
  lr01=0;
  h01=0;
  lr24=0;
  h24=0;
  _700=0;
  lastt=0;
  lasth=0;
  
  this->wbt = new list<double>();
  this->oe = new list<double>();
  this->mixing = new list<double>();
  this->virt=new list<double>();
  this->mostUnstable = new LapseRate();
  this->surfaceBased = new LapseRate();
  this->meanLayer = new LapseRate();
  this->downdraft = new LapseRate();
  this->showalter = new LapseRate();
  this->mostU500= new LapseRate();
  this->meanmostUnstable = new LapseRate();
  
  meanhum1=0;
  meand1b=0;
  
  meanhum2=0;
  meand2b=0;

  meanhum850500=0;
  meand850500=0;
  
  meanhumMIDDLE=0;
  meandMIDDLE=0;
  
  meanhum25=0;
  meand25=0;
  
  meanhum36=0;
  meand36=0;
  
  meanhum14=0;
  meand14=0;
  
  meanmxr2=0;
  meand2=0;
  
  downmr=0;
  downmrn=0;
  downo=0;
  downon=0;
  thetd=0;
  thetn=0;
  thet01d=0;
  thet01n=0;
  thet02d=0;
  thet02n=0;
  
  mthet=0;
  mthetn=0;
}
Thermodynamics::~Thermodynamics(){
  delete(this->wbt);
  delete(this->oe);
  delete(this->mixing);
  delete(this->mostUnstable);
  delete(this->surfaceBased);
  delete(this->meanLayer);
  delete(this->downdraft);
  delete(this->showalter);
  delete(this->virt);
  delete(this->mostU500);
  delete(this->meanmostUnstable);
  delete(this->MSE0);
  delete(this->MSE0_star);
  delete(this->MSE0_bar);
  delete(this->int_arg_MSE0);
}
void Thermodynamics::startConditions(int i, double p, double h, double t, double d, double a, double v, double oe)
{
  this->t0 = t;
  this->t10=t;
  this->surfaceBased->setInitialConditions(i, p, h, t, d, a, v, h0);
  this->mostUnstable->setInitialConditions(i, p, h, t, d, a, v, h0);
  this->mostU500->setInitialConditions(i, p, h, t, d, a, v, h0);
  this->lastp = p;
  maxOE = oe;
  maxOE500=-273.15;
  minTHTE = oe;
  minTHTEpos = i;
  this->lr01 = 0;
  this->lr24 = 0;
  this->h01 = 0;
  this->h24 = 0;
  
  meanhum1=ESAT(d)/ESAT(t);
  meand1b=1;
  
  meanhum2=ESAT(d)/ESAT(t);
  meand2b=1;
  
  meanhum25=ESAT(d)/ESAT(t);
  meand25=1;
   
  meanhum36=ESAT(d)/ESAT(t);
  meand36=1;
  
  meanhum14=ESAT(d)/ESAT(t);
  meand14=1;
  
  meanhumMIDDLE=ESAT(d)/ESAT(t);
  meandMIDDLE=1;

  meanhum850500=ESAT(d)/ESAT(t);
  meand850500=1;

  meanmxr2=W(d, p);
  meand2=1;
}

void Thermodynamics::putMaxTHTE(int i, double p, double h, double t, double d, double a, double v, double oe, double mr)
{
  
  if (oe > maxOE && h-h0 <= 3000){
    maxOE = oe;
    this->mostUnstable->setInitialConditions(i, p, h, t, d, a, v, h0);
  }

  if( ((fmod(abs(h-h0),100.0)==0.0)  || (h==h0)) && (meanmostUnstableUP <= 3000) ) {

    double wys = h-h0;
    
   if(wys == 0 || wys == 600 || wys == 1200 || wys == 1800 || wys == 2400 || wys == 3000){
      mumli1 = i;
      oe1 = oe;
      mh1 = h;
      mp1 = p;
      mt1 = t;
      md1 = d;
      mmr1 = mr;
      mo1 = O(t,p);
    }

   if(wys == 100 || wys == 700 || wys == 1300 || wys == 1900 || wys == 2500){
      mumli2 = i;
      oe2 = oe;
      mh2 = h;
      mp2 = p;
      mt2 = t;
      md2 = d;
      mmr2 = mr;
      mo2 = O(t,p);
    }

   if(wys == 200 || wys == 800 || wys == 1400 || wys == 2000 || wys == 2600){
      mumli3 = i;
      oe3 = oe;
      mh3 = h;
      mp3 = p;
      mt3 = t;
      md3 = d;
      mmr3 = mr;
      mo3 = O(t,p);
    }

   if(wys == 300 || wys == 900 || wys == 1500 || wys == 2100 || wys == 2700){
      mumli4 = i;
      oe4 = oe;
      mh4 = h;
      mp4 = p;
      mt4 = t;
      md4 = d;
      mmr4 = mr;
      mo4 = O(t,p);
    }

   if(wys == 400 || wys == 1000 || wys == 1600 || wys == 2200 || wys == 2800){
      mumli5 = i;
      oe5 = oe;
      mh5 = h;
      mp5 = p;
      mt5 = t;
      md5 = d;
      mmr5 = mr;
      mo5 = O(t,p);
    }

  if(wys == 500 || wys == 1100 || wys == 1700 || wys == 2300 || wys == 2900){
      mumli6 = i;
      oe6 = oe;
      mh6 = h;
      mp6 = p;
      mt6 = t;
      md6 = d;
      mmr6 = mr;
      mo6 = O(t,p);
    }

    if(wys == meanmostUnstableUP){
      oeLAST = (oe1+oe2+oe3+oe4+oe5+oe6)/6;
      mhLAST = (mh1+mh2+mh3+mh4+mh5+mh6)/6;
      mpLAST = max(max(max(max(max(mp1,mp2),mp3),mp4),mp5),mp6);
      mtLAST = (mt1+mt2+mt3+mt4+mt5+mt6)/6;
      mdLAST = (md1+md2+md3+md4+md5+md6)/6;
      mmrLAST = (mmr1+mmr2+mmr3+mmr4+mmr5+mmr6)/6;
      moLAST = (mo1+mo2+mo3+mo4+mo5+mo6)/6;
      mumliLAST = min(min(min(min(min(mumli1,mumli2),mumli3),mumli4),mumli5),mumli6); 
      meanmostUnstableUP += 100;
    }

    if(oeLAST > oeMAX){
         oeMAX = oeLAST;
         mhMAX = mhLAST;
         mpMAX = mpLAST;
         mtMAX = mtLAST;
         mdMAX = mdLAST;
         mmrMAX = mmrLAST;
         moMAX = moLAST;
         mumliMAX = mumliLAST;

    this->meanmostUnstable->setInitialConditions(mumliMAX, mpMAX, mhMAX, mtMAX, mdMAX, 0, 0, h0);
    this->meanmostUnstable->setInitialW(mmrMAX, moMAX);
    }
   }   
    if(oe>maxOE500&&h-h0<=3000 && h-h0>=500){
    maxOE500=oe;
    this->mostU500->setInitialConditions(i, p, h, t, d, a, v, h0);
  }  
}

void Thermodynamics::putMeanLayerParameters(int i, double p, double h, double t, double d, double a, double v,double mr)
{  
  if ((abs(h - h0) >= meanLayerBottom && abs(h - h0) <= meanLayerTop)  && ((fmod(abs(h-h0),100.0)==0.0)  || (h==h0)))
  {
    mh += h;
    if(p>mp)mp=p;
    mt += t;
    md += d;
    mmr += mr;
    mo += O(t,p);
    n += 1;

  }

  if((abs(h - h0) <= 1000)&&(abs(h - h0) >= 0) && (fmod(abs(h-h0),100.0)==0.0)){
    thet01d+=OE(t,d,p);
    thet01n+=1;
  }
  
  if((abs(h - h0) <= 2000)&&(abs(h - h0) >= 0) && (fmod(abs(h-h0),200.0)==0.0)){
    thet02d+=OE(t,d,p);
    thet02n+=1;
  }
  
  if((abs(h - h0) <= 5400)&&(abs(h - h0) >= 2600)&&(fmod(abs(h-h0),200.0)==0.0)){
    downmr+=W(d,p);
    downmrn+=1;
    downo+=O(t,p);
    downon+=1;
    
    thetd+=OE(t,d,p);
    thetn+=1;
  }
  
  if(t>=-40.0&&t<=-10.0){
    mthet+=OE(t,d,p);
    mthetn+=1;
  }
}
void Thermodynamics::determineDowndraft700(int i, double p, double h, double t, double d, double a, double v)
{
  int t700p = this->cache->getPressureIndex(700);
  if (t700p != _700)
  {
    _700 = t700p;
    downdraft->setInitialConditions(i, p, h, t, d, a, v, h0);
  }
  
  if (i >= _700)
  {
    downdraft->putLine(i, p, h, t, d, a, v);
  }
}

void Thermodynamics::determineDowndraftByMinTHTE(int i, double p, double h, double t, double d, double a, double v)
{
  
  if ((h-h0)==4000)
  {
    downdraft->seti700index(i);
    
  }
  
  if ((h-h0) >= 4000)
  {
    downdraft->putLine(i, p, h, t, d, a, v);
  }
}
void Thermodynamics::putShowalter(int i, double p, double h, double t, double d, double a, double v)
{
  int lindex = cache->getPressureIndex(850);
  if (i == lindex)
  {
    showalter->setInitialConditions(i, p, h, t, d, a, v, h0);
  }
  
  if (i >= lindex)
  {
    showalter->putLine(i, p, h, t, d, a, v);
  }
}
void Thermodynamics::putPWATER( int i, double p, double h, double t, double d, double a, double v)
{
  std::list<double>::iterator it = mixing->begin();
  std::advance(it, i);
  pwater += 0.5*(lastp - p) * (*(it)+*(--it));
  
  double it_eff = (*(it)+*(--it))*(ESAT(d)/ESAT(t));
  pwater_eff+=0.5 *(lastp-p)*it_eff;

}
void Thermodynamics::putLowLapseRates(int i, double p, double h, double t, double d, double a,double v)
{
  if (h <= h0 + 1000)
  {
    lr01 += t - lastt;
    h01 += h - lasth;
  }
  
  if (h >= h0 + 2000 && h<=h0+4000)
  {
    lr24 += t - lastt;
    h24 += h - lasth;
  }
}
void Thermodynamics::finishLowLapseRates()
{
  h01 /= 1000;
  h24 /= 1000;
  
  if (h01 != 0) lr01 /= h01;
  else lr01 = 0;
  if (h24 != 0) lr24 /= h24;
  else lr24 = 0;
}
void Thermodynamics::ZeroPosStartingConditions(int i, double p, double h, double t, double d, double a, double v, double wbt)
{
  wbzeropos = i;
  zeropos = i;
  mintenpos = i;
  min40pos = i;
  min25pos = i;
  zero = abs(t);
  wb0 = abs(wbt);
  minten = abs(t+10);
  min40 = abs(t+40);
  min25 = abs(t+25);
}
void Thermodynamics::putZeroPos(int i, double p, double h, double t, double d, double a, double v, double wbt)
{
  double t_zero = abs(t);
  if (t_zero < zero)
  {
    zero = t_zero;
    zeropos = i;
  }
  double t_wb0 = abs(wbt);
  
  double imten= abs(t+10);
  
  if(imten<minten){
    minten = imten;
    mintenpos = i;
  }
  
  double im40= abs(t+40);
  
  if(im40<min40){
    min40 = im40;
    min40pos = i;
  }

  double im25= abs(t+25);
  
  if(im25<min25){
    min25 = im25;
    min25pos = i;
  }
  
  if (t_wb0 < wb0)
  {
    wb0 = t_wb0;
    wbzeropos = i;
  }
}
void Thermodynamics::putSpecificLine(int i, double p, double h, double t, double d, double a, double v)
{
  double syf = 0;
  double oe=OE(t, d, p);
  double wbt=TW(t, d, p, &syf);
  double mr=W(d, p);
  double virtt = tv(t,mr);

  double MSE0_ = 0;
  double rsat_ = 0;
  double MSE0_star_ = 0;
  double qsat_ = 0; 
  double MSE0_bar_ = 0;
  double int_arg_MSE0_ = 0;
  
  if(i > 0){
      
  MSE0_ = cp * (t+kel) + xlv * (mr/1000) + g * (h-h0);
  rsat_ = compute_rsat(t+kel,p*100,0);
  qsat_ = (1 - rsat_) * rsat_;
  MSE0_star_ = cp * (t+kel) + xlv * qsat_ + g * (h-h0);
    
  if(i == 1){
     last_MSE0 = MSE0_;
  }
    
  aggregated_MSE0 += (last_MSE0 + MSE0_) * (h-lasth);  

  if(i == 1){
     lasth = h;
  }

  MSE0_bar_ = 0.5 * aggregated_MSE0 / (h-h0);
  last_MSE0 = MSE0_;

  int_arg_MSE0_ = -(g/(cp*(t+kel)))*( MSE0_bar_ - MSE0_star_);
  }
  
  this->MSE0->push_back(MSE0_);
  this->MSE0_star->push_back(MSE0_star_);
  this->MSE0_bar->push_back(MSE0_bar_);
  this->int_arg_MSE0->push_back(int_arg_MSE0_);
  this->wbt->push_back(wbt);
  this->oe->push_back(oe);
  this->mixing->push_back(mr);
  this->virt->push_back(virtt);
  if (i >= 0)
  {
    putMaxTHTE(i, p, h, t, d, a, v, oe, mr);
  }

  if (i == 0)
  {
    startConditions(i, p, h, t, d, a, v, oe);
    ZeroPosStartingConditions(i, p, h, t, d, a, v, wbt);
    putMeanLayerParameters(i, p, h, t, d, a, v, mr);	
  }
  else
  {
    putMinTHTE(i, p, h, oe);
    putMeanLayerParameters(i, p, h, t, d, a, v, mr);
    putPWATER(i, p, h, t, d, a, v);
    putLowLapseRates(i, p, h, t, d, a, v);
    putZeroPos(i, p, h, t, d, a, v, wbt);
    if((fmod(abs(h-h0),200.0)==0.0)||h==h0){
      if (abs(h - h0) <= 2000){
        meanhum2+=ESAT(d)/ESAT(t);
        meand2b+=1;
        meanmxr2+=W(d, p);
        meand2+=1;
      }
      
      if (abs(h - h0) <= 1000){
        meanhum1+=ESAT(d)/ESAT(t);
        meand1b+=1.0;
      }
      
      if (abs(h - h0) >= 2000 && abs(h - h0) <= 5000){
        meanhum25+=ESAT(d)/ESAT(t);
        meand25+=1.0;
      }
      
      if (abs(h - h0) >= 3000 && abs(h - h0) <= 6000){
        meanhum36+=ESAT(d)/ESAT(t);
        meand36+=1.0;
      }
      
      if (abs(h - h0) >= 1000 && abs(h - h0) <= 4000){
        meanhum14+=ESAT(d)/ESAT(t);
        meand14+=1.0;
      }
    }
    
    if (t<=-10&&t>=-40){
      meanhumMIDDLE+=ESAT(d)/ESAT(t);
      meandMIDDLE+=1;
    }

   if (p<=850&&p>=500){
      meanhum850500+=ESAT(d)/ESAT(t);
      meand850500+=1;
    }

    
  }
  if (abs(h - h0) <= 1000 && mr1000<mr) mr1000 = mr;
  
  surfaceBased->putLine(i, p, h, t, d, a, v);
  mostUnstable->putLine(i, p, h, t, d, a, v);
  mostU500->putLine(i, p, h, t, d, a, v);
  lastp = p;
  lastt = t;
  lasth = h;
}
void Thermodynamics::prepareMeanLayer()
{
  if (n == 0) n = 1;
  mh /= n;
  mp = mp; 
  mt /= n;
  md /= n;
  mmr /= n;
  mo /=n;
  n /= 2;
  downmr/=downmrn;
  downo/=downon;
  thetd/=thetn;
  thet01d/=thet01n;
  thet02d/=thet02n;
  mthet/=mthetn; 
  meanLayer->setInitialConditions(0, mp, mh, mt, md, 0, 0, h0);
  meanLayer->setInitialW(mmr, mo);
  downdraft->setInitialConditions(0, 0, 0, 0, 0, 0, 0, h0);
  downdraft->prepareForDCAPE();
  downdraft->setInitialW(downmr, downo);
}
void Thermodynamics::putMeanLine(int i, double p, double h, double t, double d, double a, double v)
{
  this->meanLayer->putLine(i, p, h, t, d, a, v);
  this->meanmostUnstable->putLine(i, p, h, t, d, a, v);
  determineDowndraftByMinTHTE(i, p, h, t, d, a, v);
  putShowalter(i, p, h, t, d, a, v);
}

void Thermodynamics::finish(){
  this->mostUnstable->finish();
  this->mostU500->finish();
  this->surfaceBased->finish();
  this->meanLayer->finish();
  this->meanmostUnstable->finish();
  this->downdraft->finish();
  this->finishLowLapseRates();
  pwater /= 98.1;
  pwater_eff /= 98.1;
  
  meanhum1/=meand1b;
  meanhum2/=meand2b;
  meanhum25/=meand25; 
  meanhum36/=meand36;
  meanhum14/=meand14;
  meanhumMIDDLE/=meandMIDDLE; 
  meanhum850500/=meand850500; 
  meanmxr2/=meand2;
}

class Sounding{
  friend class IndicesCollector;
private:
  
  void alloc();
  void free();
  bool checkarguments(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_);
  void insertLine(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_, int i, double dz);
  inline void insertSingleLine(double p, double h, double t, double d, Vector V);
  inline void prepareCache(double p, double h);
  inline void prepareElementaryCache(double lval, double rval, double * rarr, int rindex, list<double>* vlist, void(*pointer)(int, int, Cache*), Cache* C);
  void finish();
  void secondPhase();
  public:
    double MU_ECAPE;
    double ML_ECAPE;
    double SB_ECAPE;
    double MU_ML_ECAPE;
    double MU500_ECAPE;
    Thermodynamics *th;
    Cache* cache;
    Kinematics *ks;
    list<double>* p;
    list<double>* h;
    list<double>* t;
    list<double>* d;
    list<double>* a;
    list<double>* v;
    IndicesCollector *ic;
    Sounding(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_, int length, double dz, double* meanlayer_bottom_top, Vector storm_motion);
    ~Sounding();
    IndicesCollector * getIndicesCollectorPointer(){
      return this->ic;
    }
};



#define NULLPOINTER (void*)0
int checkCrossing(double v1, double v2, double cv){
  if(v1==cv)return 0;
  else if(v2==cv)return 1;
  else{
    if(((v1>cv)&&(v2<cv))||((v1<cv)&&(v2>cv))){
      double min1 = abs(v1-cv);
      double min2 = abs(v2-cv);
      if(min1<=min2)return 0;
      else return 1;
    }
    return -1;
  }
  return -1;
}
bool Sounding::checkarguments(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_){
  if(p_==NULLPOINTER||h_==NULLPOINTER||t_==NULLPOINTER||d_==NULLPOINTER||a_==NULLPOINTER||v_==NULLPOINTER)return false;
  return true;
}
class IndicesCollector{
private:
  Thermodynamics *th;
  Cache* cache;
  Kinematics *ks;
  Sounding *S;
public:
  IndicesCollector(Thermodynamics *t,	Cache *c, Kinematics *k, Sounding *Snd);
  
  double VMostUnstableCAPE();
  double VLLMostUnstableCAPE();
  double MUmiddlecape();
  double VMostUnstableCIN();		
  double VMostUnstableLCL();
  double VMostUnstableLFC();
  double VMostUnstableEL();
  double VMostUnstableLI();
  double VMostUnstableVmax();
  double MUELTemperature();
  double MULCLTemperature();
  double MULFCTemperature();
  double MUMRatio();			

  double VMeanMostUnstableCAPE();
  double VLLMeanMostUnstableCAPE();
  double MUMLmiddlecape();
  double VMeanMostUnstableCIN();		
  double VMeanMostUnstableLCL();
  double VMeanMostUnstableLFC();
  double VMeanMostUnstableEL();
  double VMeanMostUnstableLI();
  double VMeanMostUnstableVmax();
  double MUMLELTemperature();
  double MUMLLCLTemperature();
  double MUMLLFCTemperature();
  double MUMLMRatio();			

  double VSurfaceBasedCAPE();
  double VLLSurfaceBasedCAPE();
  double SBmiddlecape();
  double VSurfaceBasedCIN();		
  double VSurfaceBasedLCL();
  double VSurfaceBasedLFC();
  double VSurfaceBasedEL();
  double VSurfaceBasedLI();
  double VSurfaceBasedVmax();	
  double SBELTemperature();
  double SBLCLTemperature();
  double SBLFCTemperature();	
  double SBMRatio();

  double VMeanLayerCAPE();
  double VLLMeanLayerCAPE();
  double MLmiddlecape();
  double VMeanLayerCIN();		
  double VMeanLayerLCL();
  double VMeanLayerLFC();
  double VMeanLayerEL();
  double VMeanLayerLI();
  double VMeanLayerVmax();
  double MLELTemperature();
  double MLLCLTemperature();
  double MLLFCTemperature();	
  double MLMixingRatio();

  double LapseRate01();
  double lapserate03();
  double LapseRate24();	
  double LR36();
  double lapseRate500700();
  double lapseRate500800();	
  double lapseRate600800();
  double LR0500();
  double LR02();
  double LR04();
  double LR06();  
  double LR16();
  double LR26();
  double max_LR26_2km();

  double THTE_LR_LCL_to_M10();
  double THTE_LR_MU_to_M10();
  double THTE_LR14();
  double THTE_LR13();
  double THTE_LR04();
  double THTE_LR03();

  double ZeroHeight();
  double M10Height();
  double WetBulbZeroHeight();		

  double MUHeight();
  double MUMLHeight();

  double MinTHTEHeight();
  double DeltaThetaE();
  double DeltaThetaE_min04km();
  double thetae01();
  double thetae02();
  double VDCAPE();	
  double VirtualColdPoolStrength();
  double WindIndex();
  double PWATER();
  double PWATER_eff();
  double MoistureFlux();

  double RH01();
  double RH02();
  double RH25();
  double RH14();
  double RH36();
  double RHMIDDLE();		
  double RH850500();		

  double BS06_var_SD();
  double BS06_var_SI();
  double BS16_var_SD();
  double BS16_var_SI();

  double BS500();
  double BS01();
  double BS02();
  double BS03();
  double BS04();
  double BS06();
  double BS08();
  double BS36();
  double BS010();
  double BS110();
  double BS18();
  double BS16();
  double BS13();
  double BS14();
  double BS25();

  double emubs();
  double emumlbs();
  double esbbs();
  double emlbs();
  double emu500bs();

  double BulkShearSfcTen();
  double BulkShear1kmTen();
  double BulkShear2kmTen();

  double BulkShearMUMLLCLTen();
  double BulkShearMULCLTen();
  double BulkShearSBLCLTen();
  double BulkShearMLLCLTen();

  double MeanWind500();
  double MeanWind01();
  double MeanWind02();
  double MeanWind03();
  double MeanWind06();
  double MeanWind13();
  double MeanWind36();
  double MeanWind69();
  double MeanWind912();

  double SRH100RM_F();
  double SRH100LM_F();
  double SRH500RM_F();
  double SRH500LM_F();

  double SW100_RM_F();
  double SW100_LM_F();
  double SW500_RM_F();
  double SW500_LM_F();

  double SRH100RM();
  double SRH250RM();
  double SRH500RM();
  double SRH01RM();
  double SRH03RM();
  double SRH36RM();
  double SRH13_RM();
  double SRH16RM();

  double SRH100LM();
  double SRH250LM();
  double SRH500LM();
  double SRH01LM();
  double SRH03LM();
  double SRH36LM();
  double SRH13_LM();
  double SRH16LM();

  double SRH01SM();
  double SRH03SM();
  double SRH36SM();
  double SRH13_SM();
  double SRH16SM();

  double SRH01LM_eff();
  double SRH01RM_eff();
  double SRH01SM_eff();
  double SRH03LM_eff();
  double SRH03RM_eff();
  double SRH03SM_eff();

  double Bunkers_RM_A();
  double Bunkers_RM_M();
  double Bunkers_LM_A();
  double Bunkers_LM_M();
  double Bunkers_MW_A();
  double Bunkers_MW_M();

  Vector Bunkers4_RM_vector();
  Vector Bunkers4_LM_vector();
  double Bunkers4_RM_A();
  double Bunkers4_RM_M();
  double Bunkers4_LM_A();
  double Bunkers4_LM_M();

  double K_Index();
  double Showalter_Index();	
  double TotalTotals();		
  double SWEATIndex();	
  double STP();
  double STPeff();
  double SCP();
  double SCPeff();	
  double STP_LM();
  double STPeff_LM();
  double SCP_LM();
  double SCPeff_LM();
  double SHERBS3();
  double SHERBE();
  double SHERBS3_v2();
  double SHERBE_v2();
  double DEI();	 	
  double DEI_eff();
  double TIP();	 
  double SHP();
  double DCP(); 
  double DCP_eff();

  double MU_WMAXSHEAR();
  double MUML_WMAXSHEAR();
  double SB_WMAXSHEAR();
  double ML_WMAXSHEAR();
  double MU500_WMAXSHEAR();

  double MU_EFF_EWMAXSHEAR();
  double MUML_EFF_EWMAXSHEAR();
  double SB_EFF_EWMAXSHEAR();
  double ML_EFF_EWMAXSHEAR();
  double MU500_EFF_EWMAXSHEAR();

  double MU_EFF_EWMAXSHEAR_HGL();
  double MUML_EFF_EWMAXSHEAR_HGL();
  double SB_EFF_EWMAXSHEAR_HGL();
  double ML_EFF_EWMAXSHEAR_HGL();
  double MU500_EFF_EWMAXSHEAR_HGL();

  double MU_EFF_EWMAXSHEAR_3km();
  double MUML_EFF_EWMAXSHEAR_3km();
  double SB_EFF_EWMAXSHEAR_3km();
  double ML_EFF_EWMAXSHEAR_3km();
  double MU500_EFF_EWMAXSHEAR_3km();

  double Corfidi_downwind_A();
  double Corfidi_downwind_M();
  double Corfidi_upwind_A();
  double Corfidi_upwind_M();
  
  double ML_coldcape();
  double SB_coldcape();
  double MU_coldcape();
  double MUML_coldcape();
  double MU500_coldcape();
  
  double ML_coldcapeTV();
  double SB_coldcapeTV();
  double MU_coldcapeTV();
  double MUML_coldcapeTV();
  double MU500_coldcapeTV();
  double HSI();
  double HSIv2();
  
  double MSR_MW();
  double MSR_RM();
  double MSR_LM();
  double MSR_MW_HGL();
  double MSR_RM_HGL();
  double MSR_LM_HGL();
  
  double MU500_CAPE();
  double MU500_CIN();
  double MU500_LI();

  double VLLMU500CAPE();
  double MU500middlecape();

  double EHI03();
  double EHI01();
  double EHI500();
  
  double EHI03_LM();
  double EHI01_LM();
  double EHI500_LM();

  double SW100_RM();
  double SW500_RM();
  double SW01_RM();
  double SW03_RM();

  double SW100_LM();
  double SW500_LM();
  double SW01_LM();
  double SW03_LM();

  double SV_100_RM_FRA();
  double SV_500_RM_FRA();
  double SV_1000_RM_FRA();
  double SV_3000_RM_FRA();

  double SV_100_LM_FRA();
  double SV_500_LM_FRA();
  double SV_1000_LM_FRA();
  double SV_3000_LM_FRA();
  
  double MeanSR500_RM();
  double MeanSR01_RM();
  double MeanSR03_RM();
  double MeanSR16_RM();
  double MeanSR36_RM();
  
  double MeanSR500_LM();
  double MeanSR01_LM();
  double MeanSR03_LM();
  double MeanSR16_LM();
  double MeanSR36_LM();

  double MeanSR500_MW();
  double MeanSR01_MW();
  double MeanSR03_MW();
  double MeanSR36_MW();

  double MeanSR0500_RM_eff();
  double MeanSR0500_LM_eff();
  double MeanSR0500_MW_eff();

  double MeanVMSR500_RM();
  double MeanVMSR01_RM();
  double MeanVMSR03_RM();
  double MeanVMSR16_RM();
  double MeanVMSR36_RM();
  
  double MeanVMSR500_LM();
  double MeanVMSR01_LM();
  double MeanVMSR03_LM();
  double MeanVMSR16_LM();
  double MeanVMSR36_LM();

  double Peters_SR_inflow();
  double Peters_SR_inflow_eff();

  Vector Peters_vector();
  double Peters_vector_A();
  double Peters_vector_M();

  double Peters_vector_eff_A();
  double Peters_vector_eff_M();

  double VSurfaceBasedLI_M25();
  double VMeanLayerLI_M25();
  double VMostUnstableLI_M25();
  double VMeanMostUnstableLI_M25();
  double VMostU500LI_M25();	

  double MU_buoyancy();
  double MUML_buoyancy();
  double MU500_buoyancy();
  double ML_buoyancy();
  double SB_buoyancy();

  double MU_buoyancy_M10();
  double MUML_buoyancy_M10();
  double MU500_buoyancy_M10();
  double ML_buoyancy_M10();
  double SB_buoyancy_M10();

  double MU_buoyancy_3km();
  double MUML_buoyancy_3km();
  double MU500_buoyancy_3km();
  double ML_buoyancy_3km();
  double SB_buoyancy_3km();

  double MU_ebuoyancy_3km();
  double MUML_ebuoyancy_3km();
  double MU500_ebuoyancy_3km();
  double ML_ebuoyancy_3km();
  double SB_ebuoyancy_3km();

  double MU_ebuoyancy();
  double MUML_ebuoyancy();
  double MU500_ebuoyancy();
  double ML_ebuoyancy();
  double SB_ebuoyancy();

  double MU_ebuoyancy_M10();
  double MUML_ebuoyancy_M10();
  double MU500_ebuoyancy_M10();
  double ML_ebuoyancy_M10();
  double SB_ebuoyancy_M10();

  double SR_moisture_flux();
  double SR_moisture_flux_eff();
  double SR_moisture_flux_RM();
  double SR_moisture_flux_eff_RM();
  double SR_moisture_flux_LM();
  double SR_moisture_flux_eff_LM();
  double SR_moisture_flux_MW();
  double SR_moisture_flux_eff_MW();

  double SB_warm_cloud();
  double SB_cold_cloud();
  double SB_equal_layer();

  double ML_warm_cloud();
  double ML_cold_cloud();
  double ML_equal_layer();

  double MU_warm_cloud();
  double MU_cold_cloud();
  double MU_equal_layer();

  double MU_ML_warm_cloud();
  double MU_ML_cold_cloud();
  double MU_ML_equal_layer();

  double* MU_ECAPE();
  double* ML_ECAPE();
  double* MU_ML_ECAPE();
  double* SB_ECAPE();
  double* MU500_ECAPE();

  double MU_ELI();
  double ML_ELI();
  double MUML_ELI();
  double SB_ELI();
  double MU500_ELI();

  double Ventilation_16km();
  double Ventilation_36km();
  double Ventilation_69km();

  double Ventilation_16km_LM();
  double Ventilation_36km_LM();
  double Ventilation_69km_LM();

  double Ventilation_16km_RM();
  double Ventilation_36km_RM();
  double Ventilation_69km_RM();

  double Ventilation_16km_MW();
  double Ventilation_36km_MW();
  double Ventilation_69km_MW();

  double SRW_sfc_LM();
  double SRW_sfc_RM();
  double LTTP_RM();
  double LTTP_LM();
  double CA500_RM();
  double CA500_LM();
  double BS_LLmax();
  double BS_MLmax();
  double BS_ULmax();
  double WS_LLmax();
  double WS_MLmax();
  double WS_ULmax();

  double VMeanLayerCIN500();		
  double VMeanMostUnstableCIN500();		
  double VSurfaceBasedCIN500();		
  double VMostUnstableCIN500();		
  double MU500_CIN500();		

};

void Sounding::alloc(){
  this->p=new list<double>();
  this->h=new list<double>();
  this->t=new list<double>();
  this->d=new list<double>();
  this->a=new list<double>();
  this->v=new list<double>();
  this->cache = new Cache();
  this->th = new Thermodynamics();
  th->setSoundingCache(cache);
  this->ks=new Kinematics();
  ks->setSoundingCache(cache);
  this->ic=new IndicesCollector(this->th,this->cache,this->ks, this);
  ks->S=this;
}
void Sounding::free(){
  delete(this->p);
  delete(this->h);
  delete(this->t);
  delete(this->d);
  delete(this->a);
  delete(this->v);
  
  delete(this->ic);
  
  
  delete(this->th);
  delete(this->ks);
  delete(this->cache);
}
void Sounding::prepareElementaryCache(double lval, double rval, double * rarr, int rindex, list<double>* vlist, void(*pointer)(int, int, Cache*), Cache* C){
  int check1;
  int index;
  check1 = checkCrossing(lval,rval,rarr[rindex]);
  if(check1==0){
    index= (int)(vlist->size()-1);
    (pointer)(rindex,index, C);
    
  }
  else if(check1==1){
    index=(int)(vlist->size());
    (pointer)(rindex,index,C);
    
  }
}
void Sounding::prepareCache(double p, double h){
  int	plength;
  double *parray = this->cache->getArray(0,&plength);
  
  int	hlength;
  double *harray = this->cache->getArray(1,&hlength);
  
  double lp =*(--this->p->end());
  double lh =*(--this->h->end());
  
  double AGLlh=lh-this->cache->getH0();
  double AGLh=h-this->cache->getH0();
  for(int i=0;i<hlength;i++){
    if(this->h->size()>0)this->prepareElementaryCache(AGLlh,AGLh,harray,i,this->h,&(setHeightIndex),this->cache);
    if(i<plength)
      if(this->p->size()>0)this->prepareElementaryCache(lp,p,parray,i,this->p,&(setPressureIndex),this->cache);
  }
  
}
void Sounding::insertSingleLine(double p,double h, double t,double d, Vector V){
  this->prepareCache(p,h);
  int i =(this->p->size())-1;
  
  this->p->push_back(p);
  this->h->push_back(h);
  
  this->t->push_back(t);
  this->d->push_back(d);
  
  double* av = V.toAV(); 
  this->a->push_back(av[0]);
  this->v->push_back(av[1]/0.514444);
  
  this->th->putLine(i+1, p, h, t, d, av[0], av[1]/0.514444);
  this->ks->putLine(i+1, p, h, t, d, av[0], av[1]/0.514444);
  delete[] av;
}
void Sounding::insertLine(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_, int i, double dz){
  double tp = p_[i], th = h_[i], tt = t_[i], td = d_[i], ta = a_[i], tv = v_[i];
  double np = p_[i+1], nh = h_[i+1], nt = t_[i+1], nd = d_[i+1], na = a_[i+1], nv = v_[i+1]; 
  
  Vector tvec = Vector(ta,tv*0.514444);
  Vector nvec = Vector(na,nv*0.514444);
  
  double dh = nh-th;
  double dp = (np-tp)/dh;
  double dt = (nt-tt)/dh;
  double dd = (nd-td)/dh;
  Vector dVec = (nvec-tvec)/dh;
  
  dp*=dz;
  dt*=dz;
  dd*=dz;
  dVec*=dz;
  int how = (int)(abs(floor(nh-th)/dz)); 
  dh=dz;
  double _p,_h,_t,_d; Vector _V= Vector(tvec);
  for(int j=0;j<=how;j++){
    _p=tp+((double)j*dp);
    _h=th+((double)j*dh);
    _t=tt+((double)j*dt);
    _d=td+((double)j*dd);
    _V=tvec+(dVec*(double)j);
    if(_h!=nh)
      this->insertSingleLine(_p,_h,_t,_d,_V);
  }
  
}
Sounding::Sounding(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_, int length, double dz, double* meanlayer_bottom_top, Vector storm_motion){
  this->alloc();
  bool valid = this->checkarguments(p_,h_,t_,d_,a_,v_);
  if(valid){
    this->th->meanLayerBottom = meanlayer_bottom_top[0];
    this->th->meanLayerTop = meanlayer_bottom_top[1];
    int i=0;
    this->cache->setH0(h_[0]);
    
    double diff_mlb = this->th->meanLayerBottom;
    int mlb_index = 0; 
    double temp_value = 0;
    
    for (i=0; i<length-1; i++){ 
      temp_value = abs(h_[i]-h_[0]-this->th->meanLayerBottom);
      if(temp_value < diff_mlb){
        diff_mlb=temp_value;
        mlb_index = i;
      }
      this->insertLine(p_,h_,t_,d_,a_,v_,i,dz);
    }
    this->insertSingleLine(p_[length-1],h_[length-1],t_[length-1],d_[length-1],Vector(a_[length-1],v_[length-1]*0.514444));
    
    
    this->th->prepareMeanLayer();
    this->ks->finishPhase1();
    this->th->meanLayer->startIndex = mlb_index; 
    if(storm_motion.Z() != 999) this->ks->rm = storm_motion;
    if(storm_motion.Z() != 999) this->ks->lm = storm_motion;
    this->secondPhase();
    this->finish();
  }
  
  
}
Sounding::~Sounding(){
  this->free();
}
void Sounding::finish(){
  th->finish();
  ks->finish();
  
}

void Sounding::secondPhase(){
  this->ks->muheight = Get(this->h,this->th->mostUnstable->startIndex);
  this->ks->mumlheight = Get(this->h,this->th->meanmostUnstable->startIndex);
  list<double>::iterator ip;
  list<double>::iterator MSE0i = this->th->int_arg_MSE0->begin();
  list<double>::iterator ih = this->h->begin();
  list<double>::iterator it = this->t->begin();
  list<double>::iterator id = this->d->begin();
  list<double>::iterator ia = this->a->begin();
  list<double>::iterator iv = this->v->begin();
  int i=0;
  for(ip = this->p->begin(); ip!=this->p->end(); ++ip){
    double p_ = *ip;
    double h_ = *ih;
    double t_ = *it;
    double d_ = *id;
    double a_ = *ia;
    double v_ = *iv;
    this->th->putMeanLine(i, p_, h_, t_, d_, a_, v_);
    this->ks->putSecondPhaseLine(i, p_, h_, t_, d_, a_, v_);
    ++ih;++it;++id;++ia;++iv;++i;
  }
  i=0;
  
  ih = this->h->begin();
  it = this->t->begin();
  id = this->d->begin();
  ia = this->a->begin();
  iv = this->v->begin();
  double h0=*ih;
  this->th->downdraft->lasth=h0;
  
  std::list<double> vals (*this->th->downdraft->values); 
  std::list<double> virtvals (*this->th->downdraft->virtualValues); 
  
  this->th->downdraft->values->clear();
  this->th->downdraft->virtualValues->clear();
  
  for(ip = this->p->begin(); ip!=this->p->end(); ++ip){
    
    double p_ = *ip;
    double h_ = *ih;
    double t_ = *it;
    double d_ = *id;
    double a_ = *ia;
    double v_ = *iv;
    double mse0i = *MSE0i;
    this->th->mostUnstable->putNCAPE(i, mse0i);
    this->th->mostU500->putNCAPE(i, mse0i);
    this->th->surfaceBased->putNCAPE(i, mse0i);
    this->th->meanLayer->putNCAPE(i, mse0i);
    this->th->meanmostUnstable->putNCAPE(i, mse0i);

   if(h_-h0<4000)this->th->downdraft->putLine(i, p_, h_, t_, d_, a_, v_);
    ++ih;++it;++id;++ia;++iv;++i;++MSE0i;
  }
  
  list<double>::iterator vv=vals.begin(); list<double>::iterator vvv=virtvals.begin();
  while(vv!=vals.end()){
    double u = *vv;
    double w = *vvv;
    this->th->downdraft->values->push_back(u);
    this->th->downdraft->virtualValues->push_back(w);
    vv++;vvv++;
  }
  this->ks->effvec();
}

IndicesCollector::IndicesCollector(Thermodynamics *t, Cache *c, Kinematics *k,Sounding *Snd){
  this->th=t;
  this->ks=k;
  this->cache=c;
  this->S=Snd;
}


double IndicesCollector::K_Index() {
  double val = 0;
  int i500 = cache->getPressureIndex(500);
  int i700 = cache->getPressureIndex(700);
  int i850 = cache->getPressureIndex(850);
  list<double> *t = S->t;
  list<double> *d = S->d;
  val = (Get(t,i850) - Get(t,i500) + Get(d,i850) - (Get(t,i700) - Get(d,i700)));
  return val;
}

double IndicesCollector::Showalter_Index(){
  int lindex = cache->getPressureIndex(500);
  double lit = Get(S->t,lindex);
  int uindex = S->th->showalter->startIndex;
  int vindex = lindex - uindex;
  double plit = Get(S->th->showalter->virtualValues,vindex);
  double Showalter = lit - plit;
  return Showalter;
}

double IndicesCollector::TotalTotals(){
  int i500 = cache->getPressureIndex(500);
  int i850 = cache->getPressureIndex(850); 
  double d850 = Get(S->d,i850);
  double t850 = Get(S->t,i850);
  double t500 = Get(S->t,i500);
  double VerticalTotals = d850-t500;
  double CrossTotals = t850-t500;
  double result = VerticalTotals + CrossTotals;	
  return result;
}

double IndicesCollector::SWEATIndex(){
  int i500 = cache->getPressureIndex(500);
  int i850 = cache->getPressureIndex(850);
  double a500 = Get(S->a,i500);
  double a850 = Get(S->a,i850);
  double v500 = Get(S->v,i500);
  double v850 = Get(S->v,i850);
  double d850 = Get(S->d,i850);
  double angle = a500 - a850;
  double angle2 = angle / 180;
  double sinus = sin(angle2 * M_PI);
  double shear = 125 * (sinus + 0.2);
  if (a500 < 210 || a500 > 310 || a850 < 130 || a850 > 250 || a500 - a850 < 0 || v500 < 15 || v850 < 15) shear = 0;
  double term = this->TotalTotals() - 49;
  if (term < 0) term = 0;
  double x = 12 * d850;
  if (x < 0) x = 0;
  double res = x + 20 * term + 2 * v850 + v500 + shear;
  return res;
}

double* IndicesCollector::MU_ECAPE(){  
  double L = 250;
  double EL = Get(S->h, S->th->mostUnstable->vElIndex) - S->th->h0;
  double LFC = Get(S->h, S->th->mostUnstable->vLfcIndex) - S->th->h0;
  double CAPE = S->th->mostUnstable->vcape;
  double V_SR = Peters_SR_inflow_eff();
  double l = L/EL;
  double NCAPE = 0;
  double pitchfork = ksq*(alpha*alpha)*(M_PI*M_PI)*L/(4*Pr*(sigma*sigma)*EL);    
  double vsr_tilde = V_SR/sqrt(2*CAPE);
  NCAPE = S->th->mostUnstable->NCAPE / S->th->mostUnstable->NCAPE_N;  
  NCAPE *=  EL - LFC; 
  if(NCAPE < 0) NCAPE = 0;
  double N_tilde = NCAPE / CAPE;   
  double E_tilde = vsr_tilde*vsr_tilde + ( -1 - pitchfork - (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde + sqrt((((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))*((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))) + (4*(pitchfork/(vsr_tilde*vsr_tilde))*(1 - pitchfork*N_tilde) ) ) )/( 2*pitchfork/(vsr_tilde*vsr_tilde) );
  double eps = 2*ksq*L/(EL*Pr);
  double E_tilde_ = E_tilde - vsr_tilde*vsr_tilde;
  double varepsilon = 2*((1 - E_tilde_) / (E_tilde_ + N_tilde))/(EL);
  double Radius = sqrt(2*ksq*L/(Pr*varepsilon));
  double CAPE_HGL = S->th->mostUnstable->middlecape;
  double CAPE_M10 = S->th->mostUnstable->coldcape;
  double CAPE_3km = S->th->mostUnstable->vto3cape;
  if(E_tilde>10) E_tilde = 10;
  if(isnan(E_tilde)) E_tilde = 0;
  double* result = new double[7];
  result[0] = E_tilde;
  result[1] = Radius;
  result[2] = E_tilde*CAPE;
  result[3] = E_tilde*CAPE_HGL;
  result[4] = E_tilde*CAPE_M10;
  result[5] = sqrt(2*E_tilde*CAPE);
  result[6] = E_tilde*CAPE_3km;
  return result;
}

double* IndicesCollector::MU_ML_ECAPE(){  
  double L = 250;
  double EL = Get(S->h, S->th->meanmostUnstable->vElIndex) - S->th->h0;
  double LFC = Get(S->h, S->th->meanmostUnstable->vLfcIndex) - S->th->h0;
  double CAPE = S->th->meanmostUnstable->vcape;
  double V_SR = Peters_SR_inflow_eff();
  double l = L/EL;
  double NCAPE = 0;
  double pitchfork = ksq*(alpha*alpha)*(M_PI*M_PI)*L/(4*Pr*(sigma*sigma)*EL);    
  double vsr_tilde = V_SR/sqrt(2*CAPE);
  NCAPE = S->th->meanmostUnstable->NCAPE / S->th->meanmostUnstable->NCAPE_N;  
  NCAPE *=  EL - LFC; 
  if(NCAPE < 0) NCAPE = 0;
  double N_tilde = NCAPE / CAPE;   
  double E_tilde = vsr_tilde*vsr_tilde + ( -1 - pitchfork - (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde + sqrt((((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))*((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))) + (4*(pitchfork/(vsr_tilde*vsr_tilde))*(1 - pitchfork*N_tilde) ) ) )/( 2*pitchfork/(vsr_tilde*vsr_tilde) );
  double eps = 2*ksq*L/(EL*Pr);
  double E_tilde_ = E_tilde - vsr_tilde*vsr_tilde;
  double varepsilon = 2*((1 - E_tilde_) / (E_tilde_ + N_tilde))/(EL);
  double Radius = sqrt(2*ksq*L/(Pr*varepsilon));
  double CAPE_HGL = S->th->meanmostUnstable->middlecape;
  double CAPE_M10 = S->th->meanmostUnstable->coldcape;
  double CAPE_3km = S->th->meanmostUnstable->vto3cape;
  if(E_tilde>10) E_tilde = 10;
  if(isnan(E_tilde)) E_tilde = 0;
  double* result = new double[7];
  result[0] = E_tilde;
  result[1] = Radius;
  result[2] = E_tilde*CAPE;
  result[3] = E_tilde*CAPE_HGL;
  result[4] = E_tilde*CAPE_M10;
  result[5] = sqrt(2*E_tilde*CAPE);
  result[6] = E_tilde*CAPE_3km;
  return result;
}

double* IndicesCollector::SB_ECAPE(){  
  double L = 250;
  double EL = Get(S->h, S->th->surfaceBased->vElIndex) - S->th->h0;
  double LFC = Get(S->h, S->th->surfaceBased->vLfcIndex) - S->th->h0;
  double CAPE = S->th->surfaceBased->vcape;
  double V_SR = Peters_SR_inflow();
  double l = L/EL;
  double NCAPE = 0;
  double pitchfork = ksq*(alpha*alpha)*(M_PI*M_PI)*L/(4*Pr*(sigma*sigma)*EL);    
  double vsr_tilde = V_SR/sqrt(2*CAPE);
  NCAPE = S->th->surfaceBased->NCAPE / S->th->surfaceBased->NCAPE_N;  
  NCAPE *=  EL - LFC; 
  if(NCAPE < 0) NCAPE = 0;
  double N_tilde = NCAPE / CAPE;   
  double E_tilde = vsr_tilde*vsr_tilde + ( -1 - pitchfork - (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde + sqrt((((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))*((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))) + (4*(pitchfork/(vsr_tilde*vsr_tilde))*(1 - pitchfork*N_tilde) ) ) )/( 2*pitchfork/(vsr_tilde*vsr_tilde) );
  double eps = 2*ksq*L/(EL*Pr);
  double E_tilde_ = E_tilde - vsr_tilde*vsr_tilde;
  double varepsilon = 2*((1 - E_tilde_) / (E_tilde_ + N_tilde))/(EL);
  double Radius = sqrt(2*ksq*L/(Pr*varepsilon));
  double CAPE_HGL = S->th->surfaceBased->middlecape;
  double CAPE_M10 = S->th->surfaceBased->coldcape;
  double CAPE_3km = S->th->surfaceBased->vto3cape;
  if(E_tilde>10) E_tilde = 10;
  if(isnan(E_tilde)) E_tilde = 0;
  double* result = new double[7];
  result[0] = E_tilde;
  result[1] = Radius;
  result[2] = E_tilde*CAPE;
  result[3] = E_tilde*CAPE_HGL;
  result[4] = E_tilde*CAPE_M10;
  result[5] = sqrt(2*E_tilde*CAPE);
  result[6] = E_tilde*CAPE_3km;
  return result;
}

double* IndicesCollector::ML_ECAPE(){  
  double L = 250;
  double EL = Get(S->h, S->th->meanLayer->vElIndex) - S->th->h0;
  double LFC = Get(S->h, S->th->meanLayer->vLfcIndex) - S->th->h0;
  double CAPE = S->th->meanLayer->vcape;
  double V_SR = Peters_SR_inflow();
  double l = L/EL;
  double NCAPE = 0;
  double pitchfork = ksq*(alpha*alpha)*(M_PI*M_PI)*L/(4*Pr*(sigma*sigma)*EL);    
  double vsr_tilde = V_SR/sqrt(2*CAPE);
  NCAPE = S->th->meanLayer->NCAPE / S->th->meanLayer->NCAPE_N;  
  NCAPE *=  EL - LFC; 
  if(NCAPE < 0) NCAPE = 0;
  double N_tilde = NCAPE / CAPE;   
  double E_tilde = vsr_tilde*vsr_tilde + ( -1 - pitchfork - (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde + sqrt((((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))*((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))) + (4*(pitchfork/(vsr_tilde*vsr_tilde))*(1 - pitchfork*N_tilde) ) ) )/( 2*pitchfork/(vsr_tilde*vsr_tilde) );
  double eps = 2*ksq*L/(EL*Pr);
  double E_tilde_ = E_tilde - vsr_tilde*vsr_tilde;
  double varepsilon = 2*((1 - E_tilde_) / (E_tilde_ + N_tilde))/(EL);
  double Radius = sqrt(2*ksq*L/(Pr*varepsilon));
  double CAPE_HGL = S->th->meanLayer->middlecape;
  double CAPE_M10 = S->th->meanLayer->coldcape;
  double CAPE_3km = S->th->meanLayer->vto3cape;
  if(E_tilde>10) E_tilde = 10;
  if(isnan(E_tilde)) E_tilde = 0;
  double* result = new double[7];
  result[0] = E_tilde;
  result[1] = Radius;
  result[2] = E_tilde*CAPE;
  result[3] = E_tilde*CAPE_HGL;
  result[4] = E_tilde*CAPE_M10;
  result[5] = sqrt(2*E_tilde*CAPE);
  result[6] = E_tilde*CAPE_3km;
  return result;
}

double* IndicesCollector::MU500_ECAPE(){  
  double L = 250;
  double EL = Get(S->h, S->th->mostU500->vElIndex) - S->th->h0;
  double LFC = Get(S->h, S->th->mostU500->vLfcIndex) - S->th->h0;
  double CAPE = S->th->mostU500->vcape;
  double V_SR = Peters_SR_inflow_eff();
  double l = L/EL;
  double NCAPE = 0;
  double pitchfork = ksq*(alpha*alpha)*(M_PI*M_PI)*L/(4*Pr*(sigma*sigma)*EL);    
  double vsr_tilde = V_SR/sqrt(2*CAPE);
  NCAPE = S->th->mostU500->NCAPE / S->th->mostU500->NCAPE_N;  
  NCAPE *=  EL - LFC; 
  if(NCAPE < 0) NCAPE = 0;
  double N_tilde = NCAPE / CAPE;   
  double E_tilde = vsr_tilde*vsr_tilde + ( -1 - pitchfork - (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde + sqrt((((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))*((1 + pitchfork + (pitchfork/(vsr_tilde*vsr_tilde))*N_tilde))) + (4*(pitchfork/(vsr_tilde*vsr_tilde))*(1 - pitchfork*N_tilde) ) ) )/( 2*pitchfork/(vsr_tilde*vsr_tilde) );
  double eps = 2*ksq*L/(EL*Pr);
  double E_tilde_ = E_tilde - vsr_tilde*vsr_tilde;
  double varepsilon = 2*((1 - E_tilde_) / (E_tilde_ + N_tilde))/(EL);
  double Radius = sqrt(2*ksq*L/(Pr*varepsilon));
  double CAPE_HGL = S->th->mostU500->middlecape;
  double CAPE_M10 = S->th->mostU500->coldcape;
  double CAPE_3km = S->th->mostU500->vto3cape;
  if(E_tilde>10) E_tilde = 10;
  if(isnan(E_tilde)) E_tilde = 0;
  double* result = new double[7];
  result[0] = E_tilde;
  result[1] = Radius;
  result[2] = E_tilde*CAPE;
  result[3] = E_tilde*CAPE_HGL;
  result[4] = E_tilde*CAPE_M10;
  result[5] = sqrt(2*E_tilde*CAPE);
  result[6] = E_tilde*CAPE_3km;
  return result;
}

double IndicesCollector::SB_warm_cloud(){
  double LFC = Get(S->h, S->th->surfaceBased->vLfcIndex) - S->th->h0;
  double EL = Get(S->h, S->th->surfaceBased->vElIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos) - S->th->h0;
  if(FL > EL){
	  FL = EL;
  }
  double result = FL-LFC;
  if(result<0){
    result = 0;
  }
  double cape = this->VSurfaceBasedCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::ML_warm_cloud(){
  double LFC = Get(S->h, S->th->meanLayer->vLfcIndex) - S->th->h0;
  double EL = Get(S->h, S->th->meanLayer->vElIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos) - S->th->h0;
  if(FL > EL){
	  FL = EL;
  }
  double result = FL-LFC;
  if(result<0){
    result = 0;
  }
    double cape = this->VMeanLayerCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::MU_ML_warm_cloud(){
  double LFC = Get(S->h, S->th->meanmostUnstable->vLfcIndex) - S->th->h0;
  double EL = Get(S->h, S->th->meanmostUnstable->vElIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos) - S->th->h0;
  if(FL > EL){
	  FL = EL;
  }
  double result = FL-LFC;
  if(result<0){
    result = 0;
  }
    double cape = this->VMeanMostUnstableCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::MU_warm_cloud(){
  double LFC = Get(S->h, S->th->mostUnstable->vLfcIndex) - S->th->h0;
  double EL = Get(S->h, S->th->mostUnstable->vElIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos) - S->th->h0;
  if(FL > EL){
	  FL = EL;
  }  
  double result = FL-LFC;
  if(result<0){
    result = 0;
  }
    double cape = this->VMostUnstableCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::ML_cold_cloud(){
  double EL = Get(S->h, S->th->meanLayer->vElIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos) - Get(S->h,0);
  double result = EL-FL;
  if(result<0){
    result = 0;
  }
    double cape = this->VMeanLayerCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::SB_cold_cloud(){
  double EL = Get(S->h, S->th->surfaceBased->vElIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos)-Get(S->h,0);
  double result = EL-FL;
  if(result<0){
    result = 0;
  }
    double cape = this->VSurfaceBasedCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::MU_cold_cloud(){
  double EL = Get(S->h, S->th->mostUnstable->vElIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos)-Get(S->h,0);
  double result = EL-FL;
  if(result<0){
    result = 0;
  }
    double cape = this->VMostUnstableCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::MU_ML_cold_cloud(){
  double EL = Get(S->h, S->th->meanmostUnstable->vElIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos)-Get(S->h,0);
  double result = EL-FL;
  if(result<0){
    result = 0;
  }
    double cape = this->VMeanMostUnstableCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::ML_equal_layer(){
  double LFC = Get(S->h, S->th->meanLayer->vLfcIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos) - Get(S->h,0);
  double EL = Get(S->h, S->th->meanLayer->vElIndex) - S->th->h0;
  double cold = EL-FL;
  if(FL > EL){
	  FL = EL;
  } 
  double warm = FL-LFC;
  if(warm<0){
    warm = 0;
  }
  if(cold<0){
    cold = 0;
  }  
  double result = min(warm,cold);
    double cape = this->VMeanLayerCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::SB_equal_layer(){
  double LFC = Get(S->h, S->th->surfaceBased->vLfcIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos)-Get(S->h,0);
  double EL = Get(S->h, S->th->surfaceBased->vElIndex) - S->th->h0;
  double cold = EL-FL;
  if(FL > EL){
	  FL = EL;
  } 
  double warm = FL-LFC;
  if(warm<0){
    warm = 0;
  }
  if(cold<0){
    cold = 0;
  }  
  double result = min(warm,cold);
    double cape = this->VSurfaceBasedCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::MU_equal_layer(){
  double LFC = Get(S->h, S->th->mostUnstable->vLfcIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos)-Get(S->h,0);
  double EL = Get(S->h, S->th->mostUnstable->vElIndex) - S->th->h0;
  double cold = EL-FL;
  if(FL > EL){
	  FL = EL;
  } 
  double warm = FL-LFC;
  if(warm<0){
    warm = 0;
  }
  if(cold<0){
    cold = 0;
  }  
  double result = min(warm,cold);
    double cape = this->VMostUnstableCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::MU_ML_equal_layer(){
  double LFC = Get(S->h, S->th->meanmostUnstable->vLfcIndex) - S->th->h0;
  double FL = Get(S->h,S->th->mintenpos)-Get(S->h,0);
  double EL = Get(S->h, S->th->meanmostUnstable->vElIndex) - S->th->h0;
  double cold = EL-FL;
  if(FL > EL){
	  FL = EL;
  } 
  double warm = FL-LFC;
  if(warm<0){
    warm = 0;
  }
  if(cold<0){
    cold = 0;
  }  
  double result = min(warm,cold);
    double cape = this->VMeanMostUnstableCAPE();
  if(cape==0){
    result = 0; 
  }
  return result;  
}

double IndicesCollector::VMostUnstableCAPE(){  
  double result = 0;
  result = S->th->mostUnstable->vcape; 
  return result;
}

double IndicesCollector::VMeanMostUnstableCAPE(){  
  double result = 0;
  result = S->th->meanmostUnstable->vcape; 
  return result;
}

double IndicesCollector::VLLMostUnstableCAPE(){  
  double result = 0;
  result = S->th->mostUnstable->vto3cape;  
  return result;  
}

double IndicesCollector::VLLMeanMostUnstableCAPE(){  
  double result = 0;
  result = S->th->meanmostUnstable->vto3cape;  
  return result;  
}

double IndicesCollector::VMostUnstableCIN(){
  double result = 0;
  result = S->th->mostUnstable->vcin;
  double cape = this->VMostUnstableCAPE();
  if(cape==0){
    result = sqrt(-1); 
  }
  return result;
}   

double IndicesCollector::VMeanMostUnstableCIN(){
  double result = 0;
  result = S->th->meanmostUnstable->vcin;
  double cape = this->VMeanMostUnstableCAPE();
  if(cape==0){
    result = sqrt(-1); 
  }
  return result;
}   

double IndicesCollector::VMostUnstableCIN500(){
  double result = S->th->mostUnstable->vcin500;
  return result;
}   

double IndicesCollector::VMeanMostUnstableCIN500(){
  double result = S->th->meanmostUnstable->vcin500;
  return result;
}   

double IndicesCollector::VMostUnstableLCL(){  
  double result = 0;
  int index = S->th->mostUnstable->vLclIndex;
  result = Get(S->h,index)- S->th->h0;  
  return result;
}

double IndicesCollector::VMeanMostUnstableLCL(){  
  double result = 0;
  int index = S->th->meanmostUnstable->vLclIndex;
  result = Get(S->h,index)- S->th->h0;  
  return result;
}

double IndicesCollector::VMostUnstableLFC(){
  double result = 0;
  int index = S->th->mostUnstable->vLfcIndex;  
  result = Get(S->h,index)- S->th->h0;  
  double cape = this->VMostUnstableCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::VMeanMostUnstableLFC(){
  double result = 0;
  int index = S->th->meanmostUnstable->vLfcIndex;  
  result = Get(S->h,index)- S->th->h0;  
  double cape = this->VMeanMostUnstableCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::VMostUnstableEL(){
  double result = 0;
  int index = S->th->mostUnstable->vElIndex;
  result = Get(S->h,index)- S->th->h0;
  double cape = this->VMostUnstableCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
} 

double IndicesCollector::VMeanMostUnstableEL(){
  double result = 0;
  int index = S->th->meanmostUnstable->vElIndex;
  result = Get(S->h,index)- S->th->h0;
  double cape = this->VMeanMostUnstableCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
} 

double IndicesCollector::VMostUnstableLI(){
  int lindex = cache->getPressureIndex(500);  
  double lit = Get(S->t,lindex);
  int vindex = lindex - S->th->mostUnstable->startIndex;  
  double plit = Get(S->th->mostUnstable->virtualValues,vindex);
  double Showalter = lit - plit;
  return Showalter;
}

double IndicesCollector::VMeanMostUnstableLI(){
  int lindex = cache->getPressureIndex(500);  
  double lit = Get(S->t,lindex);
  int vindex = lindex - S->th->meanmostUnstable->startIndex;  
  double plit = Get(S->th->meanmostUnstable->virtualValues,vindex);
  double Showalter = lit - plit;
  return Showalter;
}

double IndicesCollector::VMostUnstableVmax(){
  return sqrt(this->VMostUnstableCAPE()*2);
}

double IndicesCollector::VMeanMostUnstableVmax(){
  return sqrt(this->VMeanMostUnstableCAPE()*2);
}

double IndicesCollector::MUELTemperature(){
  double result = Get(S->t,S->th->mostUnstable->vElIndex);
  double cape = this->VMostUnstableCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::MUMLELTemperature(){
  double result = Get(S->t,S->th->meanmostUnstable->vElIndex);
  double cape = this->VMeanMostUnstableCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::MULCLTemperature(){
  return Get(S->t,S->th->mostUnstable->vLclIndex);
}

double IndicesCollector::MUMLLCLTemperature(){
  return Get(S->t,S->th->meanmostUnstable->vLclIndex);
}

double IndicesCollector::MULFCTemperature(){
  double result = Get(S->t,S->th->mostUnstable->vLfcIndex);
  double cape = this->VMostUnstableCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::MUMLLFCTemperature(){
  double result = Get(S->t,S->th->meanmostUnstable->vLfcIndex);
  double cape = this->VMeanMostUnstableCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::VSurfaceBasedCAPE(){  
  double result = 0;
  result = S->th->surfaceBased->vcape;
  return result;
}

double IndicesCollector::VLLSurfaceBasedCAPE(){
  double result = 0;
  result = S->th->surfaceBased->vto3cape;
  return result;
}

double IndicesCollector::VSurfaceBasedCIN(){  
  double result = 0;
  result = S->th->surfaceBased->vcin;
    double cape = this->VSurfaceBasedCAPE();
  if(cape==0){
    result = sqrt(-1); 
  }
  return result;  
}   

double IndicesCollector::VSurfaceBasedCIN500(){  
  double result = S->th->surfaceBased->vcin500;
  return result;  
}   

double IndicesCollector::VSurfaceBasedLCL(){  
  double result = 0;
  int index = S->th->surfaceBased->vLclIndex;  
  result = Get(S->h,index)- S->th->h0;  
  return result;
}

double IndicesCollector::VSurfaceBasedLFC(){  
  double result = 0;
  int index = S->th->surfaceBased->vLfcIndex;  
  result = Get(S->h,index)- S->th->h0;
  double cape = this->VSurfaceBasedCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;  
}

double IndicesCollector::VSurfaceBasedEL(){  
  double result = 0;
  int index = S->th->surfaceBased->vElIndex;  
  result = Get(S->h,index)- S->th->h0;
  double cape = this->VSurfaceBasedCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;  
} 

double IndicesCollector::VSurfaceBasedLI(){
  int lindex = cache->getPressureIndex(500);  
  double lit = Get(S->t,lindex);
  int vindex = lindex - S->th->surfaceBased->startIndex;  
  double plit = Get(S->th->surfaceBased->virtualValues,vindex);  
  double Showalter = lit - plit;
  return Showalter;
}

double IndicesCollector::VSurfaceBasedVmax(){
  return sqrt(this->VSurfaceBasedCAPE()*2);
}

double IndicesCollector::SBELTemperature(){
  double result = Get(S->t,S->th->surfaceBased->vElIndex);
  double cape = this->VSurfaceBasedCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;

}

double IndicesCollector::SBLCLTemperature(){
  return Get(S->t,S->th->surfaceBased->vLclIndex);
}

double IndicesCollector::SBLFCTemperature(){
  double result = Get(S->t,S->th->surfaceBased->vLfcIndex);
  double cape = this->VSurfaceBasedCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::VMeanLayerCAPE(){
  double result = 0;
  result = S->th->meanLayer->vcape;  
  return result;
}

double IndicesCollector::VLLMeanLayerCAPE(){  
  double result = 0;
  result = S->th->meanLayer->vto3cape;  
  return result;  
}

double IndicesCollector::VMeanLayerCIN(){  
  double result = 0;
  result = S->th->meanLayer->vcin;  
    double cape = this->VMeanLayerCAPE();
  if(cape==0){
    result = sqrt(-1); 
  }
  return result;  
}        

double IndicesCollector::VMeanLayerCIN500(){  
  double result = S->th->meanLayer->vcin500;  
  return result;  
}        

double IndicesCollector::VMeanLayerLCL(){  
  double result = 0;
  int index = S->th->meanLayer->vLclIndex;  
  result = Get(S->h,index)- S->th->h0;  
  return result;
}

double IndicesCollector::VMeanLayerLFC(){  
  double result = 0;
  int index = S->th->meanLayer->vLfcIndex;  
  result = Get(S->h,index)- S->th->h0; 
  double cape = this->VMeanLayerCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::VMeanLayerEL(){  
  double result = 0;
  int index = S->th->meanLayer->vElIndex;  
  result = Get(S->h,index)- S->th->h0;  
  double cape = this->VMeanLayerCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;  
} 

double IndicesCollector::VMeanLayerLI(){
  int lindex = cache->getPressureIndex(500);  
  double lit = Get(S->t,lindex);  
  int vindex = lindex - S->th->meanLayer->startIndex;  
  double plit = Get(S->th->meanLayer->virtualValues,vindex);  
  double Showalter = lit - plit;
  return Showalter;
}

double IndicesCollector::VSurfaceBasedLI_M25(){
  int minten = S->th->min25pos;
  int mintencheck = S->th->surfaceBased->startIndex;
  if(minten<mintencheck){
    minten = mintencheck;
  }
  int vindex = minten - S->th->surfaceBased->startIndex;
  double lit = Get(S->t,minten);
  double plit = Get(S->th->surfaceBased->virtualValues,vindex);
  return lit - plit;
}

double IndicesCollector::VMostU500LI_M25(){
  int minten = S->th->min25pos;
  int mintencheck = S->th->mostU500->startIndex;
  if(minten<mintencheck){
    minten = mintencheck;
  }
  int vindex = minten - S->th->mostU500->startIndex;
  double lit = Get(S->t,minten);
  double plit = Get(S->th->mostU500->virtualValues,vindex);
  return lit - plit;
}

double IndicesCollector::VMeanLayerLI_M25(){
  int minten = S->th->min25pos;
  int mintencheck = S->th->meanLayer->startIndex;
  if(minten<mintencheck){
    minten = mintencheck;
  }
  int vindex = minten - S->th->meanLayer->startIndex;
  double lit = Get(S->t,minten);
  double plit = Get(S->th->meanLayer->virtualValues,vindex);
  return lit - plit;
}

double IndicesCollector::VMostUnstableLI_M25(){
  int minten = S->th->min25pos;
  int mintencheck = S->th->mostUnstable->startIndex;
  if(minten<mintencheck){
    minten = mintencheck;
  }
  int vindex = minten - S->th->mostUnstable->startIndex;
  double lit = Get(S->t,minten);
  double plit = Get(S->th->mostUnstable->virtualValues,vindex);
  return lit - plit;
}

double IndicesCollector::VMeanMostUnstableLI_M25(){
  int minten = S->th->min25pos;
  int mintencheck = S->th->meanmostUnstable->startIndex;
  if(minten<mintencheck){
    minten = mintencheck;
  }
  int vindex = minten - S->th->meanmostUnstable->startIndex;
  double lit = Get(S->t,minten);
  double plit = Get(S->th->meanmostUnstable->virtualValues,vindex);
  return lit - plit;
}

double IndicesCollector::VMeanLayerVmax(){
  return sqrt(this->VMeanLayerCAPE()*2);  
}

double IndicesCollector::MLELTemperature(){
  double result = Get(S->t,S->th->meanLayer->vElIndex);
  double cape = this->VMeanLayerCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::MLLCLTemperature(){
  return Get(S->t,S->th->meanLayer->vLclIndex);
}

double IndicesCollector::MLLFCTemperature(){
  double result = Get(S->t,S->th->meanLayer->vLfcIndex);
  double cape = this->VMeanLayerCAPE();
  if(cape==0)result=sqrt(-1); 
  return result;
}

double IndicesCollector::VDCAPE(){
  return S->th->downdraft->dvcape;
}

double IndicesCollector::LapseRate01(){
  return S->th->lr01;
}

double IndicesCollector::LapseRate24(){
  return S->th->lr24;
}

double IndicesCollector::lapserate03(){
  int lower = 0;
  int upper = cache->getHeightIndex(3000);  
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);  
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::LR36(){
  int lower = cache->getHeightIndex(3000);
  int upper = cache->getHeightIndex(6000);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::LR16(){
  int lower = cache->getHeightIndex(1000);
  int upper = cache->getHeightIndex(6000);  
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::LR26(){
  int lower = cache->getHeightIndex(2000);
  int upper = cache->getHeightIndex(6000);  
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::max_LR26_2km(){
  int lower = cache->getHeightIndex(2000);
  int upper = cache->getHeightIndex(4000);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  double LR12 = 1000*((tup-tlow)/(hup-hlow)); 
  int lower2 = cache->getHeightIndex(3000);
  int upper2 = cache->getHeightIndex(5000);
  double hlow2 = Get(S->h,lower2);
  double hup2 = Get(S->h,upper2);
  double tlow2 = Get(S->t,lower2);
  double tup2 = Get(S->t,upper2);
  double LR23 = 1000*((tup2-tlow2)/(hup2-hlow2)); 
  int lower3 = cache->getHeightIndex(4000);
  int upper3 = cache->getHeightIndex(6000);
  double hlow3 = Get(S->h,lower3);
  double hup3 = Get(S->h,upper3);
  double tlow3 = Get(S->t,lower3);
  double tup3 = Get(S->t,upper3);
  double LR34 = 1000*((tup3-tlow3)/(hup3-hlow3)); 
  return min(min(LR12,LR23),LR34);
}

double IndicesCollector::LR0500(){
  int lower = 0;
  int upper = cache->getHeightIndex(500);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::LR02(){
  int lower = 0;
  int upper = cache->getHeightIndex(2000);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::LR04(){
  int lower = 0;
  int upper = cache->getHeightIndex(4000);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::LR06(){
  int lower = 0;
  int upper = cache->getHeightIndex(6000);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::lapseRate500700(){
  int lower = cache->getPressureIndex(700);
  int upper = cache->getPressureIndex(500);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::lapseRate500800(){
  int lower = cache->getPressureIndex(800);
  int upper = cache->getPressureIndex(500);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::lapseRate600800(){
  int lower = cache->getPressureIndex(800);
  int upper = cache->getPressureIndex(600);
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->t,lower);
  double tup = Get(S->t,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::M10Height(){
  int zeroIndex = S->th->mintenpos;
  double h0 = Get(S->h,0);
  double hz = Get(S->h,zeroIndex);
  return hz-h0;
}

double IndicesCollector::ZeroHeight(){
  int zeroIndex = S->th->zeropos;
  double h0 = Get(S->h,0);
  double hz = Get(S->h,zeroIndex);
  return hz-h0;
}

double IndicesCollector::WetBulbZeroHeight(){
  int zeroIndex = S->th->wbzeropos;
  double h0 = Get(S->h,0);
  double hz = Get(S->h,zeroIndex);
  return hz-h0;
}

double IndicesCollector::MUHeight(){
  int zeroIndex = S->th->mostUnstable->startIndex;
  double h0 = Get(S->h,0);
  double hz = Get(S->h,zeroIndex);
  return hz-h0;
}

double IndicesCollector::MUMLHeight(){
  int zeroIndex = S->th->meanmostUnstable->startIndex;
  double h0 = Get(S->h,0);
  double hz = Get(S->h,zeroIndex);
  return hz-h0;
}

double IndicesCollector::MinTHTEHeight(){
  int zeroIndex = S->th->minTHTEpos;
  double h0 = Get(S->h,0);
  double hz = Get(S->h,zeroIndex);
  return hz-h0;
}

double IndicesCollector::THTE_LR03(){
  int lower = 0;
  int upper = cache->getHeightIndex(3000);  
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->th->oe,lower);
  double tup = Get(S->th->oe,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::THTE_LR04(){
  int lower = 0;
  int upper = cache->getHeightIndex(4000);  
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->th->oe,lower);
  double tup = Get(S->th->oe,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::THTE_LR13(){
  int lower = cache->getHeightIndex(1000);
  int upper = cache->getHeightIndex(3000);  
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->th->oe,lower);
  double tup = Get(S->th->oe,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::THTE_LR14(){
  int lower = cache->getHeightIndex(1000);
  int upper = cache->getHeightIndex(4000);  
  double hlow = Get(S->h,lower);
  double hup = Get(S->h,upper);
  double tlow = Get(S->th->oe,lower);
  double tup = Get(S->th->oe,upper);
  return 1000*((tup-tlow)/(hup-hlow));
}

double IndicesCollector::THTE_LR_LCL_to_M10(){
  int lower = S->th->meanLayer->vLclIndex;
  int upper = S->th->mintenpos;  
  double tlow = Get(S->th->oe,lower);
  double tup = Get(S->th->oe,upper);
  return tup-tlow;
}

double IndicesCollector::THTE_LR_MU_to_M10(){
  int lower = S->th->mostUnstable->startIndex;
  int upper = S->th->mintenpos;  
  double tlow = Get(S->th->oe,lower);
  double tup = Get(S->th->oe,upper);
  return tup-tlow;
}

double IndicesCollector::DeltaThetaE(){
  return Get(S->th->oe,0)-S->th->thetd;
}

double IndicesCollector::thetae01(){
  double result = 0;
  result = S->th->thet01d + 273.15;  
  return result;
}

double IndicesCollector::thetae02(){
  double result = 0;
  result = S->th->thet02d + 273.15;  
  return result;
}

double IndicesCollector::DeltaThetaE_min04km(){
  return Get(S->th->oe,0)-S->th->minTHTE;
}

double IndicesCollector::VirtualColdPoolStrength(){
  return tv(S->th->t0,Get(S->th->mixing,0))-Get(S->th->downdraft->virtualValues,0);
}

double IndicesCollector::WindIndex(){
  double _1000mr = S->th->mr1000;
  double rq = _1000mr / 12; if (rq > 1.0) rq = 1.0;
  int wbzi = S->th->wbzeropos;
  double wbzhd = this->WetBulbZeroHeight() / 1000;
  double gs = (Get(S->t,wbzi) - Get(S->t,0) ) / (wbzhd);
  double sq = wbzhd * rq * (gs * gs - 30 + _1000mr - 2 * Get(S->th->mixing,wbzi)); if (sq < 0) sq = 0;
  double wi = 5 *sqrt(sq);
  double val = wi * 0.514444;
  return val;
}

double IndicesCollector::PWATER(){
  return S->th->pwater;
}

double IndicesCollector::PWATER_eff(){
  return S->th->pwater_eff;
}

double IndicesCollector::MLMixingRatio(){
  return S->th->mmr;
}

double IndicesCollector::MUMRatio(){
  return Get(S->th->mixing,S->th->mostUnstable->startIndex);
}

double IndicesCollector::MUMLMRatio(){
  return S->th->mmrMAX;
}

double IndicesCollector::SBMRatio(){
  return Get(S->th->mixing,S->th->surfaceBased->startIndex);
}

double IndicesCollector::WS_ULmax(){
  int h1 = cache->getHeightIndex(5000); 
  int h2 = cache->getHeightIndex(5500); 
  int h3 = cache->getHeightIndex(6000); 
  int h4 = cache->getHeightIndex(6500); 
  int h5 = cache->getHeightIndex(7000); 
  int h6 = cache->getHeightIndex(7500); 
  int h7 = cache->getHeightIndex(8000); 
  int h8 = cache->getHeightIndex(8500); 
  int h9 = cache->getHeightIndex(9000); 
  double BS1 = (Get(S->ks->vw,h1)).abs();
  double BS2 = (Get(S->ks->vw,h2)).abs();
  double BS3 = (Get(S->ks->vw,h3)).abs();
  double BS4 = (Get(S->ks->vw,h4)).abs();
  double BS5 = (Get(S->ks->vw,h5)).abs();
  double BS6 = (Get(S->ks->vw,h6)).abs();
  double BS7 = (Get(S->ks->vw,h7)).abs();
  double BS8 = (Get(S->ks->vw,h8)).abs();
  double BS9 = (Get(S->ks->vw,h9)).abs();
  return max(max(max(max(max(max(max(max(BS1,BS2),BS3),BS4),BS5),BS6),BS7),BS8),BS9);
}

double IndicesCollector::WS_LLmax(){
  int h1 = cache->getHeightIndex(200); 
  int h2 = cache->getHeightIndex(400); 
  int h3 = cache->getHeightIndex(600); 
  int h4 = cache->getHeightIndex(800); 
  int h5 = cache->getHeightIndex(1000); 
  int h6 = cache->getHeightIndex(1200); 
  int h7 = cache->getHeightIndex(1400); 
  int h8 = cache->getHeightIndex(1600); 
  int h9 = cache->getHeightIndex(1800); 
  int h10 = cache->getHeightIndex(2000); 
  double BS1 = (Get(S->ks->vw,h1)).abs();
  double BS2 = (Get(S->ks->vw,h2)).abs();
  double BS3 = (Get(S->ks->vw,h3)).abs();
  double BS4 = (Get(S->ks->vw,h4)).abs();
  double BS5 = (Get(S->ks->vw,h5)).abs();
  double BS6 = (Get(S->ks->vw,h6)).abs();
  double BS7 = (Get(S->ks->vw,h7)).abs();
  double BS8 = (Get(S->ks->vw,h8)).abs();
  double BS9 = (Get(S->ks->vw,h9)).abs();
  double BS10 = (Get(S->ks->vw,h10)).abs();
  return max(max(max(max(max(max(max(max(max(BS1,BS2),BS3),BS4),BS5),BS6),BS7),BS8),BS9),BS10);
}

double IndicesCollector::WS_MLmax(){
  int h1 = cache->getHeightIndex(2000); 
  int h2 = cache->getHeightIndex(2500); 
  int h3 = cache->getHeightIndex(3000); 
  int h4 = cache->getHeightIndex(3500); 
  int h5 = cache->getHeightIndex(4000); 
  int h6 = cache->getHeightIndex(4500); 
  int h7 = cache->getHeightIndex(5000); 
  double BS1 = (Get(S->ks->vw,h1)).abs();
  double BS2 = (Get(S->ks->vw,h2)).abs();
  double BS3 = (Get(S->ks->vw,h3)).abs();
  double BS4 = (Get(S->ks->vw,h4)).abs();
  double BS5 = (Get(S->ks->vw,h5)).abs();
  double BS6 = (Get(S->ks->vw,h6)).abs();
  double BS7 = (Get(S->ks->vw,h7)).abs();
  return max(max(max(max(max(max(BS1,BS2),BS3),BS4),BS5),BS6),BS7);
}

double IndicesCollector::BS_LLmax(){
  int h1 = cache->getHeightIndex(200); 
  int h2 = cache->getHeightIndex(400); 
  int h3 = cache->getHeightIndex(600); 
  int h4 = cache->getHeightIndex(800); 
  int h5 = cache->getHeightIndex(1000); 
  int h6 = cache->getHeightIndex(1200); 
  int h7 = cache->getHeightIndex(1400); 
  int h8 = cache->getHeightIndex(1600); 
  int h9 = cache->getHeightIndex(1800); 
  int h10 = cache->getHeightIndex(2000); 
  double BS1 = (Get(S->ks->vw,h1) - Get(S->ks->vw,0)).abs();
  double BS2 = (Get(S->ks->vw,h2) - Get(S->ks->vw,0)).abs();
  double BS3 = (Get(S->ks->vw,h3) - Get(S->ks->vw,0)).abs();
  double BS4 = (Get(S->ks->vw,h4) - Get(S->ks->vw,0)).abs();
  double BS5 = (Get(S->ks->vw,h5) - Get(S->ks->vw,0)).abs();
  double BS6 = (Get(S->ks->vw,h6) - Get(S->ks->vw,0)).abs();
  double BS7 = (Get(S->ks->vw,h7) - Get(S->ks->vw,0)).abs();
  double BS8 = (Get(S->ks->vw,h8) - Get(S->ks->vw,0)).abs();
  double BS9 = (Get(S->ks->vw,h9) - Get(S->ks->vw,0)).abs();
  double BS10 = (Get(S->ks->vw,h10) - Get(S->ks->vw,0)).abs();
  return max(max(max(max(max(max(max(max(max(BS1,BS2),BS3),BS4),BS5),BS6),BS7),BS8),BS9),BS10);
}

double IndicesCollector::BS_ULmax(){
  int h1 = cache->getHeightIndex(5000); 
  int h2 = cache->getHeightIndex(5500); 
  int h3 = cache->getHeightIndex(6000); 
  int h4 = cache->getHeightIndex(6500); 
  int h5 = cache->getHeightIndex(7000); 
  int h6 = cache->getHeightIndex(7500); 
  int h7 = cache->getHeightIndex(8000); 
  int h8 = cache->getHeightIndex(8500); 
  int h9 = cache->getHeightIndex(9000); 
  double BS1 = (Get(S->ks->vw,h1) - Get(S->ks->vw,0)).abs();
  double BS2 = (Get(S->ks->vw,h2) - Get(S->ks->vw,0)).abs();
  double BS3 = (Get(S->ks->vw,h3) - Get(S->ks->vw,0)).abs();
  double BS4 = (Get(S->ks->vw,h4) - Get(S->ks->vw,0)).abs();
  double BS5 = (Get(S->ks->vw,h5) - Get(S->ks->vw,0)).abs();
  double BS6 = (Get(S->ks->vw,h6) - Get(S->ks->vw,0)).abs();
  double BS7 = (Get(S->ks->vw,h7) - Get(S->ks->vw,0)).abs();
  double BS8 = (Get(S->ks->vw,h8) - Get(S->ks->vw,0)).abs();
  double BS9 = (Get(S->ks->vw,h9) - Get(S->ks->vw,0)).abs();
  return max(max(max(max(max(max(max(max(BS1,BS2),BS3),BS4),BS5),BS6),BS7),BS8),BS9);
}

double IndicesCollector::BS_MLmax(){
  int h1 = cache->getHeightIndex(2000); 
  int h2 = cache->getHeightIndex(2500); 
  int h3 = cache->getHeightIndex(3000); 
  int h4 = cache->getHeightIndex(3500); 
  int h5 = cache->getHeightIndex(4000); 
  int h6 = cache->getHeightIndex(4500); 
  int h7 = cache->getHeightIndex(5000); 
  double BS1 = (Get(S->ks->vw,h1) - Get(S->ks->vw,0)).abs();
  double BS2 = (Get(S->ks->vw,h2) - Get(S->ks->vw,0)).abs();
  double BS3 = (Get(S->ks->vw,h3) - Get(S->ks->vw,0)).abs();
  double BS4 = (Get(S->ks->vw,h4) - Get(S->ks->vw,0)).abs();
  double BS5 = (Get(S->ks->vw,h5) - Get(S->ks->vw,0)).abs();
  double BS6 = (Get(S->ks->vw,h6) - Get(S->ks->vw,0)).abs();
  double BS7 = (Get(S->ks->vw,h7) - Get(S->ks->vw,0)).abs();
  return max(max(max(max(max(max(BS1,BS2),BS3),BS4),BS5),BS6),BS7);
}

double IndicesCollector::BS06_var_SD(){
  int H500 = cache->getHeightIndex(500); 
  int H1000 = cache->getHeightIndex(1000); 
  int H1500 = cache->getHeightIndex(1500); 
  int H2000 = cache->getHeightIndex(2000); 
  int H2500 = cache->getHeightIndex(2500); 
  int H3000 = cache->getHeightIndex(3000); 
  int H3500 = cache->getHeightIndex(3500); 
  int H4000 = cache->getHeightIndex(4000); 
  int H4500 = cache->getHeightIndex(4500); 
  int H5000 = cache->getHeightIndex(5000); 
  int H5500 = cache->getHeightIndex(5500); 
  int H6000 = cache->getHeightIndex(6000); 
  Vector V01 = Get(S->ks->vw,H500) - Get(S->ks->vw,0);
  Vector V02 = Get(S->ks->vw,H1000) - Get(S->ks->vw,H500);
  Vector V03 = Get(S->ks->vw,H1500) - Get(S->ks->vw,H1000);
  Vector V04 = Get(S->ks->vw,H2000) - Get(S->ks->vw,H1500);
  Vector V05 = Get(S->ks->vw,H2500) - Get(S->ks->vw,H2000);
  Vector V06 = Get(S->ks->vw,H3000) - Get(S->ks->vw,H2500);
  Vector V07 = Get(S->ks->vw,H3500) - Get(S->ks->vw,H3000);
  Vector V08 = Get(S->ks->vw,H4000) - Get(S->ks->vw,H3500);
  Vector V09 = Get(S->ks->vw,H4500) - Get(S->ks->vw,H4000);
  Vector V10 = Get(S->ks->vw,H5000) - Get(S->ks->vw,H4500);
  Vector V11 = Get(S->ks->vw,H5500) - Get(S->ks->vw,H5000);
  Vector V12 = Get(S->ks->vw,H6000) - Get(S->ks->vw,H5500);
  vector<double> diffs = { V01.abs() - V02.abs(),
                        V02.abs() - V03.abs(),
                        V03.abs() - V04.abs(),
                        V04.abs() - V05.abs(),
                        V05.abs() - V06.abs(),
                        V06.abs() - V07.abs(),
                        V07.abs() - V08.abs(),
                        V08.abs() - V09.abs(),
                        V09.abs() - V10.abs(),
                        V10.abs() - V11.abs(),
                        V11.abs() - V12.abs() };
  double result = computeStandardDeviation(diffs);
  return 1.0/result;
}

double IndicesCollector::BS06_var_SI(){
  int H500 = cache->getHeightIndex(500); 
  int H1000 = cache->getHeightIndex(1000); 
  int H1500 = cache->getHeightIndex(1500); 
  int H2000 = cache->getHeightIndex(2000); 
  int H2500 = cache->getHeightIndex(2500); 
  int H3000 = cache->getHeightIndex(3000); 
  int H3500 = cache->getHeightIndex(3500); 
  int H4000 = cache->getHeightIndex(4000); 
  int H4500 = cache->getHeightIndex(4500); 
  int H5000 = cache->getHeightIndex(5000); 
  int H5500 = cache->getHeightIndex(5500); 
  int H6000 = cache->getHeightIndex(6000); 
  Vector V01 = Get(S->ks->vw,H500) - Get(S->ks->vw,0);
  Vector V02 = Get(S->ks->vw,H1000) - Get(S->ks->vw,H500);
  Vector V03 = Get(S->ks->vw,H1500) - Get(S->ks->vw,H1000);
  Vector V04 = Get(S->ks->vw,H2000) - Get(S->ks->vw,H1500);
  Vector V05 = Get(S->ks->vw,H2500) - Get(S->ks->vw,H2000);
  Vector V06 = Get(S->ks->vw,H3000) - Get(S->ks->vw,H2500);
  Vector V07 = Get(S->ks->vw,H3500) - Get(S->ks->vw,H3000);
  Vector V08 = Get(S->ks->vw,H4000) - Get(S->ks->vw,H3500);
  Vector V09 = Get(S->ks->vw,H4500) - Get(S->ks->vw,H4000);
  Vector V10 = Get(S->ks->vw,H5000) - Get(S->ks->vw,H4500);
  Vector V11 = Get(S->ks->vw,H5500) - Get(S->ks->vw,H5000);
  Vector V12 = Get(S->ks->vw,H6000) - Get(S->ks->vw,H5500);  
  vector<double> diffs = { V01.abs() - V02.abs(),
                        V02.abs() - V03.abs(),
                        V03.abs() - V04.abs(),
                        V04.abs() - V05.abs(),
                        V05.abs() - V06.abs(),
                        V06.abs() - V07.abs(),
                        V07.abs() - V08.abs(),
                        V08.abs() - V09.abs(),
                        V09.abs() - V10.abs(),
                        V10.abs() - V11.abs(),
                        V11.abs() - V12.abs() };  
  double result = 1.0 / (computeStandardDeviation(diffs) / abs(computeMean(diffs)));
  return result;
}

double IndicesCollector::BS16_var_SD(){
  int H1000 = cache->getHeightIndex(1000); 
  int H1500 = cache->getHeightIndex(1500); 
  int H2000 = cache->getHeightIndex(2000); 
  int H2500 = cache->getHeightIndex(2500); 
  int H3000 = cache->getHeightIndex(3000); 
  int H3500 = cache->getHeightIndex(3500); 
  int H4000 = cache->getHeightIndex(4000); 
  int H4500 = cache->getHeightIndex(4500); 
  int H5000 = cache->getHeightIndex(5000); 
  int H5500 = cache->getHeightIndex(5500); 
  int H6000 = cache->getHeightIndex(6000); 
  Vector V03 = Get(S->ks->vw,H1500) - Get(S->ks->vw,H1000);
  Vector V04 = Get(S->ks->vw,H2000) - Get(S->ks->vw,H1500);
  Vector V05 = Get(S->ks->vw,H2500) - Get(S->ks->vw,H2000);
  Vector V06 = Get(S->ks->vw,H3000) - Get(S->ks->vw,H2500);
  Vector V07 = Get(S->ks->vw,H3500) - Get(S->ks->vw,H3000);
  Vector V08 = Get(S->ks->vw,H4000) - Get(S->ks->vw,H3500);
  Vector V09 = Get(S->ks->vw,H4500) - Get(S->ks->vw,H4000);
  Vector V10 = Get(S->ks->vw,H5000) - Get(S->ks->vw,H4500);
  Vector V11 = Get(S->ks->vw,H5500) - Get(S->ks->vw,H5000);
  Vector V12 = Get(S->ks->vw,H6000) - Get(S->ks->vw,H5500);
  vector<double> diffs = { V03.abs() - V04.abs(),
                     V04.abs() - V05.abs(),
                     V05.abs() - V06.abs(),
                     V06.abs() - V07.abs(),
                     V07.abs() - V08.abs(),
                     V08.abs() - V09.abs(),
                     V09.abs() - V10.abs(),
                     V10.abs() - V11.abs(),
                     V11.abs() - V12.abs() };
  double result = computeStandardDeviation(diffs);
  return 1.0/result;
}

double IndicesCollector::BS16_var_SI(){
  int H1000 = cache->getHeightIndex(1000); 
  int H1500 = cache->getHeightIndex(1500); 
  int H2000 = cache->getHeightIndex(2000); 
  int H2500 = cache->getHeightIndex(2500); 
  int H3000 = cache->getHeightIndex(3000); 
  int H3500 = cache->getHeightIndex(3500); 
  int H4000 = cache->getHeightIndex(4000); 
  int H4500 = cache->getHeightIndex(4500); 
  int H5000 = cache->getHeightIndex(5000); 
  int H5500 = cache->getHeightIndex(5500); 
  int H6000 = cache->getHeightIndex(6000); 
  Vector V03 = Get(S->ks->vw,H1500) - Get(S->ks->vw,H1000);
  Vector V04 = Get(S->ks->vw,H2000) - Get(S->ks->vw,H1500);
  Vector V05 = Get(S->ks->vw,H2500) - Get(S->ks->vw,H2000);
  Vector V06 = Get(S->ks->vw,H3000) - Get(S->ks->vw,H2500);
  Vector V07 = Get(S->ks->vw,H3500) - Get(S->ks->vw,H3000);
  Vector V08 = Get(S->ks->vw,H4000) - Get(S->ks->vw,H3500);
  Vector V09 = Get(S->ks->vw,H4500) - Get(S->ks->vw,H4000);
  Vector V10 = Get(S->ks->vw,H5000) - Get(S->ks->vw,H4500);
  Vector V11 = Get(S->ks->vw,H5500) - Get(S->ks->vw,H5000);
  Vector V12 = Get(S->ks->vw,H6000) - Get(S->ks->vw,H5500);
  vector<double> diffs = { V03.abs() - V04.abs(),
                     V04.abs() - V05.abs(),
                     V05.abs() - V06.abs(),
                     V06.abs() - V07.abs(),
                     V07.abs() - V08.abs(),
                     V08.abs() - V09.abs(),
                     V09.abs() - V10.abs(),
                     V10.abs() - V11.abs(),
                     V11.abs() - V12.abs() };  
  double result = 1.0 / (computeStandardDeviation(diffs) / abs(computeMean(diffs)));
  return result;
}

double IndicesCollector::BS500(){
  int tail=0;
  int head = cache->getHeightIndex(500);  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS01(){
  int tail=0;
  int head = cache->getHeightIndex(1000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS02(){
  int tail=0;
  int head = cache->getHeightIndex(2000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS03(){
  int tail=0;
  int head = cache->getHeightIndex(3000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS04(){
  int tail=0;
  int head = cache->getHeightIndex(4000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS010(){
  int tail=0;
  int head = cache->getHeightIndex(10000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS36(){
  int tail=cache->getHeightIndex(3000);
  int head = cache->getHeightIndex(6000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS06(){
  int tail=0;
  int head = cache->getHeightIndex(6000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS08(){
  int tail=0;
  int head = cache->getHeightIndex(8000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS18(){
  int tail=cache->getHeightIndex(1000);
  int head = cache->getHeightIndex(8000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS110(){
  int tail=cache->getHeightIndex(1000);
  int head = cache->getHeightIndex(10000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS14(){
  int tail=cache->getHeightIndex(1000);
  int head = cache->getHeightIndex(4000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS25(){
  int tail=cache->getHeightIndex(2000);
  int head = cache->getHeightIndex(5000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS13(){
  int tail=cache->getHeightIndex(1000);
  int head = cache->getHeightIndex(3000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BS16(){
  int tail=cache->getHeightIndex(1000);
  int head = cache->getHeightIndex(6000);
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::MeanWind500(){
  return S->ks->mean0.abs();
}

double IndicesCollector::MeanWind06(){
  return S->ks->mean06.abs();
}

double IndicesCollector::MeanWind01(){
  return S->ks->mean01.abs();
}

double IndicesCollector::MeanWind02(){
  return S->ks->mean02.abs();
}

double IndicesCollector::MeanWind03(){
  return S->ks->mean03.abs();
}

double IndicesCollector::MeanWind13(){
  return S->ks->mean13.abs();
}

double IndicesCollector::SRH100RM(){
  return S->ks->srh100rm;
}

double IndicesCollector::SRH100RM_F(){
  return S->ks->srh100rm2;
}

double IndicesCollector::SRH100LM_F(){
  return S->ks->srh100lm2;
}

double IndicesCollector::SRH250RM(){
  return S->ks->srh250rm;
}

double IndicesCollector::SRH500RM(){
  return S->ks->srh500rm;
}

double IndicesCollector::SRH500RM_F(){
  return S->ks->srh500rm2;
}

double IndicesCollector::SRH500LM_F(){
  return S->ks->srh500lm2;
}

double IndicesCollector::SRH01RM(){
  return S->ks->srh01rm;
}

double IndicesCollector::SRH03RM(){
  return S->ks->srh03rm;
}

double IndicesCollector::SRH36RM(){
  return S->ks->srh36rm - S->ks->srh03rm;
}

double IndicesCollector::SRH13_RM(){
  return S->ks->srh03rm - S->ks->srh01rm;
}

double IndicesCollector::SRH16RM(){
  return (S->ks->srh36rm + S->ks->srh03rm) - S->ks->srh01rm;
}

double IndicesCollector::SRH100LM(){
  return S->ks->srh100lm;
}

double IndicesCollector::SRH250LM(){
  return S->ks->srh250lm;
}

double IndicesCollector::SRH500LM(){
  return S->ks->srh500lm;
}

double IndicesCollector::SRH01LM(){
  return S->ks->srh01lm;
}

double IndicesCollector::SRH03LM(){
  return S->ks->srh03lm;
}

double IndicesCollector::SRH36LM(){
  return S->ks->srh36lm - S->ks->srh03lm;
}

double IndicesCollector::SRH13_LM(){
  return S->ks->srh03lm - S->ks->srh01lm;
}

double IndicesCollector::SRH16LM(){
  return (S->ks->srh36lm + S->ks->srh03lm) - S->ks->srh01lm;
}

double IndicesCollector::SRH01SM(){
  return S->ks->srh01sm;
}

double IndicesCollector::SRH03SM(){
  return S->ks->srh03sm;
}

double IndicesCollector::SRH36SM(){
  return S->ks->srh36sm - S->ks->srh03sm;
}

double IndicesCollector::SRH13_SM(){
  return S->ks->srh03sm - S->ks->srh01sm;
}

double IndicesCollector::SRH16SM(){
  return (S->ks->srh36sm + S->ks->srh03sm) - S->ks->srh01sm;
}

double IndicesCollector::SRH01LM_eff(){
  return S->ks->srh01lm_eff;
}

double IndicesCollector::SRH01RM_eff(){
  return S->ks->srh01rm_eff;
}

double IndicesCollector::SRH01SM_eff(){
  return S->ks->srh01sm_eff;
}

Vector IndicesCollector::Bunkers4_RM_vector(){
  Vector meanwind = S->ks->mean06;
  Vector tv = Vector(0, 0, 1);
  Vector dev = Vector(0, 0, 0);
  Vector tshear = S->ks->mean6 - S->ks->mean0;
  dev = Vector::vec(tshear,tv);
  dev *= 4;
  dev *= 1.0 / tshear.abs();
  Vector Bunkers4 = meanwind - dev;
  return Bunkers4;
}

double IndicesCollector::Bunkers4_RM_A(){
  Vector Bunkers4_SM = this->Bunkers4_RM_vector();
  double *tab = Bunkers4_SM.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Bunkers4_RM_M(){
  Vector Bunkers4_SM = this->Bunkers4_RM_vector();
  double *tab = Bunkers4_SM.toAV(); 
  double magnitude = tab[1];
  delete[] tab;
  return magnitude;
}

Vector IndicesCollector::Bunkers4_LM_vector(){
  Vector meanwind = S->ks->mean06;
  Vector tv = Vector(0, 0, 1);
  Vector dev = Vector(0, 0, 0);
  Vector tshear = S->ks->mean6 - S->ks->mean0;
  dev = Vector::vec(tshear,tv);
  dev *= -4;
  dev *= 1.0 / tshear.abs();
  Vector Bunkers4 = meanwind - dev;
  return Bunkers4;
}

double IndicesCollector::Bunkers4_LM_A(){
  Vector Bunkers4_SM = this->Bunkers4_LM_vector();
  double *tab = Bunkers4_SM.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Bunkers4_LM_M(){
  Vector Bunkers4_SM = this->Bunkers4_LM_vector();
  double *tab = Bunkers4_SM.toAV(); 
  double magnitude = tab[1];
  delete[] tab;
  return magnitude;
}

Vector IndicesCollector::Peters_vector(){
  double SRH_mean = S->ks->srh03sm;
  double sign_SRH = SRH_mean/abs(SRH_mean);
  double fact = 1;
  double propfac = sign_SRH*min(abs(SRH_mean)/75,fact);
  Vector meanwind = S->ks->mean06;
  Vector tv = Vector(0, 0, 1);
  Vector dev = Vector(0, 0, 0);
  Vector tshear = S->ks->mean6 - S->ks->mean0;
  dev = Vector::vec(tshear,tv);
  dev *= 7.5*propfac;
  dev *= 1.0 / tshear.abs();
  Vector Peters_SM = meanwind - dev;
  return Peters_SM;
}

double IndicesCollector::Peters_SR_inflow(){
  Vector res = S->ks->mean0 - this->Peters_vector();
  return res.abs();
}

double IndicesCollector::Peters_SR_inflow_eff(){
  Vector res = S->ks->mean01eff - this->Peters_vector();
  return res.abs();
}

double IndicesCollector::Peters_vector_A(){
  Vector Peters_SM = this->Peters_vector();
  double *tab = Peters_SM.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Peters_vector_M(){
  Vector Peters_SM = this->Peters_vector();
  double *tab = Peters_SM.toAV(); 
  double magnitude = tab[1];
  delete[] tab;
  return magnitude;
}

double IndicesCollector::Peters_vector_eff_A(){
  Vector Peters_SM = this->Peters_vector();
  double *tab = Peters_SM.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Peters_vector_eff_M(){
  Vector Peters_SM = this->Peters_vector();
  double *tab = Peters_SM.toAV(); 
  double magnitude = tab[1];
  delete[] tab;
  return magnitude;
}

double IndicesCollector::SRH03LM_eff(){
  return S->ks->srh03lm_eff;
}

double IndicesCollector::SRH03RM_eff(){
  return S->ks->srh03rm_eff;
}

double IndicesCollector::SRH03SM_eff(){
  return S->ks->srh03sm_eff;
}

double IndicesCollector::emubs(){
  Vector gVector = Get(S->ks->vw,S->th->mostUnstable->startIndex);
  int index=S->th->mostUnstable->startIndex+((S->th->mostUnstable->vElIndex-S->th->mostUnstable->startIndex)/2);
  double h0 = Get(S->h,S->th->mostUnstable->startIndex);
  double hn = Get(S->h,S->th->mostUnstable->vElIndex);  
  double middle = h0+((hn-h0)/2.0);
  double hindex = Get(S->h,index);
  int destindex =-1;  
  if(middle ==hindex )destindex = index;
  else if (middle > hindex){
    for(size_t i = index;i<S->h->size()-1;i++){
      double upper = Get(S->h,i+1);
      double lower = Get(S->h, i);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;
        break;
      }
    }
  }else{
    for(int i = index;i>1;i--){
      double upper = Get(S->h,i);
      double lower = Get(S->h, i-1);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;        
        break;
      }
    }
  }
  Vector middleVector =  Get(S->ks->vw,destindex);
  double effSHR = (middleVector - gVector).abs();
  double mucape = this->VMostUnstableCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::emumlbs(){
  Vector gVector = Get(S->ks->vw,S->th->meanmostUnstable->startIndex);
  int index=S->th->meanmostUnstable->startIndex+((S->th->meanmostUnstable->vElIndex-S->th->meanmostUnstable->startIndex)/2);
  double h0 = Get(S->h,S->th->meanmostUnstable->startIndex);
  double hn = Get(S->h,S->th->meanmostUnstable->vElIndex);  
  double middle = h0+((hn-h0)/2.0);
  double hindex = Get(S->h,index);
  int destindex =-1;  
  if(middle ==hindex )destindex = index;
  else if (middle > hindex){
    for(size_t i = index;i<S->h->size()-1;i++){
      double upper = Get(S->h,i+1);
      double lower = Get(S->h, i);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;
        break;
      }
    }
  }else{
    for(int i = index;i>1;i--){
      double upper = Get(S->h,i);
      double lower = Get(S->h, i-1);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;        
        break;
      }
    }
  }
  Vector middleVector =  Get(S->ks->vw,destindex);
  double effSHR = (middleVector - gVector).abs();
  double mucape = this->VMeanMostUnstableCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::esbbs(){
  Vector gVector = Get(S->ks->vw,S->th->surfaceBased->startIndex);
  int index=S->th->surfaceBased->startIndex+((S->th->surfaceBased->vElIndex-S->th->surfaceBased->startIndex)/2);
  double h0 = Get(S->h,S->th->surfaceBased->startIndex);
  double hn = Get(S->h,S->th->surfaceBased->vElIndex);
  double middle = h0+((hn-h0)/2.0);
  double hindex = Get(S->h,index);  
  int destindex =-1;
  if(middle ==hindex )destindex = index;
  else if (middle > hindex){
    for(size_t i = index;i<S->h->size()-1;i++){
      double upper = Get(S->h,i+1);
      double lower = Get(S->h, i);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;        
        break;
      }
    }
  }else{
    for(int i = index;i>1;i--){
      double upper = Get(S->h,i);
      double lower = Get(S->h, i-1);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;
        break;
      }
    }
  }  
  Vector middleVector =  Get(S->ks->vw,destindex);
  double effSHR = (middleVector - gVector).abs();
  double mucape = this->VSurfaceBasedCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::emu500bs(){
  Vector gVector = Get(S->ks->vw,S->th->mostU500->startIndex);
  int index=S->th->mostU500->startIndex+((S->th->mostU500->vElIndex-S->th->mostU500->startIndex)/2);
  double h0 = Get(S->h,S->th->mostU500->startIndex);
  double hn = Get(S->h,S->th->mostU500->vElIndex);  
  double middle = h0+((hn-h0)/2.0);
  double hindex = Get(S->h,index);
  int destindex =-1;  
  if(middle ==hindex )destindex = index;
  else if (middle > hindex){
    for(size_t i = index;i<S->h->size()-1;i++){
      double upper = Get(S->h,i+1);
      double lower = Get(S->h, i);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;
        break;
      }
    }
  }else{
    for(int i = index;i>1;i--){
      double upper = Get(S->h,i);
      double lower = Get(S->h, i-1);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;        
        break;
      }
    }
  }
  Vector middleVector =  Get(S->ks->vw,destindex);  
  double effSHR = (middleVector - gVector).abs();
  double mucape = this->MU500_CAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::emlbs(){
  Vector gVector = Get(S->ks->vw,S->th->meanLayer->startIndex);
  int index=S->th->meanLayer->startIndex+((S->th->meanLayer->vElIndex-S->th->meanLayer->startIndex)/2);
  double h0 = Get(S->h,S->th->meanLayer->startIndex);
  double hn = Get(S->h,S->th->meanLayer->vElIndex);  
  double middle = h0+((hn-h0)/2.0);
  double hindex = Get(S->h,index);
  int destindex =-1;  
  if(middle ==hindex )destindex = index;
  else if (middle > hindex){
    for(size_t i = index;i<S->h->size()-1;i++){
      double upper = Get(S->h,i+1);
      double lower = Get(S->h, i);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;
        break;
      }
    }
  }else{
    for(int i = index;i>1;i--){
      double upper = Get(S->h,i);
      double lower = Get(S->h, i-1);
      if(middle>=lower && middle<=upper){
        if(abs(middle-lower)>abs(upper-middle)){
          destindex = i+1;
        }else destindex=i;        
        break;
      }
    }
  }
  Vector middleVector =  Get(S->ks->vw,destindex);  
  double effSHR = (middleVector - gVector).abs();
  double mucape = this->VMeanLayerCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::MUmiddlecape(){
  return S->th->mostUnstable->middlecape;
}

double IndicesCollector::MUMLmiddlecape(){
  return S->th->meanmostUnstable->middlecape;
}

double IndicesCollector::SBmiddlecape(){
  return S->th->surfaceBased->middlecape;
}

double IndicesCollector::MLmiddlecape(){
  return S->th->meanLayer->middlecape;
}

double IndicesCollector::MU500middlecape(){
  return S->th->mostU500->middlecape;
}

double IndicesCollector::VLLMU500CAPE(){  
  double result = 0;
  result = S->th->mostU500->vto3cape;  
  return result;  
}

double IndicesCollector::LTTP_RM(){

  double CAPEml = this->VMeanLayerCAPE();
  double LCLml = this->VMeanLayerLCL();	
  double CINml = this->VMeanLayerCIN();		
  if(isnan(CINml)) CINml = 0;
	
  double tcape = max(min((CAPEml / (321.59)), (1501.74 / 321.59)), 0.0);
  double tcin = min(max((200.0 - CINml) / (150.0), 0.0), 1.0);
  double tlcl = min(max(((2100.0 - LCLml) / (700.0)), 0.0 ), 1.0);
  if(LCLml > 2100)tlcl=0;
  if(LCLml < 700)tlcl=1;
	
  double LTTPt = tcape * tlcl * tcin;

  double SRH500 = this->SRH500RM();
  double SRWSSFC = this->SRW_sfc_RM();
  double SRMLIKE = this->MeanSR500_RM() * this->MeanSR500_RM();  	
  double BWD500 = this->BS500();
  double BWD4000 = this->BS04(); 	
  double BWDULMAX = this->BS_ULmax();
  double WSLLMAX = this->WS_LLmax();
  double WSULMAX = this->WS_ULmax();
  double WS4000 = Get(S->ks->vw,cache->getHeightIndex(4000)).abs(); 	
	
  double LTTPk = max(min((SRH500 / (34.21)), (272.85 / 34.21)), 0.0) *
	         max(min((SRWSSFC / (10.96)), (24.87 / 10.96)), 0.0) *
	         max(min((SRMLIKE / (134.92)), (411.86/ 134.923)) ,0.0) * 
	         max(min((BWD500 / (4.60)), (15.48 / 4.60)), 0.0) *
	         max(min((BWD4000 / (14.22)), (30.61 / 14.22)), 0.0) *
	         max(min((BWDULMAX / (20.93)), (41.64 / 20.93)), 0.0) *
	         max(min((WSLLMAX / (11.96)),(29.81 / 11.96)), 0.0) * 
	         max(min((WS4000 / (14.97)), (34.59 / 14.97)), 0.0) * 
	         max(min((WSULMAX /(21.65)), (46.47 / 21.65)), 0.0);
	
  return sqrt(sqrt(LTTPt * LTTPk));
}

double IndicesCollector::LTTP_LM(){

  double CAPEml = this->VMeanLayerCAPE();
  double LCLml = this->VMeanLayerLCL();	
  double CINml = this->VMeanLayerCIN();	
  if(isnan(CINml)) CINml = 0;
	
  double tcape = max(min((CAPEml / (321.59)), (1501.74 / 321.59)), 0.0);
  double tcin = min(max((200.0 - CINml) / (150.0), 0.0), 1.0);
  double tlcl = min(max(((2100.0 - LCLml) / (700.0)), 0.0 ), 1.0);
  if(LCLml > 2100)tlcl=0;
  if(LCLml < 700)tlcl=1;
	
  double LTTPt = tcape * tlcl * tcin;

  double SRH500 = (this->SRH500LM())*(-1.0);
  double SRWSSFC = this->SRW_sfc_LM();
  double SRMLIKE = this->MeanSR500_LM() * this->MeanSR500_LM();  	
  double BWD500 = this->BS500();
  double BWD4000 = this->BS04(); 	
  double BWDULMAX = this->BS_ULmax();
  double WSLLMAX = this->WS_LLmax();
  double WSULMAX = this->WS_ULmax();
  double WS4000 = Get(S->ks->vw,cache->getHeightIndex(4000)).abs(); 	
	
  double LTTPk = max(min((SRH500 / (34.21)), (272.85 / 34.21)), 0.0) *
	         max(min((SRWSSFC / (10.96)), (24.87 / 10.96)), 0.0) *
	         max(min((SRMLIKE / (134.92)), (411.86/ 134.923)) ,0.0) * 
	         max(min((BWD500 / (4.60)), (15.48 / 4.60)), 0.0) *
	         max(min((BWD4000 / (14.22)), (30.61 / 14.22)), 0.0) *
	         max(min((BWDULMAX / (20.93)), (41.64 / 20.93)), 0.0) *
	         max(min((WSLLMAX / (11.96)),(29.81 / 11.96)), 0.0) * 
	         max(min((WS4000 / (14.97)), (34.59 / 14.97)), 0.0) * 
	         max(min((WSULMAX /(21.65)), (46.47 / 21.65)), 0.0);
	
  return sqrt(sqrt(LTTPt * LTTPk));
}

double IndicesCollector::STP(){
  double sbcape = this->VSurfaceBasedCAPE()/1500;
  double sblcl = this->VSurfaceBasedLCL();
  double srh1 = this->SRH01RM()/150;
  double bwd = this->BS06();
  double cin = this->VSurfaceBasedCIN();	
  if(isnan(cin)) cin = 0;
  if(sblcl<1000)sblcl=1;
  else if(sblcl>2000)sblcl=0;
  else sblcl=(2000-sblcl)/1000;  
  if(cin>-50)cin=1;
  else if(cin<-200)cin=0;
  else cin=(200+cin)/150;
  if(bwd<12.5)bwd=0.0;
  else if(bwd>30)bwd = 1.5;
  else bwd/=20;  
  return sbcape*sblcl*srh1*bwd*cin;
}

double IndicesCollector::STPeff(){
  double sbcape = this->VMeanLayerCAPE()/1500;
  double sblcl = this->VMeanLayerLCL();
  double srh1 = this->SRH500RM()/75;
  double bwd = this->emlbs();
  double cin = this->VMeanLayerCIN();	
  if(isnan(cin)) cin = 0;
  if(sblcl<1000)sblcl=1;
  else if(sblcl>2000)sblcl=0;
  else sblcl=(2000-sblcl)/1000;  
  if(cin>-50)cin=1;
  else if(cin<-200)cin=0;
  else cin=(200+cin)/150;
  if(bwd<12.5)bwd=0.0;
  else if(bwd>30)bwd = 1.5;
  else bwd/=20;  
  return sbcape*sblcl*srh1*bwd*cin;
}

double IndicesCollector::STP_LM(){
  double sbcape = this->VSurfaceBasedCAPE()/1500;
  double sblcl = this->VSurfaceBasedLCL();
  double srh1 = this->SRH01LM()/150;
  double bwd = this->BS06();
  double cin = this->VSurfaceBasedCIN();	
  if(isnan(cin)) cin = 0;
  if(sblcl<1000)sblcl=1;
  else if(sblcl>2000)sblcl=0;
  else sblcl=(2000-sblcl)/1000;  
  if(cin>-50)cin=1;
  else if(cin<-200)cin=0;
  else cin=(200+cin)/150;
  if(bwd<12.5)bwd=0.0;
  else if(bwd>30)bwd = 1.5;
  else bwd/=20;
  return sbcape*sblcl*srh1*bwd*cin;
}

double IndicesCollector::STPeff_LM(){
  double sbcape = this->VMeanLayerCAPE()/1500;
  double sblcl = this->VMeanLayerLCL();
  double srh1 = this->SRH500LM()/75;
  double bwd = this->emlbs();
  double cin = this->VMeanLayerCIN();	
  if(isnan(cin)) cin = 0;
  if(sblcl<1000)sblcl=1;
  else if(sblcl>2000)sblcl=0;
  else sblcl=(2000-sblcl)/1000;  
  if(cin>-50)cin=1;
  else if(cin<-200)cin=0;
  else cin=(200+cin)/150;
  if(bwd<12.5)bwd=0.0;
  else if(bwd>30)bwd = 1.5;
  else bwd/=20;  
  return sbcape*sblcl*srh1*bwd*cin;
}

double IndicesCollector::SHERBS3(){
  double T1 = this->BS03()/26;
  double T2 = this->lapserate03()/-(5.2);
  double T3 = this->lapseRate500700()/-(5.6);
  return T1*T2*T3;
}

double IndicesCollector::SHERBE(){
  double T1 = this->emubs()/27;
  double T2 = this->lapserate03()/-(5.2);
  double T3 = this->lapseRate500700()/-(5.6);
  return T1*T2*T3;
}

double IndicesCollector::SHERBS3_v2(){
  double T1 = this->MeanWind03()/20;
  double T2 = this->lapserate03()/-(5.2);
  double T3 = this->LR26()/-(5.6);
  return T1*T2*T3;
}

double IndicesCollector::SHERBE_v2(){
  double T1 = this->emubs()/27;
  double T2 = this->lapserate03()/-(5.2);
  double T3 = this->LR26()/-(5.6);
  return T1*T2*T3;
}

double IndicesCollector::DEI(){
  double CPS = this->VirtualColdPoolStrength();
  double WXS = this->MU_WMAXSHEAR();
  double DEI = (1560*(CPS-13)+13*WXS)/10000;
  if(DEI<(-2))DEI=(-2);
  if(WXS==0)DEI=(-2);
  return DEI;
}

double IndicesCollector::TIP(){
  double Term1 = this->VMostUnstableCAPE();
  double Term2 = this->BS06();
  double Term3 = this->PWATER();
  double Term4 = this->SRH03RM();
  if(Term2<9)Term2=9;
  return (sqrt(Term1)/32) * (Term2/18) * (Term3/25) * (1 + Term4/300);
}

double IndicesCollector::DEI_eff(){
  double CPS = this->VirtualColdPoolStrength();
  double WXS = this->MU_EFF_EWMAXSHEAR();
  double DEI = (1560*(CPS-13)+13*WXS)/10000;
  if(DEI<(-2))DEI=(-2);
  if(WXS==0)DEI=(-2);
  return DEI;
}

double IndicesCollector::SCP(){
  double mucape = this->VMostUnstableCAPE()/1000;
  double srh = this->SRH03RM()/50;
  double ewd = this->BS06();
  if(ewd<10)ewd=0;
  else if(ewd>20)ewd=1;
  else ewd/=20;
  return mucape*srh*ewd;
}

double IndicesCollector::SCPeff(){
  double mucape = this->VMostUnstableCAPE()/1000;
  double srh = this->SRH03RM_eff()/50;
  double ewd = this->emubs();
  double cin = this->VMostUnstableCIN();
  if(isnan(cin)) cin = 0;
  if(cin>-40)cin=1;
  else cin=(-40)/cin;
  if(ewd<10)ewd=0;
  else if(ewd>20)ewd=1;
  else ewd/=20;
  return mucape*srh*ewd*cin;
}

double IndicesCollector::SCP_LM(){
  double mucape = this->VMostUnstableCAPE()/1000;
  double srh = this->SRH03LM()/50;
  double ewd = this->BS06();
  if(ewd<10)ewd=0;
  else if(ewd>20)ewd=1;
  else ewd/=20;
  return mucape*srh*ewd;
}

double IndicesCollector::SCPeff_LM(){
  double mucape = this->VMostUnstableCAPE()/1000;
  double srh = this->SRH03LM_eff()/50;
  double ewd = this->emubs();
  double cin = this->VMostUnstableCIN();
  if(isnan(cin)) cin = 0;
  if(cin>-40)cin=1;
  else cin=(-40)/cin;
  if(ewd<10)ewd=0;
  else if(ewd>20)ewd=1;
  else ewd/=20;
  return mucape*srh*ewd*cin;
}

double IndicesCollector::SHP(){
  double mucape = this->VMostUnstableCAPE();
  double mumr = this->MUMRatio();
  double lr = -(this->lapseRate500700());
  double t500 = Get(S->t,cache->getPressureIndex(500));
  double dls = this->BS06();
  if(dls<7)dls=7;
  else if(dls>27)dls=27;
  if(mumr<11)mumr =11;
  else if(mumr>13.6)mumr = 13.6;
  if(t500>-5.5)t500 = -5.5;
  double ship = (mucape*mumr*lr*(-t500)*dls)/42000000;
  if(mucape<1300)ship*=(mucape/1300);
  if(lr<5.8)ship*=(lr/5.8);
  if(this->ZeroHeight()<2400)ship*=(this->ZeroHeight()/2400);
  return ship;
}

double IndicesCollector::DCP(){
  double dcape = this->VDCAPE()/980;
  double mucape = this->VMostUnstableCAPE()/2000;
  double shear = this->BS06()/(20*0.514444444 );
  double mw = this->MeanWind06()/(16*0.514444444 );
  return dcape*mucape*shear*mw;
}

double IndicesCollector::DCP_eff(){
  double dcape = this->VDCAPE()/980;
  double* CAPE_WXS = this->MU_ECAPE(); 
  double CAPE = CAPE_WXS[5];
  delete[] CAPE_WXS;
  double mucape = CAPE/2000;
  double shear = this->emubs()/(20*0.514444444);
  double mw = this->MeanWind06()/(16*0.514444444);
  return dcape*mucape*shear*mw;
}

double IndicesCollector::MU_WMAXSHEAR(){
  return this->VMostUnstableVmax()*this->BS06();
}

double IndicesCollector::MUML_WMAXSHEAR(){
  return this->VMeanMostUnstableVmax()*this->BS06();
}

double IndicesCollector::SB_WMAXSHEAR(){
  return this->VSurfaceBasedVmax()*this->BS06();
}

double IndicesCollector::ML_WMAXSHEAR(){
  return this->VMeanLayerVmax()*this->BS06();
}

double IndicesCollector::MUML_EFF_EWMAXSHEAR(){
  double* CAPE_WXS = this->MU_ML_ECAPE(); 
  double CAPE = CAPE_WXS[5];
  delete[] CAPE_WXS;
  return CAPE*this->emumlbs();
}

double IndicesCollector::MU_EFF_EWMAXSHEAR(){
  double* CAPE_WXS = this->MU_ECAPE(); 
  double CAPE = CAPE_WXS[5];
  delete[] CAPE_WXS;
  return CAPE*this->emubs();
}

double IndicesCollector::SB_EFF_EWMAXSHEAR(){
  double* CAPE_WXS = this->SB_ECAPE(); 
  double CAPE = CAPE_WXS[5];
  delete[] CAPE_WXS;
  return CAPE*this->esbbs();
}

double IndicesCollector::ML_EFF_EWMAXSHEAR(){
  double* CAPE_WXS = this->ML_ECAPE(); 
  double CAPE = CAPE_WXS[5];
  delete[] CAPE_WXS;
  return CAPE*this->emlbs();
}

double IndicesCollector::MUML_EFF_EWMAXSHEAR_HGL(){
  double* CAPE_WXS = this->MU_ML_ECAPE(); 
  double CAPE = sqrt(2*CAPE_WXS[3]);
  delete[] CAPE_WXS;
  return CAPE*this->emumlbs();
}

double IndicesCollector::MU_EFF_EWMAXSHEAR_HGL(){
  double* CAPE_WXS = this->MU_ECAPE(); 
  double CAPE = sqrt(2*CAPE_WXS[3]);
  delete[] CAPE_WXS;
  return CAPE*this->emubs();
}

double IndicesCollector::SB_EFF_EWMAXSHEAR_HGL(){
  double* CAPE_WXS = this->SB_ECAPE(); 
  double CAPE = sqrt(2*CAPE_WXS[3]);
  delete[] CAPE_WXS;
  return CAPE*this->esbbs();
}

double IndicesCollector::ML_EFF_EWMAXSHEAR_HGL(){
  double* CAPE_WXS = this->ML_ECAPE(); 
  double CAPE = sqrt(2*CAPE_WXS[3]);
  delete[] CAPE_WXS;
  return CAPE*this->emlbs();
}

double IndicesCollector::MUML_EFF_EWMAXSHEAR_3km(){
  double* CAPE_WXS = this->MU_ML_ECAPE(); 
  double CAPE = sqrt(2*CAPE_WXS[6]);
  delete[] CAPE_WXS;
  return CAPE*this->emumlbs();
}

double IndicesCollector::MU_EFF_EWMAXSHEAR_3km(){
  double* CAPE_WXS = this->MU_ECAPE(); 
  double CAPE = sqrt(2*CAPE_WXS[6]);
  delete[] CAPE_WXS;
  return CAPE*this->emubs();
}

double IndicesCollector::SB_EFF_EWMAXSHEAR_3km(){
  double* CAPE_WXS = this->SB_ECAPE(); 
  double CAPE = sqrt(2*CAPE_WXS[6]);
  delete[] CAPE_WXS;
  return CAPE*this->esbbs();
}

double IndicesCollector::ML_EFF_EWMAXSHEAR_3km(){
  double* CAPE_WXS = this->ML_ECAPE(); 
  double CAPE = sqrt(2*CAPE_WXS[6]);
  delete[] CAPE_WXS;
  return CAPE*this->emlbs();
}

double IndicesCollector::MU_ELI(){
  double* CAPE_WXS = this->MU_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->VMostUnstableLI();
  if(buoyancy<0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::SB_ELI(){
  double* CAPE_WXS = this->SB_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->VSurfaceBasedLI();
  if(buoyancy<0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::ML_ELI(){
  double* CAPE_WXS = this->ML_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->VMeanLayerLI();
  if(buoyancy<0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MUML_ELI(){
  double* CAPE_WXS = this->MU_ML_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->VMeanMostUnstableLI();
  if(buoyancy<0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MU500_ELI(){
  double* CAPE_WXS = this->MU500_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MU500_LI();
  if(buoyancy<0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MU_ebuoyancy(){
  double* CAPE_WXS = this->MU_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MU_buoyancy();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::ML_ebuoyancy(){
  double* CAPE_WXS = this->ML_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->ML_buoyancy();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MUML_ebuoyancy(){
  double* CAPE_WXS = this->MU_ML_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MUML_buoyancy();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MU500_ebuoyancy(){
  double* CAPE_WXS = this->MU500_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MU500_buoyancy();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::SB_ebuoyancy(){
  double* CAPE_WXS = this->SB_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->SB_buoyancy();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MU_ebuoyancy_M10(){
  double* CAPE_WXS = this->MU_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MU_buoyancy_M10();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::ML_ebuoyancy_M10(){
  double* CAPE_WXS = this->ML_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->ML_buoyancy_M10();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MUML_ebuoyancy_M10(){
  double* CAPE_WXS = this->MU_ML_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MUML_buoyancy_M10();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MU500_ebuoyancy_M10(){
  double* CAPE_WXS = this->MU500_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MU500_buoyancy_M10();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::SB_ebuoyancy_M10(){
  double* CAPE_WXS = this->SB_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->SB_buoyancy_M10();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::BulkShearSfcTen(){
  int tail=0;
  int head = S->th->mintenpos;
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BulkShear1kmTen(){
  int tail=cache->getHeightIndex(1000);
  int head = S->th->mintenpos;
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::BulkShear2kmTen(){
  int tail=cache->getHeightIndex(2000);
  int head = S->th->mintenpos;
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  return result.abs();
}

double IndicesCollector::MoistureFlux(){
  double result = this->MLMixingRatio() * this->MeanWind500();
  return result;
}

double IndicesCollector::SR_moisture_flux(){
  double result = this->MLMixingRatio() * this->Peters_SR_inflow();
  return result;
}

double IndicesCollector::SR_moisture_flux_eff(){
  double result = (S->th->mmrMAX) * this->Peters_SR_inflow_eff();
  return result;
}

double IndicesCollector::SR_moisture_flux_RM(){
  double result = this->MLMixingRatio() * this->MeanSR500_RM();
  return result;
}

double IndicesCollector::SR_moisture_flux_eff_RM(){
  double result = (S->th->mmrMAX) * this->MeanSR500_RM();
  return result;
}

double IndicesCollector::SR_moisture_flux_LM(){
  double result = this->MLMixingRatio() * this->MeanSR500_LM();
  return result;
}

double IndicesCollector::SR_moisture_flux_eff_LM(){
  double result = (S->th->mmrMAX) * this->MeanSR500_LM();
  return result;
}

double IndicesCollector::SR_moisture_flux_MW(){
  double result = this->MLMixingRatio() * this->MeanSR500_MW();
  return result;
}

double IndicesCollector::SR_moisture_flux_eff_MW(){
  double result = (S->th->mmrMAX) * this->MeanSR500_MW();
  return result;
}

double IndicesCollector::RH01(){
  return S->th->meanhum1;
}

double IndicesCollector::RH02(){
  return S->th->meanhum2;
}

double IndicesCollector::RH25(){
  return S->th->meanhum25;
}

double IndicesCollector::RH36(){
  return S->th->meanhum36;
}

double IndicesCollector::RH14(){
  return S->th->meanhum14;
}

double IndicesCollector::RHMIDDLE(){
  return S->th->meanhumMIDDLE;
}

double IndicesCollector::RH850500(){
  return S->th->meanhum850500;
}

double IndicesCollector::Bunkers_RM_A(){
  double *tab = S->ks->rm.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Bunkers_RM_M(){
  double *tab = S->ks->rm.toAV(); 
  double magnitude = tab[1]; 
  delete[] tab;
  return magnitude;
}

double IndicesCollector::Bunkers_LM_A(){
  double *tab = S->ks->lm.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Bunkers_LM_M(){
  double *tab = S->ks->lm.toAV(); 
  double magnitude = tab[1]; 
  delete[] tab;
  return magnitude;
}

double IndicesCollector::Bunkers_MW_A(){
  double *tab = S->ks->mean06.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Bunkers_MW_M(){
  double *tab = S->ks->mean06.toAV();
  double magnitude = tab[1]; 
  delete[] tab;
  return magnitude;
}

double IndicesCollector::Corfidi_upwind_A(){
  double *tab = S->ks->Corfidi_upwind.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Corfidi_upwind_M(){
  double *tab = S->ks->Corfidi_upwind.toAV(); 
  double magnitude = tab[1]; 
  delete[] tab;
  return magnitude;
}

double IndicesCollector::Corfidi_downwind_A(){
  double *tab = S->ks->Corfidi_downwind.toAV(); 
  double angle = tab[0];
  delete[] tab;
  return angle;
}

double IndicesCollector::Corfidi_downwind_M(){
  double *tab = S->ks->Corfidi_downwind.toAV(); 
  double magnitude = tab[1]; 
  delete[] tab;
  return magnitude;
}

double IndicesCollector::SB_coldcape(){
  double coldcape = S->th->surfaceBased->coldcape;
  return coldcape;
}

double IndicesCollector::MU_coldcape(){
  double coldcape = S->th->mostUnstable->coldcape;
  return coldcape;
}

double IndicesCollector::MUML_coldcape(){
  double coldcape = S->th->meanmostUnstable->coldcape;
  return coldcape;
}

double IndicesCollector::ML_coldcape(){
  double coldcape = S->th->meanLayer->coldcape;
  return coldcape;
}

double IndicesCollector::MU500_coldcape(){
  double coldcape = S->th->mostU500->coldcape;
  return coldcape;
}

double IndicesCollector::SB_coldcapeTV(){
  double coldcape = S->th->surfaceBased->coldcapeTV;
  return coldcape;
}

double IndicesCollector::MU_coldcapeTV(){
  double coldcape = S->th->mostUnstable->coldcapeTV;
  return coldcape;
}

double IndicesCollector::MUML_coldcapeTV(){
  double coldcape = S->th->meanmostUnstable->coldcapeTV;
  return coldcape;
}

double IndicesCollector::ML_coldcapeTV(){
  double coldcape = S->th->meanLayer->coldcapeTV;
  return coldcape;
}

double IndicesCollector::MU500_coldcapeTV(){
  double coldcape = S->th->mostU500->coldcapeTV;
  return coldcape;
}

double IndicesCollector::BulkShearMULCLTen(){
  int tail=S->th->mostUnstable->vLclIndex;
  int head = S->th->mintenpos;
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  double effSHR = result.abs();
  double mucape = this->VMostUnstableCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::BulkShearMUMLLCLTen(){
  int tail=S->th->meanmostUnstable->vLclIndex;
  int head = S->th->mintenpos;
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  double effSHR = result.abs();
  return effSHR;
}

double IndicesCollector::BulkShearMLLCLTen(){
  int tail=S->th->meanLayer->vLclIndex;
  int head = S->th->mintenpos;
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  double effSHR = result.abs(); 
  return effSHR;
}

double IndicesCollector::BulkShearSBLCLTen(){
  int tail=S->th->surfaceBased->vLclIndex;
  int head = S->th->mintenpos;
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  double effSHR = result.abs();
  return effSHR;
}

double IndicesCollector::HSI(){
  double CAPE = this->VMostUnstableCAPE();
  double BS06 = this->BS06();
  double FL = this->ZeroHeight();
  double LCL = this->VMostUnstableLCL();
  double EL = this->VMostUnstableEL();
  double LR = -(this->lapseRate500800());
  if(CAPE<201)CAPE=201;
  else if(CAPE>4000)CAPE=4000;
  if(BS06<11)BS06=11;
  else if(BS06>27)BS06=27;
  if(FL<500)FL=500;
  else if(FL>4000)FL=4000;
  if(LCL<500)LCL=500;
  else if(LCL>1500)LCL=1500;
  if(LR<5)LR=5;
  else if(LR>8)LR=8;
  double HSI = ((sqrt(10*(CAPE-200)) * (BS06-5) * (7000-FL+LCL))/194000) * sqrt(EL*(((LR-4)*(LR-4))/10000000));
  if(isnan(HSI)) HSI = 0;
  return HSI;
}

double IndicesCollector::HSIv2(){  
  double* CAPE_HSI = this->MU_ML_ECAPE(); 
  double CAPE = CAPE_HSI[2];
  delete[] CAPE_HSI;
  double BS06 = this->MSR_MW();
  double FL = this->ZeroHeight();
  double LCL = this->VMeanMostUnstableLCL();
  double EL = this->VMeanMostUnstableEL();
  double LR = -(this->LR26());
  if(CAPE<201)CAPE=201;
  else if(CAPE>2750)CAPE=2750;
  if(BS06<2)BS06=2;
  else if(BS06>16)BS06=16;
  if(FL<500)FL=500;
  else if(FL>4000)FL=4000;
  if(LCL<500)LCL=500;
  else if(LCL>1500)LCL=1500;
  if(LR<5)LR=5;
  else if(LR>8)LR=8;
  double HSI = ((sqrt(95*(CAPE-200)) * (BS06+2) * (7000-FL+LCL))/194000) * sqrt(EL*(((LR-3)*(LR-3))/10000000));
  HSI = HSI*0.48;
  if(isnan(HSI)) HSI = 0;
  //double HSI_old = this->HSI();
  //if(HSI_old > HSI) HSI = HSI_old;
  return HSI;
}

double IndicesCollector::MSR_MW(){
  Vector res = S->ks->mean02 - S->ks->mean06;
  return res.abs();
}

double IndicesCollector::MSR_RM(){
  Vector res = S->ks->mean02 - S->ks->rm;
  return res.abs();
}

double IndicesCollector::MSR_LM(){
  Vector res = S->ks->mean02 - S->ks->lm;
  return res.abs();
}

double IndicesCollector::MSR_MW_HGL(){
  Vector res = S->ks->mean020 - S->ks->mean06;
  return res.abs();
}

double IndicesCollector::MSR_RM_HGL(){
  Vector res = S->ks->mean020 - S->ks->rm;
  return res.abs();
}

double IndicesCollector::MSR_LM_HGL(){
  Vector res = S->ks->mean020 - S->ks->lm;
  return res.abs();
}

double IndicesCollector::MU500_CAPE(){
  return S->th->mostU500->vcape;
}

double IndicesCollector::MU500_CIN(){
  double result = S->th->mostU500->vcin;
  double cape = this->MU500_CAPE();
  if(cape==0){
    result = sqrt(-1); 
  }
  return result;
}

double IndicesCollector::MU500_CIN500(){
  double result = S->th->mostU500->vcin500;
  return result;
}

double IndicesCollector::MU500_LI(){
  int lindex = cache->getPressureIndex(500);  
  double lit = Get(S->t,lindex);
  int vindex = lindex - S->th->mostU500->startIndex;
  double plit = Get(S->th->mostU500->virtualValues,vindex);
  double Showalter = lit - plit;
  return Showalter;
}

double IndicesCollector::EHI03(){
  return (this->VSurfaceBasedCAPE()*this->SRH03RM())/160000;
}

double IndicesCollector::EHI01(){
  return (this->VSurfaceBasedCAPE()*this->SRH01RM())/160000;
}

double IndicesCollector::EHI500(){
  return (this->VSurfaceBasedCAPE()*this->SRH500RM())/160000;
}

double IndicesCollector::EHI03_LM(){
  return (this->VSurfaceBasedCAPE()*this->SRH03LM())/160000;
}

double IndicesCollector::EHI01_LM(){
  return (this->VSurfaceBasedCAPE()*this->SRH01LM())/160000;
}

double IndicesCollector::EHI500_LM(){
  return (this->VSurfaceBasedCAPE()*this->SRH500LM())/160000;
}

double IndicesCollector::SW100_RM_F(){
  return S->ks->sw100rm2 / 100;
}

double IndicesCollector::SW500_RM_F(){
  return S->ks->sw500rm2 / 500;
}

double IndicesCollector::SW100_LM_F(){
  return S->ks->sw100lm2 / 100;
}

double IndicesCollector::SW500_LM_F(){
  return S->ks->sw500lm2 / 500;
}

double IndicesCollector::SW100_RM(){
  return S->ks->sw100rm / 100;
}

double IndicesCollector::SW500_RM(){
  return S->ks->sw500rm / 500;
}

double IndicesCollector::SW01_RM(){
  return S->ks->sw01rm / 1000;
}

double IndicesCollector::SW03_RM(){
  return S->ks->sw03rm / 3000;
}

double IndicesCollector::SW100_LM(){
  return S->ks->sw100lm / 100;
}

double IndicesCollector::SW500_LM(){
  return S->ks->sw500lm / 500;
}

double IndicesCollector::SW01_LM(){
  return S->ks->sw01lm / 1000;
}

double IndicesCollector::SW03_LM(){
  return S->ks->sw03lm / 3000;
}

double IndicesCollector::CA500_RM(){
  double *SRW_A = (S->ks->rm - Get(S->ks->vw,0)).toAV();   
  double *BS500_A = (Get(S->ks->vw,cache->getHeightIndex(500)) - Get(S->ks->vw,0)).toAV();
  double magnitude = abs(BS500_A[0] - SRW_A[0]);
  if(magnitude>180)magnitude=abs(magnitude-360);  
  return magnitude;	
}

double IndicesCollector::CA500_LM(){
  double *SRW_A = (S->ks->lm - Get(S->ks->vw,0)).toAV();   
  double *BS500_A = (Get(S->ks->vw,cache->getHeightIndex(500)) - Get(S->ks->vw,0)).toAV();
  double magnitude = abs(BS500_A[0] - SRW_A[0]);
  if(magnitude>180)magnitude=abs(magnitude-360);  
  return magnitude;	
}

double IndicesCollector::SRW_sfc_RM(){
  Vector res = (Get(S->ks->vw,0) - S->ks->rm);   
  return res.abs();
}

double IndicesCollector::SRW_sfc_LM(){
  Vector res = (Get(S->ks->vw,0) - S->ks->lm);
  return res.abs();
}

double IndicesCollector::MeanSR500_RM(){
  Vector res = S->ks->mean0 - S->ks->rm;
  return res.abs();
}

double IndicesCollector::MeanSR01_RM(){
  Vector res = S->ks->mean01 - S->ks->rm;
  return res.abs();
}

double IndicesCollector::MeanSR03_RM(){
  Vector res = S->ks->mean03 - S->ks->rm;
  return res.abs();
}

double IndicesCollector::MeanSR36_RM(){
  Vector res = S->ks->mean36 - S->ks->rm;
  return res.abs();
}

double IndicesCollector::MeanSR16_RM(){
  Vector res = S->ks->mean16 - S->ks->rm;
  return res.abs();
}

double IndicesCollector::MeanSR500_LM(){
  Vector res = S->ks->mean0 - S->ks->lm;
  return res.abs();
}

double IndicesCollector::MeanSR01_LM(){
  Vector res = S->ks->mean01 - S->ks->lm;
  return res.abs();
}

double IndicesCollector::MeanSR03_LM(){
  Vector res = S->ks->mean03 - S->ks->lm;
  return res.abs();
}

double IndicesCollector::MeanSR36_LM(){
  Vector res = S->ks->mean36 - S->ks->lm;
  return res.abs();
}

double IndicesCollector::MeanSR36_MW(){
  Vector res = S->ks->mean36 - S->ks->mean06;
  return res.abs();
}

double IndicesCollector::MeanSR16_LM(){
  Vector res = S->ks->mean16 - S->ks->lm;
  return res.abs();
}

double IndicesCollector::MeanSR500_MW(){
  Vector res = S->ks->mean0 - S->ks->mean06;
  return res.abs();
}

double IndicesCollector::MeanSR01_MW(){
  Vector res = S->ks->mean01 - S->ks->mean06;
  return res.abs();
}

double IndicesCollector::MeanSR03_MW(){
  Vector res = S->ks->mean03 - S->ks->mean06;
  return res.abs();
}

double IndicesCollector::MeanVMSR500_RM(){
  return S->ks->SR_500_RM / S->ks->n500;
}

double IndicesCollector::MeanVMSR01_RM(){
  return S->ks->SR_1000_RM / S->ks->n1000;
}

double IndicesCollector::MeanVMSR03_RM(){
  return S->ks->SR_3000_RM / S->ks->n3000;
}

double IndicesCollector::MeanVMSR36_RM(){
  return S->ks->SR_36_RM / S->ks->n36sr;
}

double IndicesCollector::MeanVMSR16_RM(){
  return S->ks->SR_16_RM / S->ks->n16sr;
}

double IndicesCollector::MeanVMSR500_LM(){
  return S->ks->SR_500_LM / S->ks->n500;
}

double IndicesCollector::MeanVMSR01_LM(){
  return S->ks->SR_1000_LM / S->ks->n1000;
}

double IndicesCollector::MeanVMSR03_LM(){
  return S->ks->SR_3000_LM / S->ks->n3000;
}

double IndicesCollector::MeanVMSR36_LM(){
  return S->ks->SR_36_LM / S->ks->n36sr;
}

double IndicesCollector::MeanVMSR16_LM(){
  return S->ks->SR_16_LM / S->ks->n16sr;
}

double IndicesCollector::SV_100_RM_FRA(){
  return S->ks->sw100rm / S->ks->shear100m;
}

double IndicesCollector::SV_500_RM_FRA(){
  return S->ks->sw500rm / S->ks->shear500m;
}

double IndicesCollector::SV_1000_RM_FRA(){
  return S->ks->sw01rm / S->ks->shear1000m;
}

double IndicesCollector::SV_3000_RM_FRA(){
  return S->ks->sw03rm / S->ks->shear3000m;
}

double IndicesCollector::SV_100_LM_FRA(){
  return S->ks->sw100lm / S->ks->shear100m;
}

double IndicesCollector::SV_500_LM_FRA(){
  return S->ks->sw500lm / S->ks->shear500m;
}

double IndicesCollector::SV_1000_LM_FRA(){
  return S->ks->sw01lm / S->ks->shear1000m;
}

double IndicesCollector::SV_3000_LM_FRA(){
  return S->ks->sw03lm / S->ks->shear3000m;
}

double IndicesCollector::MeanSR0500_LM_eff(){
  Vector res = S->ks->mean01eff - S->ks->lm;
  return res.abs();
}

double IndicesCollector::MeanSR0500_RM_eff(){
  Vector res = S->ks->mean01eff - S->ks->rm;
  return res.abs();
}

double IndicesCollector::MeanSR0500_MW_eff(){
  Vector res = S->ks->mean01eff - S->ks->mean06;
  return res.abs();
}

double IndicesCollector::MU_buoyancy(){
  double diff = S->th->mostUnstable->peakB;
  return diff;
}

double IndicesCollector::MUML_buoyancy(){
  double diff = S->th->meanmostUnstable->peakB;
  return diff;
}

double IndicesCollector::MU500_buoyancy(){
  double diff = S->th->mostU500->peakB;
  return diff;
}

double IndicesCollector::ML_buoyancy(){
  double diff = S->th->meanLayer->peakB;
  return diff;
}

double IndicesCollector::SB_buoyancy(){
  double diff = S->th->surfaceBased->peakB;
  return diff;
}

double IndicesCollector::MU_buoyancy_3km(){
  double diff = S->th->mostUnstable->peakB_3km;
  return diff;
}

double IndicesCollector::MUML_buoyancy_3km(){
  double diff = S->th->meanmostUnstable->peakB_3km;
  return diff;
}

double IndicesCollector::MU500_buoyancy_3km(){
  double diff = S->th->mostU500->peakB_3km;
  return diff;
}

double IndicesCollector::ML_buoyancy_3km(){
  double diff = S->th->meanLayer->peakB_3km;
  return diff;
}

double IndicesCollector::SB_buoyancy_3km(){
  double diff = S->th->surfaceBased->peakB_3km;
  return diff;
}

double IndicesCollector::MU_buoyancy_M10(){
  double diff = S->th->mostUnstable->peakB_M10;
  return diff;
}

double IndicesCollector::MUML_buoyancy_M10(){
  double diff = S->th->meanmostUnstable->peakB_M10;
  return diff;
}

double IndicesCollector::MU500_buoyancy_M10(){
  double diff = S->th->mostU500->peakB_M10;
  return diff;
}

double IndicesCollector::ML_buoyancy_M10(){
  double diff = S->th->meanLayer->peakB_M10;
  return diff;
}

double IndicesCollector::SB_buoyancy_M10(){
  double diff = S->th->surfaceBased->peakB_M10;
  return diff;
}

double IndicesCollector::Ventilation_16km_RM(){
  double ventilation = distance(S->ks->mean16, S->ks->rm, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_36km_RM(){
  double ventilation = distance(S->ks->mean36, S->ks->rm, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_69km_RM(){
  double ventilation = distance(S->ks->mean69, S->ks->rm, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_16km_LM(){
  double ventilation = distance(S->ks->mean16, S->ks->lm, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_36km_LM(){
  double ventilation = distance(S->ks->mean36, S->ks->lm, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_69km_LM(){
  double ventilation = distance(S->ks->mean69, S->ks->lm, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_16km_MW(){
  double ventilation = distance(S->ks->mean16, S->ks->mean06, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_36km_MW(){
  double ventilation = distance(S->ks->mean36, S->ks->mean06, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_69km_MW(){
  double ventilation = distance(S->ks->mean69, S->ks->mean06, S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_16km(){
  double ventilation = distance(S->ks->mean16, this->Peters_vector(), S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_36km(){
  double ventilation = distance(S->ks->mean36, this->Peters_vector(), S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::Ventilation_69km(){
  double ventilation = distance(S->ks->mean69, this->Peters_vector(), S->ks->mean0);   
  return ventilation*(-1);
}

double IndicesCollector::MU_ebuoyancy_3km(){
  double* CAPE_WXS = this->MU_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MU_buoyancy_3km();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::ML_ebuoyancy_3km(){
  double* CAPE_WXS = this->ML_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->ML_buoyancy_3km();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MUML_ebuoyancy_3km(){
  double* CAPE_WXS = this->MU_ML_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MUML_buoyancy_3km();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::MU500_ebuoyancy_3km(){
  double* CAPE_WXS = this->MU500_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->MU500_buoyancy_3km();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double IndicesCollector::SB_ebuoyancy_3km(){
  double* CAPE_WXS = this->SB_ECAPE();
  double E_tilde = CAPE_WXS[0];
  delete[] CAPE_WXS;
  if(E_tilde>2)E_tilde = 2;
  double buoyancy = this->SB_buoyancy_3km();
  if(buoyancy>0)buoyancy = E_tilde * buoyancy;
  return buoyancy;
}

double * processSounding(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_, int length, double dz, Sounding **S, double* meanlayer_bottom_top, Vector storm_motion){
  *S = new Sounding(p_,h_,t_,d_,a_,v_,length, dz, meanlayer_bottom_top, storm_motion);
  double * vec = new double[315];

// vec[0] = MU_ECAPE[5]; // EWMAX
// vec[0] = ML_ECAPE[5]; // EWMAX
// vec[0] = MU_ML_ECAPE[5]; // EWMAX
// vec[0] = SB_ECAPE[5]; // EWMAX
// vec[0] = MU500_ECAPE[5]; // EWMAX
// vec[0]=(*S)->getIndicesCollectorPointer()->VMostUnstableVmax();
// vec[0]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableVmax();
// vec[0]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedVmax(); 
// vec[0]=(*S)->getIndicesCollectorPointer()->VMeanLayerVmax();
// vec[0]=(*S)->getIndicesCollectorPointer()->VMeanLayerVmax();
// vec[0]=(*S)->getIndicesCollectorPointer()->BS06_var_SD();
// vec[0]=(*S)->getIndicesCollectorPointer()->BS16_var_SD();
// vec[0]=(*S)->getIndicesCollectorPointer()->Peters_vector_eff_A();
// vec[0]=(*S)->getIndicesCollectorPointer()->Peters_vector_eff_M();
// vec[0]=(*S)->getIndicesCollectorPointer()->TIP();
//  vec[2]=(*S)->getIndicesCollectorPointer()->SB_coldcapeTV();
// vec[30] = SB_ECAPE[0]; // E_tilde
// vec[31] = SB_ECAPE[1]; // Radius
//  vec[37]=(*S)->getIndicesCollectorPointer()->ML_coldcapeTV();

// SB parcel
  vec[0]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedCAPE();
  vec[1]=(*S)->getIndicesCollectorPointer()->VLLSurfaceBasedCAPE();
  vec[2]=(*S)->getIndicesCollectorPointer()->SBmiddlecape();
  vec[3]=(*S)->getIndicesCollectorPointer()->SB_buoyancy();
  vec[4]=(*S)->getIndicesCollectorPointer()->SB_buoyancy_M10();
  vec[5]=(*S)->getIndicesCollectorPointer()->SB_buoyancy_3km();
  vec[6]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedLI();
  vec[7]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedLI_M25();
  vec[8]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedCIN();
  vec[9]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedCIN500();
  vec[10]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedLCL();
  vec[11]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedLFC();
  vec[12]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedEL();
  vec[13]=(*S)->getIndicesCollectorPointer()->SBLCLTemperature();
  vec[14]=(*S)->getIndicesCollectorPointer()->SBLFCTemperature(); 
  vec[15]=(*S)->getIndicesCollectorPointer()->SBELTemperature();
  vec[16]=(*S)->getIndicesCollectorPointer()->SB_cold_cloud(); 
  vec[17]=(*S)->getIndicesCollectorPointer()->SB_warm_cloud(); 
  vec[18]=(*S)->getIndicesCollectorPointer()->SB_equal_layer(); 
  vec[19]=(*S)->getIndicesCollectorPointer()->SBMRatio();
  vec[20]=(*S)->getIndicesCollectorPointer()->SB_WMAXSHEAR();
  double* SB_ECAPE = (*S)->getIndicesCollectorPointer()->SB_ECAPE(); 
  vec[21] = SB_ECAPE[2]; // CAPE
  vec[22] = SB_ECAPE[6]; // CAPE_3km
  vec[23] = SB_ECAPE[3]; // CAPE_HGL
  vec[24]=(*S)->getIndicesCollectorPointer()->SB_ebuoyancy();
  vec[25]=(*S)->getIndicesCollectorPointer()->SB_ebuoyancy_M10();
  vec[26]=(*S)->getIndicesCollectorPointer()->SB_ebuoyancy_3km();
  vec[27]=(*S)->getIndicesCollectorPointer()->SB_ELI();
  delete[] SB_ECAPE;
  vec[28]=(*S)->getIndicesCollectorPointer()->SB_EFF_EWMAXSHEAR();
  vec[29]=(*S)->getIndicesCollectorPointer()->SB_EFF_EWMAXSHEAR_HGL();
  vec[30]=(*S)->getIndicesCollectorPointer()->SB_EFF_EWMAXSHEAR_3km();

// ML parcel
  vec[31]=(*S)->getIndicesCollectorPointer()->VMeanLayerCAPE();
  vec[32]=(*S)->getIndicesCollectorPointer()->VLLMeanLayerCAPE();
  vec[33]=(*S)->getIndicesCollectorPointer()->MLmiddlecape();
  vec[34]=(*S)->getIndicesCollectorPointer()->ML_buoyancy();
  vec[35]=(*S)->getIndicesCollectorPointer()->ML_buoyancy_M10();
  vec[36]=(*S)->getIndicesCollectorPointer()->ML_buoyancy_3km();
  vec[37]=(*S)->getIndicesCollectorPointer()->VMeanLayerLI();
  vec[38]=(*S)->getIndicesCollectorPointer()->VMeanLayerLI_M25();
  vec[39]=(*S)->getIndicesCollectorPointer()->VMeanLayerCIN();  
  vec[40]=(*S)->getIndicesCollectorPointer()->VMeanLayerCIN500();  
  vec[41]=(*S)->getIndicesCollectorPointer()->VMeanLayerLCL();
  vec[42]=(*S)->getIndicesCollectorPointer()->VMeanLayerLFC();
  vec[43]=(*S)->getIndicesCollectorPointer()->VMeanLayerEL();
  vec[44]=(*S)->getIndicesCollectorPointer()->MLLCLTemperature();
  vec[45]=(*S)->getIndicesCollectorPointer()->MLLFCTemperature(); 
  vec[46]=(*S)->getIndicesCollectorPointer()->MLELTemperature();
  vec[47]=(*S)->getIndicesCollectorPointer()->ML_cold_cloud(); 
  vec[48]=(*S)->getIndicesCollectorPointer()->ML_warm_cloud(); 
  vec[49]=(*S)->getIndicesCollectorPointer()->ML_equal_layer(); 
  vec[50]=(*S)->getIndicesCollectorPointer()->MLMixingRatio();
  vec[51]=(*S)->getIndicesCollectorPointer()->ML_WMAXSHEAR(); 
  double* ML_ECAPE = (*S)->getIndicesCollectorPointer()->ML_ECAPE(); 
  vec[52] = ML_ECAPE[2]; // CAPE
  vec[53] = ML_ECAPE[6]; // CAPE_3km
  vec[54] = ML_ECAPE[3]; // CAPE_HGL
  vec[55]=(*S)->getIndicesCollectorPointer()->ML_ebuoyancy();
  vec[56]=(*S)->getIndicesCollectorPointer()->ML_ebuoyancy_M10();
  vec[57]=(*S)->getIndicesCollectorPointer()->ML_ebuoyancy_3km();
  vec[58]=(*S)->getIndicesCollectorPointer()->ML_ELI();
  delete[] ML_ECAPE;
  vec[59]=(*S)->getIndicesCollectorPointer()->ML_EFF_EWMAXSHEAR();
  vec[60]=(*S)->getIndicesCollectorPointer()->ML_EFF_EWMAXSHEAR_HGL();
  vec[61]=(*S)->getIndicesCollectorPointer()->ML_EFF_EWMAXSHEAR_3km();

// MU parcel
  vec[62]=(*S)->getIndicesCollectorPointer()->VMostUnstableCAPE();
  vec[63]=(*S)->getIndicesCollectorPointer()->VLLMostUnstableCAPE();
  vec[64]=(*S)->getIndicesCollectorPointer()->MUmiddlecape();
  vec[65]=(*S)->getIndicesCollectorPointer()->MU_buoyancy(); 
  vec[66]=(*S)->getIndicesCollectorPointer()->MU_buoyancy_M10();
  vec[67]=(*S)->getIndicesCollectorPointer()->MU_buoyancy_3km();
  vec[68]=(*S)->getIndicesCollectorPointer()->VMostUnstableLI();
  vec[69]=(*S)->getIndicesCollectorPointer()->VMostUnstableLI_M25();   
  vec[70]=(*S)->getIndicesCollectorPointer()->VMostUnstableCIN();  
  vec[71]=(*S)->getIndicesCollectorPointer()->VMostUnstableCIN500();  
  vec[72]=(*S)->getIndicesCollectorPointer()->VMostUnstableLCL();
  vec[73]=(*S)->getIndicesCollectorPointer()->VMostUnstableLFC();
  vec[74]=(*S)->getIndicesCollectorPointer()->VMostUnstableEL();
  vec[75]=(*S)->getIndicesCollectorPointer()->MULCLTemperature();
  vec[76]=(*S)->getIndicesCollectorPointer()->MULFCTemperature();
  vec[77]=(*S)->getIndicesCollectorPointer()->MUELTemperature();
  vec[78]=(*S)->getIndicesCollectorPointer()->MU_cold_cloud(); 
  vec[79]=(*S)->getIndicesCollectorPointer()->MU_warm_cloud(); 
  vec[80]=(*S)->getIndicesCollectorPointer()->MU_equal_layer();
  vec[81]=(*S)->getIndicesCollectorPointer()->MUMRatio();   
  vec[82]=(*S)->getIndicesCollectorPointer()->MU_WMAXSHEAR();
  double* MU_ECAPE = (*S)->getIndicesCollectorPointer()->MU_ECAPE(); 
  vec[83] = MU_ECAPE[2]; // CAPE
  vec[84] = MU_ECAPE[6]; // CAPE_3km
  vec[85] = MU_ECAPE[3]; // CAPE_HGL
  vec[86]=(*S)->getIndicesCollectorPointer()->MU_ebuoyancy();
  vec[87]=(*S)->getIndicesCollectorPointer()->MU_ebuoyancy_M10();
  vec[88]=(*S)->getIndicesCollectorPointer()->MU_ebuoyancy_3km();
  vec[89]=(*S)->getIndicesCollectorPointer()->MU_ELI();
  delete[] MU_ECAPE;
  vec[90]=(*S)->getIndicesCollectorPointer()->MU_EFF_EWMAXSHEAR();
  vec[91]=(*S)->getIndicesCollectorPointer()->MU_EFF_EWMAXSHEAR_HGL();
  vec[92]=(*S)->getIndicesCollectorPointer()->MU_EFF_EWMAXSHEAR_3km();

// MU_ML parcel
  vec[93]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableCAPE();
  vec[94]=(*S)->getIndicesCollectorPointer()->VLLMeanMostUnstableCAPE();
  vec[95]=(*S)->getIndicesCollectorPointer()->MUMLmiddlecape();
  vec[96]=(*S)->getIndicesCollectorPointer()->MUML_buoyancy();
  vec[97]=(*S)->getIndicesCollectorPointer()->MUML_buoyancy_M10();
  vec[98]=(*S)->getIndicesCollectorPointer()->MUML_buoyancy_3km();
  vec[99]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableLI();
  vec[100]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableLI_M25(); 
  vec[101]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableCIN();  
  vec[102]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableCIN500();  
  vec[103]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableLCL();
  vec[104]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableLFC();
  vec[105]=(*S)->getIndicesCollectorPointer()->VMeanMostUnstableEL();
  vec[106]=(*S)->getIndicesCollectorPointer()->MUMLLCLTemperature();
  vec[107]=(*S)->getIndicesCollectorPointer()->MUMLLFCTemperature();
  vec[108]=(*S)->getIndicesCollectorPointer()->MUMLELTemperature();
  vec[109]=(*S)->getIndicesCollectorPointer()->MU_ML_cold_cloud(); 
  vec[110]=(*S)->getIndicesCollectorPointer()->MU_ML_warm_cloud(); 
  vec[111]=(*S)->getIndicesCollectorPointer()->MU_ML_equal_layer(); 
  vec[112]=(*S)->getIndicesCollectorPointer()->MUMLMRatio();   
  vec[113]=(*S)->getIndicesCollectorPointer()->MUML_WMAXSHEAR();
  double* MU_ML_ECAPE = (*S)->getIndicesCollectorPointer()->MU_ML_ECAPE(); 
  vec[114] = MU_ML_ECAPE[2]; // CAPE
  vec[115] = MU_ML_ECAPE[6]; // CAPE_3km
  vec[116] = MU_ML_ECAPE[3]; // CAPE_HGL
  vec[117]=(*S)->getIndicesCollectorPointer()->MUML_ebuoyancy();
  vec[118]=(*S)->getIndicesCollectorPointer()->MUML_ebuoyancy_M10();
  vec[119]=(*S)->getIndicesCollectorPointer()->MUML_ebuoyancy_3km();
  vec[120]=(*S)->getIndicesCollectorPointer()->MUML_ELI();
  delete[] MU_ML_ECAPE;
  vec[121]=(*S)->getIndicesCollectorPointer()->MUML_EFF_EWMAXSHEAR();
  vec[122]=(*S)->getIndicesCollectorPointer()->MUML_EFF_EWMAXSHEAR_HGL();
  vec[123]=(*S)->getIndicesCollectorPointer()->MUML_EFF_EWMAXSHEAR_3km();
  
// MU500 parcel
  vec[124]=(*S)->getIndicesCollectorPointer()->MU500_CAPE();
  vec[125]=(*S)->getIndicesCollectorPointer()->MU500_coldcape();
  vec[126]=(*S)->getIndicesCollectorPointer()->MU500middlecape();
  vec[127]=(*S)->getIndicesCollectorPointer()->MU500_buoyancy();
  vec[128]=(*S)->getIndicesCollectorPointer()->MU500_buoyancy_M10();
  vec[129]=(*S)->getIndicesCollectorPointer()->MU500_LI();
  vec[130]=(*S)->getIndicesCollectorPointer()->VMostU500LI_M25();
  vec[131]=(*S)->getIndicesCollectorPointer()->MU500_CIN();
  vec[132]=(*S)->getIndicesCollectorPointer()->MU500_CIN500();
  double* MU500_ECAPE = (*S)->getIndicesCollectorPointer()->MU500_ECAPE(); 
  vec[133] = MU500_ECAPE[2]; // CAPE
  vec[134] = MU500_ECAPE[3]; // CAPE_HGL
  vec[135]=(*S)->getIndicesCollectorPointer()->MU500_ebuoyancy();
  vec[136]=(*S)->getIndicesCollectorPointer()->MU500_ebuoyancy_M10();
  vec[137]=(*S)->getIndicesCollectorPointer()->MU500_ELI();
  delete[] MU500_ECAPE;

  // Lapse rates
  vec[138]=(*S)->getIndicesCollectorPointer()->LR0500();
  vec[139]=(*S)->getIndicesCollectorPointer()->LapseRate01();
  vec[140]=(*S)->getIndicesCollectorPointer()->lapserate03();
  vec[141]=(*S)->getIndicesCollectorPointer()->LR04(); 
  vec[142]=(*S)->getIndicesCollectorPointer()->LR06();
  vec[143]=(*S)->getIndicesCollectorPointer()->LR16(); 
  vec[144]=(*S)->getIndicesCollectorPointer()->LapseRate24(); 
  vec[145]=(*S)->getIndicesCollectorPointer()->LR26();
  vec[146]=(*S)->getIndicesCollectorPointer()->LR36();
  vec[147]=(*S)->getIndicesCollectorPointer()->max_LR26_2km();
  vec[148]=(*S)->getIndicesCollectorPointer()->lapseRate500700();
  vec[149]=(*S)->getIndicesCollectorPointer()->lapseRate500800(); 

  // Other & theta-e
  vec[150]=(*S)->getIndicesCollectorPointer()->ZeroHeight();
  vec[151]=(*S)->getIndicesCollectorPointer()->WetBulbZeroHeight();  
  vec[152]=(*S)->getIndicesCollectorPointer()->M10Height();
  vec[153]=(*S)->getIndicesCollectorPointer()->MUHeight();
  vec[154]=(*S)->getIndicesCollectorPointer()->MUMLHeight();

  // Theta-e & downdraft 
  vec[155]=(*S)->getIndicesCollectorPointer()->DeltaThetaE();
  vec[156]=(*S)->getIndicesCollectorPointer()->DeltaThetaE_min04km();
  vec[157]=(*S)->getIndicesCollectorPointer()->THTE_LR_LCL_to_M10();
  vec[158]=(*S)->getIndicesCollectorPointer()->THTE_LR_MU_to_M10();
  vec[159]=(*S)->getIndicesCollectorPointer()->thetae01();
  vec[160]=(*S)->getIndicesCollectorPointer()->thetae02();
  vec[161]=(*S)->getIndicesCollectorPointer()->THTE_LR03();
  vec[162]=(*S)->getIndicesCollectorPointer()->THTE_LR14();
  vec[163]=(*S)->getIndicesCollectorPointer()->VDCAPE(); 
  vec[164]=(*S)->getIndicesCollectorPointer()->VirtualColdPoolStrength();

// Moisture parameters
  vec[165]=(*S)->getIndicesCollectorPointer()->RH01();
  vec[166]=(*S)->getIndicesCollectorPointer()->RH02();
  vec[167]=(*S)->getIndicesCollectorPointer()->RH14();
  vec[168]=(*S)->getIndicesCollectorPointer()->RH25();
  vec[169]=(*S)->getIndicesCollectorPointer()->RH36();
  vec[170]=(*S)->getIndicesCollectorPointer()->RHMIDDLE(); 
  vec[171]=(*S)->getIndicesCollectorPointer()->RH850500();   
  vec[172]=(*S)->getIndicesCollectorPointer()->PWATER();
  vec[173]=(*S)->getIndicesCollectorPointer()->PWATER_eff();
  vec[174]=(*S)->getIndicesCollectorPointer()->MoistureFlux(); 
  vec[175]=(*S)->getIndicesCollectorPointer()->SR_moisture_flux_MW(); 
  vec[176]=(*S)->getIndicesCollectorPointer()->SR_moisture_flux_eff_MW(); 

  // Mean wind
  vec[177]=(*S)->getIndicesCollectorPointer()->MeanWind500();
  vec[178]=(*S)->getIndicesCollectorPointer()->MeanWind01();
  vec[179]=(*S)->getIndicesCollectorPointer()->MeanWind02();
  vec[180]=(*S)->getIndicesCollectorPointer()->MeanWind03();
  vec[181]=(*S)->getIndicesCollectorPointer()->MeanWind06();
  vec[182]=(*S)->getIndicesCollectorPointer()->MeanWind13();
  vec[183]=(*S)->getIndicesCollectorPointer()->WS_LLmax();
  vec[184]=(*S)->getIndicesCollectorPointer()->WS_MLmax();
  vec[185]=(*S)->getIndicesCollectorPointer()->WS_ULmax();

  // Bulk wind shear
  vec[186]=(*S)->getIndicesCollectorPointer()->BS500();
  vec[187]=(*S)->getIndicesCollectorPointer()->BS01();
  vec[188]=(*S)->getIndicesCollectorPointer()->BS03();
  vec[189]=(*S)->getIndicesCollectorPointer()->BS06();
  vec[190]=(*S)->getIndicesCollectorPointer()->BS08();
  vec[191]=(*S)->getIndicesCollectorPointer()->BS010();
  vec[192]=(*S)->getIndicesCollectorPointer()->BS14();
  vec[193]=(*S)->getIndicesCollectorPointer()->BS16();
  vec[194]=(*S)->getIndicesCollectorPointer()->BS18();
  vec[195]=(*S)->getIndicesCollectorPointer()->BS110();
  vec[196]=(*S)->getIndicesCollectorPointer()->BS_LLmax();
  vec[197]=(*S)->getIndicesCollectorPointer()->BS_MLmax();
  vec[198]=(*S)->getIndicesCollectorPointer()->BS_ULmax();
  vec[199]=(*S)->getIndicesCollectorPointer()->esbbs();
  vec[200]=(*S)->getIndicesCollectorPointer()->emlbs();
  vec[201]=(*S)->getIndicesCollectorPointer()->emubs();
  vec[202]=(*S)->getIndicesCollectorPointer()->emumlbs();
  vec[203]=(*S)->getIndicesCollectorPointer()->emu500bs();
  vec[204]=(*S)->getIndicesCollectorPointer()->BulkShearSfcTen();
  vec[205]=(*S)->getIndicesCollectorPointer()->BulkShear1kmTen();
  vec[206]=(*S)->getIndicesCollectorPointer()->BulkShearMLLCLTen();
  vec[207]=(*S)->getIndicesCollectorPointer()->BulkShearMUMLLCLTen();
  vec[208]=(*S)->getIndicesCollectorPointer()->BS06_var_SI();

  // Storm-relative winds
  vec[209]=(*S)->getIndicesCollectorPointer()->MeanSR500_RM();
  vec[210]=(*S)->getIndicesCollectorPointer()->MeanSR500_LM();
  vec[211]=(*S)->getIndicesCollectorPointer()->MeanSR500_MW();
  vec[212]=(*S)->getIndicesCollectorPointer()->Peters_SR_inflow();
  vec[213]=(*S)->getIndicesCollectorPointer()->MeanSR01_RM();
  vec[214]=(*S)->getIndicesCollectorPointer()->MeanSR01_LM();
  vec[215]=(*S)->getIndicesCollectorPointer()->MeanSR01_MW();
  vec[216]=(*S)->getIndicesCollectorPointer()->MeanSR03_RM();
  vec[217]=(*S)->getIndicesCollectorPointer()->MeanSR03_LM();
  vec[218]=(*S)->getIndicesCollectorPointer()->MeanSR03_MW();
  vec[219]=(*S)->getIndicesCollectorPointer()->MeanSR36_RM();
  vec[220]=(*S)->getIndicesCollectorPointer()->MeanSR36_LM();
  vec[221]=(*S)->getIndicesCollectorPointer()->MeanSR36_MW();
  vec[222]=(*S)->getIndicesCollectorPointer()->MSR_RM_HGL();
  vec[223]=(*S)->getIndicesCollectorPointer()->MSR_LM_HGL();
  vec[224]=(*S)->getIndicesCollectorPointer()->MSR_MW_HGL();
  vec[225]=(*S)->getIndicesCollectorPointer()->MeanSR0500_RM_eff();
  vec[226]=(*S)->getIndicesCollectorPointer()->MeanSR0500_LM_eff();
  vec[227]=(*S)->getIndicesCollectorPointer()->MeanSR0500_MW_eff();
  vec[228]=(*S)->getIndicesCollectorPointer()->Peters_SR_inflow_eff();
  vec[229]=(*S)->getIndicesCollectorPointer()->Ventilation_16km_RM();
  vec[230]=(*S)->getIndicesCollectorPointer()->Ventilation_16km_LM();
  vec[231]=(*S)->getIndicesCollectorPointer()->Ventilation_36km_RM();
  vec[232]=(*S)->getIndicesCollectorPointer()->Ventilation_36km_LM();
  
  // Storm-relative helicity and vorticity
  vec[233]=(*S)->getIndicesCollectorPointer()->SRH100RM();
  vec[234]=(*S)->getIndicesCollectorPointer()->SRH100LM();
  vec[235]=(*S)->getIndicesCollectorPointer()->SRH100RM_F();
  vec[236]=(*S)->getIndicesCollectorPointer()->SRH100LM_F();
  vec[237]=(*S)->getIndicesCollectorPointer()->SRH500RM();
  vec[238]=(*S)->getIndicesCollectorPointer()->SRH500LM();
  vec[239]=(*S)->getIndicesCollectorPointer()->SRH500RM_F();
  vec[240]=(*S)->getIndicesCollectorPointer()->SRH500LM_F();
  vec[241]=(*S)->getIndicesCollectorPointer()->SRH01RM();
  vec[242]=(*S)->getIndicesCollectorPointer()->SRH01LM();
  vec[243]=(*S)->getIndicesCollectorPointer()->SRH03RM();
  vec[244]=(*S)->getIndicesCollectorPointer()->SRH03LM();
  vec[245]=(*S)->getIndicesCollectorPointer()->SRH16RM();
  vec[246]=(*S)->getIndicesCollectorPointer()->SRH16LM();
  vec[247]=(*S)->getIndicesCollectorPointer()->SRH01RM_eff();
  vec[248]=(*S)->getIndicesCollectorPointer()->SRH01LM_eff();
  vec[249]=(*S)->getIndicesCollectorPointer()->SRH03RM_eff();
  vec[250]=(*S)->getIndicesCollectorPointer()->SRH03LM_eff();
  vec[251]=(*S)->getIndicesCollectorPointer()->SW100_RM();
  vec[252]=(*S)->getIndicesCollectorPointer()->SW100_LM();
  vec[253]=(*S)->getIndicesCollectorPointer()->SW100_RM_F();
  vec[254]=(*S)->getIndicesCollectorPointer()->SW100_LM_F();
  vec[255]=(*S)->getIndicesCollectorPointer()->SW500_RM();
  vec[256]=(*S)->getIndicesCollectorPointer()->SW500_LM();
  vec[257]=(*S)->getIndicesCollectorPointer()->SW500_RM_F();
  vec[258]=(*S)->getIndicesCollectorPointer()->SW500_LM_F();
  vec[259]=(*S)->getIndicesCollectorPointer()->SW01_RM();
  vec[260]=(*S)->getIndicesCollectorPointer()->SW01_LM();
  vec[261]=(*S)->getIndicesCollectorPointer()->SW03_RM();  
  vec[262]=(*S)->getIndicesCollectorPointer()->SW03_LM();
  vec[263]=(*S)->getIndicesCollectorPointer()->SV_100_RM_FRA();
  vec[264]=(*S)->getIndicesCollectorPointer()->SV_100_LM_FRA();
  vec[265]=(*S)->getIndicesCollectorPointer()->SV_500_RM_FRA();
  vec[266]=(*S)->getIndicesCollectorPointer()->SV_500_LM_FRA();
  vec[267]=(*S)->getIndicesCollectorPointer()->SV_1000_RM_FRA();
  vec[268]=(*S)->getIndicesCollectorPointer()->SV_1000_LM_FRA();
  vec[269]=(*S)->getIndicesCollectorPointer()->SV_3000_RM_FRA();
  vec[270]=(*S)->getIndicesCollectorPointer()->SV_3000_LM_FRA();
	
  // Vectors
  vec[271]=(*S)->getIndicesCollectorPointer()->Bunkers_RM_A();
  vec[272]=(*S)->getIndicesCollectorPointer()->Bunkers_RM_M();
  vec[273]=(*S)->getIndicesCollectorPointer()->Bunkers_LM_A();
  vec[274]=(*S)->getIndicesCollectorPointer()->Bunkers_LM_M();
  vec[275]=(*S)->getIndicesCollectorPointer()->Bunkers_MW_A();
  vec[276]=(*S)->getIndicesCollectorPointer()->Bunkers_MW_M();
  vec[277]=(*S)->getIndicesCollectorPointer()->Bunkers4_RM_A();
  vec[278]=(*S)->getIndicesCollectorPointer()->Bunkers4_RM_M();
  vec[279]=(*S)->getIndicesCollectorPointer()->Bunkers4_LM_A();
  vec[280]=(*S)->getIndicesCollectorPointer()->Bunkers4_LM_M();
  vec[281]=(*S)->getIndicesCollectorPointer()->Peters_vector_A();
  vec[282]=(*S)->getIndicesCollectorPointer()->Peters_vector_M();
  vec[283]=(*S)->getIndicesCollectorPointer()->Corfidi_downwind_A();
  vec[284]=(*S)->getIndicesCollectorPointer()->Corfidi_downwind_M();
  vec[285]=(*S)->getIndicesCollectorPointer()->Corfidi_upwind_A();
  vec[286]=(*S)->getIndicesCollectorPointer()->Corfidi_upwind_M();

  // Composite metrics
  vec[287]=(*S)->getIndicesCollectorPointer()->K_Index();
  vec[288]=(*S)->getIndicesCollectorPointer()->TotalTotals();  
  vec[289]=(*S)->getIndicesCollectorPointer()->STP();
  vec[290]=(*S)->getIndicesCollectorPointer()->STP_LM();
  vec[291]=(*S)->getIndicesCollectorPointer()->STPeff();
  vec[292]=(*S)->getIndicesCollectorPointer()->STPeff_LM();
  vec[293]=(*S)->getIndicesCollectorPointer()->SCP();
  vec[294]=(*S)->getIndicesCollectorPointer()->SCP_LM();
  vec[295]=(*S)->getIndicesCollectorPointer()->SCPeff();
  vec[296]=(*S)->getIndicesCollectorPointer()->SCPeff_LM();
  vec[297]=(*S)->getIndicesCollectorPointer()->LTTP_RM();
  vec[298]=(*S)->getIndicesCollectorPointer()->LTTP_LM();
  vec[299]=(*S)->getIndicesCollectorPointer()->SHP();
  vec[300]=(*S)->getIndicesCollectorPointer()->HSI();
  vec[301]=(*S)->getIndicesCollectorPointer()->HSIv2();
  vec[302]=(*S)->getIndicesCollectorPointer()->DCP();
  vec[303]=(*S)->getIndicesCollectorPointer()->DCP_eff();
  vec[304]=(*S)->getIndicesCollectorPointer()->EHI500();
  vec[305]=(*S)->getIndicesCollectorPointer()->EHI500_LM();
  vec[306]=(*S)->getIndicesCollectorPointer()->EHI01();
  vec[307]=(*S)->getIndicesCollectorPointer()->EHI01_LM();
  vec[308]=(*S)->getIndicesCollectorPointer()->EHI03();
  vec[309]=(*S)->getIndicesCollectorPointer()->EHI03_LM();
  vec[310]=(*S)->getIndicesCollectorPointer()->SHERBS3();
  vec[311]=(*S)->getIndicesCollectorPointer()->SHERBE();
  vec[312]=(*S)->getIndicesCollectorPointer()->SHERBS3_v2();
  vec[313]=(*S)->getIndicesCollectorPointer()->DEI();
  vec[314]=(*S)->getIndicesCollectorPointer()->DEI_eff();
  return vec;
}

double interpol_lin(double x1,double x2, double y1,double y2, double x){
  double a =(y2-y1)/(x2-x1);
  double b = y1-a*x1;
  return a*x+b;
}

void listToArray(list<double> list, double * arr, int len){
  
  std::list<double>::iterator iter=list.begin();
  
  for(int i =0;i<len;i++){
    arr[i] =*iter;
    if(iter==list.end())break;
    iter++;
  }
}

int interpolate(double **pu, double **hu, double **tu, double **du, double **au, double **vu, int n, double *p_arr=0, double *h_arr=0, int m=0, int o=0){
  double defp[]={850,700,500}; int plen = 3;
  double defh[]={0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 15000, 16000, 17000, 18000, 19000, 20000}; int hlen = 68;
  
  double * p = *pu;
  double * h = *hu;
  double * t = *tu;
  double * d = *du;
  double * a = *au;
  double * v = *vu;
  if(p_arr==0) p_arr=defp;
  else plen=m;
  if(h_arr==0)h_arr=defh;
  else hlen=o;
  
  list<double> vp=list<double>();
  list<double> vh=list<double>();
  list<double> vt=list<double>();
  list<double> vd=list<double>();
  list<double> va=list<double>();
  list<double> vv=list<double>();
  
  int i=0;
  int j=0;
  for(i=0;i<n-1;i++){
    vp.push_back(p[i]);
    vh.push_back(h[i]);
    vt.push_back(t[i]);
    vd.push_back(d[i]);
    va.push_back(a[i]);
    vv.push_back(v[i]);
    double p1 = p[i];
    double p2= p[i+1];
    
    for(j=0;j<plen;j++){
      double pn = p_arr[j];
      
      if(pn>p2&&pn<p1){
        vp.push_back(pn);
        double h1 = h[i];
        double h2 = h[i+1];
        double hn = interpol_lin(p1,p2,h1,h2,pn);
        vh.push_back(hn);
        double t1 = t[i];
        double t2 = t[i+1];
        double tn = interpol_lin(p1,p2,t1,t2,pn);
        vt.push_back(tn);
        
        double d1 = d[i];
        double d2 = d[i+1];
        double dn = interpol_lin(p1,p2,d1,d2,pn);
        vd.push_back(dn);
        
        double a1 = a[i];
        double v1 = v[i];
        Vector w1=Vector(a1,v1);
        double a2 = a[i+1];
        double v2 = v[i+1];
        Vector w2=Vector(a2,v2);
        
        double x = interpol_lin(p1,p2,w1.X(),w2.X(),pn);
        double y = interpol_lin(p1,p2,w1.Y(),w2.Y(),pn);
        
        Vector w3 = Vector(x,y,0);
        
        double* av =w3.toAV();
        
        va.push_back(av[0]);
        vv.push_back(av[1]);
        
        delete[] av;
      }
    }
    double h1 = h[i];
    double h2= h[i+1];
    
  }
  
  vp.push_back(p[n-1]);
  vh.push_back(h[n-1]);
  vt.push_back(t[n-1]);
  vd.push_back(d[n-1]);
  va.push_back(a[n-1]);
  vv.push_back(v[n-1]);
  
  double *np, *nh, *nt, *nd, *na,*nv;
  int u = vp.size();
  
  np= new double[u];
  nh= new double[u];
  nt= new double[u];
  nd= new double[u];
  na= new double[u];
  nv= new double[u];
  
  listToArray(vp, np, u);
  
  
  listToArray(vh, nh, u);
  listToArray(vt, nt, u);
  listToArray(vd, nd, u);
  listToArray(va, na, u);
  listToArray(vv, nv, u);
  
  delete[](*pu);delete[](*hu);delete[](*tu);delete[](*du);delete[](*au);delete[](*vu);
  
  *pu=np;*hu=nh;*tu=nt;*du=nd;*au=na;*vu=nv;
  
  return u;
  
  
}
int interpolate2(double **pu, double **hu, double **tu, double **du, double **au, double **vu, int n, double *p_arr=0, double *h_arr=0, int m=0, int o=0){
  double defp[]={850,700,500}; int plen = 3;
  double defh[]={0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 15000, 16000, 17000, 18000, 19000, 20000}; int hlen = 68;
  
  double * p = *pu;
  double * h = *hu;
  double * t = *tu;
  double * d = *du;
  double * a = *au;
  double * v = *vu;
  if(p_arr==0) p_arr=defp;
  else plen=m;
  if(h_arr==0)h_arr=defh;
  else hlen=o;
  
  list<double> vp=list<double>();
  list<double> vh=list<double>();
  list<double> vt=list<double>();
  list<double> vd=list<double>();
  list<double> va=list<double>();
  list<double> vv=list<double>();
  
  int i=0;
  int j=0;
  for(i=0;i<n-1;i++){
    vp.push_back(p[i]);
    vh.push_back(h[i]);
    vt.push_back(t[i]);
    vd.push_back(d[i]);
    va.push_back(a[i]);
    vv.push_back(v[i]);
    double p1 = p[i];
    double p2= p[i+1];
    
    double h1 = h[i];
    double h2= h[i+1];
    
    for(j=0;j<hlen;j++){
      double hn = h_arr[j]+h[0];
      if(hn>h1&&hn<h2){
        double p1 = p[i];
        double p2 = p[i+1];
        double pn = interpol_lin(h1,h2,p1,p2,hn);
        vp.push_back(pn);
        vh.push_back(hn);
        
        double t1 = t[i];
        double t2 = t[i+1];
        double tn = interpol_lin(h1,h2,t1,t2,hn);
        vt.push_back(tn);
        
        double d1 = d[i];
        double d2 = d[i+1];
        double dn = interpol_lin(h1,h2,d1,d2,hn);
        vd.push_back(dn);
        
        double a1 = a[i];
        double v1 = v[i];
        Vector w1=Vector(a1,v1);
        
        double a2 = a[i+1];
        double v2 = v[i+1];
        Vector w2=Vector(a2,v2);
        
        double x = interpol_lin(h1,h2,w1.X(),w2.X(),hn);
        double y = interpol_lin(h1,h2,w1.Y(),w2.Y(),hn);
        
        Vector w3 = Vector(x,y,0);
        
        double* av =w3.toAV();
        
        va.push_back(av[0]);
        vv.push_back(av[1]);
        delete[] av;
        
      }
      
      
    }
    
  }
  
  vp.push_back(p[n-1]);
  vh.push_back(h[n-1]);
  vt.push_back(t[n-1]);
  vd.push_back(d[n-1]);
  va.push_back(a[n-1]);
  vv.push_back(v[n-1]);
  
  double *np, *nh, *nt, *nd, *na,*nv;
  int u = vp.size();
  
  np= new double[u];
  nh= new double[u];
  nt= new double[u];
  nd= new double[u];
  na= new double[u];
  nv= new double[u];
  
  listToArray(vp, np, u);
  listToArray(vh, nh, u);
  listToArray(vt, nt, u);
  listToArray(vd, nd, u);
  listToArray(va, na, u);
  listToArray(vv, nv, u);
  
  delete[](*pu);delete[](*hu);delete[](*tu);delete[](*du);delete[](*au);delete[](*vu);
  
  *pu=np;*hu=nh;*tu=nt;*du=nd;*au=na;*vu=nv;
  
  return u;
  
  
}
double * sounding_default2(double* pressure,
                           double* altitude,
                           double* temperature,
                           double* dew,
                           double* angle,
                           double* velocity,
                           int size,
                           Sounding **sret,
                           int custom_vec=1,
                           int interpolate_step=5,
                           double* meanlayer_bottom_top = NULL,
                           Vector storm_motion = Vector(999,999,999))
{
  
  
  double *p = new double[size], *h=new double[size], *d=new double[size],*t=new double[size], *a= new double[size], *v=new double[size];
  for(int i=0; i<size;i++){
    p[i]=pressure[i];
    h[i]=altitude[i];
    t[i]=temperature[i];
    d[i]=dew[i];
    a[i]=angle[i];
    v[i]=velocity[i];
  }
  
  double defh2[]={0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000, 1025, 1050, 1075, 1100, 1125, 1150, 1175, 1200, 1225, 1250, 1275, 1300, 1325, 1350, 1375, 1400, 1425, 1450, 1475, 1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725, 1750, 1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000, 2025, 2050, 2075, 2100, 2125, 2150, 2175, 2200, 2225, 2250, 2275, 2300, 2325, 2350, 2375, 2400, 2425, 2450, 2475, 2500, 2525, 2550, 2575, 2600, 2625, 2650, 2675, 2700, 2725, 2750, 2775, 2800, 2825, 2850, 2875, 2900, 2925, 2950, 2975, 3000, 3050, 3100, 3150, 3200, 3250, 3300, 3350, 3400, 3450, 3500, 3550, 3600, 3650, 3700, 3750, 3800, 3850, 3900, 3950, 4000, 4050, 4100, 4150, 4200, 4250, 4300, 4350, 4400, 4450, 4500, 4550, 4600, 4650, 4700, 4750, 4800, 4850, 4900, 4950, 5000, 5050, 5100, 5150, 5200, 5250, 5300, 5350, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000, 7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900, 8000, 8100, 8200, 8300, 8400, 8500, 8600, 8700, 8800, 8900, 9000, 9100, 9200, 9300, 9400, 9500, 9600, 9700, 9800, 9900, 10000, 10100, 10200, 10300, 10400, 10500, 10600, 10700, 10800, 10900, 11000, 11100, 11200, 11300, 11400, 11500, 11600, 11700, 11800, 11900, 12000, 12250, 12500, 12750, 13000, 13250, 13500, 13750, 14000, 14250, 14500, 14750, 15000, 15250, 15500, 15750, 16000, 16250, 16500, 16750, 17000, 17250, 17500, 17750, 18000, 18250, 18500, 18750, 19000, 19250, 19500, 19750, 20000}; int hlen = 318;
  int step = 1000000;
  double *pp=0;double *hh=0; int ch=0;
  int u = interpolate(&p, &h,&t, &d, &a, &v, size,pp, hh, 0, ch);
  if(custom_vec==2){
    hh=defh2;
    ch=hlen;
  }else if(custom_vec==3){
    
    int liczba_poziomow = 20000/interpolate_step;
    double *tmp = new double[liczba_poziomow];
    tmp[0]=0;
    for(int g=1;g<liczba_poziomow;g++){
      tmp[g]=tmp[g-1]+interpolate_step;
    }
    hh=tmp;
    ch=liczba_poziomow;
  }
  u = interpolate2(&p, &h,&t, &d, &a, &v, u,pp, hh, 0, ch);
  if(custom_vec==3)delete[] hh;
  
  double *result = processSounding(p,h,t,d,a,v,u,step,sret,meanlayer_bottom_top,storm_motion);
  
  delete[]p;delete[]h;delete[]t;delete[]d;delete[]a;delete[]v;
  return result;
  
}

//' R call to C++ function for calculating thermo- and kinematic indices derived from atmospheric profiling.
//'
//' More details in the sounding_compute() function
//' 
//' @param pressure pressure [hPa]
//' @param altitude altitude [meters]
//' @param temp temperature [degree Celsius]
//' @param dpt dew point temperature [degree Celsius]
//' @param wd wind direction [azimuth in degrees]
//' @param ws wind speed [knots]
//' @param export_profile possibility to export interpolated profile on the levels defined in accuracy configuration
//' @param accuracy accuracy of computations where 3 = high (slow), 2 = medium (recommended), 1 = low (fast)
//' @param interpolate_step interpolation step to be used for vertical interpolation. Valid only if `accuracy` is set to 3 (default is 5 m)
//' @param meanlayer_bottom_top (optional) vector of length 2 for bottom and top heights used for computing parcel starting parameters; default: 0, 500
//' @param storm_motion (optional) for moving storms only - one can define vector of length two with
//' wind speed (m/s) and wind directions (degrees) that will be used to compute adjusted SRH parameters
//' @examples 
//' pressure = c(1000, 855, 700, 500, 300, 100, 10) 
//' altitude = c(0, 1500, 2500, 6000, 8500, 12000, 25000)
//' temp = c(25, 10, 0, -15, -30, -50, -92)
//' dpt = c(20, 5, -5, -30, -55, -80, -99)
//' wd = c(0, 90, 135, 180, 270, 350, 0)
//' ws = c(5, 10, 20, 30, 40, 5, 0)
//' sounding_default(pressure, altitude, temp, dpt, wd, ws,
//'                  accuracy = 3,
//'                  export_profile = 0,
//'                  interpolate_step = 5,
//'                  storm_motion = c(999, 999),
//'                  meanlayer_bottom_top = c(0, 500))
//' @useDynLib thunder
//' @importFrom Rcpp evalCpp
//' @export
//' @return 
//' \enumerate{
//'  \item 	SB_CAPE
//'  \item 	SB_CAPE_3km
//'  \item 	SB_CAPE_HGL
//'  \item 	SB_buoy
//'  \item 	SB_buoy_HGL
//'  \item 	SB_buoy_3km
//'  \item 	SB_LI
//'  \item 	SB_LI_M25
//'  \item 	SB_CIN
//'  \item 	SB_CIN_4km
//'  \item 	SB_LCL_hgt
//'  \item 	SB_LFC_hgt
//'  \item 	SB_EL_hgt
//'  \item 	SB_LCL_tmp
//'  \item 	SB_LFC_tmp
//'  \item 	SB_EL_tmp
//'  \item 	SB_cold_cloud
//'  \item 	SB_warm_cloud
//'  \item 	SB_equal_layer
//'  \item 	SB_MIXR
//'  \item 	SB_WMAXSHEAR
//'  \item 	SB_E_CAPE
//'  \item 	SB_E_CAPE_3km
//'  \item 	SB_E_CAPE_HGL
//'  \item 	SB_E_buoy
//'  \item 	SB_E_buoy_HGL
//'  \item 	SB_E_buoy_3km
//'  \item 	SB_E_LI
//'  \item 	SB_E_WMAXSHEAR
//'  \item 	SB_E_WMAXSHEAR_HGL
//'  \item 	SB_E_WMAXSHEAR_3km
//'  \item 	ML_CAPE
//'  \item 	ML_CAPE_3km
//'  \item 	ML_CAPE_HGL
//'  \item 	ML_buoy
//'  \item 	ML_buoy_HGL
//'  \item 	ML_buoy_3km
//'  \item 	ML_LI
//'  \item 	ML_LI_M25
//'  \item 	ML_CIN
//'  \item 	ML_CIN_4km
//'  \item 	ML_LCL_hgt
//'  \item 	ML_LFC_hgt
//'  \item 	ML_EL_hgt
//'  \item 	ML_LCL_tmp
//'  \item 	ML_LFC_tmp
//'  \item 	ML_EL_tmp
//'  \item 	ML_cold_cloud
//'  \item 	ML_warm_cloud
//'  \item 	ML_equal_layer
//'  \item 	ML_MIXR
//'  \item 	ML_WMAXSHEAR
//'  \item 	ML_E_CAPE
//'  \item 	ML_E_CAPE_3km
//'  \item 	ML_E_CAPE_HGL
//'  \item 	ML_E_buoy
//'  \item 	ML_E_buoy_HGL
//'  \item 	ML_E_buoy_3km
//'  \item 	ML_E_LI
//'  \item 	ML_E_WMAXSHEAR
//'  \item 	ML_E_WMAXSHEAR_HGL
//'  \item 	ML_E_WMAXSHEAR_3km
//'  \item 	MU_CAPE
//'  \item 	MU_CAPE_3km
//'  \item 	MU_CAPE_HGL
//'  \item 	MU_buoy
//'  \item 	MU_buoy_HGL
//'  \item 	MU_buoy_3km
//'  \item 	MU_LI
//'  \item 	MU_LI_M25
//'  \item 	MU_CIN
//'  \item 	MU_CIN_4km
//'  \item 	MU_LCL_hgt
//'  \item 	MU_LFC_hgt
//'  \item 	MU_EL_hgt
//'  \item 	MU_LCL_tmp
//'  \item 	MU_LFC_tmp
//'  \item 	MU_EL_tmp
//'  \item 	MU_cold_cloud
//'  \item 	MU_warm_cloud
//'  \item 	MU_equal_layer
//'  \item 	MU_MIXR
//'  \item 	MU_WMAXSHEAR
//'  \item 	MU_E_CAPE
//'  \item 	MU_E_CAPE_3km
//'  \item 	MU_E_CAPE_HGL
//'  \item 	MU_E_buoy
//'  \item 	MU_E_buoy_HGL
//'  \item 	MU_E_buoy_3km
//'  \item 	MU_E_LI
//'  \item 	MU_E_WMAXSHEAR
//'  \item 	MU_E_WMAXSHEAR_HGL
//'  \item 	MU_E_WMAXSHEAR_3km
//'  \item 	MUML_CAPE
//'  \item 	MUML_CAPE_3km
//'  \item 	MUML_CAPE_HGL
//'  \item 	MUML_buoy
//'  \item 	MUML_buoy_HGL
//'  \item 	MUML_buoy_3km
//'  \item 	MUML_LI
//'  \item 	MUML_LI_M25
//'  \item 	MUML_CIN
//'  \item 	MUML_CIN_4km
//'  \item 	MUML_LCL_hgt
//'  \item 	MUML_LFC_hgt
//'  \item 	MUML_EL_hgt
//'  \item 	MUML_LCL_tmp
//'  \item 	MUML_LFC_tmp
//'  \item 	MUML_EL_tmp
//'  \item 	MUML_cold_cloud
//'  \item 	MUML_warm_cloud
//'  \item 	MUML_equal_layer
//'  \item 	MUML_MIXR
//'  \item 	MUML_WMAXSHEAR
//'  \item 	MUML_E_CAPE
//'  \item 	MUML_E_CAPE_3km
//'  \item 	MUML_E_CAPE_HGL
//'  \item 	MUML_E_buoy
//'  \item 	MUML_E_buoy_HGL
//'  \item 	MUML_E_buoy_3km
//'  \item 	MUML_E_LI
//'  \item 	MUML_E_WMAXSHEAR
//'  \item 	MUML_E_WMAXSHEAR_HGL
//'  \item 	MUML_E_WMAXSHEAR_3km
//'  \item 	MU5_CAPE
//'  \item 	MU5_CAPE_M10
//'  \item 	MU5_CAPE_HGL
//'  \item 	MU5_buoy
//'  \item 	MU5_buoy_HGL
//'  \item 	MU5_LI
//'  \item 	MU5_LI_M25
//'  \item 	MU5_CIN
//'  \item 	MU5_CIN_4km
//'  \item 	MU5_E_CAPE
//'  \item 	MU5_E_CAPE_HGL
//'  \item 	MU5_E_buoy
//'  \item 	MU5_E_buoy_HGL
//'  \item 	MU5_E_LI
//'  \item 	LR_0500m
//'  \item 	LR_01km
//'  \item 	LR_03km
//'  \item 	LR_04km
//'  \item 	LR_06km
//'  \item 	LR_16km
//'  \item 	LR_24km
//'  \item 	LR_26km
//'  \item 	LR_36km
//'  \item 	LR_26km_max
//'  \item 	LR_500700
//'  \item 	LR_500800
//'  \item 	HGT_ISO_0
//'  \item 	HGT_ISO_0_wetbulb
//'  \item 	HGT_ISO_M10
//'  \item 	HGT_MU
//'  \item 	HGT_MUML
//'  \item 	THETAE_delta
//'  \item 	THETAE_delta_4km
//'  \item 	THETAE_LCL_M10
//'  \item 	THETAE_MU_M10
//'  \item 	THETAE_01km
//'  \item 	THETAE_02km
//'  \item 	THETAE_LR_03km
//'  \item 	THETAE_LR_14km
//'  \item 	DCAPE
//'  \item 	Cold_Pool_Strength
//'  \item 	RH_01km
//'  \item 	RH_02km
//'  \item 	RH_14km
//'  \item 	RH_25km
//'  \item 	RH_36km
//'  \item 	RH_HGL
//'  \item 	RH_500850
//'  \item 	PRCP_WATER
//'  \item 	PRCP_WATER_eff
//'  \item 	Moisture_Flux_0500m
//'  \item 	Moisture_Flux_SR
//'  \item 	Moisture_Flux_SR_eff
//'  \item 	MW_0500m
//'  \item 	MW_01km
//'  \item 	MW_02km
//'  \item 	MW_03km
//'  \item 	MW_06km
//'  \item 	MW_13km
//'  \item 	WS_LLmax
//'  \item 	WS_MLmax
//'  \item 	WS_ULmax
//'  \item 	BS_0500m
//'  \item 	BS_01km
//'  \item 	BS_03km
//'  \item 	BS_06km
//'  \item 	BS_08km
//'  \item 	BS_010km
//'  \item 	BS_14km
//'  \item 	BS_16km
//'  \item 	BS_18km
//'  \item 	BS_110km
//'  \item 	BS_LLmax
//'  \item 	BS_MLmax
//'  \item 	BS_ULmax
//'  \item 	BS_eff_SB
//'  \item 	BS_eff_ML
//'  \item 	BS_eff_MU
//'  \item 	BS_eff_MUML
//'  \item 	BS_eff_MU5
//'  \item 	BS_0km_M10
//'  \item 	BS_1km_M10
//'  \item 	BS_ML_LCL_M10
//'  \item 	BS_MUML_LCL_M10
//'  \item 	BS_06km_smoothness
//'  \item 	SRW_0500m_RM
//'  \item 	SRW_0500m_LM
//'  \item 	SRW_0500m_MW
//'  \item 	SRW_0500m_CBV
//'  \item 	SRW_01km_RM
//'  \item 	SRW_01km_LM
//'  \item 	SRW_01km_MW
//'  \item 	SRW_03km_RM
//'  \item 	SRW_03km_LM
//'  \item 	SRW_03km_MW
//'  \item 	SRW_36km_RM
//'  \item 	SRW_36km_LM
//'  \item 	SRW_36km_MW
//'  \item 	SRW_HGL_RM
//'  \item 	SRW_HGL_LM
//'  \item 	SRW_HGL_MW
//'  \item 	SRW_eff_RM
//'  \item 	SRW_eff_LM
//'  \item 	SRW_eff_MW
//'  \item 	SRW_eff_CBV
//'  \item 	Ventilation_16km_RM
//'  \item 	Ventilation_16km_LM
//'  \item 	Ventilation_36km_RM
//'  \item 	Ventilation_36km_LM
//'  \item 	SRH_0100m_RM
//'  \item 	SRH_0100m_LM
//'  \item 	SRH_0100m_RM_G
//'  \item 	SRH_0100m_LM_G
//'  \item 	SRH_0500m_RM
//'  \item 	SRH_0500m_LM
//'  \item 	SRH_0500m_RM_G
//'  \item 	SRH_0500m_LM_G
//'  \item 	SRH_01km_RM
//'  \item 	SRH_01km_LM
//'  \item 	SRH_03km_RM
//'  \item 	SRH_03km_LM
//'  \item 	SRH_16km_RM
//'  \item 	SRH_16km_LM
//'  \item 	SRH_eff_1km_RM
//'  \item 	SRH_eff_1km_LM
//'  \item 	SRH_eff_3km_RM
//'  \item 	SRH_eff_3km_LM
//'  \item 	SV_0100m_RM
//'  \item 	SV_0100m_LM
//'  \item 	SV_0100m_RM_G
//'  \item 	SV_0100m_LM_G
//'  \item 	SV_0500m_RM
//'  \item 	SV_0500m_LM
//'  \item 	SV_0500m_RM_G
//'  \item 	SV_0500m_LM_G
//'  \item 	SV_01km_RM
//'  \item 	SV_01km_LM
//'  \item 	SV_03km_RM
//'  \item 	SV_03km_LM
//'  \item 	SV_0100m_RM_fra
//'  \item 	SV_0100m_LM_fra
//'  \item 	SV_0500m_RM_fra
//'  \item 	SV_0500m_LM_fra
//'  \item 	SV_01km_RM_fra
//'  \item 	SV_01km_LM_fra
//'  \item 	SV_03km_RM_fra
//'  \item 	SV_03km_LM_fra
//'  \item 	Bunkers_RM_A
//'  \item 	Bunkers_RM_M
//'  \item 	Bunkers_LM_A
//'  \item 	Bunkers_LM_M
//'  \item 	Bunkers_MW_A
//'  \item 	Bunkers_MW_M
//'  \item 	Bunkers_4_RM_A
//'  \item 	Bunkers_4_RM_M
//'  \item 	Bunkers_4_LM_A
//'  \item 	Bunkers_4_LM_M
//'  \item 	CBV_A
//'  \item 	CBV_M
//'  \item 	Corfidi_downwind_A
//'  \item 	Corfidi_downwind_M
//'  \item 	Corfidi_upwind_A
//'  \item 	Corfidi_upwind_M
//'  \item 	K_Index
//'  \item 	TotalTotals_Index
//'  \item 	STP_fix_RM
//'  \item 	STP_fix_LM
//'  \item 	STP_eff_RM
//'  \item 	STP_eff_LM
//'  \item 	SCP_fix_RM
//'  \item 	SCP_fix_LM
//'  \item 	SCP_eff_RM
//'  \item 	SCP_eff_LM
//'  \item 	LTTP_RM();
//'  \item 	LTTP_LM();
//'  \item 	SHIP
//'  \item 	HSI
//'  \item 	HSI_mod
//'  \item 	DCP
//'  \item 	DCP_eff
//'  \item 	EHI_0500m_RM
//'  \item 	EHI_0500m_LM
//'  \item 	EHI_01km_RM
//'  \item 	EHI_01km_LM
//'  \item 	EHI_03km_RM
//'  \item 	EHI_03km_LM
//'  \item 	SHERBS3
//'  \item 	SHERBE
//'  \item 	SHERB_mod
//'  \item 	DEI
//'  \item 	DEI_eff
//' }
 // [[Rcpp::export]]
 
 Rcpp::NumericVector sounding_default(Rcpp::NumericVector pressure,
                                      Rcpp::NumericVector altitude,
                                      Rcpp::NumericVector temp,
                                      Rcpp::NumericVector dpt,
                                      Rcpp::NumericVector wd,
                                      Rcpp::NumericVector ws,
                                      Rcpp::NumericVector export_profile,
                                      Rcpp::NumericVector accuracy,
                                      int interpolate_step,
                                      Rcpp::NumericVector meanlayer_bottom_top,
                                      Rcpp::NumericVector storm_motion)
 {
   Sounding *sret;
   int size = pressure.size();
   double * mlp = new double [2];
   mlp[0] = meanlayer_bottom_top[0];
   mlp[1] = meanlayer_bottom_top[1];
   Vector sm = Vector(storm_motion[0],storm_motion[1],storm_motion[2]);
   
   double *p = new double[size], *h=new double[size], *d=new double[size],*t=new double[size], *a= new double[size], *v=new double[size];
   for(int i=0; i<size;i++){
     p[i]=pressure[i];
     h[i]=altitude[i];
     t[i]=temp[i];
     d[i]=dpt[i];
     a[i]=wd[i];
     v[i]=ws[i];
   }
   int q = accuracy[0];
   int plen,hlen,tlen,dlen,alen,vlen,tvlen;
   int mulen,sblen,mllen,dnlen,mustart,mlstart;
   
   double *result = sounding_default2(p,h,t,d,a,v,size,&sret,q, interpolate_step, mlp, sm);
   int reslen= 315;
   int maxl=reslen;
   if(export_profile[0]==1){
     plen = sret->p->size();
     hlen = sret->h->size();
     tlen = sret->t->size();
     dlen = sret->d->size();
     alen = sret->a->size();
     vlen = sret->v->size();
     tvlen = sret->th->virt->size();
     
     mulen = sret->th->mostUnstable->getVirtualValues()->size();
     sblen = sret->th->surfaceBased->getVirtualValues()->size();
     mllen = sret->th->meanLayer->getVirtualValues()->size();
     dnlen = sret->th->downdraft->getVirtualValues()->size();
     mustart= sret->th->mostUnstable->startIndex;
     mlstart= sret->th->meanLayer->startIndex;
     maxl+=2+mulen+1+sblen+2+mllen+plen+1+hlen+1+tlen+1+dlen+1+alen+1+vlen+1+tvlen+1+dnlen+10;
   }
   
   Rcpp::NumericVector out(maxl);
   
   for(int i = 0; i < reslen; ++i) {
     out[i] = result[i];
   }
   
   if(export_profile[0]==1){
     int i= reslen;
     out[i] = mulen;i++;
     out[i] = mustart;i++;
     
     for (std::list<double>::iterator it = sret->th->mostUnstable->getVirtualValues()->begin(); it != sret->th->mostUnstable->getVirtualValues()->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     out[i]=sblen;i++;
     
     for (std::list<double>::iterator it = sret->th->surfaceBased->getVirtualValues()->begin(); it != sret->th->surfaceBased->getVirtualValues()->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     
     out[i]=mllen; i++;
     out[i] = mlstart;i++;
     
     for (std::list<double>::iterator it = sret->th->meanLayer->getVirtualValues()->begin(); it != sret->th->meanLayer->getVirtualValues()->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     
     out[i]=plen;i++;
     for (std::list<double>::iterator it = sret->p->begin(); it != sret->p->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     
     out[i]=hlen;i++;
     for (std::list<double>::iterator it = sret->h->begin(); it != sret->h->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     out[i]=tlen;i++;
     for (std::list<double>::iterator it = sret->t->begin(); it != sret->t->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     out[i]=dlen;i++;
     for (std::list<double>::iterator it = sret->d->begin(); it != sret->d->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     out[i]=alen;i++;
     for (std::list<double>::iterator it = sret->a->begin(); it != sret->a->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     out[i]=vlen;i++;
     
     for (std::list<double>::iterator it = sret->v->begin(); it != sret->v->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     
     out[i]=tvlen;i++;
     
     for (std::list<double>::iterator it = sret->th->downdraft->getVirtualValues()->begin(); it != sret->th->downdraft->getVirtualValues()->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     out[i]=dnlen;i++;
     
     for (std::list<double>::iterator it = sret->th->virt->begin(); it != sret->th->virt->end(); ++it){
       double temp = 0;
       temp= *it;
       out[i]=temp;
       i++;
     }
     
     out[i]=-999;
   }
   delete sret;
   delete[] result;
   delete[] p; delete[] h;delete[] d; delete[] t; delete[] a; delete[] v;
   
   delete[] mlp;
   
   return out;
   
 }
