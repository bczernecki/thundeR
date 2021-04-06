#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <list>
#include <cstddef>
#include <cstdlib>
#include <algorithm>
using namespace std;

const double kel = 273.15;
const double g = 9.81;
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

double SIGN(double x, double y){
  if (y < 0) return -abs(x);
  else return abs(x);
}

double ESAT(double t){
  double temp = t + kel;
  return (pow(10, 23.832241 - 5.02808 * log10(temp) - 0.00000013816 * pow(10, 11.344 - 0.0303998 * temp) + 0.0081328 * pow(10, 3.49149 - 1302.8844 / temp) - 2949.076 / temp));
}

double W(double t, double p)
{
  return 622 * ESAT(t) / (p - ESAT(t));
}

double O(double t, double p)
{
  return (t + kel) * pow(1000.0 / p, 0.288);
}

double OS(double t, double p)
{
  return (O(t, p) / exp(-2.6518986 * W(t, p) / (t + kel)));
}

double TMR(double W, double p)
{
  if (W != -622)
  {
    double zmienna = log10(W * p / (622 + W));
    return pow(10, 0.0498646455 * zmienna + 2.4082965) - 280.23475 + 38.9114 * (pow(pow(10, 0.0915 * zmienna) - 1.2035, 2));
  }
  else return -kel;
}

double TDA(double O, double p)
{
  return O * pow(p / 1000.0, 0.288) - kel;
}

double TSA(double OS, double p)
{
  double a = OS;
  double b = -2.6518986;
  double Tw = 253.16;
  for (int i = 1; i < 12; i++)
  {
    
    if ((a * exp(b * W(Tw - kel, p) / Tw) - Tw * pow(1000.0 / p, 0.288)) > 0)
      Tw += 120.0 / pow(2.0, (double)i);
    else
      Tw -= 120 / pow(2.0, (double)i);
    
  }
  return Tw - kel + ((pow(1100-p,2.5))/pow(4.18,12))-0.0035;
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
  ~Cache();
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
  delete(this->pindex);
  delete(this->p);
  
  delete(this->hindex);
  delete(this->h);
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
  double n26;
  Vector mean020;
  double n020;
  double n13;
  double n2;
  double n3;
  double n6;
  Vector mean01;
  double n1;
  Vector CorfidiA;
  Vector Corfidi_upwind;
  Vector Corfidi_downwind;
  Vector rm;
  Vector lm;
  double lasth;

  double srh100rm;
  double srh100lm;

  double srh250rm;
  double srh250lm;
  
  double srh500rm;
  double srh500lm;
  
  double srh01lm;
  double srh01rm;
  
  double srh01lmf;
  double srh01rmf;
  
  double srh13lm;
  double srh13rm;
  
  double srh13lmf;
  double srh13rmf;
  
  double srh03lm;
  double srh03rm;
  
  double srh03lmf;
  double srh03rmf;
  
  void putMandatoryVectors(int i, double p, double h, double t, double d, double a, double v, Vector v_);
  void putMeanVectors(int i, double p, double h, double t, double d, double a, double v, Vector v_);
  void putLLJ(int i, double p, double h, double t,double d,double a,double v, Vector v_);
  void putSpecificLine(int i, double p, double h, double t, double d, double a, double v);
  void prepareSupercellVectors();
  void prepareCorfidiVectors();
  void doSRH(int i, double p, double h, double t, double d, double a,double v);
  void finishMeanVectors();
  Vector shear06();
  
public:
  Kinematics();
  virtual ~Kinematics();
  void putSecondPhaseLine(int i, double p, double h, double t, double d, double a, double v)
  {
    doSRH(i, p, h, t, d, a, v);
    lasth = h;
  }
  
  void finishPhase1()
  {
    finishMeanVectors();
    prepareSupercellVectors();
    prepareCorfidiVectors();
    lasth = h0;
    
    srh100rm=0;
    srh100lm=0;
    
    srh250rm=0;
    srh250lm=0;

    srh500rm=0;
    srh500lm=0;
    
    srh01lm=0;
    srh01rm=0;
    
    srh01lmf=0;
    srh01rmf=0;
    
    srh13lm=0;
    srh13rm=0;
    
    srh13lmf=0;
    srh13rmf=0;
    
    srh03lm = 0;
    srh03rm = 0;
    
    srh03lmf = 0;
    srh03rmf = 0;
  }
  
  void finish()
  {
    
  }
};

Kinematics::Kinematics(){
  this->mean01 = Vector(0, 0, 0);
  this->mean06 = Vector(0, 0, 0);
  this->mean02 = Vector(0,0,0);
  this->mean03 = Vector(0,0,0);
  this->mean13 = Vector(0,0,0);
  this->mean0=Vector(0,0,0);
  this->mean6=Vector(0,0,0);
  
  this->mean26 = Vector(0,0,0);
  this->mean020 = Vector(0,0,0);
  
  this->lasth = h0;
  this->vw = new list<Vector>();
  this->llj = Vector(0, 0, 0);
  n2=0;
  nsurf=0;
  nsix=0;
  n6=0;
  n1=0;
  n13=0;
  srh01lm = 0;
  srh01rm = 0;
  
  srh01lmf = 0;
  srh01rmf = 0;
  
  srh13lm = 0;
  srh13rm = 0;
  
  srh13lmf = 0;
  srh13rmf = 0;
  
  srh03lm = 0;
  srh03rm = 0;
  
  srh03lmf = 0;
  srh03rmf = 0;
  n26=0;
  n020=0;
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
  if(t>=-20 && t<=0.0){
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
  if ((size_t)i < vw->size()-1&&h-h0<=3000)
  {
    std::list<Vector>::iterator it = vw->begin();
    std::advance(it, i);
    
    std::list<Vector>::iterator it2 = vw->begin();
    std::advance(it2, i+1);
    
    Vector v1 = *it;
    Vector v2 = *it2;
    
    double tmps1 = (v1.X() - rm.X()) * (v2.Y() - v1.Y()) - (v1.Y() - rm.Y()) * (v2.X() - v1.X());
    double tmps2 = (v1.X() - lm.X()) * (v2.Y() - v1.Y()) - (v1.Y() - lm.Y()) * (v2.X() - v1.X());
    
    srh03rm += tmps1;
    srh03lm += tmps2;
    if (tmps1>0)
      srh03rmf += tmps1;
    else srh03lmf -= tmps1;
    if (tmps2>0) 
      srh03lmf += tmps2;
    else srh03rmf -= tmps2;
    
    if(h-h0<500){
      srh500rm = srh03rm;
      srh500lm = srh03lm;
	}

    if(h-h0<100){
      srh100rm = srh03rm;
      srh100lm = srh03lm;
	}

    if(h-h0<250){
      srh250rm = srh03rm;
      srh250lm = srh03lm;
	}
    
    if (h - h0 <= 1000)
    {
      srh01rm = srh03rm;
      srh01lm = srh03lm;
      srh01rmf = srh03rmf;
      srh01lmf = srh03lmf;
    }
    if (h-h0>=1000.0)
    {
      srh13rm = srh03rm- srh01rm;
      srh13lm = srh03lm- srh01lm;
      srh13rmf = srh03rmf- srh01rmf;
      srh13lmf = srh03lmf-srh01lmf;
    }
  }
}

void Kinematics::finishMeanVectors()
{
  if (n6 != 0) mean06 *= 1.0 / n6;
  else mean06 = Vector(0, 0, 0);
  if (n1 != 0) mean01 *= 1.0 / n1;
  else mean01 = Vector(0, 0, 0);
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
  list<double> *values;
  list<double> *virtualValues;
  
  int lclIndex;
  int vLclIndex;
  
  int lfcIndex;
  int vLfcIndex;
  
  int elIndex;
  int vElIndex;
  
  double cape;
  double cin;
  double to3cape;
  
  double vcape;
  double vcin;
  double vto3cape;
  
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

  bool isSet;
  bool dcape_;
  
  double h0;
  
  double dcape;
  double dvcape;
  
  int i700index;
  
  void allocate();
  
  void free();
  
  void testSpecificLCL(int i,double p,double t,double tmr, double tda, int* lclInd_, int* lfcInd_, list<double>* curve_, double*os_);
  void doRest(int i, double p,double h, double t, double TSA, int* lfcInd_, double* cape_, double* to3, double* cin_, int* elInd_, list<double>* curve_);
  void putClassicLine(int i, double p, double h, double t,double d, double a, double v);
  void putVirtualLine(int i, double p, double h, double t, double d, double a, double v);
  
  
public:
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
  values = new list<double>();
  virtualValues = new list<double>();
  cape = cin = to3cape = vcape = vcin = vto3cape = os = o = w = vos=vo=vw=0;
  middlecape=0;
  coldcape=0;
  lclIndex = vLclIndex = lfcIndex = vLfcIndex = elIndex = vElIndex = -1;
  startIndex=-1;
  isSet = false;
  tcin = 0;
  tvcin = 0;
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
void LapseRate::doRest(int i, double p,double h, double t, double TSA, int* lfcInd_, double* cape_, double* to3, double* cin_, int* elInd_, list<double>* curve_){
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
    this->doRest(i, p, h, t, TSA_, &lfcIndex, &cape, &to3cape, &cin, &elIndex, values);
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
                    if (vt_parcel < t)
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
  cape = cin = to3cape = vcape = vcin = vto3cape = 0;
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
		if(t<=0&&t>=-20)middlecape+=tcap;
			if(t<=-10)coldcape+=tcap;
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
  double meanLayerZHeight;
  double n;
  double mp;
  double mh;
  double mt;
  double md;
  double mmr;
  double t0;
  double pwater;
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

  int min20pos;
  double min20;
  
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
  
  double mthet;
  double mthetn;
  
  list<double>* wbt;
  list<double>* oe;
  list<double>* mixing;
  list<double>* virt;
  
  double meanhum2;
  double meand2b;
	  
  double meanhum25;
  double meand25;
  
  double meanhum36;
  double meand36;

  double meanhumMIDDLE;
  double meandMIDDLE;

  double meanmxr2;
  double meand2;  
  
  LapseRate* mostUnstable;
  LapseRate* mostU500;
  LapseRate* surfaceBased;
  LapseRate* meanLayer;
  LapseRate* downdraft;
  LapseRate* showalter;
  
  void putMinTHTE(int i, double p, double h, double oe);
  void startConditions(int i, double p, double h, double t, double d, double a, double v, double oe);
  void putMaxTHTE(int i, double p, double h, double t, double d, double a, double v, double oe);
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
void Thermodynamics::setMlIndex(int i, double p, double h, double t, double d, double a, double v){
  this->meanLayer->setInitialConditions(i, p, h, t,d, a, v,h0);
  this->meanLayer->setInitialW(mmr, mo);
}
void Thermodynamics::putMlLine(int i, double p, double h, double t, double d, double a, double v){
  this->meanLayer->putLine(i, p, h, t, d, a, v);
}
Thermodynamics::Thermodynamics(){
  mr1000 = 0;
  meanLayerZHeight=500.0;
  n=0;
  mp=0;
  mh=0;
  mt=0;
  md=0;
  mo=0;
  mmr=0;
  t0=0;
  t10=0;
  p10=0;
  pwater=0;
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
  meanhum2=0;
  meand2b=0;

  meanhumMIDDLE=0;
  meandMIDDLE=0;
	
  meanhum25=0;
  meand25=0;
  
  meanhum36=0;
  meand36=0;

  meanmxr2=0;
  meand2=0;
  
  downmr=0;
  downmrn=0;
  downo=0;
  downon=0;
  thetd=0;
  thetn=0;
  
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
  
  meanhum2=ESAT(d)/ESAT(t);
  meand2b=1;

  meanhum25=ESAT(d)/ESAT(t);
  meand25=1;

  meanhum36=ESAT(d)/ESAT(t);
  meand36=1;

  meanhumMIDDLE=ESAT(d)/ESAT(t);
  meandMIDDLE=1;

  meanmxr2=W(d, p);
  meand2=1;
}
void Thermodynamics::putMaxTHTE(int i, double p, double h, double t, double d, double a, double v, double oe)
{
  if (oe > maxOE && h-h0 <= 3000)
  {
    maxOE = oe;
    this->mostUnstable->setInitialConditions(i, p, h, t, d, a, v, h0);
  }
  if(oe>maxOE500&&h-h0<=3000 && h-h0>=500){
	  maxOE500=oe;
	  this->mostU500->setInitialConditions(i, p, h, t, d, a, v, h0);
  }
}
void Thermodynamics::putMeanLayerParameters(int i, double p, double h, double t, double d, double a, double v,double mr)
{
  if (   ((abs(h - h0) < meanLayerZHeight)&&(fmod(abs(h-h0),100.0)==0.0))  ||(h==h0) )
  {
    mh += h;
    mp += p;
    mt += t;
    md += d;
    mmr += mr;
	mo += O(t,p);
    n += 1;
  }
  
  
  if((abs(h - h0) <= 5000)&&(abs(h - h0) >= 3000)&&(fmod(abs(h-h0),200.0)==0.0)){
	downmr+=W(d,p);
	downmrn+=1;
	downo+=O(t,p);
	downon+=1;
	
	thetd+=OE(t,d,p);
	thetn+=1;
  }
  
  if(t>=-20.0&&t<=0.0){
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
  min20pos = i;
  zero = abs(t);
  wb0 = abs(wbt);
  minten = abs(t+10);
  min20 = abs(t+20);
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

  double im20= abs(t+20);

  if(im20<min20){
    min20 = im20;
    min20pos = i;
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
  
  this->wbt->push_back(wbt);
  this->oe->push_back(oe);
  this->mixing->push_back(mr);
  this->virt->push_back(virtt);
  if (i == 0)
  {
  
    startConditions(i, p, h, t, d, a, v, oe);
    ZeroPosStartingConditions(i, p, h, t, d, a, v, wbt);
	
    
  }
  else
  {
    putMinTHTE(i, p, h, oe);
    putMaxTHTE(i, p, h, t, d, a, v, oe);
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
		   	 
                if (abs(h - h0) >= 2000 && abs(h - h0) <= 5000){
	                meanhum25+=ESAT(d)/ESAT(t);
	                meand25+=1.0;
                }

                if (abs(h - h0) >= 3000 && abs(h - h0) <= 6000){
	                meanhum36+=ESAT(d)/ESAT(t);
	                meand36+=1.0;
                }
	  }

		if (t<=0&&t>=-20){
			meanhumMIDDLE+=ESAT(d)/ESAT(t);
			meandMIDDLE+=1;
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
  mp = p0; 
  mt /= n;
  md /= n;
  mmr /= n;
  mo /=n;
  n /= 2;
  downmr/=downmrn;
  downo/=downon;
  thetd/=thetn;
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
  determineDowndraftByMinTHTE(i, p, h, t, d, a, v);
  putShowalter(i, p, h, t, d, a, v);
}
void Thermodynamics::finish(){
  this->mostUnstable->finish();
  this->mostU500->finish();
  this->surfaceBased->finish();
  this->meanLayer->finish();
  this->downdraft->finish();
  this->finishLowLapseRates();
  pwater /= 98.1;
  
  meanhum2/=meand2b;
  meanhum25/=meand25; 
  meanhum36/=meand36; 
  meanhumMIDDLE/=meandMIDDLE; 
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
public:  Thermodynamics *th;
  Cache* cache;
  Kinematics *ks;
    list<double>* p;
  list<double>* h;
  list<double>* t;
  list<double>* d;
  list<double>* a;
  list<double>* v;
  IndicesCollector *ic;
  Sounding(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_, int length, double dz=10);
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
  double ZeroHeight();
  double WetBulbZeroHeight();		
  double MUHeight();
  double MinTHTEHeight();
  double DeltaThetaE();
  double DeltaThetaE_HGL();
  double VDCAPE();	
  double VirtualColdPoolStrength();
  double WindIndex();
  double PWATER();
  double MoistureFlux();
  double RH02();
  double RH25();
  double RHMIDDLE();		
  double BS01();
  double BS02();
  double BS03();
  double BS06();
  double BS08();
  double BS36();
  double BS18();
  double emubs();
  double esbbs();
  double emlbs();
  double BulkShearSfcTen();
  double BulkShearMULFCTen();
  double BulkShearSBLFCTen();
  double BulkShearMLLFCTen();
  double MeanWind01();
  double MeanWind02();
  double MeanWind06();
  double MeanWind13();
  double SRH100RM();
  double SRH250RM();
  double SRH500RM();
  double SRH01RM();
  double SRH03RM();
  double SRH100LM();
  double SRH250LM();
  double SRH500LM();
  double SRH01LM();
  double SRH03LM();
  double Bunkers_RM_A();
  double Bunkers_RM_M();
  double Bunkers_LM_A();
  double Bunkers_LM_M();
  double Bunkers_MW_A();
  double Bunkers_MW_M();
  double K_Index();
  double Showalter_Index();	
  double TotalTotals();		
  double SWEATIndex();	
  double STP();
  double STPeff();
  double SCP();
  double SCPeff();
  double SHP();
  double DCP();
  double MU_WMAXSHEAR();
  double SB_WMAXSHEAR();
  double ML_WMAXSHEAR();
  double MU_EFF_WMAXSHEAR();
  double SB_EFF_WMAXSHEAR();
  double ML_EFF_WMAXSHEAR();

  double RH36();
  double BulkShearSfc20();
  double BulkShearSfczero();
  double BulkShear1kmzero();
  double BulkShear1km20();
  double BulkShear1kmTen();
  double BS16();
  double BulkShearMULFCzero();
  double BulkShearMULFC20();

  double BulkShear2kmzero();
  double BulkShear2km20();
  double BulkShear2kmTen();
  double BS26();
  double lapseRate600800();

  double LR0500();
  double LR02();
  double LR04();
  double LR06();
  double MeanWind03();

  double wind_3km();
  double wind_2km();
  double wind_1km();
  double wind_500m();

  double Corfidi_downwind_A();
  double Corfidi_downwind_M();
  double Corfidi_upwind_A();
  double Corfidi_upwind_M();

  double MUcoldcape_fraction();
  double SBcoldcape_fraction();
  double MLcoldcape_fraction();

  double ML_coldcape_fraction();
  double SB_coldcape_fraction();
  double MU_coldcape_fraction();
  double HSI();
  
  double MSR_MW();
  double MSR_RM();
  double MSR_LM();
  double MSR_MW_HGL();
  double MSR_RM_HGL();
  double MSR_LM_HGL();
  
  double MU500_CAPE();
  double MU500_CIN();
  double MU500_LI();
  
  double tilting();
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
Sounding::Sounding(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_, int length, double dz){
  this->alloc();
  bool valid = this->checkarguments(p_,h_,t_,d_,a_,v_);
  if(valid){
    
    int i=0;
    this->cache->setH0(h_[0]);
    for (i=0; i<length-1; i++){
      this->insertLine(p_,h_,t_,d_,a_,v_,i,dz);
    }
    this->insertSingleLine(p_[length-1],h_[length-1],t_[length-1],d_[length-1],Vector(a_[length-1],v_[length-1]*0.514444));
    
    
    this->th->prepareMeanLayer();
    this->ks->finishPhase1();
    
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
  list<double>::iterator ip;
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
  for(ip = this->p->begin(); ip!=this->p->end(); ++ip){
    double p_ = *ip;
    double h_ = *ih;
	if(h_-h0>=4000)break;
    double t_ = *it;
    double d_ = *id;
    double a_ = *ia;
    double v_ = *iv;
    this->th->downdraft->putLine(i, p_, h_, t_, d_, a_, v_);
    ++ih;++it;++id;++ia;++iv;++i;
  }
  
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
double IndicesCollector::VMostUnstableCAPE()
{
  
  double result = 0;
  result = S->th->mostUnstable->vcape;
  
  return result;
  
}

double IndicesCollector::VLLMostUnstableCAPE()
{
  
  double result = 0;
  result = S->th->mostUnstable->vto3cape;
  
  return result;
  
}

double IndicesCollector::VMostUnstableCIN()
{
  
  double result = 0;
  result = S->th->mostUnstable->vcin;
  
  return result;
  
}   
     
double IndicesCollector::VMostUnstableLCL()
{
  
  double result = 0;
  int index = S->th->mostUnstable->vLclIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
  return result;
  
}

double IndicesCollector::VMostUnstableLFC()
{
  
  double result = 0;
  int index = S->th->mostUnstable->vLfcIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
  return result;
  
}

double IndicesCollector::VMostUnstableEL()
{
  
  double result = 0;
  int index = S->th->mostUnstable->vElIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
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

double IndicesCollector::VMostUnstableVmax(){
  return sqrt(this->VMostUnstableCAPE()*2);
  
}

double IndicesCollector::MUELTemperature(){
  return Get(S->t,S->th->mostUnstable->vElIndex);
}

double IndicesCollector::MULCLTemperature(){
  return Get(S->t,S->th->mostUnstable->vLclIndex);
}

double IndicesCollector::MULFCTemperature(){
  return Get(S->t,S->th->mostUnstable->vLfcIndex);
}

double IndicesCollector::VSurfaceBasedCAPE()
{
  
  double result = 0;
  result = S->th->surfaceBased->vcape;
  
  return result;
  
}

double IndicesCollector::VLLSurfaceBasedCAPE()
{
  
  double result = 0;
  result = S->th->surfaceBased->vto3cape;
  
  return result;
  
}

double IndicesCollector::VSurfaceBasedCIN()
{
  
  double result = 0;
  result = S->th->surfaceBased->vcin;
  
  return result;
  
}   
     
double IndicesCollector::VSurfaceBasedLCL()
{
  
  double result = 0;
  int index = S->th->surfaceBased->vLclIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
  return result;
  
}

double IndicesCollector::VSurfaceBasedLFC()
{
  
  double result = 0;
  int index = S->th->surfaceBased->vLfcIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
  return result;
  
}

double IndicesCollector::VSurfaceBasedEL()
{
  
  double result = 0;
  int index = S->th->surfaceBased->vElIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
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
  return Get(S->t,S->th->surfaceBased->vElIndex);
}

double IndicesCollector::SBLCLTemperature(){
  return Get(S->t,S->th->surfaceBased->vLclIndex);
}

double IndicesCollector::SBLFCTemperature(){
  return Get(S->t,S->th->surfaceBased->vLfcIndex);
}

double IndicesCollector::VMeanLayerCAPE()
{
  
  double result = 0;
  result = S->th->meanLayer->vcape;
  
  return result;
  
}

double IndicesCollector::VLLMeanLayerCAPE()
{
  
  double result = 0;
  result = S->th->meanLayer->vto3cape;
  
  return result;
  
}

double IndicesCollector::VMeanLayerCIN()
{
  
  double result = 0;
  result = S->th->meanLayer->vcin;
  
  return result;
  
}        

double IndicesCollector::VMeanLayerLCL()
{
  
  double result = 0;
  int index = S->th->meanLayer->vLclIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
  return result;
  
}

double IndicesCollector::VMeanLayerLFC()
{
  
  double result = 0;
  int index = S->th->meanLayer->vLfcIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
  return result;
  
}

double IndicesCollector::VMeanLayerEL()
{
  
  double result = 0;
  int index = S->th->meanLayer->vElIndex;
  
  
  result = Get(S->h,index)- S->th->h0;
  
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

double IndicesCollector::VMeanLayerVmax(){
  return sqrt(this->VMeanLayerCAPE()*2);
  
}

double IndicesCollector::MLELTemperature(){
  return Get(S->t,S->th->meanLayer->vElIndex);
}

double IndicesCollector::MLLCLTemperature(){
  return Get(S->t,S->th->meanLayer->vLclIndex);
}

double IndicesCollector::MLLFCTemperature(){
  return Get(S->t,S->th->meanLayer->vLfcIndex);
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

double IndicesCollector::MinTHTEHeight(){
  int zeroIndex = S->th->minTHTEpos;
  double h0 = Get(S->h,0);
  double hz = Get(S->h,zeroIndex);
  return hz-h0;
}

double IndicesCollector::DeltaThetaE(){
  return Get(S->th->oe,0)-S->th->thetd;
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

double IndicesCollector::MLMixingRatio(){
  return S->th->mmr;
}

double IndicesCollector::MUMRatio(){
  return Get(S->th->mixing,S->th->mostUnstable->startIndex);
}

double IndicesCollector::SBMRatio(){
  return Get(S->th->mixing,S->th->surfaceBased->startIndex);
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

double IndicesCollector::BS26(){
  int tail=cache->getHeightIndex(2000);
  int head = cache->getHeightIndex(6000);
  
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

double IndicesCollector::SRH250RM(){
  return S->ks->srh250rm;
}

double IndicesCollector::SRH500RM(){
  return S->ks->srh500rm;
}

double IndicesCollector::SRH01RM(){
  return S->ks->srh01rm;
}

double IndicesCollector::SRH03RM(){
  return S->ks->srh03rm;
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

double IndicesCollector::SBmiddlecape(){
  return S->th->surfaceBased->middlecape;
}

double IndicesCollector::MLmiddlecape(){
  return S->th->meanLayer->middlecape;
}


double IndicesCollector::STP(){
  double sbcape = this->VSurfaceBasedCAPE()/1500;
  double sblcl = this->VSurfaceBasedLCL();
  double srh1 = this->SRH01RM()/150;
  double bwd = this->BS06();
  double cin = this->VSurfaceBasedCIN();	
  
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
  double bwd = this->emubs();
  double cin = this->VMeanLayerCIN();	
  
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
  double srh = this->SRH03RM()/50;
  double ewd = this->emubs();
  double cin = this->VMostUnstableCIN();

  if(cin>-40)cin=1;
  else cin=-40/cin;
	
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


double IndicesCollector::MU_WMAXSHEAR(){
  return this->VMostUnstableVmax()*this->BS06();
}

double IndicesCollector::SB_WMAXSHEAR(){
  return this->VSurfaceBasedVmax()*this->BS06();
}

double IndicesCollector::ML_WMAXSHEAR(){
  return this->VMeanLayerVmax()*this->BS06();
}

double IndicesCollector::MU_EFF_WMAXSHEAR(){
  return this->VMostUnstableVmax()*this->emubs();
}

double IndicesCollector::SB_EFF_WMAXSHEAR(){
  return this->VSurfaceBasedVmax()*this->esbbs();
}

double IndicesCollector::ML_EFF_WMAXSHEAR(){
  return this->VMeanLayerVmax()*this->emlbs();
}


double IndicesCollector::DeltaThetaE_HGL(){
  return Get(S->th->oe,0)-S->th->mthet;
}

double IndicesCollector::BulkShearSfcTen(){
  int tail=0;
  int head = S->th->mintenpos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
	
  return result.abs();
}

double IndicesCollector::BulkShearSfc20(){
  int tail=0;
  int head = S->th->min20pos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
	
  return result.abs();
}

double IndicesCollector::BulkShearSfczero(){
  int tail=0;
  int head = S->th->zeropos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
	
  return result.abs();
}

double IndicesCollector::BulkShear1kmzero(){
  int tail=cache->getHeightIndex(1000);
  int head = S->th->zeropos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  
  return result.abs();
}

double IndicesCollector::BulkShear1km20(){
  int tail=cache->getHeightIndex(1000);
  int head = S->th->min20pos;
  
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

double IndicesCollector::BulkShear2kmzero(){
  int tail=cache->getHeightIndex(2000);
  int head = S->th->zeropos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  
  return result.abs();
}

double IndicesCollector::BulkShear2km20(){
  int tail=cache->getHeightIndex(2000);
  int head = S->th->min20pos;
  
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


double IndicesCollector::BulkShearMULFCTen(){
  int tail=S->th->mostUnstable->vLfcIndex;
  int head = S->th->mintenpos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  
  double effSHR = result.abs();
  double mucape = this->VMostUnstableCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::BulkShearMULFCzero(){
  int tail=S->th->mostUnstable->vLfcIndex;
  int head = S->th->zeropos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  
  double effSHR = result.abs();
  double mucape = this->VMostUnstableCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::BulkShearMULFC20(){
  int tail = S->th->mostUnstable->vLfcIndex;
  int head = S->th->min20pos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  
  double effSHR = result.abs();
  double mucape = this->VMostUnstableCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::BulkShearMLLFCTen(){
  int tail=S->th->meanLayer->vLfcIndex;
  int head = S->th->mintenpos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  
  double effSHR = result.abs();
  double mucape = this->VMeanLayerCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::BulkShearSBLFCTen(){
  int tail=S->th->surfaceBased->vLfcIndex;
  int head = S->th->mintenpos;
  
  Vector vtail = Get(S->ks->vw,tail);
  Vector vhead = Get(S->ks->vw,head);
  Vector result = vhead-vtail;
  
  double effSHR = result.abs();
  double mucape = this->VSurfaceBasedCAPE();
  if(mucape==0)effSHR=0;  
  return effSHR;
}

double IndicesCollector::MoistureFlux(){
	double result = (S->th->meanmxr2)*(S->ks->mean02.abs());
	return result;
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

double IndicesCollector::RHMIDDLE(){
	return S->th->meanhumMIDDLE;
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

double IndicesCollector::wind_3km(){
  return Get(S->ks->v,cache->getHeightIndex(3000))* 0.514444;
}

double IndicesCollector::wind_2km(){
  return Get(S->ks->v,cache->getHeightIndex(2000))* 0.514444;
}

double IndicesCollector::wind_1km(){
  return Get(S->ks->v,cache->getHeightIndex(1000))* 0.514444;
}

double IndicesCollector::wind_500m(){
  return Get(S->ks->v,cache->getHeightIndex(500))* 0.514444;
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

double IndicesCollector::SB_coldcape_fraction(){
	double cape = this->VSurfaceBasedCAPE(); 
	double coldcape = S->th->surfaceBased->coldcape;
	double fraction=0;
	if(cape>0)fraction=coldcape/cape;
        if(cape==0)fraction=0;
	return fraction;
}

double IndicesCollector::MU_coldcape_fraction(){
	double cape = this->VMostUnstableCAPE(); 
	double coldcape = S->th->mostUnstable->coldcape;
	double fraction=0;
	if(cape>0)fraction=coldcape/cape;
        if(cape==0)fraction=0;
	return fraction;
}

double IndicesCollector::ML_coldcape_fraction(){
	double cape = this->VMeanLayerCAPE(); 
	double coldcape = S->th->meanLayer->coldcape;
	double fraction=0;
	if(cape>0)fraction=coldcape/cape;
        if(cape==0)fraction=0;
	return fraction;
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
	  S->th->mostU500->vcape;
  }
  double IndicesCollector::MU500_CIN(){
	  S->th->mostU500->vcin;
  }
  double IndicesCollector::MU500_LI(){
	   int lindex = cache->getPressureIndex(500);
  
		double lit = Get(S->t,lindex);
  
		int vindex = lindex - S->th->mostU500->startIndex;
  
  
		double plit = Get(S->th->mostU500->virtualValues,vindex);
	
		double Showalter = lit - plit;
		return Showalter;
  }
  
  double IndicesCollector::tilting(){
	  
		int m20i = S->th->min20pos;
		int lfcpos = S->th->mostUnstable->vLfcIndex;
		Vector lfc=Get(S->ks->vw,lfcpos);;
		Vector m20=Get(S->ks->vw,m20i);
		Vector shear = m20-lfc;
	  
		Vector storm = S->ks->mean06;
		Vector mean = S->ks->mean020;
		Vector sr = mean-storm;
		double pojebnia=(sr.X()*shear.X() +sr.Y()*shear.Y())/(sr.abs()*shear.abs());
		return pojebnia;
  }
  
double * processSounding(double *p_, double *h_, double *t_, double *d_, double *a_, double *v_, int length, double dz, Sounding **S){
  *S = new Sounding(p_,h_,t_,d_,a_,v_,length, dz);
  double * vec = new double[150];
  vec[0]=(*S)->getIndicesCollectorPointer()->VMostUnstableCAPE();
  vec[1]=(*S)->getIndicesCollectorPointer()->MU_coldcape_fraction();
  vec[2]=(*S)->getIndicesCollectorPointer()->VLLMostUnstableCAPE();
  vec[3]=(*S)->getIndicesCollectorPointer()->MUmiddlecape();
  vec[4]=(*S)->getIndicesCollectorPointer()->VMostUnstableCIN();  
  vec[5]=(*S)->getIndicesCollectorPointer()->VMostUnstableLCL();
  vec[6]=(*S)->getIndicesCollectorPointer()->VMostUnstableLFC();
  vec[7]=(*S)->getIndicesCollectorPointer()->VMostUnstableEL();
  vec[8]=(*S)->getIndicesCollectorPointer()->VMostUnstableLI();
  vec[9]=(*S)->getIndicesCollectorPointer()->VMostUnstableVmax();
  vec[10]=(*S)->getIndicesCollectorPointer()->MUELTemperature();
  vec[11]=(*S)->getIndicesCollectorPointer()->MULCLTemperature();
  vec[12]=(*S)->getIndicesCollectorPointer()->MULFCTemperature();
  vec[13]=(*S)->getIndicesCollectorPointer()->MUMRatio();   
  vec[14]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedCAPE();
  vec[15]=(*S)->getIndicesCollectorPointer()->SB_coldcape_fraction();
  vec[16]=(*S)->getIndicesCollectorPointer()->VLLSurfaceBasedCAPE();
  vec[17]=(*S)->getIndicesCollectorPointer()->SBmiddlecape();
  vec[18]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedCIN();  
  vec[19]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedLCL();
  vec[20]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedLFC();
  vec[21]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedEL();
  vec[22]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedLI();
  vec[23]=(*S)->getIndicesCollectorPointer()->VSurfaceBasedVmax(); 
  vec[24]=(*S)->getIndicesCollectorPointer()->SBELTemperature();
  vec[25]=(*S)->getIndicesCollectorPointer()->SBLCLTemperature();
  vec[26]=(*S)->getIndicesCollectorPointer()->SBLFCTemperature(); 
  vec[27]=(*S)->getIndicesCollectorPointer()->SBMRatio();
  vec[28]=(*S)->getIndicesCollectorPointer()->VMeanLayerCAPE();
  vec[29]=(*S)->getIndicesCollectorPointer()->ML_coldcape_fraction();
  vec[30]=(*S)->getIndicesCollectorPointer()->VLLMeanLayerCAPE();
  vec[31]=(*S)->getIndicesCollectorPointer()->MLmiddlecape();
  vec[32]=(*S)->getIndicesCollectorPointer()->VMeanLayerCIN();  
  vec[33]=(*S)->getIndicesCollectorPointer()->VMeanLayerLCL();
  vec[34]=(*S)->getIndicesCollectorPointer()->VMeanLayerLFC();
  vec[35]=(*S)->getIndicesCollectorPointer()->VMeanLayerEL();
  vec[36]=(*S)->getIndicesCollectorPointer()->VMeanLayerLI();
  vec[37]=(*S)->getIndicesCollectorPointer()->VMeanLayerVmax();
  vec[38]=(*S)->getIndicesCollectorPointer()->MLELTemperature();
  vec[39]=(*S)->getIndicesCollectorPointer()->MLLCLTemperature();
  vec[40]=(*S)->getIndicesCollectorPointer()->MLLFCTemperature(); 
  vec[41]=(*S)->getIndicesCollectorPointer()->MLMixingRatio();
  vec[42]=(*S)->getIndicesCollectorPointer()->LR0500();
  vec[43]=(*S)->getIndicesCollectorPointer()->LapseRate01();
  vec[44]=(*S)->getIndicesCollectorPointer()->LR02();
  vec[45]=(*S)->getIndicesCollectorPointer()->lapserate03();
  vec[46]=(*S)->getIndicesCollectorPointer()->LR04(); 
  vec[47]=(*S)->getIndicesCollectorPointer()->LR06();
  vec[48]=(*S)->getIndicesCollectorPointer()->LapseRate24(); 
  vec[49]=(*S)->getIndicesCollectorPointer()->LR36();
  vec[50]=(*S)->getIndicesCollectorPointer()->lapseRate500700();
  vec[51]=(*S)->getIndicesCollectorPointer()->lapseRate500800(); 
  vec[52]=(*S)->getIndicesCollectorPointer()->lapseRate600800();
  vec[53]=(*S)->getIndicesCollectorPointer()->ZeroHeight();
  vec[54]=(*S)->getIndicesCollectorPointer()->WetBulbZeroHeight();  
  vec[55]=(*S)->getIndicesCollectorPointer()->MUHeight();
  vec[56]=(*S)->getIndicesCollectorPointer()->MinTHTEHeight();
  vec[57]=(*S)->getIndicesCollectorPointer()->DeltaThetaE();
  vec[58]=(*S)->getIndicesCollectorPointer()->DeltaThetaE_HGL();
  vec[59]=(*S)->getIndicesCollectorPointer()->VDCAPE(); 
  vec[60]=(*S)->getIndicesCollectorPointer()->VirtualColdPoolStrength();
  vec[61]=(*S)->getIndicesCollectorPointer()->WindIndex();
  vec[62]=(*S)->getIndicesCollectorPointer()->PWATER();
  vec[63]=(*S)->getIndicesCollectorPointer()->MoistureFlux();
  vec[64]=(*S)->getIndicesCollectorPointer()->RH02();
  vec[65]=(*S)->getIndicesCollectorPointer()->RH25();
  vec[66]=(*S)->getIndicesCollectorPointer()->RH36();
  vec[67]=(*S)->getIndicesCollectorPointer()->RHMIDDLE(); 
  vec[68]=(*S)->getIndicesCollectorPointer()->wind_500m();
  vec[69]=(*S)->getIndicesCollectorPointer()->wind_1km();
  vec[70]=(*S)->getIndicesCollectorPointer()->wind_2km();
  vec[71]=(*S)->getIndicesCollectorPointer()->wind_3km();
  vec[72]=(*S)->getIndicesCollectorPointer()->BS01();
  vec[73]=(*S)->getIndicesCollectorPointer()->BS02();
  vec[74]=(*S)->getIndicesCollectorPointer()->BS03();
  vec[75]=(*S)->getIndicesCollectorPointer()->BS06();
  vec[76]=(*S)->getIndicesCollectorPointer()->BS08();
  vec[77]=(*S)->getIndicesCollectorPointer()->BS36();
  vec[78]=(*S)->getIndicesCollectorPointer()->BS26();
  vec[79]=(*S)->getIndicesCollectorPointer()->BS16();
  vec[80]=(*S)->getIndicesCollectorPointer()->BS18();
  vec[81]=(*S)->getIndicesCollectorPointer()->emubs();
  vec[82]=(*S)->getIndicesCollectorPointer()->esbbs();
  vec[83]=(*S)->getIndicesCollectorPointer()->emlbs();
  vec[84]=(*S)->getIndicesCollectorPointer()->BulkShearSfczero();
  vec[85]=(*S)->getIndicesCollectorPointer()->BulkShearSfcTen();	
  vec[86]=(*S)->getIndicesCollectorPointer()->BulkShearSfc20();
  vec[87]=(*S)->getIndicesCollectorPointer()->BulkShear1kmzero();
  vec[88]=(*S)->getIndicesCollectorPointer()->BulkShear1kmTen();
  vec[89]=(*S)->getIndicesCollectorPointer()->BulkShear1km20();
  vec[90]=(*S)->getIndicesCollectorPointer()->BulkShear2kmzero();
  vec[91]=(*S)->getIndicesCollectorPointer()->BulkShear2kmTen();
  vec[92]=(*S)->getIndicesCollectorPointer()->BulkShear2km20();
  vec[93]=(*S)->getIndicesCollectorPointer()->BulkShearMULFCzero();
  vec[94]=(*S)->getIndicesCollectorPointer()->BulkShearMULFCTen();
  vec[95]=(*S)->getIndicesCollectorPointer()->BulkShearMULFC20();
  vec[96]=(*S)->getIndicesCollectorPointer()->BulkShearSBLFCTen();
  vec[97]=(*S)->getIndicesCollectorPointer()->BulkShearMLLFCTen();
  vec[98]=(*S)->getIndicesCollectorPointer()->MSR_MW();
  vec[99]=(*S)->getIndicesCollectorPointer()->MSR_RM();
  vec[100]=(*S)->getIndicesCollectorPointer()->MSR_LM();
  vec[101]=(*S)->getIndicesCollectorPointer()->MSR_MW_HGL();
  vec[102]=(*S)->getIndicesCollectorPointer()->MSR_RM_HGL();
  vec[103]=(*S)->getIndicesCollectorPointer()->MSR_LM_HGL();
  vec[104]=(*S)->getIndicesCollectorPointer()->tilting();
  vec[105]=(*S)->getIndicesCollectorPointer()->MeanWind01();
  vec[106]=(*S)->getIndicesCollectorPointer()->MeanWind02();
  vec[107]=(*S)->getIndicesCollectorPointer()->MeanWind03();
  vec[108]=(*S)->getIndicesCollectorPointer()->MeanWind06();
  vec[109]=(*S)->getIndicesCollectorPointer()->MeanWind13();
  vec[110]=(*S)->getIndicesCollectorPointer()->SRH100RM();
  vec[111]=(*S)->getIndicesCollectorPointer()->SRH250RM();
  vec[112]=(*S)->getIndicesCollectorPointer()->SRH500RM();
  vec[113]=(*S)->getIndicesCollectorPointer()->SRH01RM();
  vec[114]=(*S)->getIndicesCollectorPointer()->SRH03RM();
  vec[115]=(*S)->getIndicesCollectorPointer()->SRH100LM();
  vec[116]=(*S)->getIndicesCollectorPointer()->SRH250LM();
  vec[117]=(*S)->getIndicesCollectorPointer()->SRH500LM();
  vec[118]=(*S)->getIndicesCollectorPointer()->SRH01LM();
  vec[119]=(*S)->getIndicesCollectorPointer()->SRH03LM();
  vec[120]=(*S)->getIndicesCollectorPointer()->Bunkers_RM_A();
  vec[121]=(*S)->getIndicesCollectorPointer()->Bunkers_RM_M();
  vec[122]=(*S)->getIndicesCollectorPointer()->Bunkers_LM_A();
  vec[123]=(*S)->getIndicesCollectorPointer()->Bunkers_LM_M();
  vec[124]=(*S)->getIndicesCollectorPointer()->Bunkers_MW_A();
  vec[125]=(*S)->getIndicesCollectorPointer()->Bunkers_MW_M();
  vec[126]=(*S)->getIndicesCollectorPointer()->Corfidi_downwind_A();
  vec[127]=(*S)->getIndicesCollectorPointer()->Corfidi_downwind_M();
  vec[128]=(*S)->getIndicesCollectorPointer()->Corfidi_upwind_A();
  vec[129]=(*S)->getIndicesCollectorPointer()->Corfidi_upwind_M();
  vec[130]=(*S)->getIndicesCollectorPointer()->K_Index();
  vec[131]=(*S)->getIndicesCollectorPointer()->Showalter_Index(); 
  vec[132]=(*S)->getIndicesCollectorPointer()->TotalTotals();  
  vec[133]=(*S)->getIndicesCollectorPointer()->SWEATIndex(); 
  vec[134]=(*S)->getIndicesCollectorPointer()->STP();
  vec[135]=(*S)->getIndicesCollectorPointer()->STPeff();
  vec[136]=(*S)->getIndicesCollectorPointer()->SCP();
  vec[137]=(*S)->getIndicesCollectorPointer()->SCPeff();
  vec[138]=(*S)->getIndicesCollectorPointer()->SHP();
  vec[139]=(*S)->getIndicesCollectorPointer()->HSI();
  vec[140]=(*S)->getIndicesCollectorPointer()->DCP();
  vec[141]=(*S)->getIndicesCollectorPointer()->MU_WMAXSHEAR();
  vec[142]=(*S)->getIndicesCollectorPointer()->SB_WMAXSHEAR();
  vec[143]=(*S)->getIndicesCollectorPointer()->ML_WMAXSHEAR();
  vec[144]=(*S)->getIndicesCollectorPointer()->MU_EFF_WMAXSHEAR();
  vec[145]=(*S)->getIndicesCollectorPointer()->SB_EFF_WMAXSHEAR();
  vec[146]=(*S)->getIndicesCollectorPointer()->ML_EFF_WMAXSHEAR();
  vec[147]=(*S)->getIndicesCollectorPointer()->MU500_CAPE();
  vec[148]=(*S)->getIndicesCollectorPointer()->MU500_CIN();
  vec[149]=(*S)->getIndicesCollectorPointer()->MU500_LI();
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
    double defh[]={0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 15000, 16000, 17000, 18000, 19000, 20000}; int hlen = 60;
    
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
       
        delete(*pu);delete(*hu);delete(*tu);delete(*du);delete(*au);delete(*vu);
       
        *pu=np;*hu=nh;*tu=nt;*du=nd;*au=na;*vu=nv;
       
        return u;
       
   
}
int interpolate2(double **pu, double **hu, double **tu, double **du, double **au, double **vu, int n, double *p_arr=0, double *h_arr=0, int m=0, int o=0){
    double defp[]={850,700,500}; int plen = 3;
    double defh[]={0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 15000, 16000, 17000, 18000, 19000, 20000}; int hlen = 60;
  
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
       
        delete(*pu);delete(*hu);delete(*tu);delete(*du);delete(*au);delete(*vu);
       
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
						  int custom_vec=1
                          )
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
	  
	  double *tmp = new double[4000];
	  tmp[0]=0;
	  for(int g=1;g<4000;g++){
		  tmp[g]=tmp[g-1]+5;
	  }
	  hh=tmp;
	  ch=4000;
  }
  u = interpolate2(&p, &h,&t, &d, &a, &v, u,pp, hh, 0, ch);
  if(custom_vec==3)delete[] hh;
  
  double *result = processSounding(p,h,t,d,a,v,u,step,sret);

  delete[]p;delete[]h;delete[]t;delete[]d;delete[]a;delete[]v;
  return result;
 
}

//' Generic function for calculating thermo- and kinematic indices derived from atmospheric profiling.
//' Further details given in the sounding_compute() function
//' 
//'
//' @param pressure pressure [hPa]
//' @param altitude altitude [meters]
//' @param temp temperature [degree Celsius]
//' @param dpt dew point temperature [degree Celsius]
//' @param wd wind direction [azimuth in degrees]
//' @param ws wind speed [knots]
//' @param export_profile possibility to export interpolated profile on the levels defined in accuracy configuration
//' @param accuracy accuracy of computations where 3 = high (slow), 2 = medium (recommended), 1 = low (fast)
//' @examples 
//' pressure = c(1000, 855, 700, 500, 300, 100, 10) 
//' altitude = c(0, 1500, 2500, 6000, 8500, 12000, 25000)
//' temp = c(25, 10, 0, -15, -30, -50, -92)
//' dpt = c(20, 5, -5, -30, -55, -80, -99)
//' wd = c(0, 90, 135, 180, 270, 350, 0)
//' ws = c(5, 10, 20, 30, 40, 5, 0)
//' sounding_default(pressure, altitude, temp, dpt, wd, ws, accuracy = 1, export_profile = 0)
//' @useDynLib thunder
//' @importFrom Rcpp evalCpp
//' @export
//' @return 
//' \enumerate{
//'  \item MU_CAPE 
//'  \item MU_CAPE_M10_fraction 
//'  \item MU_03km_CAPE 
//'  \item MU_HGL_CAPE 
//'  \item MU_CIN 
//'  \item MU_LCL_HGT 
//'  \item MU_LFC_HGT 
//'  \item MU_EL_HGT 
//'  \item MU_LI 
//'  \item MU_WMAX 
//'  \item MU_EL_TEMP 
//'  \item MU_LCL_TEMP 
//'  \item MU_LFC_TEMP 
//'  \item MU_MIXR 
//'  \item SB_CAPE 
//'  \item SB_CAPE_M10_fraction 
//'  \item SB_03km_CAPE 
//'  \item SB_HGL_CAPE 
//'  \item SB_CIN 
//'  \item SB_LCL_HGT 
//'  \item SB_LFC_HGT 
//'  \item SB_EL_HGT 
//'  \item SB_LI 
//'  \item SB_WMAX 
//'  \item SB_EL_TEMP 
//'  \item SB_LCL_TEMP 
//'  \item SB_LFC_TEMP 
//'  \item SB_MIXR 
//'  \item ML_CAPE 
//'  \item ML_CAPE_M10_fraction 
//'  \item ML_03km_CAPE 
//'  \item ML_HGL_CAPE 
//'  \item ML_CIN 
//'  \item ML_LCL_HGT 
//'  \item ML_LFC_HGT 
//'  \item ML_EL_HGT 
//'  \item ML_LI 
//'  \item ML_WMAX 
//'  \item ML_EL_TEMP 
//'  \item ML_LCL_TEMP 
//'  \item ML_LFC_TEMP 
//'  \item ML_MIXR 
//'  \item LR_0500m 
//'  \item LR_01km 
//'  \item LR_02km 
//'  \item LR_03km 
//'  \item LR_04km  
//'  \item LR_06km 
//'  \item LR_24km 
//'  \item LR_36km 
//'  \item LR_500700hPa 
//'  \item LR_500800hPa 
//'  \item LR_600800hPa 
//'  \item FRZG_HGT 
//'  \item FRZG_wetbulb_HGT 
//'  \item HGT_max_thetae_03km 
//'  \item HGT_min_thetae_04km 
//'  \item Delta_thetae 
//'  \item Delta_thetae_HGL 
//'  \item DCAPE 
//'  \item Cold_Pool_Strength 
//'  \item Wind_Index 
//'  \item PRCP_WATER 
//'  \item Moisture_Flux_02km 
//'  \item RH_02km 
//'  \item RH_25km 
//'  \item RH_36km 
//'  \item RH_HGL 
//'  \item Wind_500m 
//'  \item Wind_1km 
//'  \item Wind_2km 
//'  \item Wind_3km 
//'  \item BS_01km 
//'  \item BS_02km 
//'  \item BS_03km 
//'  \item BS_06km 
//'  \item BS_08km 
//'  \item BS_36km 
//'  \item BS_26km 
//'  \item BS_16km 
//'  \item BS_18km 
//'  \item BS_EFF_MU 
//'  \item BS_EFF_SB 
//'  \item BS_EFF_ML 
//'  \item BS_SFC_to_0 
//'  \item BS_SFC_to_M10 
//'  \item BS_SFC_to_M20 
//'  \item BS_1km_to_0 
//'  \item BS_1km_to_M10 
//'  \item BS_1km_to_M20 
//'  \item BS_2km_to_0 
//'  \item BS_2km_to_M10 
//'  \item BS_2km_to_M20 
//'  \item BS_MU_LFC_to_0 
//'  \item BS_MU_LFC_to_M10 
//'  \item BS_MU_LFC_to_M20 
//'  \item BS_SB_LFC_to_M10 
//'  \item BS_ML_LFC_to_M10 
//'  \item BS_MW02_to_MW06 
//'  \item BS_MW02_to_RM 
//'  \item BS_MW02_to_LM 
//'  \item BS_HGL_to_MW06 
//'  \item BS_HGL_to_RM 
//'  \item BS_HGL_to_LM 
//'  \item Tilting 
//'  \item MW_01km 
//'  \item MW_02km 
//'  \item MW_03km 
//'  \item MW_06km 
//'  \item MW_13km 
//'  \item SRH_100m_RM 
//'  \item SRH_250m_RM 
//'  \item SRH_500m_RM 
//'  \item SRH_1km_RM 
//'  \item SRH_3km_RM 
//'  \item SRH_100m_LM 
//'  \item SRH_250m_LM 
//'  \item SRH_500m_LM 
//'  \item SRH_1km_LM 
//'  \item SRH_3km_LM 
//'  \item Bunkers_RM_A 
//'  \item Bunkers_RM_M 
//'  \item Bunkers_LM_A 
//'  \item Bunkers_LM_M 
//'  \item Bunkers_MW_A 
//'  \item Bunkers_MW_M 
//'  \item Corfidi_downwind_A 
//'  \item Corfidi_downwind_M 
//'  \item Corfidi_upwind_A 
//'  \item Corfidi_upwind_M 
//'  \item K_Index 
//'  \item Showalter_Index 
//'  \item TotalTotals_Index 
//'  \item SWEAT_Index 
//'  \item STP 
//'  \item STP_new 
//'  \item SCP 
//'  \item SCP_new 
//'  \item SHIP 
//'  \item HSI 
//'  \item DCP 
//'  \item MU_WMAXSHEAR 
//'  \item SB_WMAXSHEAR 
//'  \item ML_WMAXSHEAR 
//'  \item MU_EFF_WMAXSHEAR 
//'  \item SB_EFF_WMAXSHEAR 
//'  \item ML_EFF_WMAXSHEAR 
//'  \item MU500_CAPE 
//'  \item MU500_CIN 
//'  \item MU500_LI 
//' }
// [[Rcpp::export]]

Rcpp::NumericVector sounding_default(Rcpp::NumericVector pressure,
                          Rcpp::NumericVector altitude,
                          Rcpp::NumericVector temp,
                          Rcpp::NumericVector dpt,
                          Rcpp::NumericVector wd,
                          Rcpp::NumericVector ws,
		          Rcpp::NumericVector export_profile,
			  Rcpp::NumericVector accuracy
						  )
{
  Sounding *sret;
  int size = pressure.size();
  
  double *p = new double[size], *h=new double[size], *d=new double[size],*t=new double[size], *a= new double[size], *v=new double[size];
  for(int i=0; i<size;i++){
    p[i]=pressure[i];
    h[i]=altitude[i];
    t[i]=temp[i];
    d[i]=dpt[i];
    a[i]=wd[i];
    v[i]=ws[i];
  }
  int q= accuracy[0];
  int plen,hlen,tlen,dlen,alen,vlen,tvlen;
  int mulen, sblen,mllen,mustart;

  double *result = sounding_default2(p,h,t,d,a,v,size,&sret,q);
	int reslen= 150;
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
  	mustart= sret->th->mostUnstable->startIndex;
	maxl+=2+mulen+1+sblen+1+mllen+plen+1+hlen+1+tlen+1+dlen+1+alen+1+vlen+1+tvlen+1+10;
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
	return out;
	
}
