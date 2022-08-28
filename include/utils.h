/*************************************************************************
    > File Name: utils.h
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2017 02:27:38 PM
 ************************************************************************/

#include<iostream>
#include<vector>
#include<string>
#include<assert.h>
#include <algorithm>    // std::find
#include <ctime>
#include <chrono>
#include <fstream>
#include <queue>
#include <cstring>
#include <math.h>       /* log */

using namespace std;

typedef unsigned int uint;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;
typedef std::pair<uint, uint> ipair;
typedef std::pair<double, double> dpair;
#define MP std::make_pair

#ifndef TIMES_PER_SEC
#define TIMES_PER_SEC (1.0e9)
#endif

#define SIZE(t) (int)(t.size())
#define ALL(t) (t).begin(), (t).end()
#define FOR(i, n) for(int (i)=0; (i)<((int)(n)); (i)++)

const std::string TP = "tp";
const std::string TPC = "tpc";
const std::string MC2 = "mc2";
const std::string AMC = "amc";
const std::string GEER = "geer";
const std::string SMM = "smm";
const std::string HAY = "hay";

struct Config{
    std::string strFolder;
    std::string strGraph;
    std::string strAlgo;
    double epsilon=0;
    double lambda=0;
    double delta=0;
    uint64 numwalks=0;
    uint ds=0;
    uint dt=0;
    uint lenwalk=0;
    uint tau=0;
    uint ell=0;
    uint64 eta_star_s=0;
    uint64 eta_star_t=0;
    uint64 eta_star_all=0;
    double psi_s_prime=0;
    double psi_t_prime=0;
    int edgeSeed=0;
    int numQuery=0;

    void display(){
        std::cout << "====================Configurations==================" << std::endl;
        std::cout << "data folder: " << strFolder << '\n';
        std::cout << "graph file name: " << strGraph << '\n';
        std::cout << "algorithm: " << strAlgo << '\n';
        std::cout << "absolute error epsilon: " << epsilon << '\n';
        std::cout << "failure probability: " << delta << '\n';
        std::cout << "#iterations: " << tau << '\n';
        std::cout << "use edge seeds: " << edgeSeed << '\n';
        std::cout << "#queries: " << numQuery << '\n';
        // std::cout << "#randomwalks: " << numwalks << '\n';
        //std::cout << "iteration: " << iteration << '\n';
        std::cout << "====================Configurations==================" << std::endl;
    }
    void check(){
        std::vector<std::string> Algos = {TP,GEER,AMC,SMM,TPC,MC2,HAY};
        auto f = std::find(Algos.begin(), Algos.end(), strAlgo);
        assert (f != Algos.end());
    }

    // void init(double eps, double lambda, double delta, uint tau){
    //     epsilon = eps;
    //     lambda = lambda;
    //     delta = delta;
    //     tau = tau;
    // }

    void setTP(){
        ell = uint(log(4.0/epsilon/(1-lambda))/log(1.0/lambda)-1.0);
        lenwalk = ell;
        numwalks = uint64(40*lenwalk*lenwalk*log(8*lenwalk/delta)/epsilon/epsilon);
    }

    void setTPC(){
        ell = uint(log(4.0/epsilon/(1-lambda))/log(1.0/lambda)-1.0);
        lenwalk = ell;
        // numwalks = 20000*(sqrt()+pow(lenwalk,3)/epsilon/epsilon);
    }

    void setHAY(uint64 m){
        numwalks = uint64(log(2*m/delta)/2.0/epsilon/epsilon);
    }

    void setMC2(double gamma){
        numwalks = uint64(3*log(1.0/delta)/epsilon/epsilon/gamma);
    }

    void calEll(){
        double dell = log(2.0*(1.0/ds+1.0/dt)/(1-lambda)/epsilon)/log(1.0/lambda)-1.0;
        // cout << "dell:" << dell << endl;
        if(dell<0){
            ell=0;
        }
        else{
            ell = uint(dell);
        }
        // cout << "ell:" << ell << endl;
    }

    void calEta(){
        eta_star_s = uint64(8*log(3*(tau+1)/delta)*lenwalk*lenwalk*pow(max(psi_s_prime/ds,psi_t_prime/dt),2)/epsilon/epsilon);
        eta_star_t = uint64(8*log(3*(tau+1)/delta)*lenwalk*lenwalk*psi_t_prime*psi_t_prime/epsilon/epsilon/dt/dt);
    }

    void setAMC(uint degs, uint degt){
        ds = degs;
        dt = degt;
        calEll();
        lenwalk = ell;
        psi_s_prime = 1.0;
        psi_t_prime = 1.0;
        calEta();
        cout << "lenwalk: " << lenwalk << endl;
        cout << "numwalks: " << eta_star_s << " from s and " << eta_star_s << " from t" << endl;
    }

    void setSMM(uint degs, uint degt){
        ds = degs;
        dt = degt;
        calEll();
    }

    void calEtaAll(){
        eta_star_all = uint64(8*log(tau/delta)*lenwalk*lenwalk*pow(psi_s_prime/ds+psi_t_prime/dt,2)/epsilon/epsilon);
    }

    void setGEER(double smax, double tmax, uint ellf, uint degs, uint degt){
        ds = degs;
        dt = degt;
        lenwalk = ellf;
        psi_s_prime = smax;
        psi_t_prime = tmax;
        calEta();
        calEtaAll();
    }

    double func(uint64 eta_x, double sigma, double deltaprime, double psi){
        return sqrt(2.0*sigma*log(3.0/deltaprime)/eta_x)+3*psi*log(3.0/deltaprime)/eta_x;
    }

    double calErr(uint64 eta_s, uint64 eta_t, double sigma_ss, double sigma_st, double sigma_tt){
        double err = 0;
        double deltaprime = delta/3.0/(tau+1);
        double psi_s = lenwalk*psi_s_prime/2.0;
        double psi_t = lenwalk*psi_t_prime/2.0;
        double err1 = func(eta_s,sigma_ss,deltaprime,psi_s);
        double err2 = func(eta_s,sigma_st,deltaprime,psi_t);
        double err3 = func(eta_t,sigma_tt,deltaprime,psi_t);
        cout << "sigma: " << sigma_ss << " " << sigma_st << " " << sigma_tt << endl;
        cout << "pis: " << psi_s << " " << psi_s << endl;
        cout << "eta: " << eta_s << " " << eta_t << endl;
        cout << "err: " << err1 << " " << err2 << " " << err3 << endl;
        err = err1+2*err2+err3;
        return err;
    }

    double calNewErr(uint64 eta_all, double sigma_all){
        double err_all = 0;
        double deltaprime = delta/tau;
        double psi_all = 2*lenwalk*(psi_s_prime/ds+psi_t_prime/dt);

        err_all = func(eta_all,sigma_all,deltaprime,psi_all);
        // cout << "sigma: " << sigma_all << endl;
        // cout << "pis: " << psi_all << endl;
        // cout << "eta: " << eta_all << endl;
        // cout << "err: " << err_all << endl;
        return err_all;
    }
};

void process_mem_usage(double& vm_usage, double& resident_set);
void disp_mem_usage();

uint getProcMemory();

class Timer {
public:
    static std::vector<double> timeUsed;
    static std::vector<string> timeUsedDesc;
    int id;
    std::chrono::steady_clock::time_point startTime;
    bool showOnDestroy;

    Timer(int id, string desc = "", bool showOnDestroy = false) {
        this->id = id;
        while ((int) timeUsed.size() <= id) {
            timeUsed.push_back(0);
            timeUsedDesc.push_back("");
        }
        timeUsedDesc[id] = desc;
        startTime = std::chrono::steady_clock::now();
        this->showOnDestroy = showOnDestroy;
    }

    static double used(int id) {
        return timeUsed[id] / TIMES_PER_SEC;
    }

    ~Timer() {
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - startTime).count();
        if (showOnDestroy) {
            std::cout << "time spend on " << timeUsedDesc[id] << ":" << duration / TIMES_PER_SEC << "s" << std::endl;
        }
        timeUsed[id] += duration;
    }

    static void show(bool debug = false) {
        cout << "##### Timer #####" << endl;
        for (int i = 0; i < (int) timeUsed.size(); i++) {
            if (timeUsed[i] > 0) {
                char str[100];
                sprintf(str, "%.6lf", timeUsed[i] / TIMES_PER_SEC);
                string s = str;
                if ((int) s.size() < 15) s = " " + s;
                char t[100];
                memset(t, 0, sizeof t);
                sprintf(t, "%4d %s %s", i, s.c_str(), timeUsedDesc[i].c_str());
                cout << t << endl;
            }
        }
    }

    static void reset(int id){
        if(id>=0 && id<timeUsed.size()){
            timeUsed[id] = 0;
        }
    }

    static void clearAll() {
        timeUsed.clear();
        timeUsedDesc.clear();
    }
};

const int VectorDefaultSize=20;

template <typename _T>
class iVector
{
public:
    uint m_size;
    _T* m_data;
    uint m_num;

    void free_mem()
    {
        delete[] m_data;
    }

    iVector()
    {
        //printf("%d\n",VectorDefaultSize);
        m_size = VectorDefaultSize;
        m_data = new _T[VectorDefaultSize];
        m_num = 0;
    }
    iVector( uint n )
    {
        if ( n == 0 )
        {
            n = VectorDefaultSize;
        }
//      printf("iVector allocate: %d\n",n);
        m_size = n;
        m_data = new _T[m_size];
        m_num = 0;
    }
    void push_back( _T d )
    {
        // if ( m_num == m_size )
        // {
        //     re_allocate( m_size*2 );
        // }
        m_data[m_num] = d ;
        m_num++;        
    }
    void push_back( const _T* p, uint len )
    {
        while ( m_num + len > m_size )
        {
            re_allocate( m_size*2 );
        }
        memcpy( m_data+m_num, p, sizeof(_T)*len );
        m_num += len;
    }

    void re_allocate( uint size )
    {
        if ( size < m_num )
        {
            return;
        }
        _T* tmp = new _T[size];
        memcpy( tmp, m_data, sizeof(_T)*m_num );
        m_size = size;
        delete[] m_data;
        m_data = tmp;
    }
    void Sort()
    {
        if ( m_num < 20 )
        {
            int k ;
            _T tmp;
            for ( int i = 0 ; i < m_num-1 ; ++i )
            {
                k = i ;
                for ( int j = i+1 ; j < m_num ; ++j )
                    if ( m_data[j] < m_data[k] ) k = j ;
                if ( k != i )
                {
                    tmp = m_data[i];
                    m_data[i] = m_data[k];
                    m_data[k] = tmp;
                }
            }
        }
        else sort( m_data, m_data+m_num );
    }
    void unique()
    {
        if ( m_num == 0 ) return;
        Sort();
        uint j = 0;
        for ( uint i = 0 ; i < m_num ; ++i )
            if ( !(m_data[i] == m_data[j]) )
            {
                ++j;
                if ( j != i ) m_data[j] = m_data[i];
            }
        m_num = j+1;
    }
    int BinarySearch( _T& data )
    {
        for ( int x = 0 , y = m_num-1 ; x <= y ; )
        {
            int p = (x+y)/2;
            if ( m_data[p] == data ) return p;
            if ( m_data[p] < data ) x = p+1;
            else y = p-1;
        }
        return -1;
    }
    void clean()
    {
        m_num = 0;
    }
    void assign( iVector& t )
    {
        m_num = t.m_num;
        m_size = t.m_size;
        delete[] m_data;
        m_data = t.m_data;
    }

    bool remove( _T& x )
    {
        for ( int l = 0 , r = m_num ; l < r ; )
        {
            int m = (l+r)/2;

            if ( m_data[m] == x )
            {
                m_num--;
                if ( m_num > m ) memmove( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                return true;
            }
            else if ( m_data[m] < x ) l = m+1;
            else r = m;
        }
        return false;
    }

    void sorted_insert( _T& x )
    {
        if ( m_num == 0 )
        {
            push_back( x );
            return;
        }

        if ( m_num == m_size ) re_allocate( m_size*2 );

        int l,r;

        for ( l = 0 , r = m_num ; l < r ; )
        {
            int m = (l+r)/2;
            if ( m_data[m] < x ) l = m+1;
            else r = m;
        }

        if ( l < m_num && m_data[l] == x )
        {
            //printf("Insert Duplicate....\n");
            //cout<<x<<endl;
    //      break;
        }
        else
        {
            if ( m_num > l )
            {
                memmove( m_data+l+1, m_data+l, sizeof(_T)*(m_num-l) );
            }
            m_num++;
            m_data[l] = x;
        }
    }

    bool remove_unsorted( _T& x )
    {
        for ( int m = 0 ; m < m_num ; ++m )
        {
            if ( m_data[m] == x )
            {
                m_num--;
                if ( m_num > m ) memcpy( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                return true;
            }
        }
        return false;
    }

    _T& operator[]( uint i )
    {
        //if ( i < 0 || i >= m_num ) 
        //{
        //  printf("iVector [] out of range!!!\n");
        //}
        return m_data[i];
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //close range check for [] in iVector if release

};

template <typename _T>
struct iMap
{   
    _T* m_data;
    int m_num;
    int cur;
    iVector<int> occur;
    _T nil; 
    iMap()
    {
        m_data = NULL;
        m_num = 0;
        //nil = std::make_pair((long)-9,(long)-9);
        //nil = 1073741834;
    }
    iMap(int size){
        initialize(size);
    }
    void free_mem()
    {
        delete[] m_data;
        occur.free_mem();
    }

    void initialize( int n )
    {
        occur.re_allocate(n);
        occur.clean();
        m_num = n;
        nil = -9;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = nil;
        cur = 0;
    }
    void clean()
    {
        for ( int i = 0 ; i < occur.m_num ; ++i ){
            m_data[occur[i]] = nil;
        }
        // occur.clean();
        occur.m_num = 0;
        cur = 0;
    }
    
    //init keys 0-n, value as 0
    void init_keys(int n){
        occur.re_allocate(n);
        occur.clean();
        m_num = n;
        nil = -9;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i ){
            m_data[i] = 0;
            occur.push_back( i );
            cur++;
        }
    }
    //reset all values to be zero
    void reset_zero_values(){
        // for ( int i = 0 ; i < m_num ; ++i )
            // m_data[i] = 0.0;
        memset( m_data, 0.0, m_num*sizeof(_T) );
    }

    void reset_one_values(){
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = 1.0;
        // memset( m_data, 0.0, m_num*sizeof(_T) );
    }

    _T get( int p )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap get out of range!!!\n");
        //  return -8;
        //}
        return m_data[p];
    }
    _T& operator[](  int p )
    {
        //if ( i < 0 || i >= m_num ) 
        //{
        //  printf("iVector [] out of range!!!\n");
        //}
        return m_data[p];
    }
    void erase( int p )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap get out of range!!!\n");
        //}
        m_data[p] = nil;
        cur--;
    }
    bool notexist( int p )
    {
        return m_data[p] == nil ;
    }
    bool exist( int p )
    {
        return !(m_data[p] == nil);
    }
    void insert( int p , _T d )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap insert out of range!!!\n");
        //}
        if ( m_data[p] == nil )
        {
            occur.push_back( p );
            cur++;
        }
        m_data[p] = d;
    }
    void inc( int p , _T x )
    {
        if ( m_data[p] == nil ){
            insert(p, x);
        }
        else{
            m_data[p] += x;
        }
    }
    void dec( int p )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("dec some unexisted point\n" );
        //}
        m_data[p]--;
    }
    //close range check when release!!!!!!!!!!!!!!!!!!!!    
};

template <typename _T>
struct myMap
{   
    _T* m_data;
    int* m_keys;
    long m_num;
    long cur;
    _T nil; 
    myMap(){
        m_data = NULL;
        m_num = 0;
        cur = 0;
    }
    myMap(long size){
        initialize(size);
    }
    void free_mem(){
        delete[] m_data;
        delete[] m_keys;
        m_num = 0;
        cur = 0;
    }

    void initialize( long n ){
        m_num = n;
        nil = 0;
        if ( m_data != NULL )
            delete[] m_data;
        if( m_keys != NULL )
            delete[] m_keys;
        m_data = new _T[m_num];
        m_keys = new int[m_num];
        for ( long i = 0 ; i < m_num ; ++i )
            m_data[i] = nil;
        cur = 0;
    }
    void clean(){
        for(long i=0; i<cur; i++){
            m_data[m_keys[i]]=nil;
        }
        cur = 0;
    }

    _T get( int p ){
        return m_data[p];
    }

    int getid( int p ){
        return m_keys[p];
    }

    _T& operator[](  int p ){
        return m_data[p];
    }

    bool notexist( int p ){
        return !m_data[p];
    }

    void insert( int p , _T d ){
        if(!m_data[p]){
            m_keys[cur]=p;
            cur++;
        }
        m_data[p] = d;
    }

    void inc( int p , _T x ){
        if (!m_data[p]){
            m_keys[cur]=p;
            m_data[p] = x;
            cur++;
        }
        else{
            m_data[p] += x;
        }
    }
};

static uint32_t g_seed;
void fastSrand();
uint32_t fastRand();

static uint32_t x_state;
void xorshifinit();
uint32_t xorshift32(void);
