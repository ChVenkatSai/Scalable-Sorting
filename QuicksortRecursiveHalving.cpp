#include <pthread.h>
#include <iostream>
#include <algorithm>
#include <vector>

#define CUTOFF 120

using namespace std;

vector<long long> sums;
vector<long long> sums2;// stores the sum of each sublist
vector<int> zerone;
vector<int> onezero;
vector<int>v_copy;
vector<long long> zeuss;
vector<long long> posiedon;
pthread_barrier_t barrier;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
//pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
int count1 = 0; // counter for number of threads that have broadcasted

struct ThreadArgs {
    vector<int> &v;
    int left;
    int right;
};

struct SegArgs {
    std::vector<int> &v;
    int left_sub;
    int right_sub;
    int id;
    int pivot;
    int pthread_num;
    SegArgs(std::vector<int>& vec, int left, int right, int i, int p, int num)
        : v(vec), left_sub(left), right_sub(right), id(i), pivot(p), pthread_num(num) {}
};



void *sumSublist(void *arg) {
    SegArgs *sa = (SegArgs*)arg;
    vector<int> &v = sa->v;
    int left_sub = sa->left_sub;
    int right_sub = sa->right_sub;
    int pivot = sa->pivot;
    int id = sa->id;
    int pthread_num = sa->pthread_num;
    //long long sum = 0;
    //int sublistSize = sums.size() / pthread_num;
    //int start = id * sublistSize;
    //int end = start + sublistSize;

    // calculate the sum of the sublist assigned to this thread
    //for (int i = start; i < end; i++) {
    //    sum += sums[i];
    //}
    sums[id]=0;
    sums2[id]=0;
    int i = left_sub;
    int j = right_sub;
    pthread_mutex_lock(&mutex);
    while(i<=j){
        if(v[i]<pivot){
            zerone[i]=1;
        }
        if(v[i]>pivot){
            onezero[i]=1;
        }
        i++;
    }
    pthread_mutex_unlock(&mutex);


    // store the sum in the global vector
    i = left_sub;
    j = right_sub;
    //int soma = 0;
    
    pthread_mutex_lock(&mutex);
    while (i<=j){
        sums[id]+= (long long) zerone[i];
        sums2[id]+= (long long) onezero[i];
        //soma+= zerone[i];
        zeuss[i] =  sums[id];
        posiedon[i] = sums2[id] + pivot + 1;
        i++;
    }
    //cout<<sums[id]<<"jj";
    //fflush(stdout);
    pthread_mutex_unlock(&mutex);
    //for (i = 0; i < zeuss.size(); i++) {
    //    cout << zeuss[i] << " ";
    //}
    //fflush(stdout);
    // perform tree-based sum reduction
    for (i = 1; i < pthread_num; i *= 2) {
        if (id % (2*i) == 0) {
            int partner = id + i;
            if (partner < pthread_num) {
                int left_part = left_sub + (partner - id)*(right_sub - left_sub + 1);
                int right_part = left_part + (right_sub-left_sub);
                int l = left_part;
                int p = right_part;
                pthread_mutex_lock(&mutex);
                while(l<=p){
                    zeuss[l]+=sums[id];
                    posiedon[l]+=sums2[id];
                    l++;
                    
                }
                //sums[partner] += sums[id];
                //sums2[partner] += sums2[id];
                count1++;
                pthread_mutex_unlock(&mutex);
            }
        }
        pthread_barrier_wait(&barrier);
    }

    
    
    pthread_exit(NULL);
}

void* partition(void *args) {
    ThreadArgs *ta = (ThreadArgs*)args;
    vector<int> &v = ta->v;
    int left = ta->left;
    int right = ta->right;

    // Choose pivot element
    int pivot_index = (left + right) / 2;
    int pivot = v[pivot_index];


    int pthread_num = 4; // number of threads
    pthread_t threads[pthread_num];
    //fflush(stdout);
    //std::vector<SegArgs> parts(pthread_num, SegArgs(v, 0, 0, 0, 0, 0));
    int ids[pthread_num];

    SegArgs part[4] = {
    SegArgs(v, left, right/4, 0, pivot, pthread_num),
    SegArgs(v, right/4 + 1, right/2, 1, pivot, pthread_num),
    SegArgs(v, right/2 + 1, 3*right/4, 2, pivot, pthread_num),
    SegArgs(v, 3*right/4 + 1, right, 3, pivot, pthread_num)
    };
    pthread_barrier_init(&barrier, NULL, pthread_num);
    for (int i = 0; i < pthread_num; i++) {
        ids[i] = i;
        //cout<<i<<"k";
        //fflush(stdout);
        int len = right - left + 1;
        int left_sub = left + i*((int) len/pthread_num);
        int right_sub = left_sub + ((int) len/pthread_num) - 1;
        
        //SegArgs part = {v,left_sub,right_sub,i,pivot,pthread_num};
        
        //cout<<&len<<" ";
        //fflush(stdout);
        pthread_create(&threads[i], NULL, sumSublist, (void *) &part[i]);
    }

    // wait for pthreads to finish
    for (int i = 0; i < pthread_num; i++) {

        pthread_join(threads[i], NULL);
        //cout<<i<<"kk";
        //fflush(stdout);
    }
    // Partition around pivot
    int i = left;
    int j = right;
    //cout<<i<<"kk";
    fflush(stdout);
    pthread_mutex_lock(&mutex);
    while (i<=j){
        if(v[i]<pivot){
            //cout<<zeuss[i]<<" ";
            v_copy[zeuss[i]-1] = v[i];
        }
        if(v[i]>pivot){
            //cout<<posiedon[i]<<" ";
            v_copy[posiedon[i]-1] = v[i];
        }
        i++;
    }
    v_copy[zeuss[i]] = pivot;
    i=left;
    j=right;
    while (i<=j){
        v[i] = v_copy[i];
        i++;
    }
    pthread_mutex_unlock(&mutex);
    i = left;
    j = right;
    //cout<<i<<"kk";
    //fflush(stdout);
    //std::cout<<i;
    // Recursively quicksort left and right sublists
    if (left < j && j - left > CUTOFF) {
        pthread_t left_thread;
        ThreadArgs left_args = { v, left, j };
        pthread_create(&left_thread, NULL, &partition, &left_args);
        pthread_join(left_thread, NULL);
    } else {
        sort(v.begin() + left, v.begin() + j + 1);
    }

    if (i < right && right - i > CUTOFF) {
        pthread_t right_thread;
        ThreadArgs right_args = { v, i, right };
        pthread_create(&right_thread, NULL, &partition, &right_args);
        pthread_join(right_thread, NULL);
    } else {
        sort(v.begin() + i, v.begin() + right + 1);
    }

    return NULL;
}

void quicksort(vector<int>& v) {
    ThreadArgs args = { v, 0, v.size() - 1 };
    //print(v_copy);
    pthread_t root_thread;
    pthread_create(&root_thread, NULL, &partition, &args);
    pthread_join(root_thread, NULL);
}

int main() {
    vector<int> v ;
    int n = 128; // size of the large list
    cout<<n;
    //quicksort(v);
    fflush(stdout);
    int i;
    for (i = 1; i <= n; i++)
        v.push_back(rand() % 100);

    for (i = 0; i < v.size(); i++) {
        cout << v[i] << " ";
    }    
    zeuss.resize(n);
    posiedon.resize(n);
    zerone.resize(n);
    sums.resize(n);
    sums2.resize(n);
    onezero.resize(n);
    v_copy.resize(n);
    for (i = 0; i < n; i++) {
        sums[i] = 0;
        zerone[i] = 0;
        zeuss[i] = 0;
        onezero[i] = 0;
        posiedon[i] = 0;
        sums2[i] = 0;
        v_copy[i] = 0;
    }
    //std::cout<<n;
    //fflush(stdout);
    quicksort(v);
    // initialize pthreads
    fflush(stdout);
    cout << endl;
    return 0;
}
