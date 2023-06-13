#include <pthread.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <sys/time.h>


#define CUTOFF 1000

using namespace std;

struct ThreadArgs {
    vector<int> &v;
    int left;
    int right;
};

void* partition(void *args) {
    ThreadArgs *ta = (ThreadArgs*)args;
    vector<int> &v = ta->v;
    int left = ta->left;
    int right = ta->right;

    // Choose pivot element
    int pivot_index = (left + right) / 2;
    int pivot = v[pivot_index];

    // Partition around pivot
    int i = left;
    int j = right;
    while (i <= j) {
        while (v[i] < pivot) {
            i++;
        }
        while (v[j] > pivot) {
            j--;
        }
        if (i <= j) {
            swap(v[i], v[j]);
            i++;
            j--;
        }
    }

    // Recursively quicksort left and right sublists
    if (left < j && j - left > CUTOFF) {
        pthread_t left_thread;
        ThreadArgs left_args = { v, left, j };
        pthread_create(&left_thread, NULL, partition, &left_args);
        pthread_join(left_thread, NULL);
    } else {
        sort(v.begin() + left, v.begin() + j + 1);
    }

    if (i < right && right - i > CUTOFF) {
        pthread_t right_thread;
        ThreadArgs right_args = { v, i, right };
        pthread_create(&right_thread, NULL, partition, &right_args);
        pthread_join(right_thread, NULL);
    } else {
        sort(v.begin() + i, v.begin() + right + 1);
    }

    return NULL;
}

void quicksort(vector<int>& v) {
    ThreadArgs args = { v, 0, v.size() - 1 };
    pthread_t root_thread;
    pthread_create(&root_thread, NULL, partition, &args);
    pthread_join(root_thread, NULL);
}

int partition1(int* arr, int low, int high) {
    int pivot = arr[high];
    int i = low - 1;
    for (int j = low; j < high; j++) {
        if (arr[j] <= pivot) {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    i++;
    swap(arr[i], arr[high]);
    return i;
}

void serial_quicksort(int* arr, int low, int high) {
    if (low < high) {
        int p = partition1(arr, low, high);
        serial_quicksort(arr, low,  p-1);
        serial_quicksort(arr, p+1, high);
    }
}



int main() {
    vector<int> v ;
    int n[5] = {100000,200000,500000,1000000,2000000}; // size of the large list
    //cout<<n;
    //quicksort(v);
    for (int j = 0; j <= 4; j++){

        vector<int> v ;
        timeval start, end, start1, end1;
        for (int i = 0; i <= n[j]; i++)
            v.push_back(rand() % 100);
        gettimeofday(&start, NULL);
        quicksort(v);
        //for (int i = 0; i < v.size(); i++) {
         //   cout << v[i] << " ";
        //}
        gettimeofday(&end, NULL);

        double total_time_taken = (end.tv_usec - start.tv_usec)*1e-6 + (end.tv_sec - start.tv_sec);
        vector<int> v1 ;
        for (int i = 0; i <= n[j]; i++)
            v1.push_back(rand() % 100);
        int *arr = (int*)malloc(n[j] * sizeof(int));;
     
    // Using copy() function
    // to copy elements
        copy(v1.begin(),v1.end(),arr);
        gettimeofday(&start1, NULL);
        serial_quicksort(arr, 0, n[j]-1);
        gettimeofday(&end1, NULL);
        double total_time_taken1 = (end1.tv_usec - start1.tv_usec)*1e-6 + (end1.tv_sec - start1.tv_sec);
        cout<<n[j]<<"\t\t"<<total_time_taken1/total_time_taken<<endl;
    }
    cout << endl;
    return 0;
}
