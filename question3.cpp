#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sys/time.h>
using namespace std;
const int MAXSIZE = 1000000;
// Recursive function to sort a local list of elements
int quicksort(vector<int>& local_list, int my_rank, int comm_size, MPI_Comm Comm) {
    int n = local_list.size();
    if (comm_size == 1) {

        sort(local_list.begin(), local_list.end());
        return local_list.size();
    }
    // Select pivot element as the middle element
    
    int pivot;
    if (my_rank == 0) {
        // Select pivot element as the middle element
        pivot = local_list[n / 2];
        // Broadcast pivot to all other processes
        MPI_Bcast(&pivot, 1, MPI_INT, 0, Comm);
        //cout << pivot <<"Pivot" << endl;
    }
    else {
        // Receive pivot from root process
        MPI_Bcast(&pivot, 1, MPI_INT, 0, Comm);
    }
    MPI_Barrier;
    if (comm_size % 2 != 0 && my_rank == comm_size / 2) {
        sort(local_list.begin(), local_list.end());
        return local_list.size();
    }
    vector<int> less, greater;
    // Split local list into two lists, one less than or equal to pivot, and the other greater than pivot
    for (int i = 0; i < n; i++) {
        if (local_list[i] <= pivot) {
            less.push_back(local_list[i]);
        }
        else {
            greater.push_back(local_list[i]);
        }
    }

    // Determine which processor to communicate with based on my_rank
    int partner;
    if (my_rank < comm_size / 2) {
        partner = my_rank + comm_size / 2;
    }
    else {
        partner = my_rank - comm_size / 2;
    }
    //for (int i = 0; i < n ; i++) {
    //    cout << local_list[i] << " ";
    //}
    //cout <<"My rank" << my_rank << endl;
    int color, new_rank, new_size;
    MPI_Comm new_comm;
    // Send and receive lists between processors in the lower and upper half of the ensemble
    MPI_Status status;
    MPI_Barrier;
    if (my_rank < comm_size / 2) {
        int result = MPI_Send(&greater[0], greater.size(), MPI_INT, partner, 0, Comm);
        vector<int> received(MAXSIZE);
        //cout << "neesan 1" << "  " << "partner "<<partner<<" rank "<< my_rank<<endl;
        if (result != MPI_SUCCESS) {
            char error_string[MPI_MAX_ERROR_STRING];
            int length;
            MPI_Error_string(result, error_string, &length);
            fprintf(stderr, "%s\n", error_string);
            // handle the error
        }
        result = MPI_Recv(&received[0], MAXSIZE, MPI_INT, partner, 0, Comm, &status);
        //cout << "neesan 2" <<"  " << "partner " << partner << " rank " << my_rank << endl;
        if (result != MPI_SUCCESS) {
            char error_string[MPI_MAX_ERROR_STRING];
            int length;
            MPI_Error_string(result, error_string, &length);
            fprintf(stderr, "%s\n", error_string);
            // handle the error
        }
        int received_count;
        MPI_Get_count(&status, MPI_INT, &received_count);
        less.insert(less.end(), received.begin(), received.begin() + received_count);

        local_list.assign(less.begin(), less.end());
        color = 0;
    }
    else {
        vector<int> received(MAXSIZE);
        //cout << "break 1" << "  " << "partner " << partner << " rank " << my_rank << endl;
        MPI_Recv(&received[0], MAXSIZE, MPI_INT, partner, 0, Comm, &status);
        //cout << "break 2" << "  " << "partner " << partner << " rank " << my_rank << endl;
        if (status.MPI_SOURCE == partner) {
            MPI_Send(&less[0], less.size(), MPI_INT, partner, 0, Comm);
        }
        //MPI_Recv(&received[0], less.size(), MPI_INT, partner, 0, MPI_COMM_WORLD, &status);
        //cout << "break 3" << "  " << "partner " << partner << " rank " << my_rank   << endl;
        int received_count;
        MPI_Get_count(&status, MPI_INT, &received_count);
        greater.insert(greater.end(), received.begin(), received.begin() + received_count);
        local_list.assign(greater.begin(), greater.end());
        color = 1;
    }
    // Recursively sort the less and greater lists
    //quicksort(less, my_rank, comm_size);
    //quicksort(greater, my_rank, comm_size);
    // Copy the sorted elements back into the local list
    //copy(less.begin(), less.end(), local_list.begin());
    //copy(greater.begin(), greater.end(), local_list.begin() + less.size());
    //cout << "I'm done" << my_rank << endl;
    MPI_Comm_split(Comm, color, my_rank, &new_comm);

    MPI_Comm_rank(new_comm, &new_rank);
    MPI_Comm_size(new_comm, &new_size);

    int mid = quicksort(local_list, new_rank, new_size, new_comm);

    //for (int i = 0; i < local_list.size(); i++) {
    //    cout << local_list[i] << " ";
    //}
    //cout <<"My rank" << my_rank << endl;
    return local_list.size();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int comm_size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // Generate random local list of elements
    int n = atoi(argv[1]);
    vector<int> local_list(n);
    srand(my_rank);
    for (int i = 0; i < n; i++) {
        local_list[i] = rand() % 100;
    }
    // Sort the local list using message passing quicksort algorithm
    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);
    int local = quicksort(local_list, my_rank, comm_size, MPI_COMM_WORLD);
    vector<int> recvcounts(comm_size);
    vector<int> displs(comm_size);
    // Gather the sorted local lists to the root process
    MPI_Gather(&local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    int total_recvcount = 0;
    if (my_rank == 0) {
        for (int i = 0; i < comm_size; i++) {
            displs[i] = total_recvcount;
            total_recvcount += recvcounts[i];
        }
    }
    vector<int> sorted_list(n * comm_size);
    //cout << "Im done quickly" << my_rank<<endl;
    //MPI_Gather(&local_list[0], n, MPI_INT, &sorted_list[0], n, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&local_list[0], local, MPI_INT, &sorted_list[0], &recvcounts[0], &displs[0], MPI_INT, 0, MPI_COMM_WORLD);
    // Print the sorted list in the root process
    gettimeofday(&end_time, NULL);
    if (my_rank == 0) {
        //sort(sorted_list.begin(), sorted_list.end());
        double time_taken = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

        cout << "Time taken by parallel: " << time_taken << " seconds." << endl;
        //cout << endl;
        vector<int> serial_list(n * 4);
        srand(my_rank);
        for (int i = 0; i < n * 4; i++) {
            serial_list[i] = rand() % 100;
        }
        gettimeofday(&start_time, NULL);
        sort(serial_list.begin(), serial_list.end());
        gettimeofday(&end_time, NULL);
        double time_taken2 = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

        cout << "Time taken by serial: " << time_taken2 << " seconds." << endl;
        cout << "Efficiency: " << time_taken2 / (time_taken * comm_size) << endl;
        cout << "Speed up: " << time_taken2 / time_taken << endl;
        
    }
    MPI_Finalize();
    return 0;
}

