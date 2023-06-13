# Scalable-Sorting
This repository is an implementation of the parallelized version of quicksort (using posix threads and MPI) and sample sort (using MPI).  The sorting algorithms implemented scale w.r.t the problem size and number of machines.

# Parallel Sorting Algorithms

This repository contains implementations of parallel sorting algorithms using different parallel programming models. The implemented algorithms include QuickSort with MPI and Pthreads, as well as SampleSort with MPI.

## QuickSort with MPI and Pthreads

The QuickSort algorithm is a well-known sorting algorithm that works by partitioning the array around a chosen pivot element. In this implementation, we have leveraged the power of MPI (Message Passing Interface) and Pthreads (POSIX threads) to achieve parallelism.

To compile and run the QuickSort implementation, follow these steps:

1. Ensure you have MPI and Pthreads libraries installed on your system.
2. Clone this repository and navigate to the QuickSort-MPI file(same for pthread).
3. To compile an OpenMPI program, you can use mpicc or mpic++ instead of gcc. Refer https://curc.readthedocs.io/en/latest/programming/MPI-C.html
4. Run the program with the desired number of MPI processes and threads, note that the input array size is hardcoded.

## SampleSort with MPI

SampleSort is another parallel sorting algorithm that divides the input array into equal-sized blocks and sorts them individually using MPI processes. It then selects a set of samples from each process and performs a global sample selection to determine the final split points for partitioning the array.

To compile and run the SampleSort implementation, follow these steps:

1. Ensure you have MPI library installed on your system.
2. Clone this repository and navigate to the SampleSort-MPI directory.
3. To compile an OpenMPI program, you can use mpicc or mpic++ instead of gcc. Refer https://curc.readthedocs.io/en/latest/programming/MPI-C.html
4. Run the program with the desired number of MPI processes and threads, note that the input array size is hardcoded.

## Usage and Customization

You can customize the input array or file to sort by modifying the provided code. Feel free to experiment and explore different aspects of the implementations.

Please note that these implementations serve as educational examples and may not be optimized for large-scale production environments. They are intended to showcase parallel programming concepts and algorithms.

---

Feel free to modify and enhance the README file according to your specific repository structure and implementation details.
