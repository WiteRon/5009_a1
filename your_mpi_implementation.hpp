//     [NAME]: LU JINCHENG [STUDENT ID]: 21121730
//     [EMAIL]: jlucf@connect.ust.hk

// *********************************************************************
//     NOTICE: Write your code only in this file and only submit this 
//     file to Canvas. You can add more functions, structures, classes,
//     and include other C++ standard libraries in this file if needed.
//     Do NOT change the signature of the function gen_um_lists_MPI, 
//     which will be called by the main function.
// *********************************************************************
#pragma once

#include <string>
#include <vector>
#include "utilities.hpp"
#include "mpi.h"

using namespace std;

/// @brief Generate universal minimizer lists for each read in parallel using MPI.
/// @param comm The MPI communicator.
/// @param my_rank Global rank of this process.
/// @param num_process Total number of processes.
/// @param k length of k-mer (only available in Process 0 when calling).
/// @param n number of reads (only available in Process 0 when calling).
/// @param reads The input reads in vector of strings (only available in Process 0 when calling).
/// @param reads_CSR The input reads in CSR format (only available in Process 0 when calling).
/// @param reads_CSR_offs The offsets of the input reads in CSR format (only available in Process 0 when calling).
/// @param um_lists The output lists of universal minimizers (result should be gathered to this vector in Process 0 when function ends).
void gen_um_lists_MPI (MPI_Comm comm, int my_rank, int num_process, int k, int n, vector<string>& reads, char* reads_CSR, int* reads_CSR_offs, vector<vector<kmer_t>>& um_lists) {
    int local_n;
    int process_used;
    char* local_reads_CSR;
    int* local_reads_CSR_offs;
    vector<vector<kmer_t>> local_um_lists;
    MPI_Bcast(&k,1,MPI_INT,0,comm);
    if (my_rank == 0) {
        // printf("%s\n","============");
        int ave_size = (reads_CSR_offs[n] / num_process) + 1;
        int offs = reads_CSR_offs[0];
        int offs_index = 0;
        int next_slice = ave_size;
        local_reads_CSR = reads_CSR;
        // printf("%d\n",ave_size);
        int p = 1;
        for(p = 1; offs_index < n;p++) {
            int pre_offs_index = offs_index;
            // printf("%s\n","big for loop");
            while(offs_index < n && reads_CSR_offs[offs_index] < next_slice) {
                offs_index++;
                // printf("%s\n","while loop");
            }
            local_n = offs_index - pre_offs_index;
            local_reads_CSR_offs = (int*) malloc((local_n + 1) * sizeof(int));
            for(int i = 0; i < local_n + 1; i++) {
                local_reads_CSR_offs[i] = reads_CSR_offs[pre_offs_index + i] - reads_CSR_offs[pre_offs_index];
                // printf("%s\n","inner for loop");
            }
            if(offs_index < n) {
                MPI_Send(&local_n,1,MPI_INT,p,0,comm);
                MPI_Send(local_reads_CSR_offs,local_n + 1,MPI_INT,p,0,comm);
                MPI_Send(local_reads_CSR,reads_CSR_offs[offs_index] - reads_CSR_offs[pre_offs_index],MPI_CHAR,p,0,comm);
                local_reads_CSR = reads_CSR + sizeof(char)*reads_CSR_offs[offs_index];
                next_slice = reads_CSR_offs[offs_index] + ave_size;
            }
        }
        p--;
        process_used = p;
        while(p < num_process) {
            int temp_none = -1;
            MPI_Send(&temp_none,1,MPI_INT,p,0,comm);
            p++;
        }
    } else {
        MPI_Recv(&local_n,1,MPI_INT,0,0,comm,MPI_STATUS_IGNORE);
        if(local_n != -1) {
            local_reads_CSR_offs = new int[local_n + 1];
            MPI_Recv(local_reads_CSR_offs,local_n + 1,MPI_INT,0,0,comm,MPI_STATUS_IGNORE);
            local_reads_CSR = new char[local_reads_CSR_offs[local_n]];
            MPI_Recv(local_reads_CSR,local_reads_CSR_offs[local_n],MPI_CHAR,0,0,comm,MPI_STATUS_IGNORE);
        }
    }
    for(int i = 0; i < local_n; i++) {
        // printf("%d\n",local_n);
        local_um_lists.push_back(generate_universal_minimizer_list(k, local_reads_CSR_offs[i+1] - local_reads_CSR_offs[i],local_reads_CSR));
        local_reads_CSR += sizeof(char)*(local_reads_CSR_offs[i+1] - local_reads_CSR_offs[i]);
    }

    if(my_rank == 0) {
        for(int p = 1; p < process_used; p++) {
            vector<vector<kmer_t>> temp_um_lists;
            int lists_size;
            MPI_Recv(&lists_size,1,MPI_INT,p,0,comm,MPI_STATUS_IGNORE);
            for(int i = 0; i < lists_size; i++) {
                vector<kmer_t> temp_um_list;
                int list_size;
                MPI_Recv(&list_size,1,MPI_INT,p,0,comm,MPI_STATUS_IGNORE);
                unsigned* unsigned_um_list = (unsigned*) malloc(list_size * sizeof(unsigned));
                MPI_Recv(unsigned_um_list,list_size,MPI_UNSIGNED,p,0,comm,MPI_STATUS_IGNORE);
                for(int j = 0; j < list_size; j++) {
                    kmer_t um = (kmer_t) unsigned_um_list[j];
                    temp_um_list.push_back(um);
                }
                free(unsigned_um_list);
                temp_um_lists.push_back(temp_um_list);
            }
            um_lists.insert(um_lists.end(),temp_um_lists.begin(),temp_um_lists.end());
        }
        um_lists.insert(um_lists.end(),local_um_lists.begin(),local_um_lists.end());
    } else if(local_n != -1) {
        int lists_size = local_um_lists.size();
        MPI_Send(&lists_size,1,MPI_INT,0,0,comm);
        for(vector<kmer_t> um_list : local_um_lists) {
            int list_size = um_list.size();
            MPI_Send(&list_size,1,MPI_INT,0,0,comm);
            unsigned* array = um_list.data();
            MPI_Send(array,um_list.size(),MPI_UNSIGNED,0,0,comm);
        }
    }
    
}