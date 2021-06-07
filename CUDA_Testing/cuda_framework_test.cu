//Simple test program to test to see if CUDA is up an running

#include <iostream>

using std::cout;
using std::endl;

const int ARRAY_SIZE = 1000;

__global__ void simpleKenrel(double *space){
    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if(threadId < ARRAY_SIZE){
        space[threadId] *= 2;
    }
    return;
}

int main(){

    // display GPU properties and ensure we are using the right one
    cudaDeviceProp prop;
    int numDevices;
    cudaGetDeviceCount(&numDevices);
    cout << "Number of devices detected: " << numDevices << endl;

    for(int dev = 0; dev < numDevices; dev++ ){
        cudaGetDeviceProperties(&prop, 0);
        cout << "\nDevice Number: " << dev << " \n";
        cout << "- Device name: " << prop.name << endl;
        cout << "- Device compute: " << prop.major << "." << prop.minor << endl << endl;
    }

    //Create some arrays to fill with data
    double *host_input;
    double *device_space;
    double *host_output;

    //Allocate Memory
    cudaMalloc((void**) &device_space, ARRAY_SIZE * sizeof(double));
    host_input = new double[ARRAY_SIZE];
    host_output = new double[ARRAY_SIZE];

    //Filling input: start with 3s
    for(int i = 0; i < ARRAY_SIZE; i++ ){
        host_input[i] = 3;
        host_output[i] = -1; //Should be overwritten with device_space
    }

    //Copy memory from host to device
    cudaMemcpy(device_space, host_input, ARRAY_SIZE * sizeof(double), cudaMemcpyHostToDevice);

    simpleKenrel<<<(ARRAY_SIZE+32-1)/32,32>>>(device_space);

    //Copy things off of device
    cudaMemcpy(host_output, device_space, ARRAY_SIZE * sizeof(double), cudaMemcpyDeviceToHost);

    bool check = true;
    for(int i = 0; i < ARRAY_SIZE; i++){
        //Visual inspection:
        //cout << host_output[i] << " ";
        if(host_output[i] != 6){
            check = false;
        }
    }

    if(check){
        cout << "Output is good!" << endl;
    } else {
        cout << "Output is bad!" << endl;
    }

    cout << endl;



    return 0;
}