#ifndef DEFINE_MNIST
#define DEFINE_MNIST
#include "Parameter.h"

// img_data[image index][row index][col index]
int MnistImageLoader(const std::string& img_path, double**** img_data);

void FreeMnistImgMemory(double*** img_data, int Nimg, int Nrow, int Ncol);

// label_data[image index]
void MnistLabelLoader(const std::string& img_path, unsigned int** label_data);

void FreeMnistLabelMemory(unsigned int* label_data);

// assign_label: stores given labels	 assign_data[image index][shrinked 1D image pixels]
void ExtractAssignedData(unsigned int* assign_label, int* num_list, int n_label, double*** img_data, unsigned int* label_data, int Nimg, int Nrow, int Ncol, double*** assign_data);

void FreeAssignedDataMemory(double** assign_data, int Nimg);

void CompressData2Binary(double** assign_data, int Nimg, int n_pixel, int threshold);

#endif
