#ifndef DEFINE_MNIST
#define DEFINE_MNIST
#include "Parameter.h"

void MNISTLoader(const std::string& img_path, double**** img_data);
void FreeMnistImgMemory(double*** img_data, int Nimg, int Nrow, int Ncol);

#endif
