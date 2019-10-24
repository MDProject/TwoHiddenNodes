#include "MNIST.h"
#include <fstream>



int main() {
	std::string img_path = "D:\\DeepLearning\\Application\\MNIST_DATA\\MNIST\\raw\\train-images-idx3-ubyte";
	std::string idx_path = "D:\\DeepLearning\\Application\\MNIST_DATA\\MNIST\\raw\\train-labels-idx1-ubyte";
	double*** img_data, ** assign_data;
	unsigned int* label_data;
	unsigned int assign_label[] = { 1,3,4 };
	int num_list[] = { 2000,1000,1000 };
	int Nimg = MnistImageLoader(img_path, &img_data);
	MnistLabelLoader(idx_path, &label_data);
	ExtractAssignedData(assign_label, num_list, 3, img_data, label_data, Nimg, 28, 28, &assign_data);
	FreeMnistImgMemory(img_data, 60000, 28, 28);
	FreeMnistLabelMemory(label_data);
	FreeAssignedDataMemory(assign_data, 4000);
	system("pause");
}