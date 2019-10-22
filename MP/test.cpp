#include "MNIST.h"
#include <fstream>



int main() {
	std::string img_path = "D:\\DeepLearning\\Application\\MNIST_DATA\\MNIST\\raw\\train-images-idx3-ubyte";
	double*** img_data;
	MNISTLoader(img_path, &img_data);
	FreeMnistImgMemory(img_data, 60000, 28, 28);
	system("pause");
}