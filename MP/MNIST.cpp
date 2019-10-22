#include "MNIST.h"

int ReverseInt(int i)
{
	unsigned char ch1, ch2, ch3, ch4;
	ch1 = i & 255;
	ch2 = (i >> 8) & 255;
	ch3 = (i >> 16) & 255;
	ch4 = (i >> 24) & 255;
	return((int)ch1 << 24) + ((int)ch2 << 16) + ((int)ch3 << 8) + ch4;
}

void AllocateDataMemory(double**** img_data, int Nimg, int Nrow, int Ncol) {
	(*img_data) = (double***)malloc(Nimg * sizeof(double**));
	for (int n = 0; n < Nimg; n++) {
		(*img_data)[n] = (double**)malloc(Nrow * sizeof(double*));
		for (int n_row = 0; n_row < Nrow; n_row++) {
			(*img_data)[n][n_row] = (double*)calloc(Ncol, sizeof(double));
		}
	}
}

void FreeMnistImgMemory(double*** img_data, int Nimg, int Nrow, int Ncol) {
	for (int n = 0; n < Nimg; n++) {
		for (int n_row = 0; n_row < Nrow; n_row++) {
			free(img_data[n][n_row]);
		}
		free(img_data[n]);
	}
	free(img_data);
}

void MNISTLoader(const std::string& img_path, double**** img_data) {
	std::ifstream file(img_path, std::ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		int n_rows = 0;
		int n_cols = 0;
		file.read((char*)&magic_number, sizeof(magic_number));
		magic_number = ReverseInt(magic_number);
		file.read((char*)&number_of_images, sizeof(number_of_images));
		number_of_images = ReverseInt(number_of_images);
		file.read((char*)&n_rows, sizeof(n_rows));
		n_rows = ReverseInt(n_rows);
		file.read((char*)&n_cols, sizeof(n_cols));
		n_cols = ReverseInt(n_cols);
		AllocateDataMemory(img_data, number_of_images, n_rows, n_cols);
		for (int i = 0; i < number_of_images; i++)
		{
			for (int r = 0; r < n_rows; r++)
			{
				for (int c = 0; c < n_cols; c++)
				{
					unsigned char temp = 0;
					file.read((char*)&temp, sizeof(temp));
					(*img_data)[i][r][c] = (double)temp;
				}
			}
		}
	}
	else {
		std::cout << "Fault file path input!" << std::endl;
		system("pause");
	}
}