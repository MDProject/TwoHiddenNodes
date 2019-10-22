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

int MnistImageLoader(const std::string& img_path, double**** img_data) {
	int NumOfImg = 0;
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
		NumOfImg = number_of_images;
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
	file.close();
	return NumOfImg;
}

void AllocateLabelMemory(unsigned int** label_data, int Nimg) {
	(*label_data) = (unsigned int*)malloc(Nimg * sizeof(unsigned int));
}

void FreeMnistLabelMemory(unsigned int* label_data) {
	free(label_data);
}

void MnistLabelLoader(const std::string& img_path, unsigned int** label_data) {
	std::ifstream file(img_path, std::ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		file.read((char*)&magic_number, sizeof(magic_number));
		magic_number = ReverseInt(magic_number);
		file.read((char*)&number_of_images, sizeof(number_of_images));
		number_of_images = ReverseInt(number_of_images);
		char label;
		AllocateLabelMemory(label_data, number_of_images);
		for (int i = 0; i < number_of_images; i++)
		{
			file.read(&label, sizeof(label));
			(*label_data)[i] = (unsigned int)label;
		}
	}
	else {
		std::cout << "Fault file path input!" << std::endl;
		system("pause");
	}
	file.close();
}

void AllocateAssignDataMemory(double*** assign_data, int Nimg, int Nv) {
	(*assign_data) = (double**)malloc(Nimg * sizeof(double*));
	for (int n = 0; n < Nimg; n++) {
		(*assign_data)[n] = (double*)calloc(Nv, sizeof(double));
	}
}

int ExtractAssignedData(unsigned int* assign_label, int n_label, double*** img_data, unsigned int* label_data, int Nimg, int Nrow, int Ncol, double*** assign_data) {
	int Nimg_assign = 0;
	for (int n = 0; n < Nimg; n++) {
		for (int ac = 0; ac < n_label; ac++) {
			if (label_data[n] == assign_label[ac]) {
				Nimg_assign++;
			}
		}
	}
	AllocateAssignDataMemory(assign_data, Nimg_assign, Nrow*Ncol);
	int img_assign_idx = 0;
	for (int n = 0; n < Nimg; n++) {
		for (int ac = 0; ac < n_label; ac++) {
			if (label_data[n] == assign_label[ac]) {
				for (int n_row = 0; n_row < Nrow; n_row++) {
					for (int n_col = 0; n_col < Ncol; n_col++) {
						(*assign_data)[img_assign_idx][n_row*Nrow + n_col] = img_data[n][n_row][n_col];
					}
				}
				img_assign_idx++;
			}
		}
	}
	return Nimg_assign;
}

void FreeAssignedDataMemory(double** assign_data, int Nimg) {
	for (int n = 0; n < Nimg; n++) {
		free(assign_data[n]);
	}
	free(assign_data);
}