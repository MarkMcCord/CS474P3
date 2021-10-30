#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "image.h"

using namespace std;

void readImage(char fname[], ImageType &image);
void writeImage(char fname[], ImageType &image);

ImageType padding(char fname[], int size);
void normalize(ImageType &image, vector<int> calcStore);
void salt(char fname[], int percent);

void smoothing(char fname[], int size, bool G);
void median(char fname[], int size, int percent);
void correlation(char image[], char pattern[]);
void unsharp(char fname[]);
void highboost(char fname[]);
void gradient(char fname[], bool p);
void laplacian(char fname[]);

int main(int argc, char *argv[]) {
	char boat[]    = "boat.pgm";
	char f_16[]    = "f_16.pgm";
	char lenna[]   = "lenna.pgm";
	char peppers[] = "peppers.pgm";
	char sf[]      = "sf.pgm";
	char image[]   = "Image.pgm";
	char pattern[] = "Pattern.pgm";

	// Part 1
	cout << "Part 1 start" << endl;
	correlation(image, pattern);

	// Part 2
	cout << "Part 2 start" << endl;
	smoothing(lenna, 7, 0);
	smoothing(lenna, 15, 0);
	smoothing(lenna, 7, 1);
	smoothing(lenna, 15, 1);

	smoothing(sf, 7, 0);
	smoothing(sf, 15, 0);
	smoothing(sf, 7, 1);
	smoothing(sf, 15, 1);

	// Part 3
	cout << "Part 3 start" << endl;
	salt(lenna, 50);
	char salted_l50[] = "salted_50_l.pgm";
	median(salted_l50, 7, 50);
	median(salted_l50, 15, 50);
	smoothing(salted_l50, 7, 0);
	smoothing(salted_l50, 15, 0);

	salt(boat, 50);
	char salted_b50[] = "salted_50_l.pgm";
	median(salted_b50, 7, 50);
	median(salted_b50, 15, 50);
	smoothing(salted_b50, 7, 0);
	smoothing(salted_b50, 15, 0);

	salt(lenna, 30);
	char salted_l30[] = "salted_30_l.pgm";
	median(salted_l30, 7, 30);
	median(salted_l30, 15, 30);
	smoothing(salted_l30, 7, 0);
	smoothing(salted_l30, 15, 0);

	salt(boat, 30);
	char salted_b30[] = "salted_30_l.pgm";
	median(salted_b30, 7, 30);
	median(salted_b30, 15, 30);
	smoothing(salted_b30, 7, 0);
	smoothing(salted_b30, 15, 0);

	// Part 4
	cout << "Part 4 start" << endl;
	unsharp(lenna);
	unsharp(sf);
	highboost(lenna);
	highboost(sf);

	// Part 5
	cout << "Part 5 start" << endl;
	gradient(lenna, 1);
	gradient(lenna, 0);
	gradient(sf, 1);
	gradient(sf, 0);
	laplacian(lenna);
	laplacian(sf);

	return 0;
}

ImageType padding(char fname[], int size) {
	// size is the mask size, not the number of pads you need
	ImageType image(256, 256, 255);
	readImage(fname, image);
	ImageType newImage(256 + (size - 1), 256 + (size - 1), 255);
	int temp;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			image.getPixelVal(i, j, temp);
			newImage.setPixelVal(i + (size / 2), j + (size / 2), temp);
		}
	}
	return newImage;
}

void normalize(ImageType &image, vector<int> calcStore) {
	int min   = INT_MAX;
	int max   = INT_MIN;
	int value = 0;

	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			value = calcStore[i * 256 + j];

			if (value > max) { max = value; }
			if (value < min) min = value;
		}
	}

	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			double scaled_value = 255.0 * ((calcStore[i * 256 + j] - min) / (double) (max - min));
			image.setPixelVal(i, j, scaled_value);
		}
	}
}

static int prewittY[3][3] = {{-1, -1, -1}, {0, 0, 0}, {1, 1, 1}};

static int prewittX[3][3] = {{-1, 0, 1}, {-1, 0, 1}, {-1, 0, 1}};

static int sobelY[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

static int sobelX[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};

static int lap[3][3] = {{0, 1, 0}, {1, -4, 1}, {0, 1, 0}};

void laplacian(char fname[]) {
	ImageType baseImage = padding(fname, 3);
	ImageType lapImage(256, 256, 255);

	vector<int> calcStore;
	// Laplacian calculation
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {  // for each pixel in the original image
			int temp3 = 0;
			for (int k = i; k < i + 3; k++) {
				for (int l = j; l < j + 3; l++) {  // for each weight in the mask
					int temp1;
					baseImage.getPixelVal(k, l, temp1);
					temp3 += temp1 * lap[k - i][l - j];
				}
			}
			calcStore.push_back(temp3);
		}
	}

	// Store sharpened into lapImage after normalization
	normalize(lapImage, calcStore);

	vector<int> calcStore2;
	// // Calculate values of g
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {  // for each pixel in the original image
			int temp3 = 0;
			int temp1, temp2;
			int c = -1;
			baseImage.getPixelVal(i, j, temp1);
			lapImage.getPixelVal(i, j, temp2);
			temp3 += temp1 + (c * (temp2));
			calcStore2.push_back(temp3);
		}
	}

	ImageType finalImage(256, 256, 255);
	char finImage[] = "laplacian_i.pgm";
	finImage[10]    = finImage[0];
	// Store sharpened image into finalImage after normalization
	normalize(finalImage, calcStore2);
	writeImage(finImage, finalImage);
}

void gradient(char fname[], bool p) {
	// Initialize f
	ImageType baseImage = padding(fname, 3);
	ImageType horImage(256, 256, 255);
	ImageType verImage(256, 256, 255);

	int temp;
	int sumHor = 0;
	int sumVer = 0;
	int hmin   = 1000;
	int hmax   = 0;
	int vmin   = 1000;
	int vmax   = 0;

	if (p) {
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {  // for each pixel in the original image
				for (int k = i; k < i + 3; k++) {
					for (int l = j; l < j + 3; l++) {  // for each weight in the mask
						baseImage.getPixelVal(k, l, temp);
						sumHor = sumHor + (temp * prewittX[k - i][l - j]);
						sumVer = sumVer + (temp * prewittY[k - i][l - j]);
						// cout << k << ", " << l << endl;
					}
				}
				// calculate sum
				horImage.setPixelVal(i, j, sumHor);
				verImage.setPixelVal(i, j, sumVer);
				// cout << "Pixel " << i << ", " << j << endl;

				if (sumHor > hmax) { hmax = sumHor; }
				if (sumHor < hmin) { hmin = sumHor; }
				if (sumVer > vmax) { vmax = sumVer; }
				if (sumVer < vmin) { vmin = sumVer; }

				sumHor = 0;
				sumVer = 0;
			}
		}
		char hor[] = "Hor_Prewitt_mask_i.pgm";
		hor[17]    = fname[0];
		char ver[] = "Ver_Prewitt_mask_i.pgm";
		ver[17]    = fname[0];

		double normalized;

		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {
				horImage.getPixelVal(i, j, temp);
				normalized = 255 * (temp - hmin) / (double) (hmax - hmin);
				horImage.setPixelVal(i, j, normalized);

				verImage.getPixelVal(i, j, temp);
				normalized = 255 * (temp - vmin) / (double) (vmax - vmin);
				verImage.setPixelVal(i, j, normalized);
			}
		}

		writeImage(hor, horImage);
		writeImage(ver, verImage);
	}

	if (!p) {
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {  // for each pixel in the original image
				for (int k = i; k < i + 3; k++) {
					for (int l = j; l < j + 3; l++) {  // for each weight in the mask
						baseImage.getPixelVal(k, l, temp);
						sumHor = sumHor + (temp * sobelX[k - i][l - j]);
						sumVer = sumVer + (temp * sobelY[k - i][l - j]);
						// cout << k << ", " << l << endl;
					}
				}
				// calculate sum
				horImage.setPixelVal(i, j, sumHor);
				verImage.setPixelVal(i, j, sumVer);
				// cout << "Pixel " << i << ", " << j << endl;

				if (sumHor > hmax) { hmax = sumHor; }
				if (sumHor < hmin) { hmin = sumHor; }
				if (sumVer > vmax) { vmax = sumVer; }
				if (sumVer < vmin) { vmin = sumVer; }

				sumHor = 0;
				sumVer = 0;
			}
		}
		char hor[] = "Hor_Prewitt_mask_i.pgm";
		hor[17]    = fname[0];
		char ver[] = "Ver_Prewitt_mask_i.pgm";
		ver[17]    = fname[0];

		double normalized;

		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {
				horImage.getPixelVal(i, j, temp);
				normalized = 255 * (temp - hmin) / (double) (hmax - hmin);
				horImage.setPixelVal(i, j, normalized);

				verImage.getPixelVal(i, j, temp);
				normalized = 255 * (temp - vmin) / (double) (vmax - vmin);
				verImage.setPixelVal(i, j, normalized);
			}
		}

		writeImage(hor, horImage);
		writeImage(ver, verImage);
	}

	double m;
	vector<int> calc;
	vector<int> calc2;
	int y;
	int x;
	ImageType Magnitude(256, 256, 255);

	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			horImage.getPixelVal(i, j, x);
			verImage.getPixelVal(i, j, y);
			m = x * x + y * y;
			m = sqrt(m);
			calc.push_back(m);
			Magnitude.setPixelVal(i, j, m);
		}
	}
	normalize(Magnitude, calc);
	if (p) {
		char mag[] = "Mag_Prewitt_mask_i.pgm";
		mag[17]    = fname[0];
		writeImage(mag, Magnitude);
	} else {
		char mag[] = "Mag_Sobel_mask_i.pgm";
		mag[15]    = fname[0];
		writeImage(mag, Magnitude);
	}

	ImageType Final(256, 256, 255);
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			baseImage.getPixelVal(i + 1, j + 1, x);
			Magnitude.getPixelVal(i, j, y);
			m = x + y;
			calc2.push_back(m);
			Final.setPixelVal(i, j, m);
		}
	}

	normalize(Final, calc2);
	if (p) {
		char final[] = "Prewitt_sharp_i.pgm";
		final[14]    = fname[0];
		writeImage(final, Final);
	} else {
		char final[] = "Sobel_sharp_i.pgm";
		final[12]    = fname[0];
		writeImage(final, Final);
	}
}

void unsharp(char fname[]) {
	// Initialize f
	ImageType baseImage(256, 256, 255);
	ImageType tempImage(256, 256, 255);
	readImage(fname, baseImage);

	// Initialize FLP
	smoothing(fname, 15, 1);
	char smoothedImage[] = "smooth_i_G_5.pgm";
	smoothedImage[7]     = fname[0];

	readImage(smoothedImage, tempImage);

	vector<int> calcStore;
	// Initialize and calculate Gmask
	ImageType maskImage(256, 256, 255);
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			int temp1, temp2, temp3;
			baseImage.getPixelVal(i, j, temp1);
			tempImage.getPixelVal(i, j, temp2);
			temp3 = temp1 - temp2;
			calcStore.push_back(temp3);
		}
	}
	// Create image after values have been normalized
	normalize(maskImage, calcStore);

	vector<int> calcStore2;
	int k = 1;
	// Initialize and calculate G
	ImageType finalImage(256, 256, 255);
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			int temp1, temp2, temp3;
			baseImage.getPixelVal(i, j, temp1);
			maskImage.getPixelVal(i, j, temp2);
			temp3 = (temp1 + (temp2 * k)) / 1 + k;
			calcStore2.push_back(temp3);
		}
	}

	// Create image after values have been normalized
	normalize(finalImage, calcStore2);

	char unsharpImage[] = "unsharp_i_G_5.pgm";
	unsharpImage[8]     = fname[0];
	writeImage(unsharpImage, finalImage);
}

void highboost(char fname[]) {
	// Initialize f
	ImageType baseImage(256, 256, 255);
	ImageType tempImage(256, 256, 255);
	readImage(fname, baseImage);

	// Initialize FLP
	smoothing(fname, 15, 1);
	char smoothedImage[] = "smooth_i_G_5.pgm";
	smoothedImage[7]     = fname[0];

	readImage(smoothedImage, tempImage);

	vector<int> calcStore;
	// Initialize and calculate Gmask
	ImageType maskImage(256, 256, 255);
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			int temp1, temp2, temp3;
			baseImage.getPixelVal(i, j, temp1);
			tempImage.getPixelVal(i, j, temp2);
			temp3 = temp1 - temp2;
			calcStore.push_back(temp3);
		}
	}
	// Create image after values have been normalized
	normalize(maskImage, calcStore);

	vector<int> calcStore2;
	int k = 2;
	// Initialize and calculate G
	ImageType finalImage(256, 256, 255);
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			int temp1, temp2, temp3;
			baseImage.getPixelVal(i, j, temp1);
			maskImage.getPixelVal(i, j, temp2);
			temp3 = (temp1 + (temp2 * k)) / 1 + k;
			calcStore2.push_back(temp3);
		}
	}

	// Create image after values have been normalized
	normalize(finalImage, calcStore2);

	char highboostImage[] = "highboost_i_G_5.pgm";
	highboostImage[10]    = fname[0];
	writeImage(highboostImage, finalImage);
}

void correlation(char image[], char pattern[]) {
	int rowsImage, colsImage;
	int rowsPattern, colsPattern;
	int levelsImage, levelsPattern;

	ImageType Image(442, 288, 255);
	ImageType Pattern(83, 55, 255);

	readImage(image, Image);
	readImage(pattern, Pattern);

	Image.getImageInfo(rowsImage, colsImage, levelsImage);
	Pattern.getImageInfo(rowsPattern, colsPattern, levelsPattern);

	ImageType fullImageF(442, 288, 255);
	char newfname[] = "correlatedImage.pgm";

	int temp1;
	int temp2;
	int sum = 0;
	int min = 1000000000;
	int max = 0;
	// iterate through image pixels
	for (int i = 0; i < rowsImage; i++) {
		for (int j = 0; j < colsImage; j++) {
			// iterate over mask
			for (int k = -(rowsPattern / 2); k < (rowsPattern / 2); k++) {
				for (int l = -(colsPattern / 2); l < (colsPattern / 2); l++) {
					if (i + k < 0 || i + k >= rowsImage || j + l < 0 || j + l >= colsImage) {
						sum = sum + 0;
					} else {
						Pattern.getPixelVal(k + (rowsPattern / 2), l + (colsPattern / 2), temp1);
						Image.getPixelVal(i + k, j + l, temp2);
						sum = sum + temp1 * temp2;
					}
				}
			}
			fullImageF.setPixelVal(i, j, sum);
			if (sum >= max) { max = sum; }
			if (sum <= min) { min = sum; }
			sum = 0;
		}
	}
	for (int i = 0; i < rowsImage; i++) {
		for (int j = 0; j < colsImage; j++) {
			fullImageF.getPixelVal(i, j, sum);
			double scaled = 255l * (sum - min) / (double) (max - min);
			fullImageF.setPixelVal(i, j, scaled);
		}
	}
	writeImage(newfname, fullImageF);
}

/*void correlation(char image[], char pattern[]) {
    int rowsImage, colsImage;
    int rowsPattern, colsPattern;
    int levelsImage, levelsPattern;

    ImageType fullImage(442, 288, 255);
    ImageType patternImage(83, 55, 255);

    readImage(image, fullImage);
    readImage(pattern, patternImage);

    fullImage.getImageInfo(rowsImage, colsImage, levelsImage);
    patternImage.getImageInfo(rowsPattern, colsPattern, levelsPattern);

    vector<int> corr;

    int temp;

    // iterate through image pixels
    for (int i = 0; i < rowsImage; i++) {
        for (int j = 0; j < colsImage; j++) {
            int sum   = 0;
            int temp2 = 0;
            int temp3 = 0;

            // iterate over mask
            for (int k = -rowsPattern / 2; k < rowsPattern / 2; k++) {
                for (int l = -colsPattern / 2; l < colsPattern / 2; l++) {
                    if (k < 0 || i + k >= rowsPattern || l < 0 || j + l >= colsPattern)
                        sum += 0;
                    else {
                        patternImage.getPixelVal(k, l, temp2);
                        fullImage.getPixelVal(k + i, j + l, temp3);
                        sum += temp2 * temp3;
                    }
                }
            }
            corr.push_back(sum);
        }
    }

    int temp4;

    int max = 0;
    int min = 1000000000;
    int value;

    // find min an max correlation values
    for (int i = 0; i < rowsImage; i++) {
        for (int j = 0; j < colsImage; j++) {
            value = corr[i * 288 + j];

            if (value > max) { max = value; }
            if (value < min) min = value;
        }
    }

    ImageType fullImageF(442, 288, 255);
    char newfname[] = "correlatedImage.pgm";
    for (int i = 0; i < rowsImage; i++) {
        for (int j = 0; j < colsImage; j++) { fullImageF.setPixelVal(i, j, 255 * corr[i * 288 + j] / (max - min)); }
    }

    writeImage(newfname, fullImageF);
}*/

void salt(char fname[], int percent) {
	char newfname[] = "salted_ii_i.pgm";
	newfname[10]    = fname[0];
	newfname[7]     = '0' + percent / 10;
	newfname[8]     = '0' + percent % 10;

	ImageType image(256, 256, 255);
	readImage(fname, image);
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			if (rand() % 100 + 1 <= percent) {
				if (rand() % 2 + 1 == 1) {
					image.setPixelVal(i, j, 0);
				} else {
					image.setPixelVal(i, j, 255);
				}
			}
		}
	}
	writeImage(newfname, image);
}

static int g7[7][7] = {{1, 1, 2, 2, 2, 1, 1}, {1, 2, 2, 4, 2, 2, 1}, {2, 2, 4, 8, 4, 2, 2}, {2, 4, 8, 16, 8, 4, 2},
                       {2, 2, 4, 8, 4, 2, 2}, {1, 2, 2, 4, 2, 2, 1}, {1, 1, 2, 2, 2, 1, 1}};

static int g15[15][15] = {
    {2, 2, 3, 4, 5, 5, 6, 6, 6, 5, 5, 4, 3, 2, 2},
    {2, 3, 4, 5, 7, 7, 8, 8, 8, 7, 7, 5, 4, 3, 2},
    {3, 4, 6, 7, 9, 10, 10, 11, 10, 10, 9, 7, 6, 4, 3},
    {4, 5, 7, 9, 10, 12, 13, 13, 13, 12, 10, 9, 7, 5, 4},
    {5, 7, 9, 11, 13, 14, 15, 16, 15, 14, 13, 11, 9, 7, 5},
    {5, 7, 10, 12, 14, 16, 17, 18, 17, 16, 14, 12, 10, 7, 5},
    {6, 8, 10, 13, 15, 17, 19, 19, 19, 17, 15, 13, 10, 8, 6},
    {6, 8, 11, 13, 16, 18, 19, 20, 19, 18, 16, 13, 11, 8, 6},
    {6, 8, 10, 13, 15, 17, 19, 19, 19, 17, 15, 13, 10, 8, 6},
    {5, 7, 10, 12, 14, 16, 17, 18, 17, 16, 14, 12, 10, 7, 5},
    {5, 7, 9, 11, 13, 14, 15, 16, 15, 14, 13, 11, 9, 7, 5},
    {4, 5, 7, 9, 10, 12, 13, 13, 13, 12, 10, 9, 7, 5, 4},
    {3, 4, 6, 7, 9, 10, 10, 11, 10, 10, 9, 7, 6, 4, 3},
    {2, 3, 4, 5, 7, 7, 8, 8, 8, 7, 7, 5, 4, 3, 2},
    {2, 2, 3, 4, 5, 5, 6, 6, 6, 5, 5, 4, 3, 2, 2},
};

void smoothing(char fname[], int size, bool G) {
	ImageType originalImage = padding(fname, size);
	ImageType smoothImage(256, 256, 255);

	char newfname[] = "smooth_i_i_i.pgm";
	newfname[7]     = fname[0];
	if (G) {
		newfname[9] = 'G';
	} else {
		newfname[9] = 'A';
	}
	newfname[11] = '0' + (size % 10);

	int sum    = 0;
	int factor = 0;
	int temp;
	if (!G) {
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {  // for each pixel in the original image
				for (int k = i; k < i + size; k++) {
					for (int l = j; l < j + size; l++) {  // for each weight in the mask
						originalImage.getPixelVal(k, l, temp);
						sum = sum + temp;
						// cout << k << ", " << l << endl;
					}
				}
				// calculate sum
				sum = sum / (size * size);
				smoothImage.setPixelVal(i, j, sum);
				// cout << "Pixel " << i << ", " << j << endl;
				sum = 0;
			}
		}
	} else if (size == 7) {
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {  // for each pixel in the original image
				for (int k = i; k < i + size; k++) {
					for (int l = j; l < j + size; l++) {  // for each weight in the mask
						originalImage.getPixelVal(k, l, temp);
						sum    = sum + (temp * g7[k - i][l - j]);
						factor = factor + g7[k - i][l - j];
						// cout << k << ", " << l << endl;
					}
				}
				// calculate sum
				sum = sum / factor;
				smoothImage.setPixelVal(i, j, sum);
				// cout << "Pixel " << i << ", " << j << endl;
				sum    = 0;
				factor = 0;
			}
		}
	} else if (size == 15) {
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {  // for each pixel in the original image
				for (int k = i; k < i + size; k++) {
					for (int l = j; l < j + size; l++) {  // for each weight in the mask
						originalImage.getPixelVal(k, l, temp);
						sum    = sum + (temp * g15[k - i][l - j]);
						factor = factor + g15[k - i][l - j];
						// cout << k << ", " << l << endl;
					}
				}
				// calculate sum
				sum = sum / factor;
				smoothImage.setPixelVal(i, j, sum);
				// cout << "Pixel " << i << ", " << j << endl;
				sum    = 0;
				factor = 0;
			}
		}
	}
	writeImage(newfname, smoothImage);
}

void median(char fname[], int size, int percent) {
	ImageType originalImage = padding(fname, size);
	ImageType medianImage(256, 256, 255);

	char newfname[] = "median_i_i_ii.pgm";
	newfname[7]     = fname[7];
	newfname[9]     = '0' + (size % 10);
	newfname[11]    = '0' + (percent / 10);
	newfname[12]    = '0' + (percent % 10);

	int array[size * size];
	int n = sizeof(array) / sizeof(array[0]);

	int factor = 0;
	int temp;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {  // for each pixel in the original image
			for (int k = i; k < i + size; k++) {
				for (int l = j; l < j + size; l++) {  // for each weight in the mask
					originalImage.getPixelVal(k, l, temp);
					array[((k - i) * size) + ((l - j) % size)] = temp;
				}
			}
			sort(array, array + n);
			temp = array[(size * size) / 2];
			medianImage.setPixelVal(i, j, temp);
			// cout << "Pixel " << i << ", " << j << endl;
			temp = 0;
		}
	}
	writeImage(newfname, medianImage);
}