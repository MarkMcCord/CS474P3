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

#include <math.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
# define M_PI           3.14159265358979323846  /* pi */

#include "image.h"

using namespace std;

void readImage(char fname[], ImageType &image);
void writeImage(char fname[], ImageType &image);
void fft(float data[], unsigned long nn, int isign);

void part1(bool b);

void fft2d(float **data, int size, unsigned long nn, int isign);
void imgtodata(ImageType &image, float **data);
void datatoimg(ImageType &image, float **data);

int main(int argc, char *argv[])
{

	char lenna[]   = "lenna.pgm";
/*
	//1.a 
	float test[] = {0, 2, 0, 3, 0, 4, 0, 4, 0};
	fft(test, 4, -1);
	//fft(test, 4, 1);
	for (int i = 1; i <= 8; i = i + 2){
		test[i] = test[i] / 4;
		test[i+1] = test[i+1] / 4;
		cout << test[i] << " " << test[i+1] << endl;
	}
*/
	//1.b and 1.c
	part1(true);
	part1(false);

	// Part 2
	ImageType finImage(256, 256, 255);
	readImage(lenna, finImage);
	float **data = new float *[256];

	for (int i = 0; i < 256; i++)
	{
		data[i] = new float[513];
	}

	imgtodata(finImage, data);
	//cout << "forward" << endl << endl;
	fft2d(data, 256, 256, -1);
	//cout << "inverse" << endl << endl;
	fft2d(data, 256, 256, 1);
	datatoimg(finImage, data);

	char fname[] = "part2.pgm";
	writeImage(fname, finImage);

	for (int i = 0; i < 256; ++i)
	{
		delete[] data[i];
	}

	delete[] data;
}

void imgtodata(ImageType &image, float **data)
{
	int temp = 0;
	for (int i = 0; i < 256; i++)
	{
		for (int j = 1; j < 513; j += 2)
		{
			image.getPixelVal(i, (j - 1) / 2, temp);
			data[i][j] = temp;
			data[i][j + 1] = 0;
		}
	}
}

void datatoimg(ImageType &image, float **data)
{
	int temp = 0;
	for (int i = 0; i < 256; i++)
	{
		for (int j = 1; j < 513; j += 2)
		{
			temp = data[i][j];
			image.setPixelVal(i, (j - 1) / 2, temp);
		}
	}
}

/* 

   The real part is stored in the odd locations of the array
   (data[1], data[3], data[5]. etc) and the imaginary part 
   in the even locations (data[2], data[4], data[6], etc. 

   The elements in the array data must be stored from data[1] 
   to data[2*nn] -  data[0] is not used! 
 
   nn is the length of the input which should be power of 2. 
   Warning: the program does not check this condition.

   isign: -1 Forward FFT, isign: 1  Inverse FFT (based on our definition)
   
   Warning: the FFT routine provided does not multiply by the normalization 
   factor 1/N that appears in the forward DFT equation; you should do this 
   yourself (see page 506 from the fft handout).

*/

void fft(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software 0#Y". */

void part1(bool b){
	float data [257];
	char part;
	if(b){ //1.b
		part = 'b';
		ofstream file("cosine.csv");
		for (int i = 1; i < 257; i = i + 2){
			data [i] = cos((2 * 8 * M_PI * (i / 2) / 128));
			data [i + 1] = 0;
			//cout << data[i] << " " << data[i + 1] << endl;
			file << data[i] << endl;
			data [i] = data[i]  * pow(-1, i/2);
		}
		file.close();
	}
	if(!b){ //1.c
		part = 'c';
		ofstream file("rect.csv");
		for (int i = 1; i < 65; i = i + 2){
			data[i] = 0;
			data[i + 1] = 0;
			//cout << data[i] << " " << data[i + 1] << endl;
			file << data[i] << endl;
		}
		for (int i = 65; i < 193; i = i + 2){
			data[i] = 1 * pow(-1, i/2);
			data[i + 1] = 0;
			//cout << data[i] << " " << data[i + 1] << endl;
			file << data[i] << endl;
		}
		for (int i = 193; i < 257; i = i + 2){
			data[i] = 0;
			data[i + 1] = 0;
			//cout << data[i] << " " << data[i + 1] << endl;
			file << data[i] << endl;
		}
		file.close();
	}

	char mag[] = "xmagnitude.csv";
	mag[0] = part;
	ofstream magfile(mag);
	char pha[] = "xphase.csv";
	pha[0] = part;
	ofstream phafile(pha);
	char rea[] = "xreal.csv";
	rea[0] = part;
	ofstream reafile(rea);
	char ima[] = "ximaginary.csv";
	ima[0] = part;
	ofstream imafile(ima);

	fft(data, 128, -1);

	for (int i = 1; i < 257; i = i + 2){
		data[i] = data[i] / 128;
		data[i + 1] = data[i + 1] / 128;
		reafile << data[i] << endl;
		imafile << data[i + 1] << endl;
		magfile << sqrt((data[i] * data[i]) + (data[i + 1] * data[i + 1]))<< endl;
		phafile << atan2(data[i + 1],  data[i])<< endl;
	}

	magfile.close();
	phafile.close();
	reafile.close();
	imafile.close();
}

void fft2d(float ** data, int size, unsigned long nn, int isign){
	//create dynamic 2d arrays
	float** colsdata = new float*[size];
	for(int i = 0; i < size; ++i){
    	colsdata[i] = new float[size * 2 + 1];
	}
	float** rowsdata = new float*[size];
	for(int i = 0; i < size; ++i){
    	rowsdata[i] = new float[size * 2 + 1];
	}
	//copy the image into a 2d array, flipped so the columns become the rows and vice versa
	for (int i = 0; i < size; i++){
		for (int j = 1; j < size * 2 + 1; j = j + 2){
			colsdata[i][j] = data[(j - 1) / 2][2 * i + 1];
			colsdata[i][j + 1] = data[(j - 1) / 2][2 * i + 1 + 1];
		}
	}
	//perform 1d fft on each column
	for (int i = 0; i < size; i++){
		fft(colsdata[i], size, isign);
	}
	//flip back to the way it should be so we can transform the rows
	for (int i = 0; i < size; i++){
		for (int j = 1; j < size * 2 + 1; j = j + 2){
			rowsdata[i][j] = colsdata[(j - 1) / 2][2 * i + 1];
			rowsdata[i][j + 1] = colsdata[(j - 1) / 2][2 * i + 1 + 1];
		}
	}
	//perform id fft on each row
	for (int i = 0; i < size; i++){
		fft(rowsdata[i], size, isign);
	}
	//copy to original array and do 1/N
	for (int i = 0; i < size; i++){
		for (int j = 1; j < size * 2 + 1; j = j + 2){
			data[i][j] = rowsdata[i][j] / size;
			data[i][j + 1] = rowsdata[i][j + 1] / size;
			//cout << data[i][j] << " " << data[i][j + 1] << endl;
		}
	}
	//delete the dynamic 2d arrays
	for(int i = 0; i < size; ++i){
    	delete[] colsdata[i];
	}
	delete[] colsdata;
	for(int i = 0; i < size; ++i){
    	delete[] rowsdata[i];
	}
	delete[] rowsdata;

}