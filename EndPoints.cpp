/*
 *  EndPoints.cpp
 *  
 *
 *  Created by Narendra Chaudhary on 2/27/2015.
 *  Texas A&M University.
 *
 */


#ifndef EndPoints
#include "EndPoints.h"
#endif

#include <ctime>

#ifdef WIN32
#define TIME_NORMALIZER 1000		// For Microsoft Windows Systems
#else
#define TIME_NORMALIZER 1000000		// For Mac OS X/Linux/Unix Systems
#endif


EndPoint::EndPoint()
{
	//passbool = NULL;
}

EndPoint::~EndPoint()
{
}

int EndPoint::Compression(string filename)
{ 
    
#ifndef DEBUG
	cout << filename << endl;
#else
	//cout << "File Name : " << filename << endl;
#endif

	// Start Compression ------------------------------------------------>
	unsigned int start, end;
	double elapsed;
	
	ClassEncoder original;

	start = clock();

	
	original.ReadPNG(filename);								// Read the layout image

	// Grab the image dimension
	width = original.width;
	height = original.height;
	bitdepth = original.bitdepth;
	row_bytes = original.row_bytes;

	cout << "Image Size = " << ((unsigned long long)original.width) * ((unsigned long long)original.height) << endl;
	
	original.Uncompressed_Beamlength = 0;					// store uncompressed beamlengths
	original.para_comptime = 0;

	//original.histograms = new unsigned int*[N*(N-1)];								//pointers for histograms of indiviual beams
	//unsigned int *hist = new unsigned int[Z*N*(N-1)];								// total data required for all pointers
	/*for (size_t i = 0; i < N*(N-1); i++)
	{
		unsigned int q =  i * Z;
		original.histograms[i] = (unsigned int*)hist + q;                  // setting the values of rowPtrs       
 	}*/
	/*for (int i = 0; i<N*(N-1);i++)
	{

		for(int j =0;j<=(Z-1);j++)
			original.histograms[i][j] = 0;
	}*/

	//for(int j =0;j<=255;j++)
	//	original.histogram[j] = 0;

	unsigned int sum = 0;
	unsigned long long sum_un = 0;
	unsigned int Max = 0;
	unsigned int Min = -1;
	unsigned int s_length = 0;

	s_length = original.Multibeam_split(filename);					//split function

#ifdef DEBUG
	end = clock();
	elapsed = ((double) (end-start)) / ((double) TIME_NORMALIZER);
	cout << '\t' << '\t' << "               Time taken for reading png = " << original.png_time << endl;
	cout << '\t' << '\t' << "		 Time taken for prefix compression = " << original.para_comptime/((double) TIME_NORMALIZER) << endl;
	cout << '\t' << '\t' << "	Total split and parallel encoding time = " << elapsed << endl;
	start = clock();
#else
	//original.deflate_compression(filename);
#endif

#ifdef DEBUG
	end = clock();
	elapsed = ((double) (end-start)) / ((double) TIME_NORMALIZER);
	cout << '\t' << '\t' << "                 deflate encoding time = " << elapsed << endl;
	start = clock();
#else
	//original.EntropyEncoder_AC(filename);
#endif

#ifdef DEBUG
	end = clock();
	elapsed = ((double) (end-start)) / ((double) TIME_NORMALIZER);
	cout << '\t' << '\t' << "                      AC encoding time = " << elapsed << endl;
#endif

	sum_un  = sum_un + original.Uncompressed_Beamlength;
	sum_un = sum_un*(N*(N-1));



	/*for (int i = 0; i<N*(N-1);i++)					// output histograms of indiviual beams
	{
		//cout << endl << "-----------------------------Histogram of beam: " << i << "-----------------------------" << endl; 
		
		for(int j =0;j<=(Z-1);j++)
		{
			//cout << j << '\t' << original.histograms[i][j] << endl; 
			original.histogram[j] = original.histogram[j] +  original.histograms[i][j] ;

		}
		//cout << endl << "-----------------------------Entropy of beam " << i << ": " << original.Entropy(i) << endl; 
	}*/

	/*cout << endl << "Histogram of all the beams combined: (for compressed file, ignore for uncompressed)" << endl;
	for(int j =0;j<=(Z-1);j++)
		cout << j << '\t' << original.histogram[j] << endl; 
	
	cout << endl << "Run length statistics" << endl;

	for(int j =0;j<=s_length;j++)
		cout << j << '\t' << original.RLE_histo[j] << endl; */


	//cout << endl << "-----------------------------Entropy of total data: " << original.Entropy(-1) << endl;

	cout << endl << "comp_length: " << original.comp_length << endl;
	cout <<"Total uncompressed file size (MB): " << (double(sum_un)*(10))/double(1024*1024*8) << endl;
	cout <<"Prefix compressed file size (MB): " << double(original.comp_length)/double(1024*1024*8) << endl;
	cout <<"sum_zero: " <<original.sum_zero << endl;
	cout <<"sum_un: " << sum_un << endl;
	cout << "fraction of zeros : " << double(original.sum_zero)/double(sum_un) << endl;

	cout << endl << "Prefix Compression ratio:" << (double(sum_un)*10)/double(original.comp_length) << endl;
	cout << endl << "Prefix Speedup: " << (double(sum_un)*10)/(double(original.cycles)) << endl;

	cout << "Number of empty frames: " << original.count_EF << endl;
	cout << "Number of empty columns:" << original.count_EC << endl;

	
#ifdef DEBUG
	
#else
	end = clock();
	elapsed = ((double) (end-start)) / ((double) TIME_NORMALIZER);
	cout << "\t" << "\t" << "               Time taken for reading png = " << original.png_time << endl;
	cout << '\t' << '\t' << "		 Time taken for prefix compression = " << original.para_comptime/((double) TIME_NORMALIZER) << endl;
	cout << '\t' << '\t' << "                      Total encoding time = " << elapsed << endl;
	
#endif

	// free memory
	//if (original.BeamPtrs != NULL) delete[] (unsigned char*)original.BeamPtrs;
	if (original.Beam_RLE_length != NULL) delete [] original.Beam_RLE_length;
	if (original.Total_Beam_RLE_length != NULL) delete [] original.Total_Beam_RLE_length;
//	if (hist != NULL) delete[] hist;
	if (original.RLE_histo != NULL) delete [] original.RLE_histo;
	if (original.histograms != NULL) delete [] original.histograms;
	if (original.compressed !=NULL) delete [] original.compressed;
	if (original.compressedchar !=NULL) delete [] original.compressedchar;
	//if (original.arr != NULL) delete [] original.arr;
	original.Free();

	return 0;
}

unsigned char *EndPoint::Decompression(string filename)
{
	#ifdef DEBUG
	cout << "File Name : " << filename << endl;
    #endif
	
	// Start Decompression ---------------------------------------------->
	unsigned int start, end;
	double elapsed;
	
	ClassDecoder corner;
	
	// Initialize Workspace
	corner.width = width;
	corner.height= height;
//	corner.RLE_length = RLE_Length;

	start = clock();

	//corner.EntropyDecoder_AC(filename);
	//corner.inflate_decompression(filename);
	//corner.ReadRLE(filename,0);
	

#ifdef DEBUG
	end = clock();
	elapsed = ((double) (end-start)) / ((double) TIME_NORMALIZER);
	cout << "AC Decoding Time = " << elapsed << endl;
	start = clock();
	//corner.inflate_decompression(filename);
	end = clock();
	elapsed = ((double) (end-start)) / ((double) TIME_NORMALIZER);
	cout << "inflate Decoding Time = " << elapsed << endl;
	start = clock();
#endif

	corner.WritePNG(filename);
		// Inverse Transform: Corner -> Image

	//corner.Transform_RLE_EOB_decoding();
	
	end = clock();
	// End Decompression ------------------------------------------------>
	
	// Compute Runtime
	elapsed = ((double) (end-start)) / ((double) TIME_NORMALIZER);
#ifndef DEBUG
	cout << "\t" << "\t" << "               Time taken for writing png = " << corner.png_time << endl;
	cout << "\t" << "\t" << "                      Total decoding time = " << elapsed << endl;
#else
	cout << "RLE + EOB + Inverse Transform + writing file Time = " << elapsed << endl;
#endif

	corner.Free();
	return NULL;
}
