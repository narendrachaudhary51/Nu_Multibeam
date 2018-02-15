/*
 *  ClassEncoder.cpp
 *  
 *
 *  Created by Narendra Chaudhary on 2/27/2015.
 *  Texas A&M University.
 *
 */
#ifndef ClassEncoder_H_
#include "ClassEncoder.h"
#endif

#ifndef AC_HEADER
#include "AC.h"
#endif
#include <algorithm>
#include <math.h>
#include <bitset>

#define logX(x, y) log((double) x)/log((double) y)
#define log2(x) log((double) x)/log(2.0)

ofstream csvfileout;							// stream for csv file

/* Constructor */
ClassEncoder::ClassEncoder()
{
	//MAX_SYMBOL = 0;
}

/* Destructor */
ClassEncoder::~ClassEncoder()
{
}

//----------------------------------Multibeam -------------------------------------------------
int ClassEncoder::Multibeam_split(string filename)
{
	#ifdef DEBUG
	cout << "Doing Multibeam reading and split............ " << endl;
    #endif

	//-----------------memory allocation-------
	
	if( (compressed = (char *) calloc((width*height)*5, sizeof(char))) == NULL)
	{
		cout << endl << "[Error] Cannot allocate memory for compressed ... Please check memory and try again later..." << endl;
		exit(-1);
	}

	arr = (char**) calloc(CUBELENGTH, sizeof(char*));

	for (int i = 0; i < CUBELENGTH; i++)
	{
		arr[i] = (char*) calloc(12*N*(N-1), sizeof(char));
	}

	unsigned char *buffer_front_row, *buffer_last_row;
	if( (buffer_front_row = (unsigned char *) calloc(width, sizeof(char))) == NULL)
	{
		cout << endl << "[Error] Cannot allocate memory for previous row buffer... Please check memory and try again later..." << endl;
		exit(-1);
	}
	if( (buffer_last_row = (unsigned char *) calloc(width, sizeof(char))) == NULL)
	{
		cout << endl << "[Error] Cannot allocate memory for previous row buffer... Please check memory and try again later..." << endl;
		exit(-1);
	}

	BeamPtrs = new unsigned char* [N*(N-1)];

	unsigned int buffersize  = (width + 2*N*Distance)*(ROWBUFFER+1);
	
	char *data2 = new char[buffersize];								// total data required for all pointers
	unsigned int stride1 = (buffersize/(N*(N-1))) ;  cout << "stride length for beam pointers: "<<  stride1 << endl; 
	unsigned int stripe_length = 0;

	/*if( (RLE_histo = (unsigned int *) calloc(stride1, sizeof(int))) == NULL)
	{
		cout << endl << "[Error] Cannot allocate memory for RLE_histo... Please check memory and try again later..." << endl;
		exit(-1);
	}*/
	
	for (size_t i = 0; i < N*(N-1); i++)
	{
		size_t q =  i * stride1;
		BeamPtrs[i] = (unsigned char*)data2 + q;                  // setting the values of rowPtrs       
 	}
	
	

//------------------------------------------------------loop start----------------------------------------------------

	
	unsigned int new_buffer_length = ROWBUFFER; 
	
	unsigned int count_beam = 0;
	//count_beam = new unsigned int[N*(N-1)];								//allocate memory for count_beam
	
	//for(unsigned int i=0;i < (N*(N-1));i++)
	//	count_beam[i] = 0;

	int sign = 1;															//sign for handling traversal in Y-direction 
	int temp = 0;
	int edgeflag = 0;
	div_t divresult;


#ifdef write_single_file

		csvfileout.open("compressed.csv");
		
#endif
	
	for(unsigned int row_number=0; row_number<height;)               //reading number of rows starting from row_number
	{   
		if ((row_number + ROWBUFFER) < height)
			ReadRows(row_number,ROWBUFFER);                           
		else
		{  
			new_buffer_length = (height - row_number);               // check and modify for last chunk of row buffer
            ReadRows(row_number,new_buffer_length);
		}

		//#pragma omp parallel for		
		for(int x = 0,y=0,z=(ROWBUFFER-1); (x < (width + (N-1)*Distance));)				// loop for N*(N-1) beams traversal
		{
			if(N ==2 && x==0 && y==0) x =-1;										// edge case for 2x1 case
			if(y == Distance)														// case to include the front (top) row
			{
				for(int i=0; i<(N-1) ;i++)															//loop for top row
				{
					if(((x - i*Distance) < 0) || ((x - i*Distance) >= width))
					{
						//(*(BeamPtrs[i + (N-1)*(N-1)]+count_beam[i + (N-1)*(N-1)])) = 0;
						(*(BeamPtrs[i + (N-1)*(N-1)]+count_beam)) = 0;
					}
					else
					{
						//(*(BeamPtrs[i + (N-1)*(N-1)]+count_beam[i + (N-1)*(N-1)])) = (buffer_front_row[x - i*Distance])/DIV;
						(*(BeamPtrs[i + (N-1)*(N-1)]+count_beam)) = (buffer_front_row[x - i*Distance])/DIV;
					}
				}

				for(int j = 0; j < (N-1) ;j++)						// loop for beam array rows
				{
					for(int i = 0; i<(N-1) ; i++)					// loop for beam array columns
					{
						if(((x - i*Distance) < 0) || ((x - i*Distance) >= width) || ((z - j*Distance) >= new_buffer_length))
						{
							//(*(BeamPtrs[i + j*(N-1)]+count_beam[i + j*(N-1)])) = 0;
							(*(BeamPtrs[i + j*(N-1)]+count_beam)) = 0;

						}
						else
						{

							//(*(BeamPtrs[i + j*(N-1)]+count_beam[i + j*(N-1)])) = (*(rowPtrs[z - j*Distance] + x - i*Distance))/DIV;
							(*(BeamPtrs[i + j*(N-1)]+count_beam)) = (*(rowPtrs[z - j*Distance] + x - i*Distance))/DIV;
						}
					}
				}

				if((y+1)%(Distance+1) == 0)						
					sign = (-1)*sign;
				if(y==0)
					sign = 1;

				 y = y + (sign)*1 ;						// update y by +-1  based on sign 	
				 z = z - (sign)*1 ;						// update z based on sign 

				 if(y%Distance == 0)								// update x 
					 x = x + (2*Distance) - 1;
				 else
					 x = x + (2*Distance);
				 
				 count_beam++;

			}

//----------------------------------------------Normal case---------------------------------------------------------
			else					//normal case
			{	
				//omp_set_dynamic(0);     // Explicitly disable dynamic teams
				//omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
				#pragma omp parallel for
				for(int j = 0; j < N ;j++)				// loop for beam array rows
				{
					
					for(int i = 0; i<(N-1) ; i++)				// loop for beam array columns
					{
						if(((x - i*Distance) < 0) || ((x - i*Distance) >= width) || ((z - j*Distance) >= new_buffer_length))
						{
							//(*(BeamPtrs[i + j*(N-1)]+count_beam[i + j*(N-1)])) = 0;
							(*(BeamPtrs[i + j*(N-1)]+count_beam)) = 0;
						}
						else
						{
							//(*(BeamPtrs[i + j*(N-1)]+count_beam[i + j*(N-1)])) = (*(rowPtrs[z - j*Distance] + x - i*Distance))/DIV;
							(*(BeamPtrs[i + j*(N-1)]+count_beam)) = (*(rowPtrs[z - j*Distance] + x - i*Distance))/DIV;
						}
					}
				}
				

				if((y+1)%(Distance+1) == 0)
					 sign = (-1)*sign;
				if(y==0)
					sign = 1;

				y = y + (sign)*1 ;						// update y by +-1  based on sign 	
				z = z - (sign)*1 ;						// update z based on sign 

				if(y%Distance == 0)								// update x 
					 x = x + (2*Distance) - 1;
				else
					 x = x + (2*Distance);

				count_beam++;
				
			}
		}     


		if(ROWBUFFER == new_buffer_length)
			memcpy(buffer_front_row,rowPtrs[ROWBUFFER-1],width*sizeof(char));


		/*for(int i=0;i < N*(N-1); i++)									// write all beams to seperate files
		{
			if(row_number==0)												
				WriteBeam(filename,i,count_beam,1,0);				// start of writing file
			else
				WriteBeam(filename,i,count_beam,1,1);

		}*/


		if(row_number == 0)
		{
			cout << "stripe length (count_beam): " << count_beam << endl;
			stripe_length = count_beam;
			
		}
		Uncompressed_Beamlength = Uncompressed_Beamlength + count_beam;
		count_beam = 0;

		row_number = row_number + new_buffer_length;                           //update starting row

		if(N ==2 && row_number>=height && edgeflag == 0)				//hack for last row of 2x1 case
		{
			edgeflag = 1;
			memcpy(buffer_last_row,rowPtrs[ROWBUFFER-1],width*sizeof(char));
		}
		
		//serial_prefix_compression(filename,stripe_length);
		parallel_prefix_compression(filename,stripe_length);

	}

//----------------------------------------------------------- end of long loop ----------------------------------------
	int counter = 0;
	if(N==2 && edgeflag == 1)								//hack for last column of 2 beams case
	{
		edgeflag = 0;
		for (int x = 0,y=0; (x < (width + (N-1)*Distance));)
		{
			if(x==0 && y==0)
				x = -1;

			(*(BeamPtrs[0]+counter)) = 0;
			if (y == 2)
			{
				if(x < 0 || x >= width)
				{
					(*(BeamPtrs[1]+counter)) = 0;
				}
				else
				{
					(*(BeamPtrs[1]+counter)) = buffer_last_row[x];
				
				}
			}
			else 
			{
				(*(BeamPtrs[1]+counter)) = 0;
			}

			counter++;

			if((y+1)%(Distance+1) == 0)
				sign = (-1)*sign;
			if(y==0)
				sign = 1;

			y = y + (sign)*1 ;	

			if(y%Distance == 0)								// update x 
					 x = x + (2*Distance) - 1;
				else
					 x = x + (2*Distance);

		}
		Uncompressed_Beamlength = Uncompressed_Beamlength + counter;
		
		WriteBeam(filename,0,counter,1,1);
		WriteBeam(filename,1,counter,1,1);
	}
	
	// free memory
	free(arr);
	if (BeamPtrs != NULL) 
	{
		delete[] (unsigned char*)BeamPtrs; BeamPtrs = NULL;
	}
	if (rowPtrs != NULL) 
	{
		delete [] rowPtrs; 
		rowPtrs = NULL;
	}
    if (data1 != NULL) 
	{
		delete [] data1; 
		data1 = NULL;
	}

#ifdef write_single_file
	csvfileout.close();
#endif
	return stripe_length;

}

void ClassEncoder::parallel_prefix_compression(string filename, unsigned int stripe_length)
{
	double start,stop = 0;
	start = clock();
	if (rowPtrs != NULL) 
	{
		delete [] rowPtrs; 
		rowPtrs = NULL;
	}
    if (data1 != NULL) 
	{
		delete [] data1; 
		data1 = NULL;
	}

	filename.erase(filename.end()-4,filename.end());


	unsigned int Max_comp_length = 12*N*(N-1); 
	
	unsigned int lengths[CUBELENGTH];
	

	unsigned int cubes = 0;
	unsigned int cubesize = CUBELENGTH;
	

	if(stripe_length%CUBELENGTH == 0)
		cubes = (stripe_length/CUBELENGTH);
	else
		cubes = (stripe_length/CUBELENGTH) + 1;


	unsigned int count_column = 0;
	unsigned int count_superpacket = 0;
	unsigned int count_EC_para = 0;
	unsigned int count_EF_para = 0;
	unsigned int cycles_para = 0;
	unsigned int sum_zero_para = 0;
	unsigned int sum_lengths = 0;

	for (int c = 0;c < cubes; c++)								//traversal on cubes
	{
		count_column = 0;
		count_superpacket = 0;
		count_EC_para = 0;
		count_EF_para = 0;
		cycles_para = 0;
		sum_zero_para = 0;

		for (int m = 0; m < cubesize; m++)
			lengths[m] = 0;

		if (c == cubes - 1)
			cubesize = stripe_length - c*(CUBELENGTH);


		//omp_set_dynamic(0);								// Explicitly disable dynamic teams
		//omp_set_num_threads(4);							// Use 4 threads for all consecutive parallel regions
		#pragma omp parallel for private(count_column,count_superpacket) reduction(+:count_EC_para,count_EF_para,cycles_para,sum_zero_para) 
		for ( int i=0; i < cubesize; i++)									// traversal on superpackets
		{	
			count_superpacket = 0;
			//cycles_para = cycles_para + N*(N-1);                           //v2
			cycles_para = cycles_para + N*(N-1) + 10*(N-1);					//v3
			sum_zero_para = sum_zero_para + N*(N-1);						//assume all zeros

			for (int j = 0; j < (N-1); j++)									// traversal on coloumns
			{
				count_column = 0;
			
				for(int k = 0; k < N; k++)									// traversal on rows
				{

					if((*(BeamPtrs[(N-1)*j + k] + i + c*CUBELENGTH)) == 0)
					{
						count_column++;
					
						//arr[i][lengths[i]] = 1;
						lengths[i]++;
					}
					else
					{
						sum_zero_para--;									//our assumption wrong so decrease sum_zero_para
						//arr[i][lengths[i]] = 0;
						lengths[i]++;

						//arr[i][lengths[i]] = 1;
						lengths[i]++;

						bitset<10> bits(4*(*(BeamPtrs[(N-1)*j + k] + i + c*CUBELENGTH)));
					

						//for (int l = 0; l < 10; l++)
							//arr[i][lengths[i] + l] = bits[l];
					

						lengths[i] = lengths[i] + 10;
						cycles_para = cycles_para + 11;
					}

				}
				if (count_column == N)
				{
					count_superpacket = count_superpacket + N;
					count_EC_para++;
					lengths[i] = lengths[i] - N;
					//cycles_para = cycles_para + 13;			//v1.2
					cycles_para = cycles_para + 3;				//v1.3

					//arr[i][lengths[i]] = 0;   
					lengths[i]++;

					//arr[i][lengths[i]] = 0;   
					lengths[i]++;

					//arr[i][lengths[i]] = 0;
					lengths[i]++;
				}
			}

			if(count_superpacket == N*(N-1))
			{
				count_EF_para++;
				//cycles_para = cycles_para - N*(N-1) + 3;				//v1.0
				cycles_para = cycles_para - N*(N-1) - 13*(N-1) + 3;    // v1.2

				count_EC_para = count_EC_para - (N-1);
				lengths[i] = lengths[i] - 3*(N-1);

				//arr[i][lengths[i]] = 0;   
				lengths[i]++;

				//arr[i][lengths[i]] = 0;   
				lengths[i]++;

				//arr[i][lengths[i]] = 1;
				lengths[i]++;
			}
		}

		sum_lengths = 0;
		for(int i = 0; i < cubesize; i++)
		{
			//memcpy(compressed + comp_length + sum_lengths,arr[i],lengths[i]*sizeof(char));
			sum_lengths = sum_lengths + lengths[i];
		}

		count_EF = count_EF + count_EF_para;
		count_EC = count_EC + count_EC_para;
		comp_length = comp_length + sum_lengths;
		cycles = cycles + cycles_para;
		sum_zero = sum_zero + sum_zero_para;
	}
	
	stop = clock();
	para_comptime  = para_comptime + stop - start; 
}

void ClassEncoder::serial_prefix_compression(string filename, unsigned int stripe_length)
{
	double start,stop = 0;
	start = clock();
	if (rowPtrs != NULL) 
	{
		delete [] rowPtrs; 
		rowPtrs = NULL;
	}
    if (data1 != NULL) 
	{
		delete [] data1; 
		data1 = NULL;
	}

	filename.erase(filename.end()-4,filename.end());	

	unsigned int count_column = 0;
	unsigned int count_superpacket = 0;
	unsigned int count_EC_ser = 0;
	unsigned int count_EF_ser = 0;
	unsigned long long comp_length_ser = 0;
	unsigned long long cycles_ser = 0;
	unsigned int sum_zero_ser = 0;
	
	for ( int i=0; i < stripe_length;i++)						// traversal on superpackets
	{	
		count_superpacket = 0;
		//cycles_ser = cycles_ser + N*(N-1);								// minimum 1 cycle required per symbol
		cycles_ser = cycles_ser + N*(N-1) + 10*(N-1);			//v3
		sum_zero_ser = sum_zero_ser + N*(N-1);						    // assume all are zero

		for (int j = 0; j < (N-1); j++)									// traversal on coloumns
		{
			count_column = 0;
			
			for(int k = 0; k < N; k++)									// traversal on rows
			{

				if((*(BeamPtrs[(N-1)*j + k] + i)) == 0)
				{
					count_column++;
					//compressed[comp_length + comp_length_ser] = 1;
					comp_length_ser++;

				}
				else
				{
					sum_zero_ser--;									//our assumption was wrong decrease sum_zero_ser

					//compressed[comp_length + comp_length_ser] = 0;
					comp_length_ser++;
					//compressed[comp_length + comp_length_ser] = 1;
					comp_length_ser++;

					//bitset<10> bits(4*(*(BeamPtrs[(N-1)*j + k] + i)) + (rand()%4 + 0));
					bitset<10> bits(4*(*(BeamPtrs[(N-1)*j + k] + i)));

					//for (int l = 0; l < 10; l++)
					//	compressed[comp_length + comp_length_ser + l] = bits[l];

					comp_length_ser = comp_length_ser + 10;
					cycles_ser = cycles_ser + 11;
				}

			}
			if (count_column == N)
			{
				count_superpacket = count_superpacket + N;
				count_EC_ser++;
				comp_length_ser = comp_length_ser - N;
				//cycles_ser = cycles_ser + 13;				// v1.2
				cycles_ser = cycles_ser + 3;				// v1.3

				//compressed[comp_length + comp_length_ser] = 0;   
				comp_length_ser++;
				//compressed[comp_length + comp_length_ser] = 0;   
				comp_length_ser++;
				//compressed[comp_length + comp_length_ser] = 0;
				comp_length_ser++;

			}
		}

		if(count_superpacket == N*(N-1))
		{
			count_EF_ser++;
			//cycles_ser = cycles_ser - N*(N-1) + 3;				// v 1.0
			cycles_ser = cycles_ser - N*(N-1) - 13*(N-1) + 3;		// v 1.2	
			count_EC_ser = count_EC_ser - (N-1);
			comp_length_ser = comp_length_ser - 3*(N-1);

			//compressed[comp_length + comp_length_ser] = 0;   
			comp_length_ser++;
			//compressed[comp_length + comp_length_ser] = 0;   
			comp_length_ser++;
			//compressed[comp_length + comp_length_ser] = 1;
			comp_length_ser++;

		}
	}
	
	count_EF = count_EF + count_EF_ser;
	count_EC = count_EC + count_EC_ser;
	comp_length = comp_length + comp_length_ser;
	cycles = cycles + cycles_ser;
	sum_zero = sum_zero + sum_zero_ser;
	//cout << "sum nonzero:" << N*(N-1)*stripe_length - sum_zero_ser << endl;
	//cout << "sum_zero: " << sum_zero_ser << endl;
	//cout << "comp_length: " << comp_length << endl;
	stop = clock();
	para_comptime  = para_comptime + stop - start; 
}


// Write memory to text
void ClassEncoder::WriteBeam(string filename,unsigned int Beamnumber,unsigned int length,int istext,int append)
{
	
	unsigned int start,end;
	ios_base::openmode mode;
	start = clock();
	filename.erase(filename.end()-4,filename.end());
	
	if (istext)
	{
		for (unsigned int i=0; i<length;i++)
			histograms[Beamnumber][(*(BeamPtrs[Beamnumber]+i))]++;

#ifdef write_single_file
		
		for (unsigned int i = 0; i < length;i++)
		{
			csvfileout << int(*(BeamPtrs[Beamnumber] + i)) <<",";
		}
		csvfileout << "\n";
#endif

#ifdef write_file
		char str[15];
		sprintf(str, "%d", Beamnumber);
		filename = filename + "_Beam_" + str + ".txt";
		if (append)
			mode = ios::app;
		else
			mode = ios::out;
		ofstream outTXT(filename, mode);
		
		char *Data = new char[length];
		//Data[length] = 0;
		for (unsigned int i=0; i<length;i++)
		{
			Data[i] = '0' + (*(BeamPtrs[Beamnumber]+i));
			//Data[i] = char(*(BeamPtrs[Beamnumber]+i));
			outTXT << Data[i];
		}

		//outTXT << Data;
		outTXT.close();
	
#endif
	}
	else
	{

		char str[15];
		sprintf(str, "%d", Beamnumber);
		filename = filename + "_Beam_" + str + ".txt";
		ofstream out(filename,ios::binary);
		for(unsigned int i=0;i<length;i++)
		{
			out << (*(BeamPtrs[Beamnumber]+i)) << endl;
			//cout<< RLE[i] << '\t';
		}

		out.close();
	}

	end = clock();
	//cout << "\t" << "\t" << "               Compressed file writing time = " << ((double) (end-start)) / ((double) 1000) << endl;
}

	


int ClassEncoder::EntropyEncoder_AC(string filename)						//useless right now
{
	/* Arithmetic Coding */
	// Initialize Arithmetic Coding
	unsigned int ac_bytes = 0;	// Size of the compressed file
	unsigned int ac_bits = 0;
	filename += ".enc";
	ac_encoder ace;
	ac_model acm;
	ac_encoder_init (&ace, filename.c_str());
	
	// Encode each symbol using AC
	/*ac_model_init (&acm, MAX_SYMBOL+M+K+1, NULL, 1);
	for(unsigned int j=0; j<RLE_length; j++)
	{
		ac_encode_symbol(&ace, &acm, RLE[j]);
	}*/
	ac_model_init (&acm, 256, NULL, 1);
	for(unsigned int j=0; j<charlength; j++)
	{
		ac_encode_symbol(&ace, &acm, compressedchar[j]);
	}

	ac_encode_symbol(&ace, &acm, EOF_SYMBOL);
	// Finalize Arithmetic Coder
	ac_encoder_done (&ace);
	ac_model_done (&acm);
	ac_bits = ac_encoder_bits (&ace);
	if ((ac_bits%8) !=0)
		ac_bytes = (int) ((ac_bits/8)+1);
	else
		ac_bytes = (int) (ac_bits/8);
	
#ifdef DEBUG
	cout << "                                [Done]" << endl;
	//cout << "Compressed File Size after AC Encoding in byte = " << ac_bytes << endl;
	cout << "Compressed File Size after AC Encoding in bits = " << ac_bits << endl;
#endif
	//cout << "RLE length               : " << RLE_length << endl;
	cout <<"\t" << "Arithmatic encoded bytes : " << ac_bytes << endl;
	cout.precision(15);
	double c_ratio;
	c_ratio = double(width)*double(height)/double(ac_bytes+12);
	cout << "\t" << "approx. AC coding compression ratio : " << c_ratio << endl;
	return ac_bytes;
}


int ClassEncoder::deflate_compression(string filename)
{
	//filename +=".txt";
	string filename_dest = filename + ".dft"; 
	int ret, flush;
	unsigned int i,count =0;
    unsigned have;
    z_stream strm;
    unsigned char in[32768];
    unsigned char out[32768];
	
	if( (compressedchar = (unsigned char *) calloc(((width*height)*5)/8, sizeof(char))) == NULL)
	{
		cout << endl << "[Error] Cannot allocate memory for compressedchar ... Please check memory and try again later..." << endl;
		exit(-1);
	}

	
	charlength = (comp_length/8) + 8;
	//cout << comp_length << endl;
	//cout << charlength << endl;
	for (unsigned int k = 0; k < charlength; k++)						// needs correction
	{
		for (unsigned int j = 0; j < 8; j++)
		{
			if(compressed[k*8 + j])
				compressedchar[k] |= 1 << j;
		}
	}
	
    /* allocate deflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, Z_DEFLATED);
    if (ret != Z_OK)
        return ret;
	FILE *dest;
	//FILE *source;
	//source = fopen(filename.data(),"r");
	dest = fopen(filename_dest.data(),"wb");
    /* compress until end of file */
    do {
        //strm.avail_in = fread(in, 1, 32768, source);
		count++;
		
		for (i =(count-1)*32768; i<(count*32768) && i<charlength;i++)
		//for (i =(count-1)*32768; i<(count*32768) && i<comp_length;i++)
		{
			//in[i-(count-1)*32768] = unsigned char(RLE[i]);
			in[i-(count-1)*32768] = compressedchar[i];
			//in[i-(count-1)*32768] = unsigned char(compressed[i]);
		}
		strm.avail_in = i - (count-1)*32768;
		
        /*if (ferror(source)) {
            (void)deflateEnd(&strm);
            return Z_ERRNO;
        }*/
		if((i)%(charlength)==0) {
			flush = Z_FINISH;
		}
		else flush = Z_NO_FLUSH;
       // flush = feof(source) ? Z_FINISH : Z_NO_FLUSH;

        strm.next_in = in;

        /* run deflate() on input until output buffer not full, finish
           compression if all of source has been read in */
        do {
            strm.avail_out = 32768;
            strm.next_out = out;
            ret = deflate(&strm, flush);    /* no bad return value */
            //assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            have = 32768 - strm.avail_out;
            if (fwrite(out, 1, have, dest) != have || ferror(dest)) {
                (void)deflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);

        /* done when last data in file processed */
    } while (flush != Z_FINISH);

    /* clean up and return */
    (void)deflateEnd(&strm);
	//fclose(source);
	fclose(dest);
	return 0;
}

// Entropy Encoder: RLE
int ClassEncoder::EntropyEncoder_RLE(string filename,unsigned int Beamnumber,unsigned int length,int append)
{
	filename.erase(filename.end()-4,filename.end());
	//ios_base::openmode mode;

	RLEPtrs[Beamnumber] = new unsigned char[length];

	char str[15];
	sprintf(str, "%d", Beamnumber);
	filename = filename + "_Beam_" + str + "_RLE.txt";
	unsigned int nonzero_count = 0;
	
	Beam_RLE_length[Beamnumber] = 0;							//initialize RLE_length with 0 
	for(unsigned int i=0; i<length; i++)						//loop for Data
	{	
		if(((*(BeamPtrs[Beamnumber]+i)) > 0) && ((*(BeamPtrs[Beamnumber]+i)) < MAX_SYMBOL))				// if symbol is other than zero	
		{
			nonzero_count++;
			(*(RLEPtrs[Beamnumber] + Beam_RLE_length[Beamnumber])) = (*(BeamPtrs[Beamnumber]+i));		// add symbol to particular beam of RLEPtrs
			Beam_RLE_length[Beamnumber]++;
			//histogram[(*(BeamPtrs[Beamnumber]+i))]++;   
			if (nonzero_count > 1)
				RLE_histo[0]++;							// Run length of zero
			//if(i==length)
			//	RLE_histo[length+1]++;						// if it is end of stripe
		}

		else if((*(BeamPtrs[Beamnumber]+i)) == 0)
		{
			nonzero_count = 0;
			unsigned int count = 0;
			while((*(BeamPtrs[Beamnumber]+i)) == 0)
			{
				count ++;
				i++;
			}
			//cout <<"count: " <<count <<endl;
			//cout << "RLE: " << RLE_histo[count]++ << endl;
			RLE_histo[count]++;
			if(count > 0)
			{
				i--;
				unsigned int repeat = (unsigned int) (logX(count, M));
				for(unsigned int r=0; r<=repeat; r++)
				{
					(*(RLEPtrs[Beamnumber] + Beam_RLE_length[Beamnumber])) = MAX_SYMBOL + (count % M);
					Beam_RLE_length[Beamnumber]++;
				//	histogram[MAX_SYMBOL + (count % M)]++;
//					count /= M;
					count >>= RSHIFT_M;		 // If M is 2^N, use N bit right shift instead.

				}
			}
		}
		
	}
	
		(*(RLEPtrs[Beamnumber] + Beam_RLE_length[Beamnumber])) = MAX_SYMBOL + M + 1;
		Beam_RLE_length[Beamnumber]++;
		histogram[MAX_SYMBOL + M + 1]++;
		Total_Beam_RLE_length[Beamnumber] = Total_Beam_RLE_length[Beamnumber] + Beam_RLE_length[Beamnumber] ;
        

#ifdef C_File_Write
		
		if (append)
			mode = ios::app;
		else
			mode = ios::out;
		ofstream outTXT(filename, mode);

		char *Data = new char[Beam_RLE_length[Beamnumber]];
		//Data[length] = 0;
		for (unsigned int i=0; i<Beam_RLE_length[Beamnumber];i++)
		{
			Data[i] = '0' + (*(RLEPtrs[Beamnumber]+i));
			//Data[i] = char(*(BeamPtrs[Beamnumber]+i));
			//outTXT << Data[i];
		}

		outTXT << Data;
		outTXT.close(); 
		
#endif
		if(RLEPtrs[Beamnumber] != NULL) delete [] RLEPtrs[Beamnumber];
		return 0;
}

double ClassEncoder::Entropy(int beamnumber)
{
	unsigned int *p = (unsigned int *) calloc(Z, sizeof(unsigned int));
	unsigned long sum = 0;
	double entropy = 0;
	double *prob = (double *)calloc(Z, sizeof(double));

	if(beamnumber < 0)
	{
		for (int i =0; i<=(Z-1);i++)
		{
			sum += histogram[i];
			p[i] = histogram[i];
		}
	}
	else
	{
		for (int i =0; i<=(Z-1);i++)
		{
			sum += histograms[beamnumber][i];
			p[i] = histograms[beamnumber][i];
		}
		//cout << "sum: " << sum << endl;
	}

	for (int i =0; i<=(Z-1);i++)
	{
		prob[i] = (double) p[i] / (width*height);
		if(p[i] != 0)
			entropy -= prob[i] * log2(prob[i]);
	}
	free(prob);
	free(p);
	return entropy;
}