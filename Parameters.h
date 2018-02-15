/*
 *  Parameters.h
 *  
 *
 *  Created by Narendra Chaudhary on 2/27/2015.
 *  Texas A&M University.
 *
 */

#ifndef Parameters_H_
#define Parameters_H_

#define PNGSIGSIZE 8                    // PNG signature
//#define DEBUG

//#define write_file
//#define write_single_file
//#define C_File_Write

#define N 128
#define Distance 8				   // (Distance = sqrt(N/2))
#define ROWBUFFER (N*Distance)         // based on calculation (N*Distance for zigzag scanning) and ((N-1)*Distance for raster scanning)
#define DIV 1
#define CUBELENGTH 100


#define Z 256						   // size of histogram (16 in normal case,31 in frame prediction)
	
#define WRITE_ROWBUFFER 1              // number of rows writing at a time from PNG. Adjust according to memory constraints

#define ENCODER                        // enable if compling code in encoder mode or both
//#define DECODER                      // enable if compling code in encoder mode or both
//#define PAETH						   // enables paeth filter instead of corner transform
//#define MEMORYSAVE                   // decoder runs in memory save mode if enabled(recommended to keep it enabled)

#define MAX_SYMBOL 31				   // Max symbol in array.
#define MAX_CORNER_SYMBOL 16

#define M 64					       // Base of M-ary RLE for 0s. Should be a power of 2 
#define RSHIFT_M 6					   // log2(M) for faster base M-ary RLE computation

#define K 64     					  // Base of K-ary RLE for EOBs. Should be power of 2
#define RSHIFT_K 6					  // log2(K) for faster base K-ary RLE computation

#define BlockSize 79050	              // Block Size for EOB Coding

#define EOF_SYMBOL MAX_SYMBOL+M+K     // defining end of file symbol for arithmatic encoding

#endif
