/*
 * generalpacker.h
 * 
 * Copyright 2014 Andreas Altair Redmer <altair.ibn.la.ahad.sy@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef generalpacker_h
#define generalpacker_h

#define PROGNAME "yz"
#pragma GCC diagnostic ignored "-Wreorder"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <algorithm>

//#include "vector"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <map>
#include <string>
#ifdef WIN32
#include <io.h>
#include <windows.h>
#else
#include <sys/types.h>
#include <dirent.h>
#endif

#ifndef DEFAULT_OPTION
#define DEFAULT_OPTION 4
#endif
typedef unsigned char U8;
typedef unsigned short U16;
typedef unsigned int U32;
//#ifndef min
inline int min(int a, int b)
{
	return a < b ? a : b;
}
inline int max(int a, int b)
{
	return a < b ? b : a;
}
//#endif

typedef enum
{
	DEFAULT, JPEG, EXE, BINTEXT, TEXT
} Filetype;

#define preprocFlag 1220

#define OPTION_UTF8							1
#define OPTION_USE_NGRAMS					2
#define OPTION_CAPITAL_CONVERSION			4
#define OPTION_WORD_SURROROUNDING_MODELING	8
#define OPTION_SPACE_AFTER_EOL				16
#define OPTION_EOL_CODING					32
#define OPTION_NORMAL_TEXT_FILTER			64
#define OPTION_USE_DICTIONARY				128
#define OPTION_RECORD_INTERLEAVING			256
#define OPTION_DNA_QUARTER_BYTE				512
#define OPTION_TRY_SHORTER_WORD				1024
#define OPTION_TO_LOWER_AFTER_PUNCTUATION	2048
#define OPTION_SPACELESS_WORDS				4096
#define OPTION_ADD_SYMBOLS_0_5				8192
#define OPTION_ADD_SYMBOLS_14_31			16384
#define OPTION_ADD_SYMBOLS_A_Z				32768
#define OPTION_ADD_SYMBOLS_MISC				65536
#define OPTION_SPACE_AFTER_CC_FLAG			131072
#define IF_OPTION(option) ((preprocFlag & option)!=0)

//#include "packer.h"
//int pg_main(int argc, char** argv);
//using namespace std;

//namespace yztools
//{

///////////////////////// state table ////////////////////////

// State table:
//   nex(state, 0) = next state if bit y is 0, 0 <= state < 256
//   nex(state, 1) = next state if bit y is 1
//   nex(state, 2) = number of zeros in bit history represented by state
//   nex(state, 3) = number of ones represented
//
// States represent a bit history within some context.
// State 0 is the starting state (no bits seen).
// States 1-30 represent all possible sequences of 1-4 bits.
// States 31-252 represent a pair of counts, (n0,n1), the number
//   of 0 and 1 bits respectively.  If n0+n1 < 16 then there are
//   two states for each pair, depending on if a 0 or 1 was the last
//   bit seen.
// If n0 and n1 are too large, then there is no state to represent this
// pair, so another state with about the same ratio of n0/n1 is substituted.
// Also, when a bit is observed and the count of the opposite bit is large,
// then part of this count is discarded to favor newer data over old.

static const U8 State_table[256][4] =
{
{ 1, 2, 0, 0 },
{ 3, 5, 1, 0 },
{ 4, 6, 0, 1 },
{ 7, 10, 2, 0 }, // 0-3
		{ 8, 12, 1, 1 },
		{ 9, 13, 1, 1 },
		{ 11, 14, 0, 2 },
		{ 15, 19, 3, 0 }, // 4-7
		{ 16, 23, 2, 1 },
		{ 17, 24, 2, 1 },
		{ 18, 25, 2, 1 },
		{ 20, 27, 1, 2 }, // 8-11
		{ 21, 28, 1, 2 },
		{ 22, 29, 1, 2 },
		{ 26, 30, 0, 3 },
		{ 31, 33, 4, 0 }, // 12-15
		{ 32, 35, 3, 1 },
		{ 32, 35, 3, 1 },
		{ 32, 35, 3, 1 },
		{ 32, 35, 3, 1 }, // 16-19
		{ 34, 37, 2, 2 },
		{ 34, 37, 2, 2 },
		{ 34, 37, 2, 2 },
		{ 34, 37, 2, 2 }, // 20-23
		{ 34, 37, 2, 2 },
		{ 34, 37, 2, 2 },
		{ 36, 39, 1, 3 },
		{ 36, 39, 1, 3 }, // 24-27
		{ 36, 39, 1, 3 },
		{ 36, 39, 1, 3 },
		{ 38, 40, 0, 4 },
		{ 41, 43, 5, 0 }, // 28-31
		{ 42, 45, 4, 1 },
		{ 42, 45, 4, 1 },
		{ 44, 47, 3, 2 },
		{ 44, 47, 3, 2 }, // 32-35
		{ 46, 49, 2, 3 },
		{ 46, 49, 2, 3 },
		{ 48, 51, 1, 4 },
		{ 48, 51, 1, 4 }, // 36-39
		{ 50, 52, 0, 5 },
		{ 53, 43, 6, 0 },
		{ 54, 57, 5, 1 },
		{ 54, 57, 5, 1 }, // 40-43
		{ 56, 59, 4, 2 },
		{ 56, 59, 4, 2 },
		{ 58, 61, 3, 3 },
		{ 58, 61, 3, 3 }, // 44-47
		{ 60, 63, 2, 4 },
		{ 60, 63, 2, 4 },
		{ 62, 65, 1, 5 },
		{ 62, 65, 1, 5 }, // 48-51
		{ 50, 66, 0, 6 },
		{ 67, 55, 7, 0 },
		{ 68, 57, 6, 1 },
		{ 68, 57, 6, 1 }, // 52-55
		{ 70, 73, 5, 2 },
		{ 70, 73, 5, 2 },
		{ 72, 75, 4, 3 },
		{ 72, 75, 4, 3 }, // 56-59
		{ 74, 77, 3, 4 },
		{ 74, 77, 3, 4 },
		{ 76, 79, 2, 5 },
		{ 76, 79, 2, 5 }, // 60-63
		{ 62, 81, 1, 6 },
		{ 62, 81, 1, 6 },
		{ 64, 82, 0, 7 },
		{ 83, 69, 8, 0 }, // 64-67
		{ 84, 71, 7, 1 },
		{ 84, 71, 7, 1 },
		{ 86, 73, 6, 2 },
		{ 86, 73, 6, 2 }, // 68-71
		{ 44, 59, 5, 3 },
		{ 44, 59, 5, 3 },
		{ 58, 61, 4, 4 },
		{ 58, 61, 4, 4 }, // 72-75
		{ 60, 49, 3, 5 },
		{ 60, 49, 3, 5 },
		{ 76, 89, 2, 6 },
		{ 76, 89, 2, 6 }, // 76-79
		{ 78, 91, 1, 7 },
		{ 78, 91, 1, 7 },
		{ 80, 92, 0, 8 },
		{ 93, 69, 9, 0 }, // 80-83
		{ 94, 87, 8, 1 },
		{ 94, 87, 8, 1 },
		{ 96, 45, 7, 2 },
		{ 96, 45, 7, 2 }, // 84-87
		{ 48, 99, 2, 7 },
		{ 48, 99, 2, 7 },
		{ 88, 101, 1, 8 },
		{ 88, 101, 1, 8 }, // 88-91
		{ 80, 102, 0, 9 },
		{ 103, 69, 10, 0 },
		{ 104, 87, 9, 1 },
		{ 104, 87, 9, 1 }, // 92-95
		{ 106, 57, 8, 2 },
		{ 106, 57, 8, 2 },
		{ 62, 109, 2, 8 },
		{ 62, 109, 2, 8 }, // 96-99
		{ 88, 111, 1, 9 },
		{ 88, 111, 1, 9 },
		{ 80, 112, 0, 10 },
		{ 113, 85, 11, 0 }, // 100-103
		{ 114, 87, 10, 1 },
		{ 114, 87, 10, 1 },
		{ 116, 57, 9, 2 },
		{ 116, 57, 9, 2 }, // 104-107
		{ 62, 119, 2, 9 },
		{ 62, 119, 2, 9 },
		{ 88, 121, 1, 10 },
		{ 88, 121, 1, 10 }, // 108-111
		{ 90, 122, 0, 11 },
		{ 123, 85, 12, 0 },
		{ 124, 97, 11, 1 },
		{ 124, 97, 11, 1 }, // 112-115
		{ 126, 57, 10, 2 },
		{ 126, 57, 10, 2 },
		{ 62, 129, 2, 10 },
		{ 62, 129, 2, 10 }, // 116-119
		{ 98, 131, 1, 11 },
		{ 98, 131, 1, 11 },
		{ 90, 132, 0, 12 },
		{ 133, 85, 13, 0 }, // 120-123
		{ 134, 97, 12, 1 },
		{ 134, 97, 12, 1 },
		{ 136, 57, 11, 2 },
		{ 136, 57, 11, 2 }, // 124-127
		{ 62, 139, 2, 11 },
		{ 62, 139, 2, 11 },
		{ 98, 141, 1, 12 },
		{ 98, 141, 1, 12 }, // 128-131
		{ 90, 142, 0, 13 },
		{ 143, 95, 14, 0 },
		{ 144, 97, 13, 1 },
		{ 144, 97, 13, 1 }, // 132-135
		{ 68, 57, 12, 2 },
		{ 68, 57, 12, 2 },
		{ 62, 81, 2, 12 },
		{ 62, 81, 2, 12 }, // 136-139
		{ 98, 147, 1, 13 },
		{ 98, 147, 1, 13 },
		{ 100, 148, 0, 14 },
		{ 149, 95, 15, 0 }, // 140-143
		{ 150, 107, 14, 1 },
		{ 150, 107, 14, 1 },
		{ 108, 151, 1, 14 },
		{ 108, 151, 1, 14 }, // 144-147
		{ 100, 152, 0, 15 },
		{ 153, 95, 16, 0 },
		{ 154, 107, 15, 1 },
		{ 108, 155, 1, 15 }, // 148-151
		{ 100, 156, 0, 16 },
		{ 157, 95, 17, 0 },
		{ 158, 107, 16, 1 },
		{ 108, 159, 1, 16 }, // 152-155
		{ 100, 160, 0, 17 },
		{ 161, 105, 18, 0 },
		{ 162, 107, 17, 1 },
		{ 108, 163, 1, 17 }, // 156-159
		{ 110, 164, 0, 18 },
		{ 165, 105, 19, 0 },
		{ 166, 117, 18, 1 },
		{ 118, 167, 1, 18 }, // 160-163
		{ 110, 168, 0, 19 },
		{ 169, 105, 20, 0 },
		{ 170, 117, 19, 1 },
		{ 118, 171, 1, 19 }, // 164-167
		{ 110, 172, 0, 20 },
		{ 173, 105, 21, 0 },
		{ 174, 117, 20, 1 },
		{ 118, 175, 1, 20 }, // 168-171
		{ 110, 176, 0, 21 },
		{ 177, 105, 22, 0 },
		{ 178, 117, 21, 1 },
		{ 118, 179, 1, 21 }, // 172-175
		{ 110, 180, 0, 22 },
		{ 181, 115, 23, 0 },
		{ 182, 117, 22, 1 },
		{ 118, 183, 1, 22 }, // 176-179
		{ 120, 184, 0, 23 },
		{ 185, 115, 24, 0 },
		{ 186, 127, 23, 1 },
		{ 128, 187, 1, 23 }, // 180-183
		{ 120, 188, 0, 24 },
		{ 189, 115, 25, 0 },
		{ 190, 127, 24, 1 },
		{ 128, 191, 1, 24 }, // 184-187
		{ 120, 192, 0, 25 },
		{ 193, 115, 26, 0 },
		{ 194, 127, 25, 1 },
		{ 128, 195, 1, 25 }, // 188-191
		{ 120, 196, 0, 26 },
		{ 197, 115, 27, 0 },
		{ 198, 127, 26, 1 },
		{ 128, 199, 1, 26 }, // 192-195
		{ 120, 200, 0, 27 },
		{ 201, 115, 28, 0 },
		{ 202, 127, 27, 1 },
		{ 128, 203, 1, 27 }, // 196-199
		{ 120, 204, 0, 28 },
		{ 205, 115, 29, 0 },
		{ 206, 127, 28, 1 },
		{ 128, 207, 1, 28 }, // 200-203
		{ 120, 208, 0, 29 },
		{ 209, 125, 30, 0 },
		{ 210, 127, 29, 1 },
		{ 128, 211, 1, 29 }, // 204-207
		{ 130, 212, 0, 30 },
		{ 213, 125, 31, 0 },
		{ 214, 137, 30, 1 },
		{ 138, 215, 1, 30 }, // 208-211
		{ 130, 216, 0, 31 },
		{ 217, 125, 32, 0 },
		{ 218, 137, 31, 1 },
		{ 138, 219, 1, 31 }, // 212-215
		{ 130, 220, 0, 32 },
		{ 221, 125, 33, 0 },
		{ 222, 137, 32, 1 },
		{ 138, 223, 1, 32 }, // 216-219
		{ 130, 224, 0, 33 },
		{ 225, 125, 34, 0 },
		{ 226, 137, 33, 1 },
		{ 138, 227, 1, 33 }, // 220-223
		{ 130, 228, 0, 34 },
		{ 229, 125, 35, 0 },
		{ 230, 137, 34, 1 },
		{ 138, 231, 1, 34 }, // 224-227
		{ 130, 232, 0, 35 },
		{ 233, 125, 36, 0 },
		{ 234, 137, 35, 1 },
		{ 138, 235, 1, 35 }, // 228-231
		{ 130, 236, 0, 36 },
		{ 237, 125, 37, 0 },
		{ 238, 137, 36, 1 },
		{ 138, 239, 1, 36 }, // 232-235
		{ 130, 240, 0, 37 },
		{ 241, 125, 38, 0 },
		{ 242, 137, 37, 1 },
		{ 138, 243, 1, 37 }, // 236-239
		{ 130, 244, 0, 38 },
		{ 245, 135, 39, 0 },
		{ 246, 137, 38, 1 },
		{ 138, 247, 1, 38 }, // 240-243
		{ 140, 248, 0, 39 },
		{ 249, 135, 40, 0 },
		{ 250, 69, 39, 1 },
		{ 80, 251, 1, 39 }, // 244-247
		{ 140, 252, 0, 40 },
		{ 249, 135, 41, 0 },
		{ 250, 69, 40, 1 },
		{ 80, 251, 1, 40 }, // 248-251
		{ 140, 252, 0, 41 } };  // 252, 253-255 are reserved

class GeneralPacker
{
public:
	class StateMap;
	class Mixer;
//private:
	// ilog(x) = round(log2(x) * 16), 0 <= x < 64K
	class Ilog
	{
		U8 t[65536];
	public:
		int
		operator()(U16 x) const;
		Ilog();
	};

	static class Ilog ilog;

	// Track time and memory used
	class ProgramObserver
	{
		int memused;  // bytes allocated by Array<T> now
		clock_t start_time;  // in ticks
	public:
		int maxmem;   // most bytes allocated ever
		void alloc(int n);
		ProgramObserver();
		void print() const;
	};

	static class ProgramObserver programObserver;

	//class PseudoRandomGenerator getRandomNumber;
	// 32-bit pseudo random number generator
	class PseudoRandomGenerator
	{
		U32 table[64];
		int i;
	public:
		PseudoRandomGenerator();
		U32 operator()();
	};

	static class PseudoRandomGenerator getRandomNumber;

	//////////////////////////// Stretch ///////////////////////////////

	// Inverse of squash. d = ln(p/(1-p)), d scaled by 8 bits, p by 12 bits.
	// d has range -2047 to 2047 representing -8 to 8.  p has range 0 to 4095.

	// TODO convert this inline class to method
	class Stretch
	{
		short t[4096];
	public:
		Stretch();
		int operator()(int pD) const;
	};

	static class Stretch stretch;

	//////////////////////////// FastArray ////////////////////////////

	// FastArray<T, ALIGN> a(n); creates n elements of T initialized to 0 bits.
	// Constructors for T are not called.
	// Indexing is bounds checked if assertions are on.
	// a.size() returns n.
	// a.resize(n) changes size to n, padding with 0 bits or truncating.
	// a.push_back(x) appends x and increases size by 1, reserving up to size*2.
	// a.pop_back() decreases size by 1, does not free memory.
	// Copy and assignment are not supported.
	// Memory is aligned on a ALIGN byte boundary (power of 2), default is none.

	template<class T, int ALIGN = 0> class FastArray
	{
	private:
		int n;     // user size
		int reserved;  // actual size
		char *ptr; // allocated memory, zeroed
		T* data;   // start of n elements of aligned data
		void create(int i);  // create with size i
	public:
		explicit FastArray(int i = 0)
		{
			create(i);
		}
		~FastArray();
		T& operator[](int i)
		{
			return data[i];
		}
		const T& operator[](int i) const
		{
			return data[i];
		}
		int size() const
		{
			return n;
		}
		void resize(int i);  // change size to i
		void pop_back()
		{
			if (n > 0)
				--n;
		}  // decrement size
		void push_back(const T& x);  // increment size, append x
		//private:
		FastArray(const FastArray&);  // no copy or assignment
		FastArray& operator=(const FastArray&);
	};

	////////////////////////////// FastArrayBuffer /////////////////////////////

	// FastArrayBuffer(n) buf; creates an array of n bytes (must be a power of 2).
	// buf[i] returns a reference to the i'th byte with wrap (no out of bounds).
	// buf(i) returns i'th byte back from pos (i > 0)
	// buf.size() returns n.

	class FastArrayBuffer
	{
		FastArray<U8> b;
	public:
		FastArrayBuffer(int i);
		void setsize(int i); //TODO use i only as counter
		U8& operator[](int i);
		int operator()(int i) const;
		int size() const;
	};

private:

	// Error handler: print message if any, and exit
	static inline void quitAndThrowError(const char* message = 0)
	{
		throw message;
	}

	static int fileSize, fileType;  // Set by Filter
	static long size;

	static int pos;  // Number of input bytes in buf (not wrapped)
	/////////////////////// Global context /////////////////////////

	static int compressionLevel;  // Compression level 0 to 9
#define MEM (0x10000<<compressionLevel)
	static int y;  // Last bit, 0 or 1, set by encoder

	// Global context set by Predictor and available to all models.
	static int c0; // Last 0-7 bits of the partial byte with a leading 1 bit (1-255)
	static U32 b1, b2, b3, b4, b5, b6, b7, b8, tt, c4, x4, x5, w4, w5, f4; // Last 4 whole bytes, packed.  Last byte is bits 0-7.
	static int order, bpos, cxtfl, sm_shft, sm_add, sm_add_y; // bits in c0 (0 to 7)
	static FastArrayBuffer buf;  // Rotating input queue set by Predictor

	// llog(x) accepts 32 bits
	static inline int llog(U32 x)
	{
		if (x >= 0x1000000)
			return 256 + ilog(x >> 16);
		else if (x >= 0x10000)
			return 128 + ilog(x >> 8);
		else
			return ilog(x);
	}

#define nex(state,sel) State_table[state][sel]

	///////////////////////////// Squash //////////////////////////////

	// return p = 1/(1 + exp(-d)), d scaled by 8 bits, p scaled by 12 bits
	static inline int squash(int pD)
	{
		static const int t[33] =
		{ 1, 2, 3, 6, 10, 16, 27, 45, 73, 120, 194, 310, 488, 747, 1101, 1546,
				2047, 2549, 2994, 3348, 3607, 3785, 3901, 3975, 4022, 4050,
				4068, 4079, 4085, 4089, 4092, 4093, 4094 };
		if (pD > 2047)
			return 4095;
		if (pD < -2047)
			return 0;
		int w = pD & 127;
		pD = (pD >> 7) + 16;
		return (t[pD] * (128 - w) + t[(pD + 1)] * w + 64) >> 7;
	}

	static int __attribute__ ((noinline)) dot_product(short *t, short *w,
			int n);

	static void __attribute__ ((noinline)) train(short *t, short *w, int n,
			int err);

	//////////////////////////// hash //////////////////////////////

	// Hash 2-5 ints.
	static inline U32 hash(U32 a, U32 b, U32 c = 0xffffffff)
	{
		U32 h = a * 110002499u + b * 30005491u + c * 50004239u; //+d*70004807u+e*110002499u;
		return h ^ h >> 9 ^ a >> 3 ^ b >> 3 ^ c >> 4;
	}

	// Predict to mixer m from bit history state s, using sm to map s to
	// a probability.
	static inline int mix2(Mixer& m, int s, StateMap& sm)
	{
		int p1 = sm.p(s);
		int n0 = -!nex(s, 2);
		int n1 = -!nex(s, 3);
		int st = stretch(p1);
		if (cxtfl)
		{
			m.add(st / 4);
			int p0 = 4095 - p1;
			m.add((p1 - p0) * 3 / 64);
			m.add(st * (n1 - n0) * 3 / 16);
			m.add(((p1 & n0) - (p0 & n1)) / 16);
			m.add(((p0 & n0) - (p1 & n1)) * 7 / 64);
			return s > 0;
		}
		m.add(st * 9 / 32);
		m.add(st * (n1 - n0) * 3 / 16);
		int p0 = 4095 - p1;
		m.add(((p1 & n0) - (p0 & n1)) / 16);
		m.add(((p0 & n0) - (p1 & n1)) * 7 / 64);
		return s > 0;
	}

	static int contextModel2();
	static void sparseModel(GeneralPacker::Mixer& m);
	static void recordModel(GeneralPacker::Mixer& m);
	static void wordModel(GeneralPacker::Mixer& m);
	static int matchModel(GeneralPacker::Mixer& m);

public:
	GeneralPacker();
	~GeneralPacker();
	void pack();

//////////////////////////// Mixer /////////////////////////////

// Mixer m(N, M, S=1, w=0) combines models using M neural networks with
//   N inputs each, of which up to S may be selected.  If S > 1 then
//   the outputs of these neural networks are combined using another
//   neural network (with parameters S, 1, 1).  If S = 1 then the
//   output is direct.  The weights are initially w (+-32K).
//   It is used as follows:
// m.update() trains the network where the expected output is the
//   last bit (in the global variable y).
// m.add(stretch(p)) inputs prediction from one of N models.  The
//   prediction should be positive to predict a 1 bit, negative for 0,
//   nominally +-256 to +-2K.  The maximum allowed value is +-32K but
//   using such large values may cause overflow if N is large.
// m.set(cxt, range) selects cxt as one of 'range' neural networks to
//   use.  0 <= cxt < range.  Should be called up to S times such
//   that the total of the ranges is <= M.
// m.p() returns the output prediction that the next bit is 1 as a
//   12 bit number (0 to 4095).

	class Mixer
	{
		const int N, M, S;   // max inputs, max contexts, max context sets
		FastArray<short, 16> wx; // N*M weights
		FastArray<int> cxt;  // S contexts
		int ncxt;        // number of contexts (0 to S)
		int base;        // offset of next context
		FastArray<int> pr;   // last result (scaled 12 bits)
		Mixer* mp;       // points to a Mixer to combine results
	public:
		FastArray<short, 16> tx; // N inputs from add()
		int numberOfInputsInTX;          // Number of inputs in tx, 0 to N
		Mixer(int n, int m, int s = 1, int w = 0);

		// Adjust weights to minimize coding cost of last prediction
		void update();

		void update2();

		// Input x (call up to N times)
		void add(int x);

		void mul(int x);

		// Set a context (call S times, sum of ranges <= M)
		void set(int cx, int range);

		// predict next bit
		int p();
		~Mixer();
	};

//////////////////////////// AProbabilityMap //////////////////////////////

// AProbabilityMap maps a probability and a context into a new probability
// that bit y will next be 1.  After each guess it updates
// its state to improve future guesses.  Methods:
//
// AProbabilityMap a(N) creates with N contexts, uses 66*N bytes memory.
// a.p(pr, cx, rate=8) returned adjusted probability in context cx (0 to
//   N-1).  rate determines the learning rate (smaller = faster, default 8).
//   Probabilities are scaled 16 bits (0-65535).

	class AProbabilityMap
	{
		int index;     // last p, context
//const int N;   // number of contexts
		FastArray<U16> t;  // [N][33]:  p, context -> p
	public:
		// maps p, cxt -> p initially
		AProbabilityMap(int n);
		int p(int pr = 2048, int cxt = 0, int rate = 8);
	};

//////////////////////////// StateMap //////////////////////////

// A StateMap maps a nonstationary counter state to a probability.
// After each mapping, the mapping is adjusted to improve future
// predictions.  Methods:
//
// sm.p(cx) converts state cx (0-255) to a probability (0-4095).

// Counter state -> probability * 256
	class StateMap
	{
	protected:
		int cxt;  // context
		U16 t[256]; // 256 states -> probability * 64K
	public:
		StateMap();
		int p(int cx);
	};

///////////////////////////// ByteHash ////////////////////////////////

// A ByteHash maps a 32 bit hash to an array of B bytes (checksum and B-2 values)
//
// ByteHash bh(N); creates N element table with B bytes each.
//   N must be a power of 2.  The first byte of each element is
//   reserved for a checksum to detect collisions.  The remaining
//   B-1 bytes are values, prioritized by the first value.  This
//   byte is 0 to mark an unused element.
//
// bh[i] returns a pointer to the i'th element, such that
//   bh[i][0] is a checksum of i, bh[i][1] is the priority, and
//   bh[i][2..B-1] are other values (0-255).
//   The low lg(n) bits as an index into the table.
//   If a collision is detected, up to M nearby locations in the same
//   cache line are tested and the first matching checksum or
//   empty element is returned.
//   If no match or empty element is found, then the lowest priority
//   element is replaced.

// 2 byte checksum with LRU replacement (except last 2 by priority)
	template<int B> class ByteHash
	{
		enum
		{
			M = 7
		};  // search limit
		FastArray<U8, 64> t; // elements
		U32 n; // size-1
	public:
		ByteHash(int i);
		U8* operator[](U32 i);
	};

/////////////////////////// ContextMap /////////////////////////
//
// A ContextMap maps contexts to a bit histories and makes predictions
// to a Mixer.  Methods common to all classes:
//
// ContextMap cm(M, C); creates using about M bytes of memory (a power
//   of 2) for C contexts.
// cm.set(cx);  sets the next context to cx, called up to C times
//   cx is an arbitrary 32 bit value that identifies the context.
//   It should be called before predicting the first bit of each byte.
// cm.mix(m) updates Mixer m with the next prediction.  Returns 1
//   if context cx is found, else 0.  Then it extends all the contexts with
//   global bit y.  It should be called for every bit:
//
//     if (bpos==0)
//       for (int i=0; i<C; ++i) cm.set(cxt[i]);
//     cm.mix(m);
//
// The different types are as follows:
//
// - RunContextMap.  The bit history is a count of 0-255 consecutive
//     zeros or ones.  Uses 4 bytes per whole byte context.  C=1.
//     The context should be a hash.
// - SmallStationaryContextMap.  0 <= cx < M/512.
//     The state is a 16-bit probability that is adjusted after each
//     prediction.  C=1.
// - ContextMap.  For large contexts, C >= 1.  Context need not be hashed.

// A RunContextMap maps a context into the next byte and a repeat
// count up to M.  Size should be a power of 2.  Memory usage is 3M/4.
	class RunContextMap
	{
		ByteHash<4> t;
		U8 *cp;
		int mulc;
	public:
		RunContextMap(int m, int c);
		// update count
		void set(U32 cx);
		// predict next bit
		int p();
		// return run length
		int mix(Mixer& m);
	};

// Context is looked up directly.  m=size is power of 2 in bytes.
// Context should be < m/512.  High bits are discarded.
	class SmallStationaryContextMap
	{
		FastArray<U16> t;
		int cxt, mulc;
		U16 *cp;
	public:
		SmallStationaryContextMap(int m, int c);
		void set(U32 cx);
		void mix(Mixer& m);
	};

// Context map for large contexts.  Most modeling uses this type of context
// map.  It includes a built in RunContextMap to predict the last byte seen
// in the same context, and also bit-level contexts that map to a bit
// history state.
//
// Bit histories are stored in a hash table.  The table is organized into
// 64-byte buckets alinged on cache page boundaries.  Each bucket contains
// a hash chain of 7 elements, plus a 2 element queue (packed into 1 byte)
// of the last 2 elements accessed for LRU replacement.  Each element has
// a 2 byte checksum for detecting collisions, and an array of 7 bit history
// states indexed by the last 0 to 2 bits of context.  The buckets are indexed
// by a context ending after 0, 2, or 5 bits of the current byte.  Thus, each
// byte modeled results in 3 main memory accesses per context, with all other
// accesses to cache.
//
// On bits 0, 2 and 5, the context is updated and a new bucket is selected.
// The most recently accessed element is tried first, by comparing the
// 16 bit checksum, then the 7 elements are searched linearly.  If no match
// is found, then the element with the lowest priority among the 5 elements
// not in the LRU queue is replaced.  After a replacement, the queue is
// emptied (so that consecutive misses favor a LFU replacement policy).
// In all cases, the found/replaced element is put in the front of the queue.
//
// The priority is the state number of the first element (the one with 0
// additional bits of context).  The states are sorted by increasing n0+n1
// (number of bits seen), implementing a LFU replacement policy.
//
// When the context ends on a byte boundary (bit 0), only 3 of the 7 bit
// history states are used.  The remaining 4 bytes implement a run model
// as follows: <count:7,d:1> <b1> <b2> <b3> where <b1> is the last byte
// seen, possibly repeated, and <b2> and <b3> are the two bytes seen
// before the first <b1>.  <count:7,d:1> is a 7 bit count and a 1 bit
// flag.  If d=0 then <count> = 1..127 is the number of repeats of <b1>
// and no other bytes have been seen, and <b2><b3> are not used.
// If <d> = 1 then the history is <b3>, <b2>, and <count> - 2 repeats
// of <b1>.  In this case, <b3> is valid only if <count> >= 3 and
// <b2> is valid only if <count> >= 2.
//
// As an optimization, the last two hash elements of each byte (representing
// contexts with 2-7 bits) are not updated until a context is seen for
// a second time.  This is indicated by <count,d> = <1,0>.  After update,
// <count,d> is updated to <2,0> or <2,1>.

	class ContextMap
	{
		const int C, Sz;  // max number of contexts
		class E
		{  // hash element, 64 bytes
			U16 chk[7];  // byte context checksums
			U8 last;     // last 2 accesses (0-6) in low, high nibble
		public:
			U8 bh[7][7]; // byte context, 3-bit context -> bit history state
			// bh[][0] = 1st bit, bh[][1,2] = 2nd bit, bh[][3..6] = 3rd bit
			// bh[][0] is also a replacement priority, 0 = empty
			U8* get(U16 chk, int i);  // Find element (0-6) matching checksum.
			// If not found, insert or replace lowest priority (not last).
		};
		FastArray<E, 64> t;  // bit histories for bits 0-1, 2-4, 5-7
		// For 0-1, also contains a run count in bh[][4] and value in bh[][5]
		// and pending update count in bh[7]
		FastArray<U8*> cp;   // C pointers to current bit history
		FastArray<U8*> cp0; // First element of 7 element array containing cp[i]
		FastArray<U32> cxt;  // C whole byte contexts (hashes)
		FastArray<U8*> runp; // C [0..3] = count, value, unused, unused
		StateMap *sm;    // C maps of state -> p
		int cn;          // Next context to set by set()
		void update(U32 cx, int c);  // train model that context cx predicts c
		int mix1(Mixer& m, int cc, int c1, int y1);
		// mix() with global context passed as arguments to improve speed.
	public:
		ContextMap(int m, int c = 1); // m = memory in bytes, a power of 2, C = c
		void set(U32 cx);   // set next whole byte context
		int mix(Mixer& m);
	};

//////////////////////////// Predictor /////////////////////////

// A Predictor estimates the probability that the next bit of
// uncompressed data is 1.  Methods:
// p() returns P(1) as a 12 bit number (0-4095).
// update(y) trains the predictor with the actual bit (0 or 1).

	class Predictor
	{
		int pr;  // next prediction
	public:
		Predictor();
		int predict() const;
		void update();
	};

//////////////////////////// Encoder ////////////////////////////

// An Encoder does arithmetic encoding.  Methods:
// Encoder(COMPRESS, f) creates encoder for compression to archive f, which
//   must be open past any header for writing in binary mode.
// Encoder(DECOMPRESS, f) creates encoder for decompression from archive f,
//   which must be open past any header for reading in binary mode.
// code(i) in COMPRESS mode compresses bit i (0 or 1) to file f.
// code() in DECOMPRESS mode returns the next decompressed bit from file f.
//   Global y is set to the last bit coded or decoded by code().
// compress(c) in COMPRESS mode compresses one byte.
// decompress() in DECOMPRESS mode decompresses and returns one byte.
// flush() should be called exactly once after compression is done and
//   before closing f.  It does nothing in DECOMPRESS mode.
// size() returns current length of archive
// setFile(f) sets alternate source to FILE* f for decompress() in COMPRESS
//   mode (for testing transforms).
// If level (global) is 0, then data is stored without arithmetic coding.

	typedef enum
	{
		COMPRESS, DECOMPRESS
	} Mode;
	class Encoder
	{
	private:
		Predictor predictor;
		const Mode mode;       // Compress or decompress?
		FILE* archive;         // Compressed data file
		U32 x1, x2;            // Range, initially [0, 1), scaled by 2^32
		U32 x;                 // Decompress mode: last 4 input bytes of archive
		FILE *alt;             // decompress() source in COMPRESS mode

		// Compress bit y or return decompressed bit
		int code(int i = 0);

	public:
		Encoder(Mode m, FILE* f);
//Mode getMode() const {return mode;}
//long size() const {return ftell(archive);}  // length of archive so far
//void setFile(FILE* f) {alt=f;}

		// Compress one byte
		void compress(int c);

		// Decompress and return one byte
		int decompress();

		void flush();
	};

///////////////////////////// Filter /////////////////////////////////

// A Filter is an abstract class which should be implemented to perform
// various transforms to improve compression for various file types.
// A function makeFilter(filename, filetype, encoder) will create an
// appropriate derived object either by examining filename and its
// contents (for compression) and set filetype, or as specified by
// filetype (for decompression).
//
// An implementation of Filter must provide the following 4 functions:
//
// 1. protected: void encode(FILE* f, int n) const;  f will be open for
// reading in binary mode and positioned at the beginning of the file,
// which has length n.  The function should read all n bytes of f and
// write transformed data to a temporary file, FILE* tmp, a protected
// member of Filter, which will be open in "wb+" mode (as returned
// by tmpfile()).  f and tmp should be left open.
// encode() should not modify any other data members than tmp that might
// change the behaior of decode().  It is const to prevent some errors but
// there is nothing to prevent it from modifying global variables or
// objects through pointers that might be used by decode().
//
// 2. protected: int decode();  should perform the inverse translation of
// encode() one byte at a time.  decode() will be called once for each byte
// of f.  The n'th call to decode() should return the n'th byte of f (0-255).
// decode() should get its input by calling protected member
// int read(), which returns one byte of transformed data (as stored
// in tmp).  decode() should not read tmp directly.  (tmp may not be open).
//
// 3. A public constructor taking an Encoder reference and passing it to the
// Filter constructor, e.g.
//
//   class Myfilter: public Filter {
//     public: MyFilter(Encoder& en): Filter(en){} // initialize for decode()
//
// 4. A public destructor to free any resources (memory) allocated by the
// constructor using 'new'.  A destructor is not necessary if all
// memory is allocated by creating Arrays.  Remember that a new
// Filter is created for each file in the archive.
//
// In addition an implementation should modify:
//
// 5. public: static Filter* make(const char* filename, Encoder& en);
//
// to return a pointer to a new Filter of the derived type when an
// appropriate file type is detected.  A file type might be detected by
// the filename extension or by examining the file contents.  If the
// file is opened, then it should be closed before returning.
//
//
// Filter implements the following:
//
// protected: int read() tests whether tmp is NULL.  If so, it
// decompresses a byte from the Encoder en and returns it.  Otherwise
// it reads a byte from tmp.
//
// The following are public:
//
// void decompress(FILE* f, int n);  decompresses to f (open in "wb"
// mode) by calling decode() n times and writing the bytes to f.
//
// void compare(FILE* f, int n);  decompresses n bytes by calling
// decode() n times and compares with the first n bytes of f (open
// in "rb" mode).  Prints either "identical" or a message indicating
// where the files first differ.
//
// void skip(int n);  calls decode() n times, discarding the results.
//
// void compress(FILE* f, int n);  compresses f, which is
// open in "rb" mode and has length n bytes.  It first opens tmp using
// tmpfile(), then calls encode(f) to write the transformed data to tmp.
// Then it tests the transform by rewinding tmp and f, then calling decode()
// n times and comparing with n bytes of f.  If the comparison is identical
// then filetype (global, one byte) is compressed (with filetype 0
// compression) to Encoder en, tmp is rewound again and compressed to en
// with appropriate filetype (which affects the model).  Otherwise a
// warning is written, filetype is set to 0, a 0 is compressed to en,
// f is rewound and n bytes of f are compressed to en.  Then in either case,
// tmp is closed (which deletes the file).  In addition,
// decode() must read exactly every byte of tmp and not the EOF, or else
// the transform is abandoned with a warning.
//
// A derived class defaultFilter does a 1-1 transform, equivalent
// to filetype 0.

	class Filter
	{
	private:
		Encoder* en;

	protected:
		FILE* tmp;  // temporary file
		virtual void encode(FILE* f, int n) = 0;  // user suplied
		virtual int decode() = 0;          // user supplied
		Filter(const Filter*);             // copying not allowed
		Filter& operator=(const Filter&);  // assignment not allowed
		// print progress
		static void printStatus(int n);
	public:
		Filter(Encoder* e);
		virtual ~Filter()
		{
		}

		void decompress(FILE* f, int n);
		void compare(FILE* f, int n);
		void skip(int n);
		void compress(FILE* f, int n);
		static Filter* make(const char* filename, Encoder* e);
		int read();
		int reads;  // number of calls to read()
	};

///////////////// DefaultFilter ////////////

// DefaultFilter does no translation (for testing)
	class DefaultFilter: public Filter
	{
	public:
		DefaultFilter(Encoder* e);
	protected:
		// not executed if filetype is 0
		void encode(FILE* f, int n);
		int decode();
	};

#if 0
/////////////////// ExeFilter ///////////////

// ExeFilter translates E8/E9 (call/jmp) addresses from relative to
// absolute in x86 code when the absolute address is < 16MB.
// Transform is as follows:
// 1. Write input file size as 4 bytes, MSB first
// 2. Write number of bytes to transform as 4 bytes, MSB first.
// 3. Divide input into 64KB blocks, each transformed separately.
// 4. Search right to left in block for E8/E9 xx xx xx 00/FF
// 5. Replace with absoulte address xx xx xx 00/FF + offset+5 mod 2^25,
//    in range +- 2^24, LSB first

	class ExeFilter: public Filter
	{
		enum
		{	BLOCK=0x10000};  // block size
		int offset, size, q;// decode state: file offset, input size, queue size
		int end;// where to stop coding
		U8 c[5];// queue of last 5 bytes, c[0] at front
	public:
		ExeFilter(Encoder* e);
		void encode(FILE* f, int n);
		int decode();
	};
#endif

// *************************************************************************************************
// AIR: WRT CLASS GOES IN HERE ********************************************************************
// *************************************************************************************************

#define CHAR_FIRSTUPPER		64 	// for encode lower word with first capital letter
#define CHAR_UPPERWORD		7	// for encode upper word
#define CHAR_LOWERWORD		6	// for encode lower word with a few capital letter
#define CHAR_PUNCTUATION	8	// for punctuation marks modeling
#define CHAR_NOSPACE		8   // the same as CHAR_PUNCTUATION
#define CHAR_CR_LF			14
#define CHAR_ESCAPE			12	// for encode reserved chars (CHAR_ESCAPE,CHAR_FIRSTUPPER,...)
#define NGRAM_FIRST			'A'
#define NGRAM_LAST			'Z'
#define BINARY_FIRST		128
#define BINARY_LAST			255

#define AUTO_SWITCH			8	// param for !OPTION_NORMAL_TEXT_FILTER
#define WORD_MIN_SIZE		1
#define FUNCTION_CHECK_ERRORS
#define WRT_HEADER "WRT4"

#define TOLOWER(c) ((c>='A' && c<='Z')?(c+32):((upperSet[0][c]>0)?lowerSetRev[0][upperSet[0][c]]:c))
#define TOUPPER(c) ((c>='a' && c<='z')?(c-32):((lowerSet[0][c]>0)?upperSetRev[0][lowerSet[0][c]]:c))
#define TOUPPER_SET(c) ((c>='a' && c<='z')?(c-32):((lowerSet[usedSet][c]>0)?upperSetRev[usedSet][lowerSet[usedSet][c]]:c))

#define OPTION(option) ((wrt.preprocFlag & option)!=0)
#define TURN_OFF(option) ;//{if (preprocFlag & option) preprocFlag-=option;}
#define TURN_ON(option)	 ;//{preprocFlag|=option;}
#define RESET_OPTIONS 	 ;//preprocFlag=0

#define COND_BIN_FILTER(c) (((c<32)?(c!=10 && c!=13):(0)) || (c>=BINARY_FIRST))

#define PRINT_CHARS(data) ;//printf data
#define PRINT_CODEWORDS(data) ;//printf data
#define PRINT_DICT(data) ;//printf data

#define HASH_TABLE_SIZE		(1<<21)
	static int word_hash[HASH_TABLE_SIZE];

	static bool fileCorrupted;

	typedef unsigned int uint;
	typedef unsigned char uc;

// filesize() function

	static uint flen(FILE* f)
	{
		fseek(f, 0, SEEK_END);
		uint len = ftell(f);
		fseek(f, 0, SEEK_SET);
		return len;
	}

	static Filter* WRTd_filter;

	class WRT
	{
	private:

	public:

		WRT();
		~WRT();

		enum EPreprocessType
		{
			LZ77, BWT, PPM, PAQ
		};
		enum EWordType
		{
			LOWERWORD, FIRSTUPPER, UPPERWORD
		};
		enum EEOLType
		{
			UNDEFINED, CRLF, LF
		};
		enum EUpperType
		{
			UFALSE, UTRUE, FORCE
		};
		enum ESpaceType
		{
			NONE, SPACE, EOL
		};

		bool restartEnc, initOrder, forceNormalTextFilter,
				forceWordSurroroundModeling, forceEOLcoding;
		int tryShorterBound, preprocessing, s_size, WRTd_c, WRTd_qstart,
				WRTd_qend, WRTd_type;
		int fftell, fftelld, originalFileLen, autoSwitch, WRTd_binCount;
		int bufferedChar, lastEOL, EOLcount, lastChar, llast, llbckp;
		bool swapCase, WRT_verbose, WRTd_upper;
		unsigned char WRTd_s[1024];
		unsigned char WRTd_queue[128];
		EUpperType upperWord;
		EEOLType EOLType;
		ESpaceType spaceBefore;
		EPreprocessType preprocType;

#define DECODE_GETC(c,file)\
{\
	if (fftelld<originalFileLen) \
	{ \
		c=WRTd_filter->read(); \
		fftelld++; \
	} \
	else \
		c=EOF; \
}

#define ENCODE_PUTC(c,file)\
{ \
	putc(c,file); \
}

#define MAX_FREQ_ORDER1		2520
#define ORDER1_STEP	4

		int mZero[MAX_FREQ_ORDER1];
		int mOne[MAX_FREQ_ORDER1];

#define	UPDATE_ORDER1(prev,value)	UpdateOrder1(prev,value,ORDER1_STEP)
#define	ENCODE_ORDER1(prev,value)	EncodeOrder1(prev,value)
#define	DECODE_ORDER1(prev)		DecodeOrder1(prev)
#define	INIT_ORDER1			InitOrder1(MAX_FREQ_ORDER1)

#define DICTNAME_EXT ".dic"
#define DICTNAME "temp_HKCC_dict"
#define SHORT_DICTNAME "temp_HKCC_dict_sh"
#define WRT_DICT_DIR "./"

#define HASH_DOUBLE_MULT	29
#define HASH_MULT		23

		int sizeDict;
		unsigned char** dict;
		unsigned char* dictlen;
		unsigned char* dictmem;

		int ngram_hash[256][256];

#define CHARSET_COUNT		6

		int lowerSet[CHARSET_COUNT][256];
		int upperSet[CHARSET_COUNT][256];
		int lowerSetRev[CHARSET_COUNT][256];
		int upperSetRev[CHARSET_COUNT][256];
		int freeUpper[CHARSET_COUNT], freeLower[CHARSET_COUNT];
		int usedSet;

		int reservedSet[256];
		int addSymbols[256];
		int sym2codeword[256];
		int codeword2sym[256];
		int value[256];

		int dictionary, dict1size, dict2size, dict3size, dict4size,
				dict1plus2plus3, dict1plus2;
		int bound4, bound3, dict123size, dict12size;

// convert upper string to lower
		inline void toLower(unsigned char* s, int s_size)
		{
			for (int i = 0; i < s_size; i++)
				s[i] = TOLOWER(s[i]);
		}

// convert lower string to upper
		inline void toUpper(unsigned char* s, int s_size)
		{
			for (int i = 0; i < s_size; i++)
				s[i] = TOUPPER(s[i]);
		}

#define ORIGINAL_CHARSET(c)\
{\
	if (usedSet>0)\
	{\
		if (lowerSet[0][c]>0)\
			c=lowerSetRev[usedSet][lowerSet[0][c]];\
		else\
		if (upperSet[0][c]>0)\
			c=upperSetRev[usedSet][upperSet[0][c]];\
	}\
}

// make hash from string
		inline unsigned int stringHash(const unsigned char *ptr, int len)
		{
			unsigned int hash;
			for (hash = 0; len > 0; len--, ptr++)
				hash = HASH_MULT * hash + *ptr;

			return hash & (HASH_TABLE_SIZE - 1);
		}

// check if word "s" does exist in the dictionary using hash "h"
		inline int checkHashExactly(const unsigned char* s, int s_size, int h)
		{
			int i;

			i = word_hash[h];
			if (i > 0)
			{
				if (dictlen[i] != s_size || memcmp(dict[i], s, s_size) != 0)
				{
					i = word_hash[(h + s_size * HASH_DOUBLE_MULT)
							& (HASH_TABLE_SIZE - 1)];
					if (i > 0)
					{
						if (dictlen[i] != s_size
								|| memcmp(dict[i], s, s_size) != 0)
						{
							i = word_hash[(h
									+ s_size * HASH_DOUBLE_MULT
											* HASH_DOUBLE_MULT)
									& (HASH_TABLE_SIZE - 1)];
							if (i > 0)
							{
								if (dictlen[i] != s_size
										|| memcmp(dict[i], s, s_size) != 0)
									i = -1;
							}
							else
								i = -1;
						}
					}
					else
						i = -1;
				}
			}
			else
				i = -1;

			if (i > dictionary)
				i = -1;

			return i;
		}

// check if word "s" (prefix of original word) does exist in the dictionary using hash "h"
		inline int checkHash(const unsigned char* s, int s_size, int h)
		{
			int i;

			i = word_hash[h];
			if (i > 0)
			{
				if (dictlen[i] > s_size || memcmp(dict[i], s, s_size) != 0)
				{
					i = word_hash[(h + s_size * HASH_DOUBLE_MULT)
							& (HASH_TABLE_SIZE - 1)];
					if (i > 0)
					{
						if (dictlen[i] > s_size
								|| memcmp(dict[i], s, s_size) != 0)
						{
							i = word_hash[(h
									+ s_size * HASH_DOUBLE_MULT
											* HASH_DOUBLE_MULT)
									& (HASH_TABLE_SIZE - 1)];
							if (i > 0)
							{
								if (dictlen[i] > s_size
										|| memcmp(dict[i], s, s_size) != 0)
									i = -1;
							}
							else
								i = -1;
						}
					}
					else
						i = -1;
				}
			}
			else
				i = -1;

			if (i > dictionary)
				i = -1;

			return i;
		}

// check if word "s" or prefix of word "s" does exist in the dictionary using hash "h"
		inline int findShorterWord(const unsigned char* s, int s_size)
		{
			int ret, i, best;
			unsigned int hash;

			hash = 0;
			for (i = 0; i < WORD_MIN_SIZE + tryShorterBound; i++)
				hash = HASH_MULT * hash + s[i];

			best = -1;
			for (; i < s_size; i++)
			{
				ret = checkHash(s, i, hash & (HASH_TABLE_SIZE - 1));
				if (ret >= 0)
					best = ret;
				hash = HASH_MULT * hash + s[i];
			}

			return best;
		}

		inline int findShorterWordRev(const unsigned char* s, int s_size)
		{
			int ret, i;

			for (i = s_size - 1; i >= WORD_MIN_SIZE + tryShorterBound; i--)
			{
				ret = checkHash(s + s_size - i, i,
						stringHash(s + s_size - i, i));
				if (ret >= 0)
					return ret;
			}

			return -1;
		}

// encode word (should be lower case) using n-gram array (when word doesn't exist in the dictionary)
		inline void encodeAsText(unsigned char* s, int s_size, FILE* fileout)
		{
			int i, ngram;

			if (spaceBefore != NONE)
			{
				if (spaceBefore == SPACE)
					ENCODE_PUTC(' ', fileout);
				spaceBefore = NONE;
			}

			if (usedSet > 0)
			{
				for (i = 0; i < s_size; i++)
				{
					ORIGINAL_CHARSET(s[i]);
					ENCODE_PUTC(s[i], fileout);
				}
			}
			else
			{
				ngram = 0;
				for (i = 0; i < s_size;)
				{
					if (IF_OPTION(OPTION_USE_NGRAMS))
						ngram = ngram_hash[s[i]][s[i + 1]];

					if (ngram > 0 && ngram < dict1size)	///// && preprocType!=LZ77)
					{
						encodeCodeWord(ngram, fileout);
						i += 2;
					}
					else
					{
						ENCODE_PUTC(s[i], fileout);
						i++;
					}
				}
			}
		}

		inline void encodeCodeWord(int i, FILE* fileout)
		{
			int first, second, third;

			first = i - 1;

			if (first >= 80 * 49) //bound3)
			{
				first -= 80 * 49; //bound3;

				third = first / dict12size;
				first = first % dict12size;
				second = first / dict1size;
				first = first % dict1size;

				ENCODE_PUTC(sym2codeword[dict1plus2 + third], fileout);PRINT_CODEWORDS(
						("1st=%d(%d) ",sym2codeword[dict1plus2+third],third));

				ENCODE_PUTC(sym2codeword[dict1size + second], fileout);PRINT_CODEWORDS(
						("2nd=%d(%d) ",sym2codeword[dict1size+second],second));

				ENCODE_PUTC(sym2codeword[first], fileout);PRINT_CODEWORDS(
						("3rd=%d(%d) ",sym2codeword[first],first));
			}
			else if (first >= dict1size)
			{
				first -= dict1size;

				second = first / dict1size;
				first = first % dict1size;

				ENCODE_PUTC(sym2codeword[dict1size + second], fileout);PRINT_CODEWORDS(
						("1st=%d ",sym2codeword[dict1size+second]));

				ENCODE_PUTC(sym2codeword[first], fileout);PRINT_CODEWORDS(
						("2nd=%d ",sym2codeword[first]));
			}
			else
			{
				ENCODE_PUTC(sym2codeword[first], fileout);PRINT_CODEWORDS(
						("1st=%d ",sym2codeword[first]));
			}

			PRINT_CODEWORDS((" no=%d %s\n", no-1,dict[no]));
		}

// encode word "s" using dictionary
		inline void encodeWord(FILE* fileout, unsigned char* s, int s_size,
				EWordType wordType)
		{
			int i, j, d, e;
			int size = 0;
			int flagToEncode = -1;

			if (s_size < 1)
			{
				if (spaceBefore != NONE)
				{
					if (spaceBefore == SPACE)
						ENCODE_PUTC(' ', fileout);
					spaceBefore = NONE;
				}
				return;
			}

			s[s_size] = 0;

			if (wordType != LOWERWORD)
			{
				if (IF_OPTION(OPTION_CAPITAL_CONVERSION))
				{
					if (wordType == FIRSTUPPER)
					{
						flagToEncode = CHAR_FIRSTUPPER;
						s[0] = TOLOWER(s[0]);
					}
					else // wordType==UPPERWORD
					{
						flagToEncode = CHAR_UPPERWORD;
						toLower(s, s_size);
					}
				}
				else
					wordType = LOWERWORD;
			}

			if (IF_OPTION(OPTION_USE_DICTIONARY) && s_size >= WORD_MIN_SIZE)
			{
				i = checkHashExactly(s, s_size, stringHash(s, s_size));
				PRINT_CODEWORDS(("checkHashExactly i=%d %d=%s\n",i,s_size,s));

				if (i < 0 && IF_OPTION(OPTION_TRY_SHORTER_WORD))
				{
					// try to find shorter version of word in dictionary
					i = findShorterWord(s, s_size);
					j = findShorterWordRev(s, s_size);
					PRINT_CODEWORDS(("findShorterWord i=%d\n",i));

					d = e = 0;
					if (i >= 0)
						d = dictlen[i] - (i > 80) - (i > 3920) - 1;
					if (j >= 0)
						e = dictlen[j] - (j > 80) - (j > 3920) - 1;
					if (d >= e)
					{
						if (d > 0)
							size = dictlen[i];
					}
					else if (!IF_OPTION(OPTION_SPACELESS_WORDS))
					{
						i = j;
						PRINT_CODEWORDS(("findShorterWordRev i=%d\n",i));
						if (e > 0)
						{
							if (wordType != LOWERWORD)
							{
								ENCODE_PUTC(flagToEncode, fileout);
								if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
									ENCODE_PUTC(' ', fileout);
								wordType = LOWERWORD;
							}

							s[s_size - dictlen[i]] = 0;
							encodeAsText(s, s_size - dictlen[i], fileout);
							ENCODE_PUTC(CHAR_NOSPACE, fileout);
							s += dictlen[i];
							s_size -= dictlen[i];
						}
					}
				}
			}
			else
				i = -1;

			if (i >= 0)
			{
				if (wordType != LOWERWORD)
				{
					ENCODE_PUTC(flagToEncode, fileout);
					if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
						ENCODE_PUTC(' ', fileout);
				}

				if (spaceBefore == NONE)
				{
					if (IF_OPTION(OPTION_SPACELESS_WORDS))
						ENCODE_PUTC(CHAR_NOSPACE, fileout);
				}
				else
					spaceBefore = NONE;

				encodeCodeWord(i, fileout);

				if (size > 0)
					encodeAsText(s + size, s_size - size, fileout);
			}
			else
			{
				if (wordType != LOWERWORD)
				{
					if (spaceBefore != NONE)
					{
						if (spaceBefore == SPACE)
							ENCODE_PUTC(' ', fileout);
						spaceBefore = NONE;
					}

					ENCODE_PUTC(flagToEncode, fileout);
					if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
						ENCODE_PUTC(' ', fileout);
				}

				encodeAsText(s, s_size, fileout);
			}
		}

// decode word using dictionary
#define DECODE_WORD(dictNo,i)\
{\
		switch (dictNo)\
		{\
			case 4:\
				i+=bound4;\
				break;\
			case 3:\
				i+=80*49; /*bound3;*/\
				break;\
			case 2:\
				i+=dict1size;\
				break;\
		}\
\
		if (i>=0 && i<sizeDict)\
		{\
			PRINT_CODEWORDS(("i=%d ",i)); \
			i++;\
			s_size=dictlen[i];\
			memcpy(s,dict[i],s_size+1);\
			PRINT_CODEWORDS(("%s\n",dict[i])); \
		}\
		else\
		{\
			printf("File is corrupted!\n");\
			fileCorrupted=true;\
		}\
}

		inline int decodeCodeWord(FILE* file, unsigned char* s, int& c)
		{
			#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
			int i, dictNo, s_size;
			#pragma GCC diagnostic warning "-Wmaybe-uninitialized"

			if (codeword2sym[c] < dict1size)
			{
				i = codeword2sym[c];
				dictNo = 1;
				DECODE_WORD(dictNo, i);
				return s_size;
			}
			i = dict1size * (codeword2sym[c] - dict1size);

			DECODE_GETC(c, file);PRINT_CODEWORDS(("DC1 c=%d i=%d\n",c,i));

			if (codeword2sym[c] < dict1size)
			{
				i += codeword2sym[c];
				dictNo = 2;
				DECODE_WORD(dictNo, i);
				return s_size;
			}
			{
				i = (i - dict12size) * dict2size;
				PRINT_CODEWORDS(("DC2b c=%d\n",codeword2sym[c]-dict1size));
				i += dict1size * (codeword2sym[c] - dict1size);
			}

			DECODE_GETC(c, file);PRINT_CODEWORDS(("DC2 c=%d i=%d\n",c,i));

			{
				PRINT_CODEWORDS(("DC3b c=%d\n",codeword2sym[c]));
				i += codeword2sym[c];
				dictNo = 3;
				DECODE_WORD(dictNo, i);
				return s_size;
			}

		}

		unsigned char*
		loadDictionary(FILE* file, unsigned char* mem, int word_count);

		int loadCharset(FILE* file, int& freeChar, int* charset,
				int* charsetRev, bool *joinCharsets = NULL);

		void initializeCodeWords();

// read dictionary from files to arrays
		bool initialize(unsigned char* dictName, unsigned char* shortDictName,
				bool encoding);

		void deinitialize();

#define MAX_DICT_NUMBER		255
#define SAMPLE_WORDS_COUNT	250
#define SAMPLE_WORDS_COUNT_MAX	(SAMPLE_WORDS_COUNT*CHARSET_COUNT)

		int lang[MAX_DICT_NUMBER];
		int langCount, langSum;
		unsigned char* langName[MAX_DICT_NUMBER];
		int longDictLen, shortDictLen, longDict, shortDict, lastShortDict;
		bool joinCharsets[256];

		std::multimap<std::string, int> map;
		std::multimap<std::string, int>::iterator it;

		inline void checkWord(unsigned char* s, int s_size)
		{
			if (s_size < WORD_MIN_SIZE)
				return;

			std::string str;
			str.append((char*) s, s_size);

			it = map.find(str);
			if (it == map.end())
				return;

			do
			{
				lang[it->second / SAMPLE_WORDS_COUNT_MAX]++;
				langSum++;

				it++;
			} while (it != map.end() && it->first == str);

			return;
		}

		int detectFileType(FILE* file, int part_length, int parts,
				int& recordLen);

		void set_options(char c, char c2);

		void get_options(int& c, int& c2);

		int defaultSettings(int argc, char* argv[]);

		inline bool addWord(std::string s, int& sizeFullDict)
		{
			std::pair<std::string, int> pair(s, sizeFullDict++);

			map.insert(pair);

			return true;
		}

		bool readDicts(char* pattern, char* dictPath, int dictPathLen);

		void freeNames();

		int getSourcePath(char* buf, int buf_size);

		int getFileType(FILE* file, int& recordLen);

#define SWAP_CASE(c) \
{\
\
	if (swapCase)\
	{\
		if (c!=' ' && c!='\r' && c!='\n' && c!='\'' && c!='"' && c!='\t')\
			swapCase=false;\
\
		if (c>='a' && c<='z')\
		{\
			c-=32; /* toupper(c); */ \
		}\
		else\
		if (c>='A' && c<='Z')\
		{\
			c+=32; /* tolower(c); */ \
		}\
	}\
	else\
	if (c=='.' || c=='!' || c=='?')\
		swapCase=true;\
}

#define ENCODE_GETC(c,file)\
{\
	if (llbckp!=0) \
	{ \
		if (llbckp==32*32) \
		{ \
			c=32; \
			llbckp=32; \
		} \
		else \
		{ \
			c=llbckp; \
			llbckp=0; \
		} \
	} \
	else \
	{ \
		c=getc(file); \
 \
		if (IF_OPTION(OPTION_SPACE_AFTER_EOL) && llast==10) \
		{ \
			if (c==32) \
			{ \
				c=getc(file); \
 \
				if (c==32) \
					llbckp=32*32; \
			} \
			else \
			{ \
				llbckp=c; \
				c=32; \
			} \
		} \
 \
		llast=c; \
	} \
 \
 	fftell++; \
 \
	if (c>127) \
	{ \
		if (IF_OPTION(OPTION_UTF8)) \
		{ \
			if (c!=194 && c!=195) \
			{ \
				TURN_OFF(OPTION_UTF8); \
				restartEnc=true; \
				return; \
			} \
			else \
			{ \
				int c2=fgetc(file); \
				if (c2<128 || c2>191) \
				{ \
					TURN_OFF(OPTION_UTF8); \
					restartEnc=true; \
					return; \
				} \
				else \
					c=c2+(c-194)*64; \
			} \
		} \
	} \
 \
	if (IF_OPTION(OPTION_TO_LOWER_AFTER_PUNCTUATION))\
		SWAP_CASE(c);\
}

#define DECODE_QUEUE(c)\
{\
	if (IF_OPTION(OPTION_SPACE_AFTER_EOL) && llast==10) \
	{ \
		if ((c)!=32) \
		{ \
			WRTd_queue[WRTd_qend++]=32;\
			WRTd_queue[WRTd_qend++]=c;\
		} \
 \
	} \
 	else \
	{ \
		WRTd_queue[WRTd_qend++]=c;\
	} \
}

#define DECODE_PUTC(c)\
{\
	if (IF_OPTION(OPTION_TO_LOWER_AFTER_PUNCTUATION)) \
		SWAP_CASE(c);\
 \
 \
	if (c>127) \
	{ \
		if (IF_OPTION(OPTION_UTF8)) \
		{ \
			DECODE_QUEUE((c >> 6) | 0xc0); \
			c = (c & 0x3f) | 0x80; \
		} \
	} \
 \
	DECODE_QUEUE(c); \
 \
	llast=c; \
 	fftell++; \
}

// preprocess the file
		void encode(FILE* file, FILE* fileout, int fileLen);

		inline void hook_putc(int c)
		{
			if (bufferedChar < 0 && c == ' ')
			{
				bufferedChar = c;
				return;
			}

			if (bufferedChar >= 0)
			{
				DECODE_PUTC(bufferedChar);

				if (c == ' ')
				{
					lastChar = bufferedChar;
					bufferedChar = c;
					return;
				}

				bufferedChar = -1;
			}

			if (c == 10)
				lastEOL = fftell;

			lastChar = c;

			if (c == EOF)
				return;

			DECODE_PUTC(c);
			return;
		}

		int readEOLstream(FILE* file);

		void writeEOLstream(FILE* fileout);

		inline void decode(FILE* file)
		{
			PRINT_CHARS(("c=%d (%c)\n",WRTd_c,WRTd_c));

			if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && preprocessing > 0)
			{
				DECODE_PUTC(WRTd_c);

				if (!IF_OPTION(
						OPTION_NORMAL_TEXT_FILTER) && COND_BIN_FILTER(WRTd_c))
				{
					preprocessing = autoSwitch;
					WRTd_binCount++;
				}
				else
				{
					preprocessing--;
					PRINT_CHARS(
							("preprocessing=%d c=%c(%d)\n",preprocessing,WRTd_c,WRTd_c));
					if (preprocessing == 0)
					{
						initOrder = true;
						if (WRTd_binCount * 100 / (fftell + 5000) > 25)
						{
							autoSwitch = AUTO_SWITCH * 16;
							preprocessing = autoSwitch;
						}
					}
				}

				DECODE_GETC(WRTd_c, file);
				return;
			}

			if (addSymbols[WRTd_c] && IF_OPTION(OPTION_USE_DICTIONARY))
			{
				PRINT_CHARS(
						("addSymbols[c] && IF_OPTION(OPTION_USE_DICTIONARY) upperWord=%d\n",upperWord));

				if (spaceBefore == SPACE)
				{
					if (upperWord == FORCE)
						upperWord = UTRUE;
					else
						upperWord = UFALSE;
				}

				s_size = decodeCodeWord(file, WRTd_s, WRTd_c);

				if (WRTd_upper)
				{
					WRTd_upper = false;
					WRTd_s[0] = TOUPPER(WRTd_s[0]);
				}

				if (upperWord != UFALSE)
					toUpper(WRTd_s, s_size);

				if (IF_OPTION(OPTION_SPACELESS_WORDS))
				{
					if (spaceBefore == SPACE)
						hook_putc(' ');
					else
						spaceBefore = SPACE;
				}

				int i;
				PRINT_CHARS(("word="));
				for (i = 0; i < s_size; i++)
					PRINT_CHARS(("%c",WRTd_s[i]));PRINT_CHARS(
						(" upperWord=%d\n",upperWord));

				for (i = 0; i < s_size; i++)
				{
					ORIGINAL_CHARSET(WRTd_s[i]);
					hook_putc(WRTd_s[i]);
				}

				DECODE_GETC(WRTd_c, file);
				return;
			}

			if (reservedSet[WRTd_c])
			{
				PRINT_CHARS(
						("reservedSet[%d] OPTION_SPACELESS_WORDS=%d\n",WRTd_c,IF_OPTION(OPTION_SPACELESS_WORDS)));

				if (WRTd_c == CHAR_ESCAPE)
				{
					WRTd_upper = false;
					upperWord = UFALSE;

					DECODE_GETC(WRTd_c, file);PRINT_CHARS(
							("c==CHAR_ESCAPE, next=%x\n",WRTd_c));
					hook_putc(WRTd_c);

					if (!IF_OPTION(
							OPTION_NORMAL_TEXT_FILTER) && COND_BIN_FILTER(WRTd_c))
						preprocessing = autoSwitch;

					DECODE_GETC(WRTd_c, file);
					return;
				}

				if (WRTd_c == CHAR_NOSPACE)
				{
					PRINT_CHARS(("c==CHAR_NOSPACE\n"));

					if (upperWord == FORCE)
						upperWord = UTRUE;

					DECODE_GETC(WRTd_c, file);
					spaceBefore = NONE;
					return;
				}

				if (WRTd_c == CHAR_CR_LF)
				{
					PRINT_CHARS(("c==CHAR_CR_LF\n"));

					hook_putc('\r');
					WRTd_c = '\n';
				}

				if (WRTd_c == CHAR_FIRSTUPPER)
				{
					PRINT_CHARS(("c==CHAR_FIRSTUPPER\n"));

					if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
					{
						DECODE_GETC(WRTd_c, file); // skip space
					}
					WRTd_upper = true;
					upperWord = UFALSE;
					DECODE_GETC(WRTd_c, file);PRINT_CHARS(
							("c==CHAR_FIRSTUPPER WRTd_c=%d\n",WRTd_c));
					return;
				}

				if (WRTd_c == CHAR_UPPERWORD)
				{
					PRINT_CHARS(("c==CHAR_UPPERWORD\n"));

					if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
					{
						DECODE_GETC(WRTd_c, file); // skip space
					}
					upperWord = FORCE;
					DECODE_GETC(WRTd_c, file);
					return;
				}

				if (WRTd_c == CHAR_LOWERWORD)
				{
					PRINT_CHARS(("c==CHAR_LOWERWORD\n"));

					WRTd_upper = false;
					upperWord = UFALSE;
					DECODE_GETC(WRTd_c, file);
					if (WRTd_c == ' ') // skip space
					{
						DECODE_GETC(WRTd_c, file);
					}
					return;
				}
			}

			PRINT_CHARS(
					("other c=%d (%d %d)\n",WRTd_c,fftell+((bufferedChar>=0)?1:0),originalFileLen));

			if (upperWord != UFALSE)
			{
				if (upperWord == FORCE)
					upperWord = UTRUE;

				if ((WRTd_c >= 'a' && WRTd_c <= 'z')
						|| lowerSet[usedSet][WRTd_c] > 0)
					WRTd_c = TOUPPER_SET(WRTd_c);
				else
					upperWord = UFALSE;
			}
			else if (WRTd_upper == true)
			{
				WRTd_upper = false;
				WRTd_c = TOUPPER_SET(WRTd_c);
			}

			hook_putc(WRTd_c);

			DECODE_GETC(WRTd_c, file);
			return;
		}

		void start_encoding(FILE* file, FILE* fileout, unsigned int fileLen,
				bool type_detected);

		void start_decoding(FILE* file, FILE* fileout, int header);

		void prepare_decoding();

		int decode_char(FILE* file, FILE* fileout, int header);

	};
// end class

// **************************************************************************************************************
// **************************************************************************************************************
// AIR: include end
// **************************************************************************************************************
// **************************************************************************************************************

	class TextFilter: public Filter
	{
	private:
		FILE* dtmp;
		bool first;
	public:
		TextFilter(Encoder* e);
		~TextFilter();
		void encode(FILE* f, int n);
		int decode();
		void reset();
	};

	static int pg_main(int argc, char** argv);
	static void compress(FILE *inputFile, FILE *outputFile, int level);

};
// end: class

//}

#endif
