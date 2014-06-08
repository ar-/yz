/*
 * generalpacker.cxx
 * 
 * Copyright 2014 Andreas Altair Redmer <altair.ibn.la.ahad.sy@gmail.com>
 *
 * Original algorithm and reference implementation by Matt Mahoney
 * and P.Skibinski <inikep@o2.pl> released in 2006 (GPL2+). Some parts of the code
 * are from there (same licence) but it was too 'optimised' to reuse it,
 * or build it in here. This is a reimplementation, which contains bugfixes,
 * major refactorings and large improvements in readablitiy and
 * extensibility (for the price of minor performance costs).
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

#include "generalpacker.h"
#include "packer.h"

#pragma GCC diagnostic ignored "-Wparentheses"

// using namespace yztools;

//namespace yztools
//{

using namespace std;
#define NDEBUG  // remove for debugging (turns on Array bound checks)
#include <assert.h>

int GeneralPacker::y = 0;
int GeneralPacker::fileSize, GeneralPacker::fileType = DEFAULT;
long GeneralPacker::size;
int GeneralPacker::c0 = 1;
U32 GeneralPacker::b1 = 0, GeneralPacker::b2 = 0, GeneralPacker::b3 = 0,
		GeneralPacker::b4 = 0, GeneralPacker::b5 = 0, GeneralPacker::b6 = 0,
		GeneralPacker::b7 = 0, GeneralPacker::b8 = 0, GeneralPacker::tt = 0,
		GeneralPacker::c4 = 0, GeneralPacker::x4 = 0, GeneralPacker::x5 = 0,
		GeneralPacker::w4 = 0, GeneralPacker::w5 = 0, GeneralPacker::f4 = 0;
int GeneralPacker::order, GeneralPacker::bpos = 0, GeneralPacker::cxtfl = 3,
		GeneralPacker::sm_shft = 7, GeneralPacker::sm_add = 65535 + 127,
		GeneralPacker::sm_add_y = 0;

bool GeneralPacker::fileCorrupted = false;
int GeneralPacker::compressionLevel = DEFAULT_OPTION;

GeneralPacker::Filter* GeneralPacker::WRTd_filter = NULL;

// the remaining ones also initalizes to make the linker happy
int GeneralPacker::pos = 0;
GeneralPacker::Stretch GeneralPacker::stretch; //= GeneralPacker::Stretch();
GeneralPacker::ProgramObserver GeneralPacker::programObserver;
GeneralPacker::PseudoRandomGenerator GeneralPacker::getRandomNumber;
GeneralPacker::Ilog GeneralPacker::ilog;

int GeneralPacker::word_hash[HASH_TABLE_SIZE];
GeneralPacker::FastArrayBuffer GeneralPacker::buf =
		GeneralPacker::FastArrayBuffer(0);

//GeneralPacker::GeneralPacker()
//
//{
////	fileCorrupted = false;
////	WRTd_filter = NULL;
//
//}
GeneralPacker::PseudoRandomGenerator::PseudoRandomGenerator()
{
	table[0] = 123456789;
	table[1] = 987654321;
	for (int j = 0; j < 62; j++)
		table[j + 2] = table[j + 1] * 11 + table[j] * 23 / 16;
	i = 0;
}

U32 GeneralPacker::PseudoRandomGenerator::operator()()
{
	return ++i, table[i & 63] = table[i - 24 & 63] ^ table[i - 55 & 63];
}

void GeneralPacker::ProgramObserver::alloc(int n)
{  // report memory allocated, may be negative
	memused += n;
	if (memused > maxmem)
		maxmem = memused;
}
GeneralPacker::ProgramObserver::ProgramObserver() :
		memused(0), maxmem(0)
{
	start_time = clock();
	assert(sizeof(U8)==1);
	assert(sizeof(U16)==2);
	assert(sizeof(U32)==4);
	assert(sizeof(short)==2);
	assert(sizeof(int)==4);
}
void GeneralPacker::ProgramObserver::print() const
{  // print time and memory used
	printf("Time %1.2f sec, used %d bytes of memory\n",
			double(clock() - start_time) / CLOCKS_PER_SEC, maxmem);
}

template<class T, int ALIGN>
void GeneralPacker::FastArray<T, ALIGN>::resize(int i)
{
	if (i <= reserved)
	{
		n = i;
		return;
	}
	char *saveptr = ptr;
	T *savedata = data;
	int saven = n;
	create(i);
	if (savedata && saveptr)
	{
		memcpy(data, savedata, sizeof(T) * min(i, saven));
		programObserver.alloc(-ALIGN - n * sizeof(T));
		free(saveptr);
	}
}

template<class T, int ALIGN>
void GeneralPacker::FastArray<T, ALIGN>::create(int i)
{
	n = reserved = i;
	if (i <= 0)
	{
		data = 0;
		ptr = 0;
		return;
	}
	const int sz = ALIGN + n * sizeof(T);
	programObserver.alloc(sz);
	ptr = (char*) calloc(sz, 1);
	if (!ptr)
		quitAndThrowError("Out of memory");
	data =
			(ALIGN ?
					(T*) (ptr + ALIGN - (((long) ptr) & (ALIGN - 1))) : (T*) ptr);
	assert((char*)data>=ptr && (char*)data<=ptr+ALIGN);
}

template<class T, int ALIGN>
GeneralPacker::FastArray<T, ALIGN>::~FastArray()
{
	programObserver.alloc(-ALIGN - n * sizeof(T));
	free(ptr);
}

template<class T, int ALIGN>
void GeneralPacker::FastArray<T, ALIGN>::push_back(const T& x)
{
	if (n == reserved)
	{
		int saven = n;
		resize(max(1, n * 2));
		n = saven;
	}
	data[n++] = x;
}

GeneralPacker::FastArrayBuffer::FastArrayBuffer(int i = 0) :
		b(i)
{
}
void GeneralPacker::FastArrayBuffer::setsize(int i)
{
	if (!i)
		return;
	assert(i>0 && (i&(i-1))==0);
	b.resize(i);
}
U8&
GeneralPacker::FastArrayBuffer::operator[](int i)
{
	return b[i & b.size() - 1];
}
int GeneralPacker::FastArrayBuffer::operator()(int i) const
{
	assert(i>0);
	return b[pos - i & b.size() - 1];
}
int GeneralPacker::FastArrayBuffer::size() const
{
	return b.size();
}

void GeneralPacker::pack()
{
	// run it here

}

int GeneralPacker::Ilog::operator()(U16 x) const
{
	return t[x];
}
// Compute lookup table by numerical integration of 1/x
GeneralPacker::Ilog::Ilog()
{
	U32 x = 14155776;
	for (int i = 2; i < 65536; ++i)
	{
		x += 774541002 / (i * 2 - 1);  // numerator is 2^29/ln 2
		t[i] = x >> 24;
	}
}

int GeneralPacker::Stretch::operator()(int pD) const
{
	assert(pD>=0 && pD<4096); //TODO check if this is not assert(p>=0 && p<4096); What is the other p?
	return t[pD];
}

GeneralPacker::Stretch::Stretch()
{
	int pi = 0;
	for (int x = -2047; x <= 2047; ++x)
	{  // invert squash()
		int i = GeneralPacker::squash(x);
		for (int j = pi; j <= i; ++j)
			t[j] = x;
		pi = i + 1;
	}
	t[4095] = 2047;
}

#if __x86_64__ // this is GCC specific
/* 64-bit */
#else
/* 32-bit */
#endif
// noinline because jump labels will be redefined otherwise
// dot_product returns dot product t*w of n elements.  n is rounded
// up to a multiple of 8.  Result is scaled down by 8 bits.
#pragma GCC diagnostic ignored "-Wreturn-type"
int __attribute__ ((noinline))
GeneralPacker::dot_product(short *t, short *w, int n)
{
	//asm source is in intel syntax
	asm(".intel_syntax noprefix");
	asm(
			"mov rcx, rdx\n\t"
			"mov rax, rdi\n\t"
			"mov rdx, rsi\n\t"
			"add rcx, 7\n\t"
			"and rcx, -8\n\t"
			"jz .done\n\t"
			"sub rax, 16\n\t"
			"sub rdx, 16\n\t"
			"pxor xmm0, xmm0\n\t"

			".loop:\n\t"

			"movdqa xmm1, [rax+rcx*2]\n\t"
			"pmaddwd xmm1, [rdx+rcx*2]\n\t"
			"psrad xmm1, 8\n\t"
			"paddd xmm0, xmm1\n\t"
			"sub rcx, 8\n\t"
			"ja .loop\n\t"
			"movdqa xmm1, xmm0\n\t"
			"psrldq xmm1, 8\n\t"
			"paddd xmm0, xmm1\n\t"
			"movdqa xmm1, xmm0\n\t"
			"psrldq xmm1, 4\n\t"
			"paddd xmm0, xmm1\n\t"
			"movd rax, xmm0\n\t"

			".done:\n\t"
			"ret"
	);

	// now switch back to normal syntax
	asm(".att_syntax prefix");
}
#pragma GCC diagnostic warning "-Wreturn-type"

// Train neural network weights w[n] given inputs t[n] and err.
// w[i] += t[i]*err, i=0..n-1.  t, w, err are signed 16 bits (+- 32K).
// err is scaled 16 bits (representing +- 1/2).  w[i] is clamped to +- 32K
// and rounded.  n is rounded up to a multiple of 8.

void __attribute__ ((noinline))
GeneralPacker::train(short *t, short *w, int n, int err)
{
	asm(".intel_syntax noprefix");
	asm(
			"mov rax, rcx\n\t"
			"and rax, 0xffff\n\t"
			"movd xmm0, rax\n\t"
			"movd xmm1, rax\n\t"
			"pslldq xmm1, 2\n\t"
			"por xmm0, xmm1\n\t"
			"movdqa xmm1, xmm0\n\t"
			"pslldq xmm1, 4\n\t"
			"por xmm0, xmm1\n\t"
			"movdqa xmm1, xmm0\n\t"
			"pslldq xmm1, 8\n\t"
			"por xmm0, xmm1;\n\t"
			"pcmpeqb xmm1, xmm1\n\t"
			"psrlw xmm1, 15\n\t"
			"mov rcx, rdx\n\t"
			"mov rax, rdi\n\t"
			"mov rdx, rsi\n\t"
			"add rcx, 7\n\t"
			"and rcx, -8\n\t"
			"sub rax, 16\n\t"
			"sub rdx, 16\n\t"
			"jz .done\n\t"
			".loop1:\n\t"
			"movdqa xmm2, [rdx+rcx*2]\n\t"
			"movdqa xmm3, [rax+rcx*2]\n\t"
			"paddsw xmm3, xmm3\n\t"
			"pmulhw xmm3, xmm0 \n\t"
			"paddsw xmm3, xmm1\n\t"
			"psraw xmm3, 1\n\t"
			"paddsw xmm2, xmm3\n\t"
			"movdqa [rdx+rcx*2], xmm2\n\t"
			"sub rcx, 8\n\t"
			"ja .loop\n\t"

			".done1:\n\t"
			"ret"
	);

	// now switch back to normal syntax
	asm(".att_syntax prefix");
}

GeneralPacker::Mixer::Mixer(int n, int m, int s, int w) :
		N((n + 7) & -8), M(m), S(s), tx(N), wx(N * M), cxt(S), ncxt(0), base(0), numberOfInputsInTX(
				0), pr(S), mp(0)
{
	assert(n>0 && N>0 && (N&7)==0 && M>0);
	int i;
	for (i = 0; i < S; ++i)
		pr[i] = 2048;
	for (i = 0; i < N * M; ++i)
		wx[i] = w;
	if (S > 1)
		mp = new Mixer(S, 1, 1, 0x7fff);
}

void GeneralPacker::Mixer::update()
{
	for (int i = 0; i < ncxt; ++i)
	{
		int err = ((y << 12) - pr[i]) * 7;
		assert(err>=-32768 && err<32768);
		train(&tx[0], &wx[cxt[i] * N], numberOfInputsInTX, err);
	}
	numberOfInputsInTX = base = ncxt = 0;
}

void GeneralPacker::Mixer::update2()
{
	train(&tx[0], &wx[0], numberOfInputsInTX, ((y << 12) - base) * 3 / 2);
	numberOfInputsInTX = 0;
}

void GeneralPacker::Mixer::add(int x)
{
	assert(numberOfInputsInTX<N);
	tx[numberOfInputsInTX++] = x;
}

void GeneralPacker::Mixer::mul(int x)
{
	int z = tx[numberOfInputsInTX];
	z = z * x / 4;
	tx[numberOfInputsInTX++] = z;
}

void GeneralPacker::Mixer::set(int cx, int range)
{
	assert(range>=0);
	assert(ncxt<S);
	assert(cx>=0);
	assert(base+cx<M);
	cxt[ncxt++] = base + cx;
	base += range;
}

int GeneralPacker::Mixer::p()
{
	while (numberOfInputsInTX & 7)
		tx[numberOfInputsInTX++] = 0;  // pad
	if (mp)
	{  // combine outputs
		mp->update2();
		for (int i = 0; i < ncxt; ++i)
		{
			int dp = dot_product(&tx[0], &wx[cxt[i] * N], numberOfInputsInTX);
			dp = (dp * 9) >> 9;
			pr[i] = squash(dp);
			mp->add(dp);
		}
		return mp->p();
	}
	else
	{  // S=1 context
		int z = dot_product(&tx[0], &wx[0], numberOfInputsInTX);
		base = squash((z * 15) >> 13);
		return squash(z >> 9);
	}
}
GeneralPacker::Mixer::~Mixer()
{
	delete mp;
}

GeneralPacker::AProbabilityMap::AProbabilityMap(int n) :
		index(0), t(n * 33)
{
	for (int j = 0; j < 33; ++j)
		t[j] = squash((j - 16) * 128) * 16;
	for (int i = 33; i < n * 33; ++i)
		t[i] = t[i - 33];
}
int GeneralPacker::AProbabilityMap::p(int pr, int cxt, int rate)
{
	assert(pr>=0 && pr<4096 && cxt>=0 && cxt<N && rate>0 && rate<32);
	pr = stretch(pr);
	int g = (y << 16) + (y << rate) - y * 2;
	t[index] += g - t[index] >> rate;
	t[index + 1] += g - t[index + 1] >> rate;
	const int w = pr & 127;  // interpolation weight (33 points)
	index = (pr + 2048 >> 7) + cxt * 33;
	return t[index] * (128 - w) + t[index + 1] * w >> 11;
}

GeneralPacker::StateMap::StateMap() :
		cxt(0)
{
	for (int i = 0; i < 256; ++i)
	{
		int n0 = nex(i, 2);
		int n1 = nex(i, 3);
		if (n0 == 0)
			n1 *= 128;
		if (n1 == 0)
			n0 *= 128;
		t[i] = 65536 * (n1 + 1) / (n0 + n1 + 2);
	}
}
int GeneralPacker::StateMap::p(int cx)
{
	assert(cx>=0 && cx<t.size());
	int q = t[cxt];
	t[cxt] = q + (sm_add_y - q >> sm_shft);
	return t[cxt = cx] >> 4;
}

template<int B>
GeneralPacker::ByteHash<B>::ByteHash(int i) :
		t(i * B), n(i - 1)
{
	assert(B>=2 && i>0 && (i&(i-1))==0); // size a power of 2?
}

template<int B>
inline U8*
GeneralPacker::ByteHash<B>::operator[](U32 i)
{
	int chk = (i >> 16 ^ i) & 0xffff;
	i = i * M & n;
	U8 *p = &t[i * B - B];
	int j;
	for (j = 0; j < M; ++j)
	{
		p += B;
		if (p[2] == 0)
		{
			*(U16*) p = chk;
			break;
		}
		if (*(U16*) p == chk)
			break;  // found
	}
	if (j == 0)
		return p;  // front

	if (j == M)
	{
		--j;
		if ( /*M>2&&*/p[2] > p[-2])
			--j;
	}
	else
		chk = *(int*) p;
	p = &t[i * 4];
	memmove(p + 4, p, j * 4);
	*(int*) p = chk;
	return p;
}

GeneralPacker::RunContextMap::RunContextMap(int m, int c) :
		t(m / 4), mulc(c)
{
	cp = t[0] + 2;
}
void GeneralPacker::RunContextMap::set(U32 cx)
{
	if (cp[0] == 0 || cp[1] != b1)
		cp[0] = 1, cp[1] = b1;
	else if (cp[0] < 255)
		++cp[0];
	cp = t[cx] + 2;
}
int GeneralPacker::RunContextMap::p()
{
	if (cp[1] + 256 >> 8 - bpos == c0)
		return ((cp[1] >> 7 - bpos & 1) * 2 - 1) * ilog(cp[0] + 1) * mulc;
	else
		return 0;
}
int GeneralPacker::RunContextMap::mix(Mixer& m)
{
	m.add(p());
	return cp[0] != 0;
}

GeneralPacker::SmallStationaryContextMap::SmallStationaryContextMap(int m,
		int c) :
		t(m / 2), cxt(0), mulc(c)
{
	assert((m/2&m/2-1)==0); // power of 2?
	for (int i = 0; i < t.size(); ++i)
		t[i] = 32768;
	cp = &t[0];
}
void GeneralPacker::SmallStationaryContextMap::set(U32 cx)
{
	cxt = cx * 256 & t.size() - 256;
}
void GeneralPacker::SmallStationaryContextMap::mix(Mixer& m)
{
	if (pos < 4000000)
		*cp += (y << 16) - *cp + (1 << 8) >> 9;
	else
		*cp += (y << 16) - *cp + (1 << 9) >> 10;
	cp = &t[cxt + c0];
	m.add(stretch(*cp >> 4) * mulc / 32);
}

// Find or create hash element matching checksum ch
inline U8*
GeneralPacker::ContextMap::E::get(U16 ch, int j)
{
	ch += j;
	if (chk[last & 15] == ch)
		return &bh[last & 15][0];
	int b = 0xffff, bi = 0;
	for (int i = 0; i < 7; ++i)
	{
		if (chk[i] == ch)
			return last = last << 4 | i, &bh[i][0];
		int pri = bh[i][0];
		if ((last & 15) != i && last >> 4 != i && pri < b)
			b = pri, bi = i;
	}
	return last = 0xf0 | bi, chk[bi] = ch, (U8*) memset(&bh[bi][0], 0, 7);
}

// Construct using m bytes of memory for c contexts
GeneralPacker::ContextMap::ContextMap(int m, int c) :
		C(c), t(m >> 6), Sz((m >> 6) - 1), cp(c), cp0(c), cxt(c), runp(c), cn(0)
{
	assert(m>=64 && (m&m-1)==0);  // power of 2?
	assert(sizeof(E)==64);
	sm = new StateMap[C];
	for (int i = 0; i < C; ++i)
	{
		cp0[i] = cp[i] = &t[0].bh[0][0];
		runp[i] = cp[i] + 3;
	}
}

// Set the i'th context to cx
inline void GeneralPacker::ContextMap::set(U32 cx)
{
	int i = cn++;
	assert(i>=0 && i<C);
	cx = cx * 123456791 + i; // permute (don't hash) cx to spread the distribution
	cx = cx << 16 | cx >> 16;
	cxt[i] = cx * 987654323 + i;
}

// Update the model with bit y1, and predict next bit to mixer m.
// Context: cc=c0, bp=bpos, c1=buf(1), y1=y.
int GeneralPacker::ContextMap::mix1(Mixer& m, int cc, int c1, int y1)
{

	// Update model with y
	int result = 0;
	for (int i = 0; i < cn; ++i)
	{
		U8 *cpi = cp[i];
		if (cpi)
		{
			assert(cpi>=&t[0].bh[0][0] && cpi<=&t[Sz].bh[6][6]);
			assert((long(cpi)&63)>=15);
			int ns = nex(*cpi, y1);
			if (ns >= 204 && getRandomNumber() << (452 - ns >> 3))
				ns -= 4;  // probabilistic increment
			*cpi = ns;
		}

		// Update context pointers
		if (bpos > 1 && runp[i][0] == 0)
			cpi = 0;
		else if (bpos == 1 || bpos == 3 || bpos == 6)
			cpi = cp0[i] + 1 + (cc & 1);
		else if (bpos == 4 || bpos == 7)
			cpi = cp0[i] + 3 + (cc & 3);
		else
		{
			cp0[i] = cpi = t[cxt[i] + cc & Sz].get(cxt[i] >> 16, i);

			// Update pending bit histories for bits 2-7
			if (bpos == 0)
			{
				if (cpi[3] == 2)
				{
					const int c = cpi[4] + 256;
					U8 *p = t[cxt[i] + (c >> 6) & Sz].get(cxt[i] >> 16, i);
					p[0] = 1 + ((c >> 5) & 1);
					p[p[0]] = 1 + ((c >> 4) & 1);
					p[3 + ((c >> 4) & 3)] = 1 + ((c >> 3) & 1);
					p = t[cxt[i] + (c >> 3) & Sz].get(cxt[i] >> 16, i);
					p[0] = 1 + ((c >> 2) & 1);
					p[p[0]] = 1 + ((c >> 1) & 1);
					p[3 + ((c >> 1) & 3)] = 1 + (c & 1);
					cpi[6] = 0;
				}

				U8 c0 = runp[i][0];
				// Update run count of previous context
				if (c0 == 0)  // new context
					c0 = 2, runp[i][1] = c1;
				else if (runp[i][1] != c1)  // different byte in context
					c0 = 1, runp[i][1] = c1;
				else if (c0 < 254)  // same byte in context
					c0 += 2;
				runp[i][0] = c0;
				runp[i] = cpi + 3;
			}
		}

		// predict from last byte in context
		int rc = runp[i][0];  // count*2, +1 if 2 different bytes seen
		if (runp[i][1] + 256 >> 8 - bpos == cc)
		{
			int b = (runp[i][1] >> 7 - bpos & 1) * 2 - 1; // predicted bit + for 1, - for 0
			int c = ilog(rc + 1);
			if (rc & 1)
				c = (c * 15) / 4;
			else
				c *= 13;
			m.add(b * c);
		}
		else
			m.add(0);

		// predict from bit context
		result += mix2(m, cpi ? *cpi : 0, sm[i]);
		cp[i] = cpi;
	}
	if (bpos == 7)
		cn = 0;
	return result;
}

int GeneralPacker::ContextMap::mix(Mixer& m)
{
	return mix1(m, c0, b1, y);
}

//////////////////////////// Models //////////////////////////////

// All of the models below take a Mixer as a parameter and write
// predictions to it.

//////////////////////////// matchModel ///////////////////////////

// matchModel() finds the longest matching context and returns its length
int GeneralPacker::matchModel(GeneralPacker::Mixer& m)
{
	const int MAXLEN = 2047;  // longest allowed match + 1
	static GeneralPacker::FastArray<int> t(MEM); // hash table of pointers to contexts
	static int h = 0;  // hash of last 7 bytes
	static int ptr = 0;  // points to next byte of match if any
	static int len = 0;  // length of match, or 0 if no match
	static int result = 0;

	if (!bpos)
	{
		h = h * 887 * 8 + b1 + 1 & t.size() - 1;  // update context hash
		if (len)
			++len, ++ptr;
		else
		{  // find match
			ptr = t[h];
			if (ptr && pos - ptr < buf.size())
				while (buf(len + 1) == buf[ptr - len - 1] && len < MAXLEN)
					++len;
		}
		t[h] = pos;  // update hash table
		result = len;
//    if (result>0 && !(result&0xfff)) printf("pos=%d len=%d ptr=%d\n", pos, len, ptr);
	}

	// predict
	if (len > MAXLEN)
		len = MAXLEN;
	int sgn;
	if (len && b1 == buf[ptr - 1] && c0 == buf[ptr] + 256 >> 8 - bpos)
	{
		if (buf[ptr] >> 7 - bpos & 1)
			sgn = 8;
		else
			sgn = -8;
	}
	else
		sgn = len = 0;
	m.add(sgn * ilog(len));
	m.add(sgn * 8 * min(len, 32));
	return result;
}

#if 0
//////////////////////////// picModel //////////////////////////

// Model a 1728 by 2376 2-color CCITT bitmap image, left to right scan,
// MSB first (216 bytes per row, 513216 bytes total).  Insert predictions
// into m.

void GeneralPacker::picModel(GeneralPacker::Mixer& m)
{
	static U32 r0, r1, r2, r3;  // last 5 rows, bit 8 is over current pixel
	static GeneralPacker::FastArray<U8> t(0x10200);// model: cxt -> state
	const int N=3;// number of contexts
	static int cxt[N];// contexts
	static GeneralPacker::StateMap sm[N];
	int i;

	// update the model
	for (i=0; i<N; ++i)
	t[cxt[i]]=nex(t[cxt[i]],y);

	// update the contexts (pixels surrounding the predicted one)
	r0+=r0+y;
	r1+=r1+((buf(215)>>(7-bpos))&1);
	r2+=r2+((buf(431)>>(7-bpos))&1);
	r3+=r3+((buf(647)>>(7-bpos))&1);
	cxt[0]=r0&0x7|r1>>4&0x38|r2>>3&0xc0;
	cxt[1]=0x100+(r0&1|r1>>4&0x3e|r2>>2&0x40|r3>>1&0x80);
	cxt[2]=0x200+(r0&0x3f^r1&0x3ffe^r2<<2&0x7f00^r3<<5&0xf800);

	// predict
	for (i=0; i<N; ++i)
	m.add(stretch(sm[i].predict(t[cxt[i]])));
}
#endif

static U32 col, frstchar = 0, spafdo = 0, spaces = 0, spacecount = 0, words = 0,
		wordcount = 0, fails = 0, failz = 0, failcount = 0;

//////////////////////////// wordModel /////////////////////////

// Model English text (words and columns/end of line)
void GeneralPacker::wordModel(GeneralPacker::Mixer& m)
{
	static U32 word0 = 0, word1 = 0, word2 = 0, word3 = 0, word4 = 0;  // hashes
	static GeneralPacker::ContextMap cm(MEM * 64, 46);
	static int nl1 = -3, nl = -2;  // previous, current newline position
	static U32 t1[256];
	static U16 t2[0x10000];

	// Update word hashes
	if (bpos == 0)
	{
		U32 c = b1, f = 0;

		if (spaces & 0x80000000)
			--spacecount;
		if (words & 0x80000000)
			--wordcount;
		spaces = spaces * 2;
		words = words * 2;

		if ((c - 'a') <= ('z' - 'a') || c == 8 || c == 6
				|| (c > 127 && b2 != 12))
		{
			++words, ++wordcount;
			word0 = word0 * 263 * 8 + c;
		}
		else
		{
			if (c == 32 || c == 10)
			{
				++spaces, ++spacecount;
				if (c == 10)
					nl1 = nl, nl = pos - 1;
			}
			if (word0)
			{
				word4 = word3 * 43;
				word3 = word2 * 47;
				word2 = word1 * 53;
				word1 = word0 * 83;
				word0 = 0;
				if (c == '.' || c == 'O' || c == ('}' - '{' + 'P'))
					f = 1, spafdo = 0;
				else
				{
					++spafdo;
					spafdo = min(63, spafdo);
				}
			}
		}

		U32 h = word0 * 271 + c;
		cm.set(word0);
		cm.set(h + word1);
		cm.set(word0 * 91 + word1 * 89);
		cm.set(h + word1 * 79 + word2 * 71);

		cm.set(h + word2);
		cm.set(h + word3);
		cm.set(h + word4);
		cm.set(h + word1 * 73 + word3 * 61);
		cm.set(h + word2 * 67 + word3 * 59);

		if (f)
		{
			word4 = word3 * 31;
			word3 = word2 * 37;
			word2 = word1 * 41;
			word1 = '.';
		}

		cm.set(b3 | b4 << 8);
		cm.set(spafdo * 8 * ((w4 & 3) == 1));

		col = min(31, pos - nl);
		if (col <= 2)
		{
			if (col == 2)
				frstchar = min(c, 96);
			else
				frstchar = 0;
		}
		if (frstchar == '[' && c == 32)
		{
			if (b3 == ']' || b4 == ']')
				frstchar = 96;
		}
		cm.set(frstchar << 11 | c);

		int above = buf[nl1 + col]; // text column context

		// Text column models
		cm.set(col << 16 | c << 8 | above);
		cm.set(col << 8 | c);
		cm.set(col * (c == 32));

		h = wordcount * 64 + spacecount;
		cm.set(spaces & 0x7fff);
		cm.set(frstchar << 7);
		cm.set(spaces & 0xff);
		cm.set(c * 64 + spacecount / 2);
		cm.set((c << 13) + h);
		cm.set(h);

		U32 d = c4 & 0xffff;
		h = w4 << 6;
		cm.set(c + (h & 0xffffff00));
		cm.set(c + (h & 0x00ffff00));
		cm.set(c + (h & 0x0000ff00));
		h <<= 6;
		cm.set(d + (h & 0xffff0000));
		cm.set(d + (h & 0x00ff0000));
		h <<= 6, f = c4 & 0xffffff;
		cm.set(f + (h & 0xff000000));

		U16& r2 = t2[f >> 8];
		r2 = r2 << 8 | c;
		U32& r1 = t1[d >> 8];
		r1 = r1 << 8 | c;
		U32 t = c | t1[c] << 8;
		cm.set(t & 0xffff);
		cm.set(t & 0xffffff);
		cm.set(t);
		cm.set(t & 0xff00);
		t = d | t2[d] << 16;
		cm.set(t & 0xffffff);
		cm.set(t);

		cm.set(x4 & 0x00ff00ff);
		cm.set(x4 & 0xff0000ff);
		cm.set(x4 & 0x00ffff00);
		cm.set(c4 & 0xff00ff00);
		cm.set(c + b5 * 256 + (1 << 17));
		cm.set(c + b6 * 256 + (2 << 17));
		cm.set(b4 + b8 * 256 + (4 << 17));

		cm.set(d);
		cm.set(w4 & 15);
		cm.set(f4);
		cm.set((w4 & 63) * 128 + (5 << 17));
		cm.set(d << 9 | frstchar);
		cm.set((f4 & 0xffff) << 11 | frstchar);
	}
	cm.mix(m);
}

//////////////////////////// recordModel ///////////////////////

// Model 2-D data with fixed record length.  Also order 1-2 models
// that include the distance to the last match.

void GeneralPacker::recordModel(GeneralPacker::Mixer& m)
{
	static int cpos1[256]; //, cpos2[256], cpos3[256], cpos4[256]; //buf(1)->last 3 pos
	static int wpos1[0x10000]; // buf(1..2) -> last position
///  static int rlen=2, rlen1=3, rlen2=4;  // run length and 2 candidates
///  static int rcount1=0, rcount2=0;  // candidate counts
	static GeneralPacker::ContextMap cm(32768 / 4, 2), cn(32768 / 2, 5), co(
			32768, 4), cp(32768 * 2, 3), cq(32768 * 4, 3);

	// Find record length
	if (!bpos)
	{
		int c = b1, w = (b2 << 8) + c, d = w & 0xf0ff, e = c4 & 0xffffff;
#if 0
		int r=pos-cpos1[c];
		if (r>1 && r==cpos1[c]-cpos2[c]
				&& r==cpos2[c]-cpos3[c] && r==cpos3[c]-cpos4[c]
				&& (r>15 || (c==buf(r*5+1)) && c==buf(r*6+1)))
		{
			if (r==rlen1) ++rcount1;
			else if (r==rlen2) ++rcount2;
			else if (rcount1>rcount2) rlen2=r, rcount2=1;
			else rlen1=r, rcount1=1;
		}
		if (rcount1>15 && rlen!=rlen1) rlen=rlen1, rcount1=rcount2=0;
		if (rcount2>15 && rlen!=rlen2) rlen=rlen2, rcount1=rcount2=0;

		// Set 2 dimensional contexts
		assert(rlen>0);
#endif
		cm.set(c << 8 | (min(255, pos - cpos1[c]) / 4));
		cm.set(w << 9 | llog(pos - wpos1[w]) >> 2);
///    cm.set(rlen|buf(rlen)<<10|buf(rlen*2)<<18);
		cn.set(w);
		cn.set(d << 8);
		cn.set(c << 16);
		cn.set((f4 & 0xffff) << 3);
		int col = pos & 3;
		cn.set(col | 2 << 12);

		co.set(c);
		co.set(w << 8);
		co.set(w5 & 0x3ffff);
		co.set(e << 3);

		cp.set(d);
		cp.set(c << 8);
		cp.set(w << 16);

		cq.set(w << 3);
		cq.set(c << 19);
		cq.set(e);

		// update last context positions
///    cpos4[c]=cpos3[c];
///    cpos3[c]=cpos2[c];
///    cpos2[c]=cpos1[c];
		cpos1[c] = pos;
		wpos1[w] = pos;
	}
	co.mix(m);
	cp.mix(m);
	cxtfl = 0;
	cm.mix(m);
	cn.mix(m);
	cq.mix(m);
	cxtfl = 3;
}

//////////////////////////// sparseModel ///////////////////////

// Model order 1-2 contexts with gaps.
void GeneralPacker::sparseModel(GeneralPacker::Mixer& m)
{
	static GeneralPacker::ContextMap cn(MEM * 2, 5);
	static GeneralPacker::SmallStationaryContextMap scm1(0x20000, 17), scm2(
			0x20000, 12), scm3(0x20000, 12), scm4(0x20000, 13), scm5(0x10000,
			12), scm6(0x20000, 12), scm7(0x2000, 12), scm8(0x8000, 13), scm9(
			0x1000, 12), scma(0x10000, 16);

	if (bpos == 0)
	{
		cn.set(words & 0x1ffff);
		cn.set((f4 & 0x000fffff) * 7);
		cn.set((x4 & 0xf8f8f8f8) + 3);
		cn.set((tt & 0x00000fff) * 9);
		cn.set((x4 & 0x80f0f0ff) + 6);
		scm1.set(b1);
		scm2.set(b2);
		scm3.set(b3);
		scm4.set(b4);
		scm5.set(words & 127);
		scm6.set((words & 12) * 16 + (w4 & 12) * 4 + (b1 >> 4));
		scm7.set(w4 & 15);
		scm8.set(spafdo * ((w4 & 3) == 1));
		scm9.set(col * (b1 == 32));
		scma.set(frstchar);
	}
	cn.mix(m);
	scm1.mix(m);
	scm2.mix(m);
	scm3.mix(m);
	scm4.mix(m);
	scm5.mix(m);
	scm6.mix(m);
	scm7.mix(m);
	scm8.mix(m);
	scm9.mix(m);
	scma.mix(m);
}

int primes[] =
{ 0, 257, 251, 241, 239, 233, 229, 227, 223, 211, 199, 197, 193, 191 };
static U32 WRT_mpw[16] =
{ 3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 }, tri[4] =
{ 0, 4, 3, 7 }, trj[4] =
{ 0, 6, 6, 12 };
static U32 WRT_mtt[16] =
{ 0, 0, 1, 2, 3, 4, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7 };

//////////////////////////// contextModel //////////////////////

// file types (order is important: the last will be sorted by filetype detection as the first)

// This combines all the context models with a Mixer.

int GeneralPacker::contextModel2()
{
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

	static GeneralPacker::ContextMap cm(MEM * 32, 7);
	static GeneralPacker::RunContextMap rcm7(MEM / 4, 14), rcm9(MEM / 4, 18),
			rcm10(MEM / 2, 20);
	static GeneralPacker::Mixer m(456, 128 * (16 + 14 + 14 + 12 + 14 + 16), 6,
			512);
	static U32 cxt[16];  // order 0-11 contexts
	static Filetype filetype = DEFAULT;
	static int size = 0;  // bytes remaining in block
//  static const char* typenames[4]={"", "jpeg ", "exe ", "text "};

	// Parse filetype and size
	if (bpos == 0)
	{
		--size;
		if (size == -1)
			filetype = (Filetype) b1;
		if (size == -5)
		{
			size = c4;
//      if (filetype<=3) printf("(%s%d)", typenames[filetype], size);
/////      if (filetype==EXE) size+=8;
		}
	}

	m.update();
	m.add(64);

	// Test for special file types
/////  int isjpeg=jpegModel(m);  // 1 if JPEG is detected, else 0
	int ismatch = matchModel(m);  // Length of longest matching context
/////  int isbmp=bmpModel(m);  // Image width (bytes) if BMP or TIFF detected, or 0

	//if (ismatch>1024) {   // Model long matches directly
	//  m.set(0, 8);
	//  return m.p();
	//}
/////  else if (isjpeg) {
/////    m.set(1, 8);
/////    m.set(c0, 256);
/////    m.set(buf(1), 256);
/////    return m.p();
/////  }
/////  else if (isbmp>0) {
/////    static int col=0;
/////    if (++col>=24) col=0;
/////    m.set(2, 8);
/////    m.set(col, 24);
/////    m.set(buf(isbmp)+buf(3)>>4, 32);
/////    m.set(c0, 256);
/////    return m.p();
/////  }

	// Normal model
	if (bpos == 0)
	{
		int i = 0, f2 = buf(2);

		if (f2 == '.' || f2 == 'O' || f2 == 'M' || f2 == '!' || f2 == ')'
				|| f2 == ('}' - '{' + 'P'))
		{
#pragma GCC diagnostic ignored "-Wsign-compare"
			if (b1 != f2 && buf(3) != f2)
#pragma GCC diagnostic warning "-Wsign-compare"
				i = 13, x4 = x4 * 256 + f2;
		}

		for (; i > 0; --i)  // update order 0-11 context hashes
			cxt[i] = cxt[i - 1] * primes[i];

		for (i = 13; i > 0; --i)  // update order 0-11 context hashes
			cxt[i] = cxt[i - 1] * primes[i] + b1;

		cm.set(cxt[3]);
		cm.set(cxt[4]);
		cm.set(cxt[5]);
		cm.set(cxt[6]);
		cm.set(cxt[8]);
		cm.set(cxt[13]);
		cm.set(0);

		rcm7.set(cxt[7]);
		rcm9.set(cxt[9]);
		rcm10.set(cxt[11]);

		x4 = x4 * 256 + b1;
	}
	rcm7.mix(m);
	rcm9.mix(m);
	rcm10.mix(m);
	int qq = m.numberOfInputsInTX;
	order = cm.mix(m) - 1;
	if (order < 0)
		order = 0;
	int zz = (m.numberOfInputsInTX - qq) / 7;

	m.numberOfInputsInTX = qq + zz * 3;
	for (qq = zz * 2; qq != 0; --qq)
		m.mul(5);
	for (qq = zz; qq != 0; --qq)
		m.mul(6);
	for (qq = zz; qq != 0; --qq)
		m.mul(9);

	if (compressionLevel >= 4)
	{
		wordModel(m);
		sparseModel(m);
		recordModel(m);
		//indirectModel(m);
/////    if (filetype==EXE) exeModel(m);
/////    if (fsize==513216) picModel(m);
	}

//m.tx[420]=0;
//m.tx[425]=0;
//m.tx[430]=0;
	U32 c1 = b1, c2 = b2, c;
	if (c1 == 9 || c1 == 10 || c1 == 32)
		c1 = 16;
	if (c2 == 9 || c2 == 10 || c2 == 32)
		c2 = 16;

	m.set(256 * order + (w4 & 240) + (c2 >> 4), 256 * 7);

	c = (words >> 1) & 63;
	m.set((w4 & 3) * 64 + c + order * 256, 256 * 7);

	c = (w4 & 255) + 256 * bpos;
	m.set(c, 256 * 8);

	if (bpos)
	{
		c = c0 << (8 - bpos);
		if (bpos == 1)
			c += b3 / 2;
		c = (min(bpos, 5)) * 256 + (tt & 63) + (c & 192);
	}
	else
		c = (words & 12) * 16 + (tt & 63);
	m.set(c, 1536);

	c = bpos;
	c2 = (c0 << (8 - bpos)) | (c1 >> bpos);
	m.set(order * 256 + c + (c2 & 248), 256 * 7);

	c = c * 256 + ((c0 << (8 - bpos)) & 255);
	c1 = (words << bpos) & 255;
	m.set(c + (c1 >> bpos), 2048);

	return m.p();
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

}

GeneralPacker::Predictor::Predictor() :
		pr(2048)
{
}
int GeneralPacker::Predictor::predict() const
{
	assert(pr>=0 && pr<4096);
	return pr;
}

void GeneralPacker::Predictor::update()
{
	static AProbabilityMap a1(256), a2(0x8000), a3(0x8000), a4(0x20000), a5(
			0x10000), a6(0x10000);

	// Update global context: pos, bpos, c0, c4, buf
	c0 += c0 + y;
	if (c0 >= 256)
	{
		buf[pos++] = c0;
		c0 -= 256;
		if (pos <= 1024 * 1024)
		{
			if (pos == 1024 * 1024)
				sm_shft = 9, sm_add = 65535 + 511;
			if (pos == 512 * 1024)
				sm_shft = 8, sm_add = 65535 + 255;
			sm_add_y = sm_add & (-y);
		}
		int i = WRT_mpw[c0 >> 4];
		w4 = w4 * 4 + i;
		if (b1 == 12)
			i = 2;
		w5 = w5 * 4 + i;
		b8 = b7, b7 = b6, b6 = b5, b5 = b4, b4 = b3;
		b3 = b2;
		b2 = b1;
		b1 = c0;
		if (c0 == '.' || c0 == 'O' || c0 == 'M' || c0 == '!' || c0 == ')'
				|| c0 == ('}' - '{' + 'P'))
		{
			w5 = (w5 << 8) | 0x3ff, x5 = (x5 << 8) + c0, f4 = (f4 & 0xfffffff0)
					+ 2;
			if (c0 != '!' && c0 != 'O')
				w4 |= 12;
			if (c0 != '!')
				b2 = '.', tt = (tt & 0xfffffff8) + 1;
		}
		c4 = (c4 << 8) + c0;
		x5 = (x5 << 8) + c0;
		if (c0 == 32)
			--c0;
		f4 = f4 * 16 + (c0 >> 4);
		tt = tt * 8 + WRT_mtt[c0 >> 4];
		c0 = 1;
	}
	bpos = (bpos + 1) & 7;

	if (fails & 0x00000080)
		--failcount;
	fails = fails * 2;
	failz = failz * 2;
	if (y)
		pr ^= 4095;
	if (pr >= 1820)
		++fails, ++failcount;
	if (pr >= 848)
		++failz;

	// Filter the context model with APMs
	pr = contextModel2();

	int rate = 6 + (pos > 14 * 256 * 1024) + (pos > 28 * 512 * 1024);
	int pt, pu = a1.p(pr, c0, 3) + 7 * pr + 4 >> 3, pv, pz = failcount + 1;
	pz += tri[(fails >> 5) & 3];
	pz += trj[(fails >> 3) & 3];
	pz += trj[(fails >> 1) & 3];
	if (fails & 1)
		pz += 8;
	pz = pz / 2;

	pu = a4.p(pu,
			(c0 * 2) ^ hash(b1, (x5 >> 8) & 255, (x5 >> 16) & 0x80ff) & 0x1ffff,
			rate);
	pv = a2.p(pr, (c0 * 8) ^ hash(29, failz & 2047) & 0x7fff, rate + 1);
	pv = a5.p(pv, hash(c0, w5 & 0xfffff) & 0xffff, rate);
	pt = a3.p(pr, (c0 * 32) ^ hash(19, x5 & 0x80ffff) & 0x7fff, rate);
	pz = a6.p(pu, (c0 * 4) ^ hash(min(9, pz), x5 & 0x80ff) & 0xffff, rate);

	if (fails & 255)
		pr = pt * 6 + pu + pv * 11 + pz * 14 + 16 >> 5;
	else
		pr = pt * 4 + pu * 5 + pv * 12 + pz * 11 + 16 >> 5;
}

int GeneralPacker::Encoder::code(int i)
{
	int p = predictor.predict();
	assert(predict>=0 && predict<4096);
	p += p < 2048;
	U32 xmid = x1 + (x2 - x1 >> 12) * p + ((x2 - x1 & 0xfff) * p >> 12);
	assert(xmid>=x1 && xmid<x2);
	if (mode == DECOMPRESS)
		y = x <= xmid;
	else
		y = i;
	y ? (x2 = xmid, p = sm_add) : (x1 = xmid + 1, p = 0);
	sm_add_y = p;
	predictor.update();
	while (((x1 ^ x2) & 0xff000000) == 0)
	{  // pass equal leading bytes of range
		if (mode == COMPRESS)
			putc(x2 >> 24, archive);
		x1 <<= 8;
		x2 = (x2 << 8) + 255;
		if (mode == DECOMPRESS)
			x = (x << 8) + (getc(archive) & 255);  // EOF is OK
	}
	return y;
}

GeneralPacker::Encoder::Encoder(Mode m, FILE* f) :
		mode(m), archive(f), x1(0), x2(0xffffffff), x(0), alt(0)
{
	if (compressionLevel > 0 && mode == DECOMPRESS)
	{  // x = first 4 bytes of archive
		for (int i = 4; i != 0; --i)
			x = (x << 8) + (getc(archive) & 255);
	}
}

void GeneralPacker::Encoder::compress(int c)
{
	assert(mode==COMPRESS);
	if (compressionLevel == 0)
	{
		putc(c, archive);
		return;
	}
	if (c >= '{' && c < 127)
		c += 'P' - '{';
	else if (c >= 'P' && c < 'T')
		c -= 'P' - '{';
	else if ((c >= ':' && c <= '?') || (c >= 'J' && c <= 'O'))
		c ^= 0x70;
	if (c == 'X' || c == '`')
		c ^= 'X' ^ '`';
	for (int i = 7; i >= 0; --i)
		code((c >> i) & 1);
}

int GeneralPacker::Encoder::decompress()
{
	if (mode == COMPRESS)
	{
		assert(alt);
		return getc(alt);
	}
	else if (compressionLevel == 0)
		return getc(archive);
	else
	{
		int c = 0;
		for (int i = 8; i != 0; --i)
			c += c + code();
		if (c >= '{' && c < 127)
			c += 'P' - '{';
		else if (c >= 'P' && c < 'T')
			c -= 'P' - '{';
		else if ((c >= ':' && c <= '?') || (c >= 'J' && c <= 'O'))
			c ^= 0x70;
		if (c == 'X' || c == '`')
			c ^= 'X' ^ '`';
		return c;
	}
}

void GeneralPacker::Encoder::flush()
{
	if (mode == COMPRESS && compressionLevel > 0)
		putc(x1 >> 24, archive);  // Flush first unequal byte of range
}

int GeneralPacker::Filter::read()
{
	return ++reads, tmp ? getc(tmp) : en->decompress();
}

void GeneralPacker::Filter::printStatus(int n)
{
	if (n > 0 && !(n & 0x3fff))
		printf("%10d \b\b\b\b\b\b\b\b\b\b\b", n), fflush(stdout);
}

GeneralPacker::Filter::Filter(Encoder* e) :
		en(e), reads(0), tmp(0)
{
}

void GeneralPacker::Filter::decompress(FILE* f, int n)
{
	for (int i = 0; i < n; ++i)
	{
		putc(decode(), f);
		printStatus(i);
	}
	printf("extracted  \n");
}

void GeneralPacker::Filter::compare(FILE* f, int n)
{
	bool found = false;  // mismatch?
	for (int i = 0; i < n; ++i)
	{
		printStatus(i);
		int c1 = found ? EOF : getc(f);
		int c2 = decode();
		if (c1 != c2 && !found)
		{
			printf("differ at %d: file=%d archive=%d\n", i, c1, c2);
			found = true;
		}
	}
	if (!found && getc(f) != EOF)
		printf("file is longer\n");
	else if (!found)
		printf("identical  \n");
}

void GeneralPacker::Filter::skip(int n)
{
	for (int i = 0; i < n; ++i)
	{
		decode();
		printStatus(i);
	}
	printf("skipped    \n");
}

// Compress n bytes of f. If filetype > 0 then transform
void GeneralPacker::Filter::compress(FILE* f, int n)
{

	// No transform
	if (!fileType)
	{
		en->compress(0);
		for (int i = 0; i < n; ++i)
		{
			printStatus(i);
			en->compress(getc(f));
		}
		return;
	}

	// Try transform
	tmp = tmpfile();
	if (!tmp)
		perror("tmpfile"), exit(1);
	encode(f, n);

	// Test transform.  decode() should produce copy of input and read
	// exactly all of tmp but not EOF.
	int tmpn = ftell(tmp);
	rewind(tmp);
	rewind(f);

	bool found = false;  // mismatch?
#ifndef NO_UNWRT_CHECK
	for (int i = 0; i < n; ++i)
	{
		int c1 = getc(f);
		int c2 = decode();
		if (c1 != c2 && !found)
		{
			found = true;
			printf("filter %d failed at %d: %d -> %d, skipping\n", fileType, i,
					c1, c2);
			break;
		}
	}

	if (!found && tmpn != reads)
		printf("filter %d reads %d/%d bytes, skipping\n", fileType, reads,
				tmpn);
	if (found || tmpn != reads)
	{  // bug in Filter, code untransformed
		rewind(f);
		fileType = 0;
		en->compress(0);
		for (int i = 0; i < n; ++i)
		{
			printStatus(i);
			en->compress(getc(f));
		}
	}
	else
#endif
	{  // transformed
		rewind(tmp);
		///if (filetype==EXE) printf("-> %d (x86) ", tmpn);
		///else
		if (fileType == TEXT)
			printf("-> %d (text) ", tmpn);
		else if (fileType == BINTEXT)
			printf("-> %d (binary+text) ", tmpn);

		int t = fileType;
		fileType = 0;
		en->compress(t);
		fileType = t;
		for (int i = 0; i < tmpn; ++i)
		{
			printStatus(i);
			en->compress(getc(tmp));
		}
	}
	fclose(tmp);
}

GeneralPacker::DefaultFilter::DefaultFilter(Encoder* e) :
		Filter(e)
{
}
void GeneralPacker::DefaultFilter::encode(FILE* f, int n)
{
	while (n--)
		putc(getc(f), tmp);
}
int GeneralPacker::DefaultFilter::decode()
{
	return read();
}

#if 0

GeneralPacker::ExeFilter::ExeFilter(Encoder* e): Filter(e), offset(-8), size(0), q(0), end(0)
{}

void GeneralPacker::ExeFilter::encode(FILE* f, int n)
{
	FastArray<U8> buf(BLOCK);
	fprintf(tmp, "%c%c%c%c", n>>24, n>>16, n>>8, n);  // size, MSB first

	// Scan for jpeg and mark end
	U32 buf1=0, buf2=0;
	int end;
	for (end=0; end<n; ++end)
	{
		buf2=buf2<<8|buf1>>24;
		buf1=buf1<<8|getc(f);
		if (buf2==0xffe00010 && buf1==0x4a464946)
		{  // APP0 16 "JFIF"
			end-=8;
			break;
		}
	}
	fprintf(tmp, "%c%c%c%c", end>>24, end>>16, end>>8, end);
	rewind(f);

	// Transform
	for (int offset=0; offset<n; offset+=BLOCK)
	{
		int bytesRead=fread(&buf[0], 1, BLOCK, f);
		for (int i=bytesRead-1; i>=4; --i)
		{
			if ((buf[i-4]==0xe8||buf[i-4]==0xe9) && (buf[i]==0||buf[i]==0xff)
					&& offset+i+1<end)
			{
				int a=(buf[i-3]|buf[i-2]<<8|buf[i-1]<<16|buf[i]<<24)+offset+i+1;
				a<<=7;
				a>>=7;
				buf[i]=a>>24;
				buf[i-1]=a>>16;
				buf[i-2]=a>>8;
				buf[i-3]=a;
			}
		}
		fwrite(&buf[0], 1, bytesRead, tmp);
	}
}

int GeneralPacker::ExeFilter::decode()
{

	// Read file size and end from first 8 bytes, MSB first
	while (offset<-4)
	size=size*256+read(), ++offset;
	while (offset<0)
	end=end*256+read(), ++offset;

	// Fill queue
	while (offset<size && q<5)
	{
		memmove(c+1, c, 4);
		c[0]=read();
		++q;
		++offset;
	}

	// E8E9 transform: E8/E9 xx xx xx 00/FF -> subtract location from x
	if (q==5 && (c[4]==0xe8||c[4]==0xe9) && (c[0]==0||c[0]==0xff)
			&& ((offset-1^offset-5)&-BLOCK)==0// not crossing block boundary
			&& offset<end)
	{
		int a=(c[3]|c[2]<<8|c[1]<<16|c[0]<<24)-offset;
		a<<=7;
		a>>=7;
		c[3]=a;
		c[2]=a>>8;
		c[1]=a>>16;
		c[0]=a>>24;
	}

	// return oldest byte in queue
	return c[--q];
}
#endif

//#include "textfilter.hpp"

// **************************************************************************************************************
// **************************************************************************************************************
// AIR: include start
// **************************************************************************************************************
// **************************************************************************************************************

GeneralPacker::WRT::WRT() :
		WRT_verbose(false), preprocType(PAQ), dict(NULL), dictlen(NULL), dictmem(
		NULL), langCount(0), lastShortDict(-1), restartEnc(false)
{
}
;
GeneralPacker::WRT::~WRT()
{
	deinitialize();
	freeNames();
}

unsigned char*
GeneralPacker::WRT::loadDictionary(FILE* file, unsigned char* mem,
		int word_count)
{
	unsigned char* word;
	int i, j, c, collision, bound;

	collision = 0;
	bound = sizeDict + word_count;

	while (!feof(file))
	{
		word = mem;
		do
		{
			c = getc(file);
			if (usedSet == CHARSET_COUNT - 1)
			{
				if (lowerSet[0][c] > 0)
					c = lowerSetRev[CHARSET_COUNT - 1][lowerSet[0][c]];
			}
			word[0] = c;
			word++;
		} while (c > 32);

		if (c == EOF)
			break;
		if (c == '\r')
			c = getc(file);

		word[-1] = 0;
		i = word - mem - 1;

		dictlen[sizeDict] = i;
		dict[sizeDict] = mem;

		j = stringHash(mem, i);
		mem += (i / 4 + 1) * 4;

		if (word_hash[j] != 0)
		{
			if (dictlen[sizeDict] != dictlen[word_hash[j]]
					|| memcmp(dict[sizeDict], dict[word_hash[j]],
							dictlen[sizeDict]) != 0)
			{
				c = (j + i * HASH_DOUBLE_MULT) & (HASH_TABLE_SIZE - 1);
				if (word_hash[c] != 0)
				{
					if (dictlen[sizeDict] != dictlen[word_hash[c]]
							|| memcmp(dict[sizeDict], dict[word_hash[c]],
									dictlen[sizeDict]) != 0)
					{
						c = (j + i * HASH_DOUBLE_MULT * HASH_DOUBLE_MULT)
								& (HASH_TABLE_SIZE - 1);
						if (word_hash[c] != 0)
						{
							collision++;
						}
						else
						{
							word_hash[c] = sizeDict++;
						}
					}
				}
				else
				{
					word_hash[c] = sizeDict++;
				}
			}
		}
		else
		{
			word_hash[j] = sizeDict++;

		}

		if (sizeDict > dictionary || sizeDict >= bound)
		{
			sizeDict--;
			break;
		}

	}

	return mem;
}

int GeneralPacker::WRT::loadCharset(FILE* file, int& freeChar, int* charset,
		int* charsetRev, bool *joinCharsets)
{
	int c, res, mult;

	res = 0;
	c = getc(file);
	mult = 100;
	while (c > 32)
	{
		{
			if (joinCharsets)
				joinCharsets[c] = true;

			charsetRev[freeChar] = c;
			charset[c] = freeChar++;

		}

		res += mult * value[c];
		mult--;

		c = getc(file);
	}

	if (c == 13)
		c = getc(file); // skip CR+LF or LF
	return res;
}

void GeneralPacker::WRT::initializeCodeWords()
{
	int c, charsUsed;

	for (c = 0; c < 256; c++)
	{
		addSymbols[c] = 0;
		codeword2sym[c] = 0;
		sym2codeword[c] = 0;
		reservedSet[c] = 0;
	}

	for (c = BINARY_FIRST; c <= BINARY_LAST; c++)
		addSymbols[c] = 1;

	if (IF_OPTION(OPTION_ADD_SYMBOLS_MISC))
	{
		addSymbols[35] = 1;
		addSymbols[38] = 1;
		addSymbols[60] = 1;
		//	addSymbols[61]=1;
		addSymbols[62] = 1;
		addSymbols[64] = 1;
		addSymbols[94] = 1;
		addSymbols[96] = 1;
		addSymbols[123] = 1;
		addSymbols[124] = 1;
		addSymbols[125] = 1;
		addSymbols[126] = 1;
	}

	if (IF_OPTION(OPTION_ADD_SYMBOLS_0_5))
		for (c = 0; c <= 5; c++)
			addSymbols[c] = 1;

	if (IF_OPTION(OPTION_ADD_SYMBOLS_14_31))
		for (c = 14; c <= 31; c++)
			addSymbols[c] = 1;

	for (c = 0; c < 256; c++)
	{
		if (IF_OPTION(OPTION_USE_DICTIONARY) && addSymbols[c]
				&& (lowerSet[usedSet][c] != 0))
			addSymbols[c] = 0;

		if (IF_OPTION(
				OPTION_USE_DICTIONARY) && addSymbols[c]
				&& (upperSet[usedSet][c]!=0) && !IF_OPTION(OPTION_CAPITAL_CONVERSION))
			addSymbols[c] = 0;

		if ((IF_OPTION(OPTION_USE_DICTIONARY) && addSymbols[c])
				|| c == CHAR_ESCAPE || c == CHAR_LOWERWORD
				|| c == CHAR_FIRSTUPPER || c == CHAR_UPPERWORD
				|| c == CHAR_CR_LF || c == CHAR_NOSPACE
				|| (IF_OPTION(OPTION_WORD_SURROROUNDING_MODELING)
						&& c == CHAR_PUNCTUATION))
			reservedSet[c] = 1;
		else
			reservedSet[c] = 0;
	}

	if (IF_OPTION(
			OPTION_ADD_SYMBOLS_A_Z) && IF_OPTION(OPTION_CAPITAL_CONVERSION))
		for (c = NGRAM_FIRST; c <= NGRAM_LAST; c++)
			addSymbols[c] = 1;

	if (IF_OPTION(OPTION_USE_DICTIONARY))
	{
		charsUsed = 0;
		for (c = 0; c < 256; c++)
		{
			if (addSymbols[c])
			{
				codeword2sym[c] = charsUsed;
				sym2codeword[charsUsed] = c;
				charsUsed++;

			}
		}

		dict1size = 80;
		dict2size = 32;
		dict3size = 16;
		dict4size = 0;

		//dictionary=(dict1size*dict2size*dict3size*dict4size+dict1size*dict2size*dict3size+dict1size*dict2size+dict1size);
		bound4 = dict1size * dict2size * dict3size + dict1size * dict2size
				+ dict1size;
		bound3 = dict1size * dict2size + dict1size;
		dict123size = dict1size * dict2size * dict3size;
		dict12size = dict1size * dict2size;
		dict1plus2 = dict1size + dict2size;
		dict1plus2plus3 = dict1size + dict2size + dict3size;
		dictionary = dict123size + dict1size * (dict2size + dict3size + 1);

		dict = (unsigned char**) calloc(
				sizeof(unsigned char*) * (dictionary + 1), 1);
		dictlen = (unsigned char*) calloc(
				sizeof(unsigned char) * (dictionary + 1), 1);

		PRINT_DICT(
				("usedSet=%d preprocType=%d %d %d %d %d(%d) charsUsed=%d sizeDict=%d\n",usedSet,preprocType,dict1size,dict2size,dict3size,dict4size,dictionary,charsUsed,sizeDict));
	}
}

// read dictionary from files to arrays
bool GeneralPacker::WRT::initialize(unsigned char* dictName,
		unsigned char* shortDictName, bool encoding)
{
	PRINT_DICT(("dictName=%s shortDictName=%s\n",dictName,shortDictName));

	int i, j, c, set[CHARSET_COUNT], fileLen;
	FILE* file, *file2;
	unsigned char* mem;

	deinitialize();
	sizeDict = 0;

	memset(&word_hash[0], 0, HASH_TABLE_SIZE * sizeof(word_hash[0]));
	memset(lowerSet, 0, sizeof(lowerSet));
	memset(upperSet, 0, sizeof(upperSet));
	memset(lowerSetRev, 0, sizeof(lowerSetRev));
	memset(upperSetRev, 0, sizeof(upperSetRev));
	memset(ngram_hash, 0, sizeof(ngram_hash));

	if (dictName == NULL && shortDictName == NULL)
	{
		initializeCodeWords();
		return true;
	}

	if (dictName == NULL && shortDictName != NULL)
	{
		dictName = shortDictName;
		shortDictName = NULL;
	}

	if (IF_OPTION(OPTION_USE_DICTIONARY))
	{
		/////if (WRT_verbose)
		/////	printf("- loading dictionary %s\n",dictName);

		file = fopen((const char*) dictName, "rb");
		if (file == NULL)
		{
			printf("Can't open dictionary %s\n", dictName);
			return false;
		}

		fileLen = flen(file);

		#pragma GCC diagnostic ignored "-Wunused-result"
		fscanf(file, "%d", &dict123size);
		#pragma GCC diagnostic warning "-Wunused-result"

		do
		{
			c = getc(file);
		} while (c >= 32);
		if (c == 13)
			c = getc(file); // skip CR+LF or LF

		for (i = 0; i < CHARSET_COUNT; i++)
		{
			freeUpper[i] = 1;
			set[i] = loadCharset(file, freeUpper[i], upperSet[i],
					upperSetRev[i]);
			freeLower[i] = 1;
			set[i] += loadCharset(file, freeLower[i], lowerSet[i],
					lowerSetRev[i]);
		}

		if (encoding)
		{
			j = 0;
			for (i = 1; i < CHARSET_COUNT - 1; i++)
			{
				if (set[i] > set[j])
					j = i;
			}

			usedSet = j;

			if (set[usedSet] == 0 && freeLower[CHARSET_COUNT - 1] > 1)
				usedSet = CHARSET_COUNT - 1;
		}

		if (freeUpper[0] != freeLower[0])
		{
			/////printf("Bad the first charset in %s (the 2nd and the 3rd line length are different)\n",dictName);
			fclose(file);
			return false;
		}

		if (freeUpper[usedSet] != freeLower[usedSet]
				|| freeLower[usedSet] != freeLower[0])
		{
			/////printf("Bad the second charset in %s (the %dth and the %dth line length are different)\n",dictName,usedSet*2-2,usedSet*2-1);
			fclose(file);
			return false;
		}

		dictmem = (unsigned char*) calloc(fileLen * 2, 1);
		mem = dictmem;
		sizeDict = 1;

		if (!dictmem)
		{
			initializeCodeWords();
			return true;
		}

		if (shortDictName)
		{
			file2 = fopen((const char*) shortDictName, "rb");
			if (file2 == NULL)
			{
				printf("Can't open dictionary %s\n", shortDictName);
				return false;
			}

			do
				c = getc(file2);
			while (c >= 32);
			if (c == 13)
				c = getc(file2);
			for (i = 0; i < CHARSET_COUNT; i++)
			{
				loadCharset(file2, freeUpper[i], upperSet[i], upperSetRev[i]);
				loadCharset(file2, freeLower[i], lowerSet[i], lowerSetRev[i]);
			}

			if (freeUpper[0] != freeLower[0])
			{
				/////printf("Bad the first charset in %s (the 2nd and the 3rd line length are different)\n",shortDictName);
				fclose(file2);
				fclose(file);
				return false;
			}

			if (freeUpper[usedSet] != freeLower[usedSet]
					|| freeLower[usedSet] != freeLower[0])
			{
				/////printf("Bad the second charset in %s (the %dth and the %dth line length are different)\n",dictName,usedSet*2-2,usedSet*2-1);
				fclose(file2);
				fclose(file);
				return false;
			}

			initializeCodeWords();

			if (dict == NULL || dictlen == NULL)
				return false;

			mem = loadDictionary(file2, mem, dictionary);
			fclose(file2);

		}
		else
		{
			initializeCodeWords();

			if (dict == NULL || dictlen == NULL)
				return false;
		}

		mem = loadDictionary(file, mem, dictionary);

		if (encoding && usedSet == CHARSET_COUNT - 1)
		{
			memset(lowerSet, 0, sizeof(lowerSet));
			memset(upperSet, 0, sizeof(upperSet));
			memset(lowerSetRev, 0, sizeof(lowerSetRev));
			memset(upperSetRev, 0, sizeof(upperSetRev));
		}

		fclose(file);

		/////if (WRT_verbose)
		/////	printf("- loaded dictionary %d/%d words\n",sizeDict,dictionary);
	}
	else
	{
		initializeCodeWords();
	}

	return true;
}

void GeneralPacker::WRT::deinitialize()
{
	if (dict)
	{
		free(dict);
		dict = NULL;
	}
	if (dictlen)
	{
		free(dictlen);
		dictlen = NULL;
	}
	if (dictmem)
	{
		free(dictmem);
		dictmem = NULL;
	}

	sizeDict = 0;
}

int GeneralPacker::WRT::detectFileType(FILE* file, int part_length, int parts,
		int& recordLen)
{
	unsigned char s[1024];
	int skip_first, d, i, last_c, c, flen, length;
	int XML_MULT, xml, s_size, shortLangSum, binCount, EOLCount, EOLCountBad,
			punctCount, punctCountBad, fftell, spaceAfterLF;
	float max, current;
	int lastPos[256];
	int fc[1024], fc_max;
	/////int quarterByte=0;
	bool nonlatin;

	dictionary = 1 << 30;
	memset(&lang[0], 0, MAX_DICT_NUMBER * sizeof(lang[0]));
	memset(value, 0, sizeof(value));

	longDictLen = 1;
	longDict = shortDictLen = 0;
	shortDict = lastShortDict = -1;

	return preprocFlag;
}

void GeneralPacker::WRT::set_options(char c, char c2)
{
}

void GeneralPacker::WRT::get_options(int& c, int& c2)
{
	if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
		TURN_OFF(OPTION_CAPITAL_CONVERSION);

	if (!IF_OPTION(OPTION_USE_DICTIONARY))
		TURN_OFF(OPTION_USE_NGRAMS);

	if ((IF_OPTION(OPTION_WORD_SURROROUNDING_MODELING))
			&& (IF_OPTION(OPTION_SPACELESS_WORDS)))
	{
		/////printf("warning: OPTION_WORD_SURROROUNDING_MODELING and OPTION_SPACELESS_WORDS collision\n");
		TURN_OFF(OPTION_WORD_SURROROUNDING_MODELING);
	}

	c = c2 = 0;
	if (IF_OPTION(OPTION_USE_NGRAMS))
		c = c + 128;
	if (IF_OPTION(OPTION_USE_DICTIONARY))
		c = c + 64;
	if (IF_OPTION(OPTION_WORD_SURROROUNDING_MODELING))
		c = c + 16;
	if (IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
		c = c + 8;
	if (IF_OPTION(OPTION_SPACE_AFTER_EOL))
		c = c + 4;

	c += preprocType;

	if (IF_OPTION(OPTION_CAPITAL_CONVERSION))
		c2 = c2 + 32;
	if (IF_OPTION(OPTION_UTF8))
		c2 = c2 + 16;
}

int GeneralPacker::WRT::defaultSettings(int argc, char* argv[])
{
	forceNormalTextFilter = false;
	forceWordSurroroundModeling = false;

	RESET_OPTIONS;

	TURN_ON(OPTION_TRY_SHORTER_WORD);TURN_ON(OPTION_USE_DICTIONARY);TURN_ON(
			OPTION_NORMAL_TEXT_FILTER);TURN_ON(OPTION_CAPITAL_CONVERSION);TURN_ON(
			OPTION_UTF8);

	tryShorterBound = 6;
	TURN_ON(OPTION_WORD_SURROROUNDING_MODELING);

	int optCount = 0;

	return optCount;
}

bool GeneralPacker::WRT::readDicts(char* pattern, char* dictPath,
		int dictPathLen)
{
	cout << "reading dicts ..." << endl;

	FILE* file;
	bool nonlatin;
	int c, i, sizeDict, sizeFullDict = 0;

	memset(joinCharsets, 0, sizeof(joinCharsets));

	map.clear();

	langSum = 0;
	langCount = 0;
	sizeDict = 0;
	dictPath[dictPathLen] = 0;
	strcat(dictPath, pattern);
#ifdef WIN32
	struct _finddata_t c_file;
	long hFile=_findfirst( dictPath, &c_file );

	if (hFile != -1L)
	do
	{
		char* filename=c_file.name;
#else
	{
#pragma GCC diagnostic ignored "-Wwrite-strings"
		char* filename = "temp_HKCC_dict1.dic";
#pragma GCC diagnostic warning  "-Wwrite-strings"
#endif

		dictPath[dictPathLen] = 0;
		strcat(dictPath, filename);

		file = fopen((const char*) dictPath, "rb");
		if (file == NULL)
			return false; //continue;

		i = strlen(filename);

		//toLower((unsigned char*)filename,i);
		langName[langCount] = (unsigned char*) malloc(i + 1);
		memcpy(langName[langCount], (const char*) filename, i + 1);

		memset(lowerSet, 0, sizeof(lowerSet));
		memset(lowerSetRev, 0, sizeof(lowerSetRev));

		do
			c = getc(file);
		while (c >= 32);
		if (c == 13)
			c = getc(file);

		for (i = 0; i < CHARSET_COUNT; i++)
		{
			do
				c = getc(file);
			while (c >= 32);
			if (c == 13)
				c = getc(file);
			freeLower[i] = 1;
			loadCharset(file, freeLower[i], lowerSet[i], lowerSetRev[i],
					joinCharsets);
		}

		sizeFullDict = sizeDict;
		std::string s;

		while (!feof(file))
		{
			s.erase();
			nonlatin = false;
			while (true)
			{
				c = getc(file);
				if (c < 32)
					break;
				s.append(1, (char) c);
				if (lowerSet[0][c] > 0)
					nonlatin = true;
			}
			if (c == EOF)
				break;

			if (c == 13)
				c = getc(file); // skip CR+LF or LF

			if (addWord(s, sizeFullDict))
			{
				sizeDict++;
				if (sizeDict % SAMPLE_WORDS_COUNT == 0)
					break;
			}
			else
				continue;

			if (nonlatin)
			{
				std::string t;

				for (i = 1; i < CHARSET_COUNT - 1; i++)
				{
					if (freeLower[i] > 1)
					{
						t.erase();
#pragma GCC diagnostic ignored "-Wsign-compare"
						for (c = 0; c < s.size(); c++)
#pragma GCC diagnostic warning "-Wsign-compare"
						{
							unsigned char uc = s[c];
							if (lowerSet[0][uc] > 0)
								t.append(1,
										(char) lowerSetRev[i][lowerSet[0][uc]]);
							else
								t.append(1, (char) uc);
						}

						if (addWord(t, sizeFullDict))
						{
							if (sizeDict % SAMPLE_WORDS_COUNT == 0)
								break;
						}
						else
							continue;
					}
				}
			} // end if (nonlatin)

			if (sizeDict % SAMPLE_WORDS_COUNT == 0)
				break;
		}

		if (sizeDict % SAMPLE_WORDS_COUNT_MAX != 0)
			sizeDict = ((sizeDict / SAMPLE_WORDS_COUNT_MAX) + 1)
					* SAMPLE_WORDS_COUNT_MAX;

		langCount++;

		fclose(file);
	}
#ifdef WIN32
	while (_findnext( hFile, &c_file ) == 0);

	_findclose( hFile );

#endif

	return true;
}

void GeneralPacker::WRT::freeNames()
{
	for (int i = 0; i < langCount; i++)
		free(langName[i]);
}

int GeneralPacker::WRT::getSourcePath(char* buf, int buf_size)
{
#ifdef WIN32
	int pos;

	pos=GetModuleFileName(NULL,buf,buf_size);

	if (pos>0)
	{
		for (int i=pos-1; i>=0; i--)
		if (buf[i]=='\\')
		{
			buf[i+1]=0;
			pos=i+1;
			break;
		}
	}
	else
	buf[0]=0;

	return pos;
#else
	buf[0] = 0;
	return 0;
#endif
}

int GeneralPacker::WRT::getFileType(FILE* file, int& recordLen)
{
	int dictPathLen;
	unsigned char dictPath[256];

	if (map.size() == 0)
	{
		getSourcePath((char*) dictPath, sizeof(dictPath));
		strcat((char*) dictPath, WRT_DICT_DIR);
		dictPathLen = strlen((char*) dictPath);

#pragma GCC diagnostic ignored "-Wwrite-strings"
		readDicts(DICTNAME "*" DICTNAME_EXT, (char*) dictPath, dictPathLen);
#pragma GCC diagnostic ignored "-Wwrite-strings"
	}

	return detectFileType(file, 10240, 5, recordLen);
}

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
void GeneralPacker::WRT::encode(FILE* file, FILE* fileout, int fileLen) //,int fileType)
{
#pragma GCC diagnostic ignored "-Wsign-compare"
	unsigned char s[1024];
	EWordType wordType;
	int last_c, c, next_c;
	int binCount = 0;

	preprocessing = 0;
	s_size = 0;
	last_c = 0;
	lastEOL = 0;
	EOLType = UNDEFINED;
	wordType = LOWERWORD;
	spaceBefore = NONE;
	initOrder = true;
	lastEOL = 0;

	if (!IF_OPTION(
			OPTION_NORMAL_TEXT_FILTER) && !IF_OPTION(OPTION_USE_DICTIONARY))
	{
		autoSwitch = 1 << 31 - 1; // MaxSignedInt
		preprocessing = autoSwitch;
	}
	else if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
		autoSwitch = AUTO_SWITCH * 4;
	else
		autoSwitch = AUTO_SWITCH;

	ENCODE_GETC(c, file);
	fftell = 0;

	while (!feof(file))
	{
		if (restartEnc)
			return;

		if (fileCorrupted)
			return;

		PRINT_CHARS(("c=%d (%c)\n",c,c));

		if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && preprocessing > 0)
		{
			if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && COND_BIN_FILTER(c))
			{
				binCount++;
				preprocessing = autoSwitch;
				PRINT_CHARS(("preprocessing=%d c=%c(%d)\n",preprocessing,c,c));
			}
			else
			{
				preprocessing--;
				PRINT_CHARS(("preprocessing=%d c=%c(%d)\n",preprocessing,c,c));
				if (preprocessing == 0)
				{
					initOrder = true;
					if (binCount * 100 / (fftell + 5000) > 25)
					{
						autoSwitch = AUTO_SWITCH * 16;
						preprocessing = autoSwitch;
					}
				}
			}

			ENCODE_PUTC(c, fileout);
			ENCODE_GETC(c, file);
			continue;
		}

		if (reservedSet[c])
		{
			PRINT_CHARS(("reservedSet[c] c=%d (%c)\n",c,c));

			encodeWord(fileout, s, s_size, wordType);
			s_size = 0;
			ENCODE_PUTC(CHAR_ESCAPE, fileout);
			ENCODE_PUTC(c, fileout);

			if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && COND_BIN_FILTER(c))
				preprocessing = autoSwitch;

			last_c = c;
			ENCODE_GETC(c, file);
			continue;
		}

		if (c >= 'a' && c <= 'z' || lowerSet[usedSet][c] > 0)
		{
			if (usedSet > 0 && lowerSet[usedSet][c] > 0)
				c = lowerSetRev[0][lowerSet[usedSet][c]];

			PRINT_CHARS(("a-z c=%d (%c)\n",c,c));

			if (s_size == 0)
			{
				wordType = LOWERWORD;
			}
			else
			{
				if (wordType == UPPERWORD)
				{
					encodeWord(fileout, s, s_size, wordType);
					if (IF_OPTION(OPTION_CAPITAL_CONVERSION))
						ENCODE_PUTC(CHAR_LOWERWORD, fileout);

					wordType = LOWERWORD;
					s_size = 1;
					s[0] = c;
					last_c = c;
					ENCODE_GETC(c, file);
					continue;
				}
			}

			s[s_size++] = c;
			if (s_size >= sizeof(s) - 1)
			{
				encodeWord(fileout, s, s_size, wordType);
				s_size = 0;
			}
			last_c = c;
			ENCODE_GETC(c, file);
			continue;
		}

		if (c <= 'Z' && c >= 'A' || upperSet[usedSet][c] > 0)
		{
			if (usedSet > 0 && upperSet[usedSet][c] > 0)
				c = upperSetRev[0][upperSet[usedSet][c]];

			PRINT_CHARS(("A-Z c=%d (%c)\n",c,c));

			if (s_size == 0)
			{
				wordType = FIRSTUPPER;
			}
			else
			{
				if (wordType == FIRSTUPPER)
				{
					if (s_size == 1)
						wordType = UPPERWORD;
					else
					{
						encodeWord(fileout, s, s_size, wordType);

						wordType = FIRSTUPPER;
						s_size = 1;
						s[0] = c;
						last_c = c;
						ENCODE_GETC(c, file);
						continue;
					}
				}
				else if (wordType == LOWERWORD)
				{
					encodeWord(fileout, s, s_size, wordType);

					wordType = FIRSTUPPER;
					s_size = 1;
					s[0] = c;
					last_c = c;
					ENCODE_GETC(c, file);
					continue;
				}
			}

			s[s_size++] = c;
			if (s_size >= sizeof(s) - 1)
			{
				encodeWord(fileout, s, s_size, wordType);
				s_size = 0;
			}

			last_c = c;
			ENCODE_GETC(c, file);
			continue;
		}

		encodeWord(fileout, s, s_size, wordType);

		s_size = 0;

		PRINT_CHARS(("other c=%d\n",c));

		ENCODE_GETC(next_c, file);

		{
			if (c != EOF)
			{
				if ((IF_OPTION(OPTION_WORD_SURROROUNDING_MODELING)) && c > ' ')
				{
					if ((last_c >= 'a' && last_c <= 'z')
							|| (last_c >= 'A' && last_c <= 'Z'))
					{
						if (fftell < fileLen)
						{
							ENCODE_PUTC(' ', fileout);
							ENCODE_PUTC(CHAR_PUNCTUATION, fileout);
							ENCODE_PUTC(c, fileout);
						}
						else
							ENCODE_PUTC(c, fileout);
					}
					else if (next_c >= 'a' && next_c <= 'z')
					{
						ENCODE_PUTC(c, fileout);
						ENCODE_PUTC(CHAR_LOWERWORD, fileout);
						ENCODE_PUTC(' ', fileout);
					}
					else
						ENCODE_PUTC(c, fileout);
				}
				else if ((IF_OPTION(OPTION_SPACELESS_WORDS)) && c == ' ')
					spaceBefore = SPACE;
				else
					ENCODE_PUTC(c, fileout);
			}
		}

		last_c = c;
		c = next_c;
	}

	encodeWord(fileout, s, s_size, wordType);
	s_size = 0;
#pragma GCC diagnostic warning "-Wsign-compare"
}

int GeneralPacker::WRT::readEOLstream(FILE* file)
{
	unsigned int EOLstream_len = 0;
	unsigned int fileLen;

	fseek(file, -4, SEEK_END);
	fileLen = ftell(file) + 4;

	for (int i = 0; i < 4; i++)
		EOLstream_len = EOLstream_len * 256 + fgetc(file);

	return EOLstream_len;
}

void GeneralPacker::WRT::writeEOLstream(FILE* fileout)
{
	unsigned int EOLstream_len;

	EOLstream_len = 0;

	fprintf(fileout, "%c%c%c%c", EOLstream_len >> 24, EOLstream_len >> 16,
			EOLstream_len >> 8, EOLstream_len);
}

void GeneralPacker::WRT::start_encoding(FILE* file, FILE* fileout,
		unsigned int fileLen, bool type_detected)
{
	int i, c, c2, recordLen = 0, dictPathLen;
	unsigned char s[256];
	unsigned char t[256];
	unsigned char dictPath[256];
	s[0] = 0;
	t[0] = 0;

	llbckp = 0;
	swapCase = false;
	usedSet = 0;

	getSourcePath((char*) dictPath, sizeof(dictPath));
	strcat((char*) dictPath, WRT_DICT_DIR);
	dictPathLen = strlen((char*) dictPath);

	if (!type_detected)
		getFileType(file, recordLen);

	if (IF_OPTION(OPTION_USE_DICTIONARY) && longDictLen > 0)
		memcpy(s, langName[longDict],
				strlen((const char*) langName[longDict]) + 1);
	else
		longDictLen = 0;

	if (IF_OPTION(OPTION_USE_DICTIONARY) && shortDictLen > 0)
		memcpy(t, langName[shortDict],
				strlen((const char*) langName[shortDict]) + 1);
	else
		shortDictLen = 0;

	if (dictPathLen > 0)
	{
		dictPath[dictPathLen] = 0;
		strcat((char*) dictPath, (char*) s);
		strcpy((char*) s, (char*) dictPath);

		dictPath[dictPathLen] = 0;
		strcat((char*) dictPath, (char*) t);
		strcpy((char*) t, (char*) dictPath);
	}

	restart:

	int pos = ftell(fileout);
	fprintf(fileout, "%c%c%c%c", 0, 0, 0, 0);

	get_options(c, c2); // before initialize
	putc(c, fileout);
	putc(c2, fileout);

	/////WRT_print_options();

	deinitialize();

	if (!initialize((longDictLen <= 0) ? NULL : s,
			(shortDictLen <= 0) ? NULL : t, true))
		return;

	if (IF_OPTION(OPTION_USE_DICTIONARY))
	{
		putc(usedSet, fileout);

		putc(shortDictLen, fileout);
		c = strlen(SHORT_DICTNAME);
		for (i = c; i < shortDictLen + c; i++)
			fputc(langName[shortDict][i], fileout);

		putc(longDictLen, fileout);
		c = strlen(DICTNAME);
		for (i = c; i < longDictLen + c; i++)
			fputc(langName[longDict][i], fileout);
	}

	{
		encode(file, fileout, fileLen);
		if (restartEnc)
		{
			restartEnc = false;
			fseek(fileout, pos, SEEK_SET);
			fseek(file, 0, SEEK_SET);
			llbckp = 0;
			swapCase = false;
			goto restart;
		}

		unsigned int fileLen = ftell(fileout) + 4;
		fseek(fileout, pos, SEEK_SET);
		fprintf(fileout, "%c%c%c%c", fileLen >> 24, fileLen >> 16, fileLen >> 8,
				fileLen);
		fseek(fileout, fileLen - 4, SEEK_SET);
	}
}

void GeneralPacker::WRT::start_decoding(FILE* file, FILE* fileout, int header)
{
	int i, c, c2, recordLen = 0, dictPathLen;
	unsigned char s[256];
	unsigned char t[256];
	unsigned char dictPath[256];
	unsigned int fileLen;
	s[0] = 0;
	t[0] = 0;

	llbckp = 0;
	swapCase = false;
	usedSet = 0;
	WRTd_binCount = 0;

	for (i = 0, fileLen = 0; i < 4; i++)
		fileLen = fileLen * 256 + fgetc(file);

	i = 0;
	c = getc(file);
	c2 = getc(file);

	preprocType = (EPreprocessType) (c % 4); // { LZ77, BWT, PPM, PAQ };

	defaultSettings(0, NULL); // after setting preprocType

	set_options(c, c2);

	if (IF_OPTION(OPTION_USE_DICTIONARY))
	{
		usedSet = getc(file);
		shortDictLen = getc(file);
		c = strlen(SHORT_DICTNAME);
		memcpy(t, SHORT_DICTNAME, c);
		for (i = c; i < shortDictLen + c; i++)
			t[i] = getc(file);
		c = strlen(DICTNAME_EXT);
		memcpy(t + i, DICTNAME_EXT, c);
		i += c;
		t[i] = 0;

		longDictLen = getc(file);
		c = strlen(DICTNAME);
		memcpy(s, DICTNAME, c);
		for (i = c; i < longDictLen + c; i++)
			s[i] = getc(file);
		c = strlen(DICTNAME_EXT);
		memcpy(s + i, DICTNAME_EXT, c);
		i += c;
		s[i] = 0;

		i = longDictLen + 1 + shortDictLen + 1
				+ (IF_OPTION(OPTION_USE_DICTIONARY) ? 1 : 0); // usedSet
	}
	else
	{
		longDictLen = 0;
		shortDictLen = 0;
	}
	header += 4;
	i += 2 + header; // WRT4

	getSourcePath((char*) dictPath, sizeof(dictPath));
	strcat((char*) dictPath, WRT_DICT_DIR);
	dictPathLen = strlen((char*) dictPath);

	if (dictPathLen > 0)
	{
		dictPath[dictPathLen] = 0;
		strcat((char*) dictPath, (char*) s);
		strcpy((char*) s, (char*) dictPath);

		dictPath[dictPathLen] = 0;
		strcat((char*) dictPath, (char*) t);
		strcpy((char*) t, (char*) dictPath);
	}

	deinitialize();

	if (!initialize((longDictLen <= 0) ? NULL : s,
			(shortDictLen <= 0) ? NULL : t, false))
		return;

	{
		int EOLlen = 4;

		fseek(file, i, SEEK_SET); // skip "WRTx" header
		EOLlen += i;	// header + fileLen

		WRTd_filter->reads += i;

		originalFileLen = fileLen - EOLlen;
		bufferedChar = -1;
		lastChar = 0;
		fftell = 0;
		fftelld = 0;
		WRTd_upper = false;
		upperWord = UFALSE;
		preprocessing = 0;
		s_size = 0;
		initOrder = true;
		lastEOL = -1;
		EOLType = UNDEFINED;

		if (!IF_OPTION(
				OPTION_NORMAL_TEXT_FILTER) && !IF_OPTION(OPTION_USE_DICTIONARY))
		{
			autoSwitch = 1 << 31 - 1; // MaxSignedInt
			preprocessing = autoSwitch;
		}
		else if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
			autoSwitch = AUTO_SWITCH * 4;
		else
			autoSwitch = AUTO_SWITCH;

		if (IF_OPTION(OPTION_SPACELESS_WORDS))
			spaceBefore = SPACE;
		else
			spaceBefore = NONE;

		DECODE_GETC(WRTd_c, file);PRINT_CHARS(
				("WRT_start_decoding WRTd_c=%d ftell=%d\n",WRTd_c,ftell(file)));
	}
}

void GeneralPacker::WRT::prepare_decoding()
{
	WRTd_type = 0;
}

int GeneralPacker::WRT::decode_char(FILE* file, FILE* fileout, int header)
{
	switch (WRTd_type)
	{
	default:
	case 0:
		start_decoding(file, fileout, header);
		WRTd_qstart = WRTd_qend = 0;
		WRTd_type = 1;
		/////if (IF_OPTION(OPTION_DNA_QUARTER_BYTE) || IF_OPTION(OPTION_RECORD_INTERLEAVING))
		/////	return EOF;
	case 1:
		if (WRTd_c != EOF)
		{
			while (WRTd_qstart >= WRTd_qend && WRTd_c != EOF)
			{
				WRTd_qstart = WRTd_qend = 0;
				decode(file);
				if (fileCorrupted)
					WRTd_type = 2;
			}

			if (WRTd_qstart < WRTd_qend)
				return WRTd_queue[WRTd_qstart++];
		}
		hook_putc(EOF);
		WRTd_type = 2;
	case 2:
		if (WRTd_qstart < WRTd_qend)
			return WRTd_queue[WRTd_qstart++];
		else
			return -1;
	}
}

// **************************************************************************************************************
// **************************************************************************************************************
// AIR: include end
// **************************************************************************************************************
// **************************************************************************************************************

GeneralPacker::WRT wrt;

GeneralPacker::TextFilter::TextFilter(Encoder* e) :
		Filter(e), dtmp(NULL)
{
	reset();
	WRTd_filter = this;
}
GeneralPacker::TextFilter::~TextFilter()
{
	reset();
	tmp = NULL;
}

void GeneralPacker::TextFilter::encode(FILE* f, int n)
{
//	wrt.defaultSettings(0,NULL);
	wrt.start_encoding(f, tmp, n, true); // wrt.WRT_getFileType() called in make()
}

int GeneralPacker::TextFilter::decode()
{
	if (first)
	{
		first = false;
		if (!tmp)
		{
			if (dtmp)
				fclose(dtmp);

			dtmp = tmpfile();
			if (!dtmp)
				perror("WRT tmpfile"), exit(1);

			unsigned int size = 0, i;
			for (i = 4; i != 0; --i)
			{
				int c = read();
				size = size * 256 + c;
				putc(c, dtmp);
			}

			size -= 4;
			for (; i < size; i++)
			{
				putc(read(), dtmp);
				printStatus(i);
			}

			fseek(dtmp, 0, SEEK_SET);
			tmp = dtmp;
		}
	}

	return wrt.decode_char(tmp, NULL, 0);
}

void GeneralPacker::TextFilter::reset()
{
	first = true;
	wrt.prepare_decoding();
	if (dtmp)
		fclose(dtmp);
}
;

////////////////// Filter::make ////////////

// Create a new Filter of an appropriate type, either by examining
// filename (and maybe its contents) and setting filetype (for compression)
// or if filename is NULL (decompression) then use the supplied value
// of filetype.

GeneralPacker::Filter*
GeneralPacker::Filter::make(const char* filename, Encoder* e)
{
	if (filename)
	{
		fileType = 0;
		const char* ext = strrchr(filename, '.');

		FILE* file = fopen(filename, "rb");
		if (file)
		{
			if (fgetc(file) == 'M' && fgetc(file) == 'Z')
				fileType = EXE;
			else
			///if (!ext || (!equals(ext, ".dbf") && !equals(ext, ".mdb") && !equals(ext, ".tar")
			///	&& !equals(ext, ".c") && !equals(ext, ".cpp") && !equals(ext, ".h")
			///	&& !equals(ext, ".hpp") && !equals(ext, ".ps") && !equals(ext, ".hlp")
			///	&& !equals(ext, ".ini") && !equals(ext, ".inf") ))
			{
				// fseek(file, 0, SEEK_SET ); <- unnecessary for WRT
				wrt.defaultSettings(0, NULL);
				int recordLen; // unused in PAQ
				if (wrt.getFileType(file, recordLen) > 0) // 0 = binary or not known
				{
					if (IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
						fileType = TEXT;
					else
						fileType = BINTEXT;
				}
			}
			fclose(file);
		}
	}
	if (e)
	{
		///if (filetype==EXE)
		///return new ExeFilter(e);
		///else
		if (fileType == TEXT || fileType == BINTEXT)
			return new TextFilter(e);
		else
			return new DefaultFilter(e);
	}
	return NULL;
}

//////////////////////////// main program ////////////////////////////

// Read one line, return NULL at EOF or ^Z.  f may be opened ascii or binary.
// Trailing \r\n is dropped.  Lines longer than MAXLINE-1=511 are truncated.

char *
getline(FILE *f = stdin)
{
	const int MAXLINE = 512;
	static char s[MAXLINE];
	int len = 0, c;
	while ((c = getc(f)) != EOF && c != 26 && c != '\n' && len < MAXLINE - 1)
	{
		if (c != '\r')
			s[len++] = c;
	}
	s[len] = 0;
	if (c == EOF || c == 26)
		return 0;
	return s;
}

// Test if files exist and get their sizes, store in archive header
void store_in_header(FILE* f, char* filename, long& total_size)
{
	FILE *fi = fopen(filename, "rb");
	if (fi)
	{
		fseek(fi, 0, SEEK_END);  // get size
		long size = ftell(fi);
		total_size += size;
		if ((size & ~0x7fffffffL) || (total_size & ~0x7fffffffL))
		{
			fprintf(stderr, "File sizes must total less than 2 gigabytes\n");
			fprintf(f, "-1\tError: over 2 GB\r\n");
			exit(1);
		}
		fclose(fi);
		if (size != -1)
		{
			fprintf(f, "%ld\t%s\r\n", size, filename);
			return;
		}
	}
	perror(filename);
}

// Compress/decompress files.  Usage: yz archive files...
// If archive does not exist, it is created and the named files are
// compressed.  If there are no file name arguments after the archive,
// then file names are read from input up to a blank line or EOF.
// If archive already exists, then the files in the archive are either
// extracted, or compared if the file already exists.  The files
// are extracted to the names listed on the command line in the
// order they were stored, defaulting to the stored names.

// TODO add DIC files as gcc dependencies
int GeneralPacker::pg_main(int argc, char** argv)
{
	clock_t start_time = clock();  // start timer
	long total_size = 0;  // uncompressed size of all files
	FILE *f;

	if (argc < 2)
	{
		//TODO throw ex
		return 0;
	}

	// Get option
	int option = '4';
	if (argv[1][0] == '-')
		option = argv[1][1], ++argv, --argc;
#ifndef SHORTEN_CODE
	if (option < 32)
		option = 32;
	if (option < '0' || option > '9')
		fprintf(stderr, "Bad option -%c (use -0 to -9)\n", option), exit(1);
#endif

	// Test for archive.  If none, create one and write a header.
	// The first line is PROGNAME.  This is followed by a list of
	// file sizes (as decimal numbers) and names, separated by a tab
	// and ending with \r\n.  The last entry is followed by ^Z
	Mode mode = DECOMPRESS;

	f = fopen(argv[1], "rb");
	if (!f)
	{
		mode = COMPRESS;

		f = fopen(argv[1], "wb");
		if (!f)
			perror(argv[1]), exit(1);
		fprintf(f, "%s -%c\r\n", PROGNAME, option);

		// Get filenames
#ifndef SHORTEN_CODE
		if (argc == 2)
			printf(
					"Enter names of files to compress, followed by blank line\n");
#endif
		int i = 2;
		std::multimap<int, std::string> filetypes;
		std::multimap<int, std::string>::iterator it;
		std::vector<std::string> filenames;
		char *filename;

		while (true)
		{

#ifndef SHORTEN_CODE
			if (argc == 2)
			{
				filename = getline();
				if (!filename || !filename[0])
					break;
			}
			else
#endif
			{
				if (i == argc)
					break;
				filename = argv[i++];
			}
			filenames.push_back(filename);
		} // end while

#pragma GCC diagnostic ignored "-Wsign-compare"
		for (i = 0; i < filenames.size(); i++)
#pragma GCC diagnostic warning "-Wsign-compare"
		{
			Filter::make(filenames[i].c_str(), NULL);
			std::pair<int, std::string> p(
					(1 << 31 - 1)
							- (fileType * 8 * 16 * 2048 + preprocFlag
									+ (wrt.longDict + 1) * 16 * 2048
									+ (wrt.shortDict + 1) * 2048),
					filenames[i]);
			PRINT_DICT((" %s\n",filenames[i].c_str()));
			filetypes.insert(p);
		}

		for (it = filetypes.begin(); it != filetypes.end(); it++) // multimap is sorted by filetype
			store_in_header(f, (char*) it->second.c_str(), total_size);

		fputc(26, f);  // EOF
		fclose(f);
		f = fopen(argv[1], "r+b");
		if (!f)
			perror(argv[1]), exit(1);
	}

	// Read existing archive. Two pointers (header and body) track the
	// current filename and current position in the compressed data.
	long header, body;  // file positions in header, body
	char *filename = getline(f);  // check header
	if (!filename || strncmp(filename, PROGNAME " -", strlen(PROGNAME) + 2))
		fprintf(stderr, "%s: not a " PROGNAME " file\n", argv[1]), exit(1);
	option = filename[strlen(filename) - 1];
	compressionLevel = option - '0';
	if (compressionLevel < 0 || compressionLevel > 9)
		compressionLevel = DEFAULT_OPTION;

	{
		buf.setsize(MEM * 8);
		FILE *dictfile = fopen("./to_train_models.dic", "rb"), *tmpfi = fopen(
				"./tmp_tmp_tmp_tmp.dic", "wb");
		fileType = 0;
		// AIR
//		Encoder en(COMPRESS, tmpfi);
//		en.compress(0);
//		for (int i = 0; i < 465211; ++i)
//			en.compress(getc(dictfile));
//		en.flush();
//		fclose(tmpfi);
	}
	header = ftell(f);

	// Initialize encoder at end of header
	if (mode == COMPRESS)
		fseek(f, 0, SEEK_END);
	else
	{  // body starts after ^Z in file
		int c;
		while ((c = getc(f)) != EOF && c != 26)
			;
		if (c != 26)
			fprintf(stderr, "Archive %s is incomplete\n", argv[1]), exit(1);
	}
	Encoder en(mode, f);
	body = ftell(f);

	// Compress/decompress files listed on command line, or header if absent.
	int filenum = 1;  // command line index
	total_size = 0;
	while (1)
	{
		fseek(f, header, SEEK_SET);
		if ((filename = getline(f)) == 0)
			break;
		size = atol(filename);  // parse size and filename, separated by tab
		total_size += size;
		while (*filename && *filename != '\t')
			++filename;
		if (*filename == '\t')
			++filename;
		printf("%ld\t%s: ", size, filename);
		fileSize = (int) size;
		/*   if (mode==DECOMPRESS && ++filenum<argc  // doesn't work with sorting depend on type of file
		 && strcmp(argv[filenum], filename)) {
		 printf(" -> %s", argv[filenum]);
		 filename=argv[filenum];
		 } */
		if (size < 0 || total_size < 0)
			break;
		header = ftell(f);
		fseek(f, body, SEEK_SET);

		// If file exists in COMPRESS mode, compare, else compress/decompress
		FILE *fi = fopen(filename, "rb");
		fileType = 0;
		if (mode == COMPRESS)
		{
			if (!fi)
				perror(filename), exit(1);
			Filter* fp = Filter::make(filename, &en);   // sets filetype
			fp->compress(fi, size);
			printf(" -> %4ld   \n", ftell(f) - body);
			delete fp;
		}
		else
		{  // DECOMPRESS, first byte determines filter type
			fileType = en.decompress();
			Filter* fp = Filter::make(0, &en);
			if (fi)
				fp->compare(fi, size);
			else
			{  // extract
				fi = fopen(filename, "wb");
				if (fi)
					fp->decompress(fi, size);
				else
				{
					perror(filename);
					fp->skip(size);
				}
			}
			delete fp;
		}
		if (fi)
			fclose(fi);
		body = ftell(f);
	}
	fseek(f, body, SEEK_SET);
	en.flush();

	// Print stats
	if (f)
	{
		if (mode == DECOMPRESS && filenum < argc - 1)
			printf("No more files to extract\n");
		long compressed_size = ftell(f);
		if (mode == COMPRESS)
			printf("%ld -> %ld", total_size, compressed_size);
		else
			printf("%ld -> %ld", compressed_size, total_size);
		start_time = clock() - start_time;
		if (compressed_size > 0 && total_size > 0 && start_time > 0)
		{
			printf(" (%1.4f bpc) in %1.2f sec (%1.3f KB/sec), %d Kb\n",
					8.0 * compressed_size / total_size,
					(double) start_time / CLOCKS_PER_SEC,
					0.001 * total_size * CLOCKS_PER_SEC / start_time,
					programObserver.maxmem / 1024);
		}
	}
	return 0;
}

//}
