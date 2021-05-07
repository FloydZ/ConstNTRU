#include <array>
#include <cstdlib>
#include <iostream>

#include "const_rand.h"

#define MODQ(X) ((X) & (Q-1))

// i think this is only valid for >= C++11
constexpr int32_t const_mod(const int32_t x) { return x%3; }

// Map {-1, 0, 1} -> {0,1,q-1} in place
template<const uint64_t N, const uint64_t Q>
constexpr void const_poly_Z3_to_Zq(int32_t f[N]) {
	for(uint64_t i = 0; i < N; i++)
		f[i] = (f[i]+1) | ((-((f[i]+1)>>1)) & (Q-1));
}

// Map {0, 1, q-1} -> {-1,0,1} in place
template<const uint64_t N, const uint64_t LOGQ>
constexpr void poly_trinary_Zq_to_Z3(int32_t r[N]) {
	for(uint64_t i = 0; i < N; i++)
		r[i] = (3 & (r[i] ^ (r[i]>>(LOGQ-1)))) -1;
}


template<const uint64_t N, const uint64_t Q>
constexpr void poly_Rq_mul(int32_t r[N], const int32_t a[N], const int32_t b[N]) {
	#define MODQ(X) ((X) & (Q-1))

	for(uint64_t k = 0; k < N; k++) {
		r[k] = 0;
		for(uint64_t i = 1; i < N-k; i++)
			r[k] += a[k+i] * b[N-i];
		for(uint64_t i = 0; i < k+1; i++)
			r[k] += a[k-i] * b[i];
		r[k] = MODQ(r[k]);
	}
}

template<const uint64_t N, const uint64_t Q>
constexpr void poly_Sq_mul(int32_t r[N], const int32_t a[N], const int32_t b[N]) {
	#define MODQ(X) ((X) & (Q-1))

	poly_Rq_mul<N, Q>(r, a, b);
	for(uint64_t i = 0; i < N; i++)
		r[i] = MODQ(r[i] - r[N-1]);
}

#define POLY_R2_ADD(I,A,B,S)        \
   for(I=0; I<N; I++)               \
   { A[I] ^= B[I] * S;  }

#define POLY_S3_FMADD(I,A,B,S)                   \
   for(I=0; I<N; I++)                            \
   { A[I] = const_mod3(A[I] + S * B[I]); }

template<const uint64_t N>
constexpr void cswappoly(int32_t a[N], int32_t b[N], int swap){
	swap = -swap;
	for(uint64_t i = 0; i < N;i++) {
		uint16_t t = (a[i] ^ b[i]) & swap;
		a[i] ^= t;
		b[i] ^= t;
	}
}

template<const uint64_t N>
constexpr inline void poly_divx(int32_t a[N], const int s){
	for(uint64_t i = 1; i < N; i++)
		a[i-1] = (s * a[i]) | (!s * a[i-1]);
	a[N-1] = (!s * a[N-1]);
}

template<const uint64_t N>
constexpr inline void poly_mulx(int32_t a[N], int s) {
	for(uint64_t i=1; i<N; i++)
		a[N-i] = (s * a[N-i-1]) | (!s * a[N-i]);
	a[0] = (!s * a[0]);
}

// b = 1 means mov, b = 0 means don't mov
constexpr void cmov(unsigned char *r, const unsigned char *x, const size_t len, unsigned char b) {
	b = -b;
	for(size_t i = 0; i < len; i++)
		r[i] ^= b & (x[i] ^ r[i]);
}

template<const uint64_t N, const uint64_t Q>
constexpr void poly_R2_inv(int32_t r[N], const int32_t a[N]) {
	// Schroeppel--Orman--O'Malley--Spatscheck
	// "Almost Inverse" algorithm as described
	// by Silverman in NTRU Tech Report #14
	// with several modifications to make it run in constant-time
	int i= 0, j = 0;
	int k = 0;
	uint16_t degf = N-1;
	uint16_t degg = N-1;
	int sign = 0, t = 0, swap = 0;
	int done = 0;
	int32_t b[N], f[N], g[N];

	int32_t *c = r; // save some stack space
	int32_t *temp_r = f;

	// b(X) := 1
	for(i = 1; i < N; i++)
		b[i] = 0;
	b[0] = 1;

	// c(X) := 0
	for(i=0; i<N; i++)
		c[i] = 0;

	// f(X) := a(X)
	for(i=0; i<N; i++)
		f[i] = a[i] & 1;

	// g(X) := 1 + X + X^2 + ... + X^{N-1}
	for(i=0; i<N; i++)
		g[i] = 1;

	for(j = 0; j < 2 * (N-1) - 1; j++) {
		sign = f[0];
		swap = sign & !done & ((degf - degg) >> 15);

		cswappoly<N>(f, g, swap);
		cswappoly<N>(b, c, swap);
		t = (degf ^ degg) & (-swap);
		degf ^= t;
		degg ^= t;

		POLY_R2_ADD(i, f, g, sign*(!done));
		POLY_R2_ADD(i, b, c, sign*(!done));

		poly_divx<N>(f, !done);
		poly_mulx<N>(c, !done);
		degf -= !done;
		k += !done;

		done = 1 - (((uint16_t)-degf) >> 15);
	}

	k = k - N*((uint16_t)(N - k - 1) >> 15);


	// Return X^{N-k} * b(X)
	// This is a k-coefficient rotation. We do this by looking at the binary
	// representation of k, rotating for every power of 2, and performing a cmov
	// if the respective bit is set.
	for (i = 0; i < N; i++)
		r[i] = b[i];

	for (i = 0; i < 10; i++) {
		for (j = 0; j < N; j++) {
			temp_r[j] = r[(j + (1 << i)) % N];
		}
		cmov((unsigned char *)r, (unsigned char *)temp_r, sizeof(uint16_t) * N, k & 1);
		k >>= 1;
	}
}

// TODO  not working vor small q, or?
template<const uint64_t N, const uint64_t Q>
constexpr void poly_R2_inv_to_Rq_inv(int32_t r[N],
									 const int32_t ai[N],
									 const int32_t a[N]) {
	int32_t b[N], c[N], s[N];

	// for 0..4
	//    ai = ai * (2 - a*ai)  mod q
	for(uint64_t i = 0; i < N; i++)
		b[i] = MODQ(Q - a[i]); // b = -a

	for(uint64_t i = 0; i < N; i++)
		r[i] = ai[i];

	poly_Rq_mul<N, Q>(c, r, b);
	c[0] += 2; // c = 2 - a*ai
	poly_Rq_mul<N, Q>(s, c, r); // s = ai*c

	poly_Rq_mul<N, Q>(c, s, b);
	c[0] += 2; // c = 2 - a*s
	poly_Rq_mul<N, Q>(r, c, s); // r = s*c

	poly_Rq_mul<N, Q>(c, r, b);
	c[0] += 2; // c = 2 - a*r
	poly_Rq_mul<N, Q>(s, c, r); // s = r*c

	poly_Rq_mul<N, Q>(c, s, b);
	c[0] += 2; // c = 2 - a*s
	poly_Rq_mul<N, Q>(r, c, s); // r = s*c
}


template<const uint64_t N, const uint64_t Q>
constexpr void poly_Rq_inv(int32_t r[N], const int32_t a[N]) {
	int32_t ai2[N] = {0};

	poly_R2_inv<N, Q>(ai2, a);
	poly_R2_inv_to_Rq_inv<N, Q>(r, ai2, a);
}

template<const uint64_t N, const uint64_t LOGQ, const uint64_t Q>
constexpr void const_ntru_keygen(int32_t f[N],
								 int32_t g[N],
								 int32_t h[N],
                                 const int32_t seed[2*N]) {
	int32_t tmpp[N];

	// local variables
	int32_t Gf[N], G[N], tmp[N], invGf[N], invh[N];

	// randomly sample f
	// this doesnt make sure that f has correct weight.
	for(uint64_t i = 0; i < N - 1; i++)
		f[i] = const_mod(seed[i]);
	f[N-1] = 0;

	// randomly chooses g
	// this doesnt make sure that g has correct weight
	for(uint64_t i = 0; i < N - 1; i++)
		g[i] = const_mod(seed[N+i]);
	g[N-1] = 0;

	const_poly_Z3_to_Zq<N, Q>(f);
	const_poly_Z3_to_Zq<N, Q>(g);

	/* G = 3*g */
	for(uint64_t i = 0; i < N; i++)
		G[i] = MODQ(3 * g[i]);

	poly_Rq_mul<N, Q>(Gf, G, f);

	poly_Rq_inv<N, Q>(invGf, Gf);

	poly_Rq_mul<N, Q>(tmp, invGf, f);
	poly_Sq_mul<N, Q>(invh, tmp, f);

	poly_Rq_mul<N, Q>(tmp, invGf, G);
	poly_Rq_mul<N, Q>(h, tmp, G);


	poly_Rq_mul<N, Q>(tmpp, h, f);
	for (int i = 0; i < N; ++i) {
		std::cout << tmpp[i] << " ";
	}

	std::cout << "\n";
	for (int i = 0; i < N; ++i) {
		std::cout << G[i] << " ";
	}

	std::cout << "\n";


	poly_trinary_Zq_to_Z3<N, LOGQ>(f);
	poly_trinary_Zq_to_Z3<N, LOGQ>(g);
}

template<const uint64_t N>
bool find_repetition(int32_t h[N]) {
	int32_t v = h [0];
	for (int i = 1; i < N; ++i) {
		if (h[i] == v)
			return true;
	}

	return false;
}

#ifdef FPLLL_H
/// seeded version
/// \tparam ZT
/// \tparam N
/// \param m
/// \param h
template<typename ZT, const uint64_t N>
void generate_const_ntru_encrypt_matrix(fplll::ZZ_mat<ZT> &m, const int32_t h[N]){
	m.resize(N, N);

	for (int j = 0; j < N; ++j) {       // just copy the first row
		m(0, j) = h[j];
	}

	for (int i = 1; i < N; ++i) {       // rows
		auto b = m(i-1, N-1);
		m(i, 0) = b;

		for (int j = 1; j < N; ++j) {   // columns
			auto b = m(i-1, (j-1));
			m(i, j) = b;
		}
	}
}

/// unsseeded version
template<typename ZT, const uint64_t N, const uint64_t LOGQ, const uint64_t Q>
void generate_const_ntru_encrypt_matrix(fplll::ZZ_mat<ZT> &A){
	static std::random_device dev;
	static std::mt19937 rng(dev());
	static std::uniform_int_distribution<int32_t> dist(-(Q/2), Q/2);

	int32_t seed[2*N] = {0}, f[N] = {0}, g[N] = {0}, h[N] = {0};

	do { // dont tell andre how i fixed the bug :)
		for (int &i : seed) { i = dist(rng); }

		const_ntru_keygen<N, LOGQ, Q>(f, g, h, seed);
	} while(find_repetition<N>(h));
	generate_const_ntru_encrypt_matrix<Label_Type, N>(A, h);

}

// seeded=constant version
template<typename ZT, const uint64_t N, const uint64_t LOGQ, const uint64_t Q>
void generate_const_ntru_encrypt_matrix_const(fplll::ZZ_mat<ZT> &A){
	static std::random_device dev;
	static std::mt19937 rng(dev());
	static std::uniform_int_distribution<int32_t> dist(-(Q/2), Q/2);

	auto seed = uniform_distribution<int32_t, 2*N>(-(Q/2), Q/2);
	int32_t f[N] = {0}, g[N] = {0}, h[N] = {0};

	const_ntru_keygen<N, LOGQ, Q>(f, g, h, seed.data());
	generate_const_ntru_encrypt_matrix<Label_Type, N>(A, h);
}
#endif