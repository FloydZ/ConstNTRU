#include "const_ntru.h"
#include "const_rand.h"

// NTRU-HPS, paramters
#ifdef NTRUHPS2048509
//  -ntruhps2048509
#define N       509u
#define LOGQ    11u
#define Q       int64_t(uint64_t (1) << LOGQ)
#endif

#ifdef NTRUHPS2048677
//  -ntruhps2048677
#define N       677u
#define LOGQ    11u
#define Q       int64_t(uint64_t (1) << LOGQ)
#endif

#ifdef NTRUHPS4096821
//  -ntruhps4096821
#define N       821u
#define LOGQ    12u
#define Q       int64_t(uint64_t (1) << LOGQ)
#endif

#ifdef NTRUHRSS701
// -ntruhrss701
// These parameters are correct for any choice of√q > 8 2(n − 1).
// The q recommended here is the smallest power of two that provides correctness
#define N       701u
#define LOGQ    13u
#define Q       int64_t(uint64_t (1) << LOGQ)
#endif

#define N       509u
#define LOGQ    11u
#define Q       int64_t(uint64_t(1) << LOGQ)

int main() {
	auto seed = const_uniform_distribution_array<int32_t, 2*N>
	        (-(Q/2), Q/2);

	for (auto &el: seed)
		std::cout << signed(el) << " ";
	std::cout << "\n";

	int32_t f[N] = {0};
	int32_t g[N] = {0};
	int32_t h[N] = {0};

	const_ntru_keygen<N, LOGQ, Q>(f, g, h, seed.data());

	std::cout << "f:";
	for (auto &el: f)
		std::cout << el << " ";
	std::cout << "\n";

	std::cout << "g:";
	for (auto &el: g)
		std::cout << el << " ";
	std::cout << "\n";

	std::cout << "h:";
	for (auto &el: h)
		std::cout << el << " ";
	std::cout << "\n";
}