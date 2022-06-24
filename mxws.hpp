
#include <numeric>
#include <random>

template <typename RN>
	requires
std::same_as<RN, uint32_t> ||
std::same_as<RN, uint64_t>
class mxws
{

private:

	std::random_device r;

public:

	uint64_t x, w, x1, x2, w1, w2;

	typedef RN result_type;

	void seed()
	{
		init();
		x1 = x2 = 1;
	}

	mxws(const std::seed_seq& seq)
	{
		if (seq.size() == 2)
		{
			std::vector<uint32_t> seeds(seq.size());
			seq.param(seeds.rbegin());
			w = (uint64_t(seeds[1]) << 32) | seeds[0];
			x = 1;

			w1 = w;
			w2 = w1 + 1;
			x1 = x2 = 1;
		}

		else init();
	}

	mxws()
	{
		init();
	}

	mxws(const uint64_t& seed)
	{
		init();
	}

	void init()
	{
		w = (uint64_t(r()) << 32) | r();
		x = 1;
		w1 = w;
		w2 = w1 + 1;
		x1 = x2 = 1;
	}

	void init(const uint64_t& seed)
	{
		w = seed;
		x = 1;
	}

	virtual ~mxws() = default;

	static constexpr RN min() { return std::numeric_limits<RN>::min(); }
	static constexpr RN max() { return std::numeric_limits<RN>::max(); }

	inline RN operator()()
		requires
	std::same_as<RN, uint32_t>
	{
		w += x = std::rotr(x *= w, 32);
		return RN(x);
	}

	inline RN operator()()
		requires
	std::same_as<RN, uint64_t>
	{
		x1 *= w1;
		x1 = std::rotr(x1, 32);
		w1 += x1;

		x2 *= w2;
		x2 = std::rotr(x2, 32);
		w2 += x2;

		return (x1 << 32) | uint32_t(x2);
	}

	template <typename T>
		requires std::floating_point<T>&&
	std::same_as<RN, uint64_t>
		inline T operator()(const T& f)
	{
		return ((*this)() >> 11) / T(9007199254740992) * f;
	}

	template <typename T>
		requires std::floating_point<T>&&
	std::same_as<RN, uint32_t>
		inline T operator()(const T& f)
	{
		return (*this)() / T(4294967296) * f;
	}

	template <typename T>
		requires std::floating_point<T>
	inline T operator()(const T& min, const T& max)
	{
		return (*this)(1.0) * (max - min) + min;
	}

	template <typename T, typename U>
		requires std::integral<T>&& std::floating_point<U>
	inline U operator()(const T& min, const U& max)
	{
		return (*this)(1.0) * (max - min) + min;
	}

	template <typename T>
		requires std::integral<T>
	inline T operator()(const T& max)
	{
		return (*this)() % (max + 1);
	}

	template <typename T>
		requires std::integral<T>
	inline T operator()(const T& min, const T& max)
	{
		return min + ((*this)() % (max - min + 1));
	}

	template <typename T>
		requires std::floating_point<T>
	T round_to_half(T in)
	{
		in *= 2; in = round(in); return in /= 2;
	}

	template <typename T>
		requires std::floating_point<T>
	inline T to_float(const T& in)
	{
		return in / T(4294967296);
	}

	template <typename T>
		requires std::integral<T>
	inline uint32_t to_int(const T& n)
	{
		return uint32_t(n >> 32);
	}

	template <typename T>
		requires std::floating_point<T>
	T inline normalRandom(const T& mean, const T& sigma)
	{
		// return a normally distributed random value
		T v1 = (*this)(1.0);
		T v2 = (*this)(1.0);
		return std::cos(2 * std::numbers::pi * v2) * std::sqrt(-2 * std::log(v1)) * sigma + mean;
	}

	template <typename T, typename L>
		requires std::floating_point<T>&&
	std::integral<L>
		T inline error_function(const T& x, const L& N)
	{
		const auto f = [](const auto& t) {return exp(-pow(t, 2)); };

		const auto ymin = f(x);

		L	t = 0;

		for (auto i = 0; i < N; i++)
			if ((ymin + (*this)(1. - ymin)) <= f((*this)(x)))
				t++;

		T k = T(t) / N;
		k *= x * (1. - ymin);
		k += x * ymin;
		k *= 2 / sqrt(std::numbers::pi);

		std::cout << "error function(" << x << ") = " << k << std::endl;
		std::cout << " std::function(" << x << ") = "
			<< std::erf(x) << std::endl;

		return k;
	}

	template <typename R, typename I, typename L>
		requires
	std::same_as<R, double>&&
		std::integral<I>&&
		std::same_as<L, std::uint64_t>
		std::tuple<R, I> inline Probability_Wave(const I& cycle_SIZE,
			std::vector<I>& cycle, const I& N_cycles, const L& TRIALS) {

		const I cycle_size = I(round(log(cycle_SIZE * 6) * pow(tan(36 / 5.), 2)));
		const I rn_range = I(floor(cycle_SIZE / sqrt(log2(cycle_SIZE))));

		//const I cycle_size = cycle_SIZE;
		//const I rn_range = I(round(sqrt(cycle_size) + log(cycle_size / 4)));

		L random_walk = 0;

		for (L i = 0; i < TRIALS; i++, random_walk = 0)
		{
			for (I j = 0; j < cycle_size; j++)
				random_walk += (*this)();

			cycle[to_int(random_walk * rn_range) % cycle_SIZE]++;
		}

		return std::make_tuple(rn_range, cycle_size);
	}
};