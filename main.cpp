#include <stdio.h>
#include <cmath>
#include <random>
#include <vector>

// -1 means non deterministic
#define SEED() 0

#define PRINT_NUMBERS() true

static const float c_goldenRatio = (1.0f + std::sqrt(5.0f)) / 2.0f;
static const float c_goldenRatioFract = std::fmod(c_goldenRatio, 1.0f); // Same value as 1.0f / golden ratio aka golden ratio conjugate

static const unsigned int c_itemCounts[] = { 10, 37, 64, 100 };// , 1000, 1337, 1546 };

template <typename T>
T Clamp(T value, T themin, T themax)
{
	if (value <= themin)
		return themin;
	else if (value >= themax)
		return themax;
	else
		return value;
}

// From https://blog.demofox.org/2015/01/24/programmatically-calculating-gcd-and-lcm/
unsigned int CalculateGCD(unsigned int smaller, unsigned int larger)
{
	// make sure A <= B before starting
	if (larger < smaller)
		std::swap(smaller, larger);

	// loop
	while (1)
	{
		// if the remainder of larger / smaller is 0, they are the same
		// so return smaller as the GCD
		unsigned int remainder = larger % smaller;
		if (remainder == 0)
			return smaller;

		// otherwise, the new larger number is the old smaller number, and
		// the new smaller number is the remainder
		larger = smaller;
		smaller = remainder;
	}
}

bool IsCoprime(unsigned int A, unsigned int B)
{
	return CalculateGCD(A, B) == 1;
}

template <typename LAMBDA>
void DoTest_Single(unsigned int itemCount, const LAMBDA& lambda)
{
	std::vector<int> visitCount(itemCount, 0);
	int errorCount = 0;
	for (unsigned int index = 0; index < itemCount; ++index)
	{
		unsigned int shuffleItem = lambda();

		if(PRINT_NUMBERS())
			printf("%u ", shuffleItem);
		visitCount[shuffleItem]++;
		if (visitCount[shuffleItem] > 1)
			errorCount++;
	}

	if (PRINT_NUMBERS())
	{
		printf("\n");
		if (errorCount > 0)
		{
			printf("Duplicates:");
			for (unsigned int index = 0; index < itemCount; ++index)
			{
				if (visitCount[index] > 1)
					printf(" %u", index);
			}
			printf("\n");

			printf("Missing:");
			for (unsigned int index = 0; index < itemCount; ++index)
			{
				if (visitCount[index] == 0)
					printf(" %u", index);
			}
			printf("\n");
		}
	}

	printf("%i Errors (%0.2f%%)\n\n", errorCount, 100.0f * float(errorCount) / float(itemCount));
}

void DoTest(unsigned int itemCount, std::mt19937& rng)
{
	printf("========== %i items ==========\n\n", itemCount);

	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float startValue = dist(rng);

	// Golden Ratio LDS
	{
		printf("GRLDS:");

		if (PRINT_NUMBERS())
			printf("\n");
		else
			printf(" ");

		float valueF = startValue;
		DoTest_Single(itemCount,
			[&valueF, itemCount] ()
			{
				unsigned int ret = Clamp<unsigned int>((unsigned int)(valueF * float(itemCount)), 0, itemCount - 1);
				valueF = std::fmod(valueF + c_goldenRatioFract, 1.0f);
				return ret;
			}
		);
	}

	// Integer solution
	{
		unsigned int coprime = 0;

		// Get the coprime to itemCount which is the closest to the golden ratio index
		unsigned int goldenRatioIndex = unsigned int(c_goldenRatioFract * float(itemCount) + 0.5f);
		{
			unsigned int offset = 0;
			while (1)
			{
				if (offset < goldenRatioIndex)
				{
					coprime = goldenRatioIndex - offset;
					if (IsCoprime(coprime, itemCount))
						break;
				}

				coprime = goldenRatioIndex + offset + 1;
				if (coprime < itemCount && IsCoprime(coprime, itemCount))
					break;

				offset++;
			}
		}

		printf("GRINT (%u vs %0.2f):", coprime, c_goldenRatioFract * float(itemCount));

		if (PRINT_NUMBERS())
			printf("\n");
		else
			printf(" ");

		unsigned int valueI = Clamp<unsigned int>((unsigned int)(startValue * float(itemCount)), 0, itemCount - 1);

		DoTest_Single(itemCount,
			[&valueI, itemCount, coprime]()
			{
				unsigned int ret = valueI;
				valueI = (valueI + coprime) % itemCount;
				return ret;
			}
		);
	}
}

int main(int argc, char** argv)
{
	// initialize RNG
	std::random_device rd;
	unsigned int seed = (SEED() == -1) ? rd() : SEED();
	printf("Seed = %u\n\n", seed);
	std::mt19937 rng(seed);

	for (unsigned int itemCount : c_itemCounts)
	{
		DoTest(itemCount, rng);
	}

	return 0;
}

/*
TODO:
- work out 2D stuff? should we show some kind of image demo for it?

Blog:
- talk about how to extend it to 2d with R2?
 - need to think through how to make this work without repeats though. maybe co-irrational to each other, as well as the edge lengths? idk.
- GCD: https://blog.demofox.org/2015/01/24/programmatically-calculating-gcd-and-lcm/
- the coprime version of the golden ratio is real close to the actual thing. but drift.
- show where the sequnces diverge! lil bit of drift due to rounding to integer, but we want that cause it's a coprime generator.
- mention that you are not going to be analzying the quality of the sequence, but you can see it's very, very similar to the pure golden ratio one

find marc reynolds post about golden ratio integer LDS thing.
and link this too: ttps://mastodon.gamedev.place/@demofox/112452866510781822

*/