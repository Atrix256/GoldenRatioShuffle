#define _CRT_SECURE_NO_WARNINGS // for stb

#include <stdio.h>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <direct.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

typedef unsigned int uint;

typedef std::array<uint, 2> uint2;

typedef std::array<float, 2> float2;
typedef std::array<float, 3> float3;

// -1 means non deterministic
#define SEED() 0

#define DO_TEST_COPRIMES() true
// When looking for coprimes near the golden ratio integer, this is the range looked at
static const uint c_coprimeTestStart = 2;
static const uint c_coprimeTestEnd = 1000;

#define DO_TEST_1D() true
#define PRINT_NUMBERS_1D() false
static const uint c_itemCounts1D[] = { 10, 37, 64, 100 , 1000, 1337, 1546 };

#define DO_TEST_2D() true
#define PRINT_NUMBERS_2D() false
#define PRINT_FRAMES_2D() false
#define MAKE_VIDEOS() (PRINT_FRAMES_2D() && true)
static const uint2 c_imageSizes2D[] = { {16, 32} };// { 64, 64 }, { 256, 256 }, { 640, 480 } };

template <typename T>
T Frac(T f)
{
	return std::fmod(f, T(1.0));
}

static const float c_goldenRatio = (1.0f + std::sqrt(5.0f)) / 2.0f;
static const float c_goldenRatioFract = Frac(c_goldenRatio); // Same value as 1.0f / golden ratio aka golden ratio conjugate

// A generalization of the golden ratio to 2D, from http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
static const float c_R2_g = 1.32471795724474602596f;
static const float c_R2_a1 = 1 / c_R2_g;
static const float c_R2_a2 = 1 / (c_R2_g * c_R2_g);

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

// R2 is from http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
// R2 Low discrepancy sequence
float2 R2LDS(int index)
{
	return float2{ Frac(float(index) * c_R2_a1), Frac(float(index) * c_R2_a2) };
}

// From https://blog.demofox.org/2015/01/24/programmatically-calculating-gcd-and-lcm/
uint CalculateGCD(uint smaller, uint larger)
{
	// make sure A <= B before starting
	if (larger < smaller)
		std::swap(smaller, larger);

	// loop
	while (1)
	{
		// if the remainder of larger / smaller is 0, they are the same
		// so return smaller as the GCD
		uint remainder = larger % smaller;
		if (remainder == 0)
			return smaller;

		// otherwise, the new larger number is the old smaller number, and
		// the new smaller number is the remainder
		larger = smaller;
		smaller = remainder;
	}
}

bool IsCoprime(uint A, uint B)
{
	return CalculateGCD(A, B) == 1;
}

// Find the coprime of n that is nearest to target
uint GetNearestCoprime(uint target, uint n)
{
	uint coprime = 0;

	uint offset = 0;
	while (1)
	{
		if (offset < target)
		{
			coprime = target - offset;
			if (IsCoprime(coprime, n))
				break;
		}

		coprime = target + offset + 1;
		if (coprime < n && IsCoprime(coprime, n))
			break;

		offset++;
	}

	return coprime;
}

// Find the coprime of m and n that is nearest to target
uint GetNearestCoprime(uint target, uint m, uint n)
{
	uint coprime = 0;

	uint offset = 0;
	while (1)
	{
		if (offset < target)
		{
			coprime = target - offset;
			if (IsCoprime(coprime, m) && IsCoprime(coprime, n))
				break;
		}

		coprime = target + offset + 1;
		if (coprime < n && IsCoprime(coprime, m) && IsCoprime(coprime, n))
			break;

		offset++;
	}

	return coprime;
}

void MakeVideo(const char* label, const uint2& dims)
{
	if (!MAKE_VIDEOS())
		return;

	char buffer[1024];
	sprintf_s(buffer, "ffmpeg -framerate 60 -i out/%s_%u_%u_%%d.png -c:v libx264 -r 60 %s_%u_%u.mp4", label, dims[0], dims[1], label, dims[0], dims[1]);
	system(buffer);
}

template <typename LAMBDA>
void DoTest1D_Single(uint itemCount, const LAMBDA& lambda)
{
	std::vector<int> visitCount(itemCount, 0);
	int errorCount = 0;
	for (uint index = 0; index < itemCount; ++index)
	{
		uint shuffleItem = lambda();

		if(PRINT_NUMBERS_1D())
			printf("%u ", shuffleItem);
		visitCount[shuffleItem]++;
		if (visitCount[shuffleItem] > 1)
			errorCount++;
	}

	if (PRINT_NUMBERS_1D())
	{
		printf("\n");
		if (errorCount > 0)
		{
			printf("Duplicates:");
			for (uint index = 0; index < itemCount; ++index)
			{
				if (visitCount[index] > 1)
					printf(" %u", index);
			}
			printf("\n");

			printf("Missing:");
			for (uint index = 0; index < itemCount; ++index)
			{
				if (visitCount[index] == 0)
					printf(" %u", index);
			}
			printf("\n");
		}
	}

	printf("%i Errors (%0.2f%%)\n\n", errorCount, 100.0f * float(errorCount) / float(itemCount));
}

void DoTest1D(uint itemCount, std::mt19937& rng)
{
	printf("========== 1D: %i items ==========\n\n", itemCount);

	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float startValue = dist(rng);

	// Golden Ratio LDS
	{
		printf("GR LDS:");

		if (PRINT_NUMBERS_1D())
			printf("\n");
		else
			printf(" ");

		float valueF = startValue;
		DoTest1D_Single(itemCount,
			[&valueF, itemCount] ()
			{
				uint ret = Clamp<uint>((uint)(valueF * float(itemCount)), 0, itemCount - 1);
				valueF = Frac(valueF + c_goldenRatioFract);
				return ret;
			}
		);
	}

	// Integer solution
	{
		uint coprime = GetNearestCoprime(uint(c_goldenRatioFract * float(itemCount) + 0.5f), itemCount);

		printf("GR INT (%u vs %0.2f):", coprime, c_goldenRatioFract * float(itemCount));

		if (PRINT_NUMBERS_1D())
			printf("\n");
		else
			printf(" ");

		uint valueI = Clamp<uint>((uint)(startValue * float(itemCount)), 0, itemCount - 1);

		DoTest1D_Single(itemCount,
			[&valueI, itemCount, coprime]()
			{
				uint ret = valueI;
				valueI = (valueI + coprime) % itemCount;
				return ret;
			}
		);
	}
}

template <typename LAMBDA>
void DoTest2D_Single(const uint2& dims, const LAMBDA& lambda)
{
	std::vector<int> visitCount(dims[0]*dims[1], 0);

	int errorCount = 0;
	for (uint index = 0; index < dims[0] * dims[1]; ++index)
	{
		float percent = float(index) / float(dims[0] * dims[1] - 1);
		uint2 shuffleItem = lambda(index, percent);

		if (PRINT_NUMBERS_2D())
			printf("(%u,%u) ", shuffleItem[0], shuffleItem[1]);

		uint shuffleItemFlat = shuffleItem[1] * dims[0] + shuffleItem[0];

		visitCount[shuffleItemFlat]++;
		if (visitCount[shuffleItemFlat] > 1)
			errorCount++;
	}

	if (PRINT_NUMBERS_2D())
	{
		printf("\n");
		if (errorCount > 0)
		{
			printf("Duplicates:");
			for (uint index = 0; index < dims[0] * dims[1]; ++index)
			{
				if (visitCount[index] > 1)
					printf(" (%u,%u)", index / dims[0], index % dims[0]);
			}
			printf("\n");

			printf("Missing:");
			for (uint index = 0; index < dims[0] * dims[1]; ++index)
			{
				if (visitCount[index] == 0)
					printf(" (%u,%u)", index / dims[0], index % dims[0]);
			}
			printf("\n");
		}
	}

	printf("%i Errors (%0.2f%%)\n\n", errorCount, 100.0f * float(errorCount) / float(dims[0] * dims[1]));
}

void DoTest2D(const uint2& dims, std::mt19937& rng)
{
	printf("========== 2D: %u x %u ==========\n\n", dims[0], dims[1]);

	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float3 startValue = { dist(rng), dist(rng), dist(rng) };

	// R2 LDS
	if(true)
	{
		std::vector<unsigned char> image(dims[0] * dims[1], 0);
		std::vector<unsigned char> image2(dims[0] * dims[1], 0);

		printf("R2 LDS:");

		if (PRINT_NUMBERS_2D())
			printf("\n");
		else
			printf(" ");

		float2 valueF = float2{ startValue[0], startValue[1] };
		DoTest2D_Single(dims,
			[&valueF, &dims, &image, &image2](uint index, float percent)
			{
				uint2 ret;
				ret[0] = Clamp<uint>((uint)(valueF[0] * float(dims[0])), 0, dims[0] - 1);
				ret[1] = Clamp<uint>((uint)(valueF[1] * float(dims[1])), 0, dims[1] - 1);

				valueF[0] = Frac(valueF[0] + c_R2_a1);
				valueF[1] = Frac(valueF[1] + c_R2_a2);

				image[ret[1] * dims[0] + ret[0]] = (unsigned char)Clamp(percent * 255.0f, 0.0f, 255.0f);
				image2[ret[1] * dims[0] + ret[0]] = 255;

				if (PRINT_FRAMES_2D())
				{
					char fileName[256];
					sprintf_s(fileName, "out/R2_%u_%u_%u.png", dims[0], dims[1], index + 1);
					stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image2.data(), 0);
				}

				return ret;
			}
		);

		char fileName[256];
		sprintf_s(fileName, "out/R2_%u_%u.png", dims[0], dims[1]);
		stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image.data(), 0);

		MakeVideo("R2", dims);
	}

	// Integer offset solution
	if (false)
	{
		// Note: the idea was that the first dims[1] points of R2Int are good. so, we could do them, then do them again at a horizontal offset.
		// We could LDS shuffle those offsets, and possibly (LDS?) shuffle the order they are done in each section.
		// It turns out though it isn't dims[1], it's the "a1 integer" which is 48 for an image with width 64.
		// Going to pursue a different idea

		std::vector<unsigned char> image(dims[0] * dims[1], 0);
		std::vector<unsigned char> image2(dims[0] * dims[1], 0);

		uint2 R2Offset = {
			uint(c_R2_a1 * float(dims[0]) + 0.5f),
			uint(c_R2_a2 * float(dims[1]) + 0.5f)
		};

		uint R2OffsetFlat = R2Offset[1] * dims[0] + R2Offset[0];

		uint coprime = GetNearestCoprime(R2OffsetFlat, dims[0] * dims[1]);

		uint2 coprime2D = { coprime % dims[0], coprime / dims[0] };

		printf("R2 Int off (%u, %u vs %0.2f, %0.2f):", coprime2D[0], coprime2D[1], c_R2_a1 * float(dims[0]), c_R2_a2 * float(dims[1]));

		// NOTE: The first "next to" point happens when the 48th point occurs (frame 47).  The a1 value is 48.31.  The a2 value is 36.47.

		if (PRINT_NUMBERS_2D())
			printf("\n");
		else
			printf(" ");

		uint2 valueI2D;
		valueI2D[0] = Clamp<uint>((uint)(startValue[0] * float(dims[0])), 0, dims[0] - 1);
		valueI2D[1] = Clamp<uint>((uint)(startValue[1] * float(dims[1])), 0, dims[1] - 1);

		uint valueI = valueI2D[1] * dims[0] + valueI2D[0];

		DoTest2D_Single(dims,
			[&valueI, &valueI2D, &dims, coprime, &image, &image2](uint index, float percent)
			{
				/*
				if ((index % dims[1]) == 0)
				{
					valueI = valueI2D[1] * dims[0] + valueI2D[0];
				}
				*/

				// TODO: golden ratio shuffled horizontal offset

				uint2 ret = { valueI % dims[0], valueI / dims[0] };
				valueI = (valueI + coprime) % (dims[0] * dims[1]);

				image[ret[1] * dims[0] + ret[0]] = (unsigned char)Clamp(percent * 255.0f, 0.0f, 255.0f);
				image2[ret[1] * dims[0] + ret[0]] = 255;

				if (PRINT_FRAMES_2D())
				{
					char fileName[256];
					sprintf_s(fileName, "out/R2IntOff_%u_%u_%u.png", dims[0], dims[1], index + 1);
					stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image2.data(), 0);
				}

				return ret;
			}
		);

		char fileName[256];
		sprintf_s(fileName, "out/R2IntOff_%u_%u.png", dims[0], dims[1]);
		stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image.data(), 0);

		MakeVideo("R2IntOff", dims);
	}

	// Integer solution
	if (false)
	{
		std::vector<unsigned char> image(dims[0] * dims[1], 0);
		std::vector<unsigned char> image2(dims[0] * dims[1], 0);

		uint2 R2Offset = {
			uint(c_R2_a1 * float(dims[0]) + 0.5f),
			uint(c_R2_a2 * float(dims[1]) + 0.5f)
		};

		uint R2OffsetFlat = R2Offset[1] * dims[0] + R2Offset[0];

		uint coprime = GetNearestCoprime(R2OffsetFlat, dims[0] * dims[1]);

		uint2 coprime2D = { coprime % dims[0], coprime / dims[0] };

		printf("R2 Int (%u, %u vs %0.2f, %0.2f):", coprime2D[0], coprime2D[1], c_R2_a1 * float(dims[0]), c_R2_a2 * float(dims[1]));

		if (PRINT_NUMBERS_2D())
			printf("\n");
		else
			printf(" ");

		uint2 valueI2D;
		valueI2D[0] = Clamp<uint>((uint)(startValue[0] * float(dims[0])), 0, dims[0] - 1);
		valueI2D[1] = Clamp<uint>((uint)(startValue[1] * float(dims[1])), 0, dims[1] - 1);

		uint valueI = valueI2D[1] * dims[0] + valueI2D[0];

		DoTest2D_Single(dims,
			[&valueI, &dims, coprime, &image, &image2](uint index, float percent)
			{
				uint2 ret = { valueI % dims[0], valueI / dims[0] };
				valueI = (valueI + coprime) % (dims[0] * dims[1]);

				image[ret[1] * dims[0] + ret[0]] = (unsigned char)Clamp(percent * 255.0f, 0.0f, 255.0f);
				image2[ret[1] * dims[0] + ret[0]] = 255;

				if (PRINT_FRAMES_2D())
				{
					char fileName[256];
					sprintf_s(fileName, "out/R2Int_%u_%u_%u.png", dims[0], dims[1], index + 1);
					stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image2.data(), 0);
				}

				return ret;
			}
		);

		char fileName[256];
		sprintf_s(fileName, "out/R2Int_%u_%u.png", dims[0], dims[1]);
		stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image.data(), 0);

		MakeVideo("R2Int", dims);
	}

	// 1D Integer solution
	if (false)
	{
		std::vector<unsigned char> image(dims[0] * dims[1], 0);
		std::vector<unsigned char> image2(dims[0] * dims[1], 0);

		uint coprime = GetNearestCoprime(uint(c_goldenRatioFract * float(dims[0] * dims[1]) + 0.5f), dims[0] * dims[1]);

		uint2 coprime2D = { coprime % dims[0], coprime / dims[0] };

		printf("GR Int (%u, %u):", coprime2D[0], coprime2D[1]);

		if (PRINT_NUMBERS_2D())
			printf("\n");
		else
			printf(" ");

		uint2 valueI2D;
		valueI2D[0] = Clamp<uint>((uint)(startValue[0] * float(dims[0])), 0, dims[0] - 1);
		valueI2D[1] = Clamp<uint>((uint)(startValue[1] * float(dims[1])), 0, dims[1] - 1);

		uint valueI = valueI2D[1] * dims[0] + valueI2D[0];

		DoTest2D_Single(dims,
			[&valueI, &dims, coprime, &image, &image2](uint index, float percent)
			{
				uint2 ret = { valueI % dims[0], valueI / dims[0] };
				valueI = (valueI + coprime) % (dims[0] * dims[1]);

				image[ret[1] * dims[0] + ret[0]] = (unsigned char)Clamp(percent * 255.0f, 0.0f, 255.0f);
				image2[ret[1] * dims[0] + ret[0]] = 255;

				if (PRINT_FRAMES_2D())
				{
					char fileName[256];
					sprintf_s(fileName, "out/GRInt_%u_%u_%u.png", dims[0], dims[1], index + 1);
					stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image2.data(), 0);
				}

				return ret;
			}
		);

		char fileName[256];
		sprintf_s(fileName, "out/GRInt_%u_%u.png", dims[0], dims[1]);
		stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image.data(), 0);

		MakeVideo("GRInt", dims);
	}

	// Multi Shuffle
	if (true)
	{
		// Note: the idea here is we could go y=0 to height and do a 1D LDS shuffle to get each x.
		// That happens in vertical order so we could (LDS?) shuffle the order it happens in.
		// After that we only have 1 pixel set per row.
		// We can do the process again with an (LDS?) shuffled x axis offset.

		std::vector<unsigned char> image(dims[0] * dims[1], 0);
		std::vector<unsigned char> image2(dims[0] * dims[1], 0);

		// This irrational controls the LDS shuffle as we go vertically down the screen
		static const float c_irrational1 = Frac(c_goldenRatio);
		static const float c_irrational2 = Frac(std::sqrt(2.0f));
		static const float c_irrational3 = Frac(std::sqrt(3.0f));

		uint coprime1 = GetNearestCoprime(uint(c_irrational1 * float(dims[0]) + 0.5f), dims[0]);
		uint coprime2 = GetNearestCoprime(uint(c_irrational2 * float(dims[1]) + 0.5f), dims[1]);
		uint coprime3 = GetNearestCoprime(uint(c_irrational3 * float(dims[0]) + 0.5f), dims[0]);

		printf("Multi Shuffle (%u,%u,%u vs %0.2f,%0.2f,%0.2f):", coprime1, coprime2, coprime3, c_irrational1 * float(dims[0]), c_irrational2 * float(dims[1]), c_irrational3 * float(dims[0]));

		if (PRINT_NUMBERS_2D())
			printf("\n");
		else
			printf(" ");

		uint valueI1 = Clamp<uint>((uint)(startValue[0] * float(dims[0])), 0, dims[0] - 1);
		uint valueI2 = Clamp<uint>((uint)(startValue[1] * float(dims[1])), 0, dims[1] - 1);
		uint valueI3 = Clamp<uint>((uint)(startValue[2] * float(dims[1])), 0, dims[1] - 1);

		DoTest2D_Single(dims,
			[&valueI1, &valueI2, &valueI3, coprime1, coprime2, coprime3, &image, &image2, &dims](uint index, float percent)
			{
				if ((index % dims[0]) == 0)
					valueI3 = (valueI3 + coprime3) % dims[0];

				uint2 ret;
				ret[0] = (valueI1 + valueI3) % dims[0];
				ret[1] = valueI2;

				valueI1 = (valueI1 + coprime1) % dims[0];
				valueI2 = (valueI2 + coprime2) % dims[1];

				image[ret[1] * dims[0] + ret[0]] = (unsigned char)Clamp(percent * 255.0f, 0.0f, 255.0f);
				image2[ret[1] * dims[0] + ret[0]] = 255;

				if (PRINT_FRAMES_2D())
				{
					char fileName[256];
					sprintf_s(fileName, "out/MS_%u_%u_%u.png", dims[0], dims[1], index + 1);
					stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image2.data(), 0);
				}

				return ret;
			}
		);

		char fileName[256];
		sprintf_s(fileName, "out/MS_%u_%u.png", dims[0], dims[1]);
		stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image.data(), 0);

		MakeVideo("MS", dims);
	}
}

int main(int argc, char** argv)
{
	_mkdir("out");

	// initialize RNG
	std::random_device rd;
	uint seed = (SEED() == -1) ? rd() : SEED();
	printf("Seed = %u\n\n", seed);

	// Coprime test
	if (DO_TEST_COPRIMES())
	{
		FILE* file = nullptr;
		fopen_s(&file, "out/coprimes.csv", "wb");

		fprintf(file, "\"itemCount\",\"goldenIndexF\",\"goldenIndexI\",\"coprime\",\"diffF\",\"diffI\"\n");

		uint minDiff = 0;
		uint maxDiff = 0;

		for (uint itemCount = c_coprimeTestStart; itemCount < c_coprimeTestEnd; ++itemCount)
		{
			uint target = uint(c_goldenRatioFract * float(itemCount) + 0.5f);
			uint coprime = GetNearestCoprime(target, itemCount);

			uint diff = (target >= coprime) ? target - coprime : coprime - target;

			if (itemCount == c_coprimeTestStart)
			{
				minDiff = diff;
				maxDiff = diff;
			}
			else
			{
				minDiff = std::min(minDiff, diff);
				maxDiff = std::max(maxDiff, diff);
			}

			float diffF = std::abs(c_goldenRatioFract * float(itemCount) - float(coprime));

			fprintf(file, "\"%u\",\"%f\",\"%u\",\"%u\",\"%f\",\"%u\"\n", itemCount, c_goldenRatioFract * float(itemCount), uint(c_goldenRatioFract * float(itemCount) + 0.5f), coprime, diffF, diff);
		}

		printf(
			"When looking for coprimes near the golden ratio index, the diffs were:\n"
			"  min: %u\n"
			"  max: %u\n"
			"\n",
			minDiff, maxDiff
		);

		fclose(file);
	}

	// 1D Tests
	if (DO_TEST_1D())
	{
		for (uint itemCount : c_itemCounts1D)
		{
			std::mt19937 rng(seed);
			DoTest1D(itemCount, rng);
		}
	}

	// 2D Tests
	if (DO_TEST_2D())
	{
		std::mt19937 rng(seed);
		for (const uint2& size : c_imageSizes2D)
			DoTest2D(size, rng);
	}

	return 0;
}

// TODO: at 640x480, MS had 50% error rate! maybe you mixed up w/h somewhere? could try a smaller non square texture 

/*
TODO:


Clean this up?
think about how this should go out

- ffmpeg -framerate 30 -i out/GRInt_64_64_%d.png -c:v libx264 -r 30 GRInt_64_64.mp4
 - https://shotstack.io/learn/use-ffmpeg-to-convert-images-to-video/

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