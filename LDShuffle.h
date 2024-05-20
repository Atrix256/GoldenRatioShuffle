// A stateless / constant time shuffle iterator which shuffles the values into a low discrepancy sequence.
// Supports random access and inversion.
// By Alan Wolfe.  May, 2024
// https://github.com/Atrix256/GoldenRatioShuffle
// https://blog.demofox.org/
// MIT licensed. Enjoy!

#pragma once

#include <cmath>
#include <algorithm>

class LDShuffle
{
public:
	typedef unsigned int uint;

	static inline const float c_goldenRatio = (1.0f + std::sqrt(5.0f)) / 2.0f;

	LDShuffle(uint numItems, uint seed, float irrational = c_goldenRatio)
		: m_numItems(numItems)
		, m_seed(seed)
		, m_currentItem(seed)
	{
		CalculateCoprime(irrational);

		// Calculate "steps to unity" which is needed for inversion
		uint s;
		uint GCD = ExtendedEuclidianAlgorithm(m_coprime, m_numItems, s, m_stepsToUnity);
	}

	uint Next()
	{
		uint ret = m_currentItem;
		m_currentItem = (m_currentItem + m_coprime) % m_numItems;
		return ret;
	}

	uint GetValueAtIndex(uint index)
	{
		return ((index % m_numItems) * m_coprime + m_seed) % m_numItems;
	}

	uint GetIndexOfValue(uint value)
	{
		uint stepsToValue = (m_stepsToUnity * value) % m_numItems;
		uint stepsToSeed = (m_stepsToUnity * m_seed) % m_numItems;

		return (stepsToValue + m_numItems - stepsToSeed) % m_numItems;
	}

	uint GetCoprime() const { return m_coprime; }

private:
	void CalculateCoprime(float irrational)
	{
		uint target = uint(std::fmod(irrational, 1.0f) * float(m_numItems) + 0.5f);

		uint offset = 0;
		while (1)
		{
			if (offset < target)
			{
				m_coprime = target - offset;
				if (IsCoprime(m_coprime, m_numItems))
					break;
			}

			m_coprime = target + offset + 1;
			if (m_coprime < m_numItems && IsCoprime(m_coprime, m_numItems))
				break;

			offset++;
		}
	}

	static inline bool IsCoprime(uint A, uint B)
	{
		return CalculateGCD(A, B) == 1;
	}

	// From https://blog.demofox.org/2015/01/24/programmatically-calculating-gcd-and-lcm/
	static inline uint CalculateGCD(uint smaller, uint larger)
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

	// From https://blog.demofox.org/2015/09/10/modular-multiplicative-inverse/
	static inline uint ExtendedEuclidianAlgorithm(uint smaller, uint larger, uint& s, uint& t)
	{
		// make sure A <= B before starting
		bool swapped = false;
		if (larger < smaller)
		{
			swapped = true;
			std::swap(smaller, larger);
		}

		// set up our storage for the loop.  We only need the last two values so will
		// just use a 2 entry circular buffer for each data item
		uint remainders[2] = {larger, smaller};
		uint ss[2] = {1, 0};
		uint ts[2] = {0, 1};
		uint indexNeg2 = 0;
		uint indexNeg1 = 1;

		// loop
		while (1)
		{
			// calculate our new quotient and remainder
			uint newQuotient = remainders[indexNeg2] / remainders[indexNeg1];
			uint newRemainder = remainders[indexNeg2] - newQuotient * remainders[indexNeg1];

			// if our remainder is zero we are done.
			if (newRemainder == 0)
			{
				// return our s and t values as well as the quotient as the GCD
				s = ss[indexNeg1];
				t = ts[indexNeg1];
				if (swapped)
					std::swap(s, t);
				return remainders[indexNeg1];
			}

			// calculate this round's s and t
			int newS = ss[indexNeg2] - newQuotient * ss[indexNeg1];
			int newT = ts[indexNeg2] - newQuotient * ts[indexNeg1];

			// store our values for the next iteration
			remainders[indexNeg2] = newRemainder;
			ss[indexNeg2] = newS;
			ts[indexNeg2] = newT;

			// move to the next iteration
			std::swap(indexNeg1, indexNeg2);
		}
	}

private:
	uint m_numItems = 0;
	uint m_seed = 0;

	uint m_coprime = 0;
	uint m_stepsToUnity = 0;

	uint m_currentItem = 0;
};