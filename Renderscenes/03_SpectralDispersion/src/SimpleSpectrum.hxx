/*
* Converts RGB to a suitable SPD.
* The algorithm is based on a paper by Brian Smits: http://www.cs.utah.edu/~bes/papers/color/paper.pdf
* Code is partly derived from a C# project which aims to visually illustrate the algorithm: https://github.com/wip-/RgbToSpectrum 
*/

#include <vector>
#define BinsCount 21

enum Primary : int { R, G, B, C, M, Y, W, Count };

class SimpleSpectrum
{
public:
	static std::vector<double> Lambdas;

	static std::vector<double> Rspectrum;
	static std::vector<double> Gspectrum;
	static std::vector<double> Bspectrum;

	static std::vector<double> Cspectrum;
	static std::vector<double> Mspectrum;
	static std::vector<double> Yspectrum;

	static std::vector<double> Wspectrum;

	double static SamplePrimarySpectrum(Primary primary, int index)
	{
		std::vector<std::vector<double>> Spectra { Rspectrum, Gspectrum, Bspectrum, Cspectrum, Mspectrum, Yspectrum, Wspectrum };
		return Spectra[(int)primary][index];
	}

	double values[BinsCount];

	//static double LambdaMin;
	static double getLambdaMin() { return Lambdas[0]; }
	//static double LambdaMax;
	static double getLambdaMax() { return Lambdas[BinsCount - 1]; }
	//static double LambdaStep;
	static double getLambdaStep() { return Lambdas[1] - Lambdas[0]; }

		// R,  G,  B must be between [0,  1]
		SimpleSpectrum(double r, double g, double b)
	{
		double Rweight = 0;
		double Gweight = 0;
		double Bweight = 0;
		double Cweight = 0;
		double Mweight = 0;
		double Yweight = 0;
		double Wweight = 0;

		if (r <= g && g <= b)
		{
			Wweight = r;
			Cweight = g - r;
			Bweight = b - g;
		}
		else if (r <= b && b <= g)
		{
			Wweight = r;
			Cweight = b - r;
			Gweight = g - b;
		}
		else if (g <= r && r <= b)
		{
			Wweight = g;
			Mweight = r - g;
			Bweight = b - r;
		}
		else if (g <= b && b <= r)
		{
			Wweight = g;
			Mweight = b - g;
			Rweight = r - b;
		}
		else if (b <= r && r <= g)
		{
			Wweight = b;
			Yweight = r - b;
			Gweight = g - r;
		}
		else if (b <= g && g <= r)
		{
			Wweight = b;
			Yweight = g - b;
			Rweight = r - g;
		}
		else
		{
			//Debugger.Break();
		}

		for (int i = 0; i<BinsCount; ++i)
		{
			values[i] =
				Wweight * Wspectrum[i] +
				Rweight * Rspectrum[i] +
				Gweight * Gspectrum[i] +
				Bweight * Bspectrum[i] +
				Cweight * Cspectrum[i] +
				Mweight * Mspectrum[i] +
				Yweight * Yspectrum[i];
		}
	}

	double Sample(double lambda)
	{
		// find closest smaller lambda
		int i;
		for (i = 0; i < BinsCount - 1; ++i)
		{
			if (Lambdas[i + 1] > lambda)
				break;
		}
		if (lambda >= Lambdas[BinsCount - 1])
			i = BinsCount - 1;

		return values[i];
	}
};

//initialization of static variables
// values read from Brian Smits' paper curve images pixels - approximative for the least
std::vector<double> SimpleSpectrum::Lambdas = { 361.856, 383.505, 405.155, 426.804, 448.454, 470.103, 491.753, 513.402, 535.052, 556.701, 578.351, 600.000, 621.649, 643.299, 664.948, 686.598, 708.247, 729.897, 751.546, 773.196, 794.845 };
std::vector<double> SimpleSpectrum::Rspectrum = { 0.09559, 0.09559, 0.08824, 0.07353, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.08824, 0.69853, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000 };
std::vector<double> SimpleSpectrum::Gspectrum = { 0.00000, 0.00000, 0.00000, 0.00000, 0.03676, 0.38971, 0.78676, 1.00000, 1.00000, 0.77941, 0.31618, 0.00735, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000 };
std::vector<double> SimpleSpectrum::Bspectrum = { 0.99265, 1.00000, 1.00000, 1.00000, 0.86765, 0.61029, 0.30882, 0.08088, 0.00000, 0.00000, 0.00000, 0.02941, 0.05147, 0.06618, 0.06618, 0.06618, 0.06618, 0.06618, 0.06618, 0.06618, 0.06618 };

std::vector<double> SimpleSpectrum::Cspectrum = { 0.99265, 0.99265, 0.98529, 0.92647, 0.98529, 1.00000, 1.00000, 1.00000, 1.00000, 0.91912, 0.30147, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00735, 0.00735, 0.00735, 0.00735, 0.00735 };
std::vector<double> SimpleSpectrum::Mspectrum = { 1.00000, 1.00000, 1.00000, 1.00000, 0.97794, 0.61765, 0.19853, 0.00000, 0.00000, 0.24265, 0.68382, 0.98529, 1.00000, 1.00000, 1.00000, 0.99265, 0.97794, 0.97794, 0.98529, 0.98529, 0.98529 };
std::vector<double> SimpleSpectrum::Yspectrum = { 0.00735, 0.00735, 0.00000, 0.00000, 0.30147, 0.39706, 0.69853, 0.92647, 1.00000, 1.00000, 1.00000, 0.97059, 0.95588, 0.95588, 0.95588, 0.96324, 0.97059, 0.98529, 0.99265, 1.00000, 1.00000 };

//std::vector<double> SimpleSpectrum::Wspectrum = { 0.98750, 0.99554, 1.05714, 1.07589, 0.89643, 0.99554, 1.11607, 1.07857, 1.00089, 0.84286, 1.01964, 1.04643, 1.05179, 1.06250, 1.06518, 1.06518, 1.06250, 1.06518, 1.06786, 1.06786, 1.06518 };
std::vector<double> SimpleSpectrum::Wspectrum = { 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000 };