//Jacob Zahn, ECE1305 sec 004, Professor Ian Scott-Flemming, May 3rd 2017
//This program creates pgms of the night sky as veiwed from a satilite with an orbital path around the moon, Radius 200k km, closest approach 15k km, focal point in front of the moon on R
//This program uses the file ecliptic_stars_16_deg.txt, which should be located in the same folder as the .cpp file, no user input is needed.
//Partner: Dylan Scott



/*Attempt 2 and 3 are the same with the only difference being that 2 was commented, attempt 1 was incomplete*/
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>


using namespace std;

void fillArray(double *&p, double *&q, double *&r, int size);//p=ra,q=dec,r=mag
void inFOV(double *ra, vector<int> &stars, double hour);//checks if stars are in fov
void makePGM(int x, double*ra, double*dec, double*mag);//makes pgm and calls make pix
void makePix(int x, double*ra, double*dec, double*mag, int size, unsigned char *&pix);//makes an array that is written into the pgm
double makeHour(int x);//calculates RA hour we're looking at
double getGamma(int x);//gets the angle our radius is relative to RA(0hrs)
double getC(double a, double b, double gamma);//gets our distance from focal point for use in calculations
double getAlpha(double gamma, double a, double c);// gets the angle we're looking at above the RA in degrees
double getDPF(int x);//gets degrees that alpha changes per frame
double getAhDegree(double gamma, double alpha);//gets the degree we're looking at above the RA

											   //the program outputs 7200 pgm images of stars and a white moon
											   /*I was able to get as far as having an orbit with a varrying viewpoint of the moon that got close and then further away, never got to map the texture of the moon. The stars fade in a bell curve. Ra of view is
											   modeled correctly with inverse quadratics to calculate its change. moon is drawn and grows and rotates from right of view to left.*/

int main()
{
	double *RA = NULL;
	double *DEC = NULL;
	double *MAG = NULL;
	vector<int> *stars = new vector<int>{ 0 };
	cout << "calling fillArray" << endl;
	fillArray(RA, DEC, MAG, 1098);
	cout << "called fillArray" << endl;
	/*for (int x = 0; x < 1098; x++)
	{
	cout << x << " " << RA[x] << " " << DEC[x] << " " << MAG[x] << endl;
	}*/
	cout << "running 7200 for loop" << endl;
	for (int x = 0; x < 7200; x++) //7200 frames for 60fps over 2mins
	{
		cout << "calling make pgm " << x << endl;
		makePGM(x, RA, DEC, MAG);
		//cout << "called make pgm" << endl;
	}
	cout << "end of main loop" << endl;










	return 0;
}

void fillArray(double *&p, double *&q, double *&r, int size)
{
	ifstream ifs;
	ifs.open("ecliptic_stars_16_deg.txt");//opens data

	p = new double[size];
	q = new double[size];
	r = new double[size];
	for (int g = 0; g < size; g++)//fills corresponding data
	{
		ifs.ignore(4);
		ifs >> p[g];
		ifs >> q[g];
		ifs >> r[g];
		ifs.ignore(1);
	}
	ifs.close();
}

void inFOV(double *ra, vector<int> &stars, double hour)//checks if we can see stars +or- .8 hours from where we're looking
{

	double min = (-hour - 0.8);
	double max = (-hour + 0.8);
	//cout << min << "\t\t" << max << endl;
	if (min <= 0.0)
		min = min + 24.0;

	if (min <= 0)
		min = min + 24.0;

	if (max >= 24.0)
		max = max - 24.0;
	else if (max <= 0)
		max = max + 24.0;
	//cout << min << "\t\t\t" << max << endl;
	for (int x = 0; x < 1098; x++)
	{
		if (min < max)
		{
			if ((ra[x] > min) && (ra[x] < max))
				stars.push_back(x);
		}

		else if (max < min)
		{
			if (((ra[x] < max) || (ra[x] > min)))
				stars.push_back(x);
		}

	}
	//for (int i = 0; i < stars.size(); i++)
	//	cout << stars[i] << " ";
}

void makePGM(int x, double*ra, double*dec, double*mag)//instanciates variables that are neccesary to continue, then calls make pix, after pix is filled it fills the pgm binarily
{

	int size = 1200 * 800;
	unsigned char *pix = NULL;

	//cout << "calling make pix" << endl;
	makePix(x, ra, dec, mag, size, pix);
	if (pix == NULL)
		return;

	//cout << "called make pix" << endl;
	//cout << "creating file" << endl;
	string frameName = "image_frame_" + to_string(x) + ".pgm";
	ofstream ofs(frameName, ios::binary);
	//cout << "writing header" << endl;
	ofs << "P5\n" << "1200" << " " << "800" << "\n" << "255" << "\n";
	//cout << "filling pgm with pix" << endl;
	//for (int i = 0; i < size; i++)
	//ofs << pix[i] << " ";
	ofs.write((char *)(pix), size * sizeof(char));
	//cout << "filled pgm with pix" << endl;
	ofs.close();
	delete[] pix;
	pix = NULL;
	//deletes pix for memory management
}

void makePix(int x, double * ra, double * dec, double * mag, int size, unsigned char *&pix)
{
	vector<int> stars;//makes a variably sized container for stars
	double hour = makeHour(x);//calculates hour


							  //cout << "calling inFOV" << endl;
	inFOV(ra, stars, hour);//get a list of stars that are in the fov
						   //cout << "called inFOV" << endl;
	pix = new unsigned char[size];
	//cout << "made pix" << endl;
	for (int r = 0; r < 800; r++)
	{
		for (int c = 0; c<1200; c++)
			pix[((r) * 1200) + c] = 0;
	}
	//cout << "filled pix with black" << endl;

	//cout << "running first i for loop "<< stars.size() << endl;
	/*for (int i = 0; i < stars.size(); i++)
	cout << stars[i] << " ";*/
	//int count = 0;
	double *RA = NULL;
	double *DEC = NULL;
	double *MAG = NULL;
	RA = new double[stars.size()];
	DEC = new double[stars.size()];
	MAG = new double[stars.size()];


	for (int i = 0; i < stars.size(); i++)
	{
		RA[i] = ra[stars[i]];
		DEC[i] = dec[stars[i]];
		MAG[i] = mag[stars[i]];
	}

	for (int i = 0; i < stars.size(); i++)//for length of stars recieved
	{
		//cout << "running second r for loop "<< i << endl;
		double xPos = ((-RA[i] - (hour - .8))*(1200.0 / 1.6));//gets x y and magnitude per star
		double yPos = (((8.0) - DEC[i])*(800.0 / 16.0));
		double MagI = MAG[i];
		for (int r = 0; r < 800; r++)
		{
			//cout << "asking first if "<< r << endl;
			if ((r >= (yPos - 14.0)) && (r <= (yPos + 14.0)))//if star center is in radius of 14
			{
				//cout << "running third c for loop" << endl;
				for (int c = 0; c < 1200; c++)
				{
					//cout << "asking second if" << endl;
					if (xPos>-18000.0)//for stars out of pixel range
					{

						//cout << "passed first if " << ++count << endl;
						if ((c >= (xPos + 18000.0 - 14.0)) && (c <= (xPos + 18000.0 + 14.0)))
						{//            sqrt(         (x1   -            x2)                                                       ^2  +          ( y1  -                  y2                         )    ^2 )
							double d = sqrt((pow(((double)c - round(xPos + 18000.0)), 2.0) + pow((double)r - round(yPos), 2.0)));// calculation for distance from center of star
							double b = exp(-pow(d, 2.0) / (2.0*pow(((12.0 - MagI) / 2.354), 2.0))) *255.0;//calculation for brightness based of gausian curve
																										  //cout << "passed first if "<<++count << endl;
																										  //cout << "filling pix value" << endl;
							if (d < 13.0 - MagI)//checks for within magnitude
							{
								if ((int)pix[(1200 * r) + c] + b + 1 >= 255)//sets to limit if stars pass it
									pix[(1200 * r) + c] = 255;
								else if ((int)pix[(1200 * r) + c] + b >= 0)//combines the brightness
									pix[(1200 * r) + c] = ((int)pix[(1200 * r) + c] + b);
							}


						}
					}//( pow (((double)c - (((-ra[stars[i]] - (hour - .8))*(1200.0 / 1.6)) + 18000.0)), 2.0) + pow ((double)r - (((8.0) - dec[stars[i]])*(800.0 / 16.0)), 2.0))
					if (xPos > -18000.0 * 2)//for stars really out of pixel range
					{

						//cout << "passed first if " << ++count << endl;
						if ((c >= (xPos + 18000.0*2.0 - 14.0)) && (c <= (xPos + 18000.0*2.0 + 14.0)))
						{//            sqrt(         (x1   -            x2)                                                       ^2  +          ( y1  -                  y2                         )    ^2 )
							double d = sqrt((pow(((double)c - round(xPos + 18000.0*2.0)), 2.0) + pow((double)r - round(yPos), 2.0)));// calculation for distance from center of star
							double b = exp(-pow(d, 2.0) / (2.0*pow(((12.0 - MagI) / 2.354), 2.0))) *255.0;//calculation for brightness based of gausian curve
																										  //cout << "passed first if "<<++count << endl;
																										  //cout << "filling pix value" << endl;
							if (d < 13.0 - MagI)
							{
								if ((int)pix[(1200 * r) + c] + b + 1 >= 255)
									pix[(1200 * r) + c] = 255;
								else if ((int)pix[(1200 * r) + c] + b >= 0)
									pix[(1200 * r) + c] = ((int)pix[(1200 * r) + c] + b);
							}

						}
					}
					if ((c >= xPos - 14.0) && (c <= xPos + 14.0))
					{
						double d = sqrt((pow(((double)c - round(xPos)), 2.0) + pow((double)r - round(yPos), 2.0)));// calculation for distance from center of star
						double b = exp(-pow(d, 2.0) / (2.0*pow(((12.0 - MagI) / 2.354), 2.0))) *255.0;//calculation for brightness based of gausian curve
																									  //cout << "passed second if " << ++count << endl;
																									  //cout << "filling pix value" << endl;
						if (d < 13.0 - MagI)
						{
							if ((int)pix[(1200 * r) + c] + b + 1 >= 255)
								pix[(1200 * r) + c] = 255;
							else if ((int)pix[(1200 * r) + c] + b >= 0)
								pix[(1200 * r) + c] = ((int)pix[(1200 * r) + c] + b);
						}


					}
				}
			}
		}
	}

	delete[] RA;
	RA = NULL;
	delete[] DEC;
	DEC = NULL;
	delete[] MAG;
	MAG = NULL;
	//deletes and resets arrays for memory management
	/*for (int r = 0; r < 800; r++)
	{
	for (int c = 0; c < 1200; c++)
	{
	cout << pix[1200 * r + c] << " ";
	}
	cout << endl;
	}


	*/
	// Creates needed values to calculate the moon's offset in radians
	double moonr = 1737.0;
	double a = 181526;
	double b = 200000;
	double gamma = getGamma(x);
	double c = getC(a, b, gamma);
	double mangle = (180.0 / M_PI)*atan(moonr / (sqrt(((a + moonr) * (a + moonr)) + b * b - 2 * (a + moonr) * b * cos((M_PI / 180.0) * gamma))));
	double angle = (M_PI / 180.0)*getAhDegree(gamma, getAlpha(gamma, a, c));

	double mult = 1.0;
	double theta;
	if (angle < M_PI) theta = M_PI - angle;
	else if (angle >= M_PI) {
		theta = angle - M_PI;
		mult = -1;
	}


	double centdist = sqrt((c * c) + (moonr * moonr) - (2 * c * moonr * cos(theta)));
	double offset = -((asin((moonr * sin(theta)) / centdist)) * (180 / M_PI) * mult * 50); // offset in pixels

	double pixrad = mangle * 50; // Multiplies degrees of moon radius by pixels per degree
								 //cout << "Pixrad: " << pixrad << "\t\tOffset: " << offset << endl;
	for (int r = 0; r < 800; r++) {
		if (r < 400 + pixrad && r > 400 - pixrad) {
			for (int p = 0; p < 1200; p++) {
				if (p > 600 - offset - pixrad && p < 600 - offset + pixrad) {
					if (sqrt(pow(p - 600 + offset, 2.0) + pow(r - 400, 2.0)) < pixrad) {
						pix[(1200 * r) + p] = 255;
					}
				}
			}
		}
	}

}


double makeHour(int x)//distance from center of circle r=b=200000km, distance from center of orbit to moon:  a=181526km
{
	double gamma = getGamma(x);//gets the degree we're at in our orbit
	double ahDegree = getAhDegree(gamma, getAlpha(gamma, 181526.0, getC(181526.0, 200000.0, gamma)));//gets the degree were looking at
	double dpf = getDPF(x);//gets degrees per frame alpha (degree were looking at relative to gamma a b and c) changes
	double dph = 15.0;//degrees per hour
					  //cout << dpf << endl;

	double hour = -(((ahDegree)+dpf * (double)x) / dph) - ((dpf / dph)*(double)x);
	return hour;
}

double getGamma(int x)//gamma in terms of frame, const rate of change
{
	double gamma = -90.0 + ((double)x*(180.0 / (7200.0)));
	return gamma;
}

double getC(double a, double b, double gamma)//distance from focal point dependent on gamma
{
	double c;

	if (gamma != 0.0)
		c = sqrt(pow(a, 2.0) + pow(b, 2.0) - (2.0*a*b*(cos((M_PI / 180.0)*gamma))));
	else
		c = 18474.0;
	return c;
}

double getAlpha(double gamma, double a, double c)//where were looking relative and depending on gamma and c
{
	double alpha;
	if (c != 0.0)

		alpha = -180.0 / M_PI * asin(a*(sin(gamma*M_PI / 180.0)) / c);

	else
		alpha = 0.0;
	return alpha;
}

double getDPF(int x)//degrees per frame, really a polar equation but can be modeled with two seperate quadratics
{
	double dpf;
	double ahDegree = getAhDegree(getGamma(x), getAlpha(getGamma(x), 181526.0, getC(181526.0, 200000.0, getGamma(x))));
	//cout <<ahDegree<<" : " << getAlpha(getGamma(x), 181526.0, getC(181526.0, 200000.0, getGamma(x))) <<" : "<< getGamma(x)<<" : "<< getC(181526.0, 200000.0, getGamma(x))<<  endl;
	if (ahDegree >= 42.0&&ahDegree <= 180.0)
		dpf = (-0.001665 + 8.86*pow(10.0, -9.0)*(getC(181526.0, 200000.0, getGamma(x)) - 18474.0)) / 7200.0;
	else
		dpf = (0.001665 - 8.86*pow(10.0, -9.0)*(getC(181526.0, 200000.0, getGamma(x)) - 18474.0)) / 7200.0;
	return dpf;
}

double getAhDegree(double gamma, double alpha)//gets degree we're looking at relative to RA, dependent on gamma, equation changes at halfway point
{
	double AHD;
	if (alpha >= 0.0)
		AHD = (90.0 + ((90.0 - fabs(gamma)) - alpha));
	else
		AHD = (270.0 - ((90.0 - gamma) + alpha));
	return AHD;
}





