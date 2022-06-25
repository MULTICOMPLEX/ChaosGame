
#include "nebulabrot.hpp"

constexpr auto WIDTH = 1000;//1000;
constexpr auto HEIGHT = 900;//900;

mxws<uint32_t> RNG;

std::vector<Pixel> pix;

template <typename T>
std::vector<Pixel> polygonPoints(int numPoints, const T& WIDTH, const T& HEIGHT) {
	std::vector<Pixel> coordinates;
	int radius = HEIGHT / 2;
	double theta = -1.5 * (((numPoints - 2) * std::numbers::pi) / numPoints);
	int centerX = WIDTH / 2;
	int centerY = HEIGHT / 2;
	for (int i = 0; i < numPoints; i++)
	{
		auto X = radius * cos(2 * std::numbers::pi * i / numPoints + theta) + centerX;
		auto Y = radius * sin(2 * std::numbers::pi * i / numPoints + theta) + centerY;
		coordinates.push_back({ int(X),int(Y) });
	}
	return coordinates;
}

template <typename T>
std::vector<Pixel> setcorners_triangle(const T& WIDTH, const T& HEIGHT) {
	std::vector<Pixel> pix(3);
	pix[0].x = WIDTH / 2;
	pix[0].y = 0;
	pix[1].x = 0;
	pix[1].y = HEIGHT;
	pix[2].x = WIDTH;
	pix[2].y = HEIGHT;
	return pix;
}

template <typename T>
std::vector<Pixel> setcorners_square(const T& WIDTH, const T& HEIGHT) {
	std::vector<Pixel> pix(4);
	pix[0].x = 0;
	pix[0].y = 0;
	pix[1].x = WIDTH;
	pix[1].y = 0;
	pix[2].x = 0;
	pix[2].y = HEIGHT;
	pix[3].x = WIDTH;
	pix[3].y = HEIGHT;
	return pix;
}

template <typename T>
std::vector<Pixel> Vicsek_fractal(const T& WIDTH, const T& HEIGHT) {
	std::vector<Pixel> pix(5);
	pix[0].x = 0;
	pix[0].y = 0;
	pix[1].x = WIDTH;
	pix[1].y = 0;
	pix[2].x = 0;
	pix[2].y = HEIGHT;
	pix[3].x = WIDTH;
	pix[3].y = HEIGHT;
	pix[4].x = WIDTH / 2;
	pix[4].y = HEIGHT / 2;
	return pix;
}

template <typename T>
std::vector<Pixel> setcorners_ngon(const T& vertices, const T& WIDTH, const T& HEIGHT) {

	return polygonPoints(vertices, WIDTH, HEIGHT);
}

constexpr auto ITER = 1000000;

std::vector<unsigned char> image((WIDTH + 1)* (HEIGHT + 1) * 4);
template <typename im, typename I, typename col>
void putpixel(im& image, const I& x, const I& y, const col& color);
void encodeOneStep(const char* filename, const std::vector<unsigned char>& image,
	unsigned width, unsigned height);

int main(int argc, char* argv[])
{
	//Buddhabrot buddhabrot(WIDTH, HEIGHT);
	//buddhabrot.gen_fractal();
	//encodeOneStep("Buddhabrot.png", buddhabrot.image, WIDTH, HEIGHT);
	//Test_Buddhabrot();
	//return 0;

	Surface surf(WIDTH, HEIGHT);
	Color color;

	color.r = 100, color.g = 100, color.b = 255, color.a = 255;

	Pixel last = { WIDTH, HEIGHT };
	Pixel point = {};

	const int fractal_select = 2;

	if (fractal_select == 0)
		pix = setcorners_triangle(WIDTH, HEIGHT);
	else if (fractal_select == 1)
		pix = setcorners_square(WIDTH, HEIGHT);
	else if (fractal_select == 2)
		pix = setcorners_ngon(5, WIDTH, HEIGHT);
	else if (fractal_select == 3)
		pix = Vicsek_fractal(WIDTH, HEIGHT);

	int choice = 0, choice2 = 0;

	for (auto i = 0; i < ITER; i++) {

		surf.RectFill(point.x, point.y, 1, 1, color);
		putpixel(image, point.x, point.y, color);

		if (fractal_select == 0)
			choice = RNG(2);

		else if (fractal_select == 3) 
			choice = RNG(int(pix.size() - 1));

		else {
			while (choice2 == choice)
				choice = RNG(int(pix.size() - 1));
		}

		choice2 = choice;

		auto dx = pix[choice].x - last.x;
		auto dy = pix[choice].y - last.y;
		if (fractal_select == 3) {
			point.x = int(round(last.x + dx * (2./3)));
			point.y = int(round(last.y + dy * (2./3)));
		}
		else 
		{
			point.x = int(round(last.x + dx / 2));
			point.y = int(round(last.y + dy / 2));
		}
		last = point;

	}

	surf.Save("Chaos_Game.bmp");
	encodeOneStep("Chaos_Game.png", image, WIDTH, HEIGHT);

	return 0;
}

void encodeOneStep(const char* filename, const std::vector<unsigned char>& image,
	unsigned width, unsigned height) {
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

template <typename im, typename I, typename col>
void putpixel(im& image, const I& x, const I& y, const col& color)
{
	image[4ull * WIDTH * y + 4ull * x + 0] = color.r;
	image[4ull * WIDTH * y + 4ull * x + 1] = color.g;
	image[4ull * WIDTH * y + 4ull * x + 2] = color.b;
	image[4ull * WIDTH * y + 4ull * x + 3] = color.a;
}
