
#include <iostream>
#include <numbers>
#include <vector>
#include "mxws.hpp"
#include "surface.hpp"
#include "lodepng.h"

constexpr auto WIDTH = 1000;
constexpr auto HEIGHT = 900;

mxws<uint32_t> RNG;

std::vector<Pixel> pix(3);

constexpr auto ITER = 1000000;

std::vector<unsigned char> image((WIDTH+4)* (HEIGHT+4) * 4);
template <typename im, typename I, typename col>
void putpixel(im& image, const I& x, const I& y, const col& color);
void encodeOneStep(const char* filename, const std::vector<unsigned char>& image, 
	unsigned width, unsigned height);

int main(int argc, char* argv[])
{

	Surface surf(WIDTH, HEIGHT);
	Color color;

	color.r = 100, color.g = 100, color.b = 255, color.a = 255;

	Pixel last = { WIDTH, HEIGHT };
	Pixel point = {};

	pix[0].x = WIDTH / 2;
	pix[0].y = 0;
	pix[1].x = 0;
	pix[1].y = HEIGHT;
	pix[2].x = WIDTH;
	pix[2].y = HEIGHT;

	for (auto i = 0; i < ITER; i++) {

		surf.RectFill(point.x, point.y, 1, 1, color);
		putpixel(image, point.x, point.y, color);

		int choice = RNG(2);
		auto dx = pix[choice].x - last.x;
		auto dy = pix[choice].y - last.y;

		point.x = int(round(last.x + dx / 2.));
		point.y = int(round(last.y + dy / 2.));

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