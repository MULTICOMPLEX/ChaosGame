
#include <iostream>
#include <numbers>
#include <vector>
#include "mxws.hpp"
#include "surface.hpp"

constexpr auto WIDTH = 1000;
constexpr auto HEIGHT = 900;

mxws<uint32_t> RNG;

std::vector<Pixel> pix(3);

constexpr auto ITER = 1000000;

int main(int argc, char* argv[])
{

	Surface surf(WIDTH, HEIGHT);
	Color color;

	color.r = 100, color.g = 100, color.b = 255;

	Pixel last = {WIDTH, HEIGHT};
	Pixel point = {};

	pix[0].x = WIDTH / 2;
	pix[0].y = 0;
	pix[1].x = 0;
	pix[1].y = HEIGHT;
	pix[2].x = WIDTH;
	pix[2].y = HEIGHT;

	for (auto i = 0; i < ITER; i++) {

		surf.RectFill(point.x, point.y, 1, 1, color);

		int choice = RNG(2);
		auto dx = pix[choice].x - last.x;
		auto dy = pix[choice].y - last.y;

		point.x = int(round(last.x + dx / 2.));
		point.y = int(round(last.y + dy / 2.));

		last = point;

	}

	surf.Save("Chaos_Game.bmp");

	return 0;
}



