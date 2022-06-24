#pragma once

#include <fstream>

struct Color {
	union {
		struct { unsigned char b, g, r, a; };
		unsigned char bytes[4];
		unsigned int ref;
	};
	Color(unsigned int ref = 0) : ref(ref) {}
};

class Surface {
	int width, height;
	std::vector<Color> pixels;
public:
	void Fill(Color color) { std::fill(pixels.begin(), pixels.end(), color); }
	void HLine(int left, int y, int len, Color color) {
		if (y < 0 || y >= height || left >= width) { return; }
		if (left < 0) { len += left; left = 0; }
		if (left + len > width) { len = width - left; }
		int offset = y * width + left;
		std::fill(pixels.begin() + offset, pixels.begin() + offset + len, color);
	}
	void RectFill(int x, int y, int w, int h, Color color) {
		for (int i = 0; i < h; ++i) { HLine(x, y + i, w, color); }
	}
	Surface(int width, int height) :
		width(width),
		height(height),
		pixels(width* height, Color())
	{}
	template <typename I>       Color& operator () (const I& x, const I& y) { return pixels[y * width + x]; }
	template <typename I> const Color& operator () (const I& x, const I& y) const { return pixels[y * width + x]; }

	class Writer {
		std::ofstream ofs;
	public:
		Writer(const char* filename) : ofs(filename, std::ios_base::out | std::ios_base::binary) {}
		void operator () (const void* pbuf, int size) { ofs.write(static_cast<const char*>(pbuf), size); }
		template <typename T> void operator () (const T& obj) { operator () (&obj, sizeof(obj)); }
	};

	struct BIH {
		unsigned int   sz;
		int            width, height;
		unsigned short planes;
		short          bits;
		unsigned int   compress, szimage;
		int            xppm, yppm;
		unsigned int   clrused, clrimp;
	};

	void Save(const char* filename) const {
		Writer w(filename);;
		w("BM", 2);
		BIH bih = { sizeof(bih) };
		bih.width = width;
		bih.height = -height;
		bih.planes = 1;
		bih.bits = 32;
		const unsigned int headersize = sizeof(bih) + 14;
		const int szbuf = int(sizeof(Color) * pixels.size());
		const unsigned int filesize = static_cast<unsigned int>(headersize + szbuf);
		w(filesize);
		const unsigned short z = 0;
		w(z);
		w(z);
		w(headersize);
		w(bih);
		w(pixels.data(), szbuf);
	}
};

typedef struct Pixel
{
	int x;
	int y;
} Pixel;
