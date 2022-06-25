#pragma once

#include <iostream>
#include <numbers>
#include <vector>
#include <assert.h>
#include "mxws.hpp"
#include "surface.hpp"
#include "lodepng.h"
#include "buddhabrot.h"

const int MAXITER = 1100;
mxws<uint32_t> rng;

bool IsMandelbrot(double c_r, double c_i) {
  double z_r = 0, z_i = 0, z_temp;
  double q = pow(c_r - .25, 2.0) + pow(c_i, 2.0);

  //Checks if point is in main cardioid
  if (q * (q + (c_r - .25)) < .25 * pow(c_i, 2.0)) {
    return true;
  }

  //Check if point is in period 2 bulb
  if (pow(c_r + 1, 2.0) + pow(c_i, 2.0) < .0625) {
    return true;
  }

  for (int iter = 0; iter < MAXITER; iter++) {
    z_temp = pow(z_r, 2.0) - pow(z_i, 2.0) + c_r;
    z_i = 2 * z_r * z_i + c_i;
    z_r = z_temp;

    if (pow(z_r, 2.0) + pow(z_i, 2.0) > 4) {
      return false;
    }
  }
  return true;
}

unsigned int Buddhabrot::get_width() const
{
  return m_width;
}

unsigned int Buddhabrot::get_height() const
{
  return m_height;
}

void Buddhabrot::gen_fractal()
{
  const int width = get_width(), height = get_height();
  const int NUM_PIXELS = width * height;
  std::vector<double> buckets(NUM_PIXELS);
  int max_value = 0;

  for (int i = 0; i < NUM_PIXELS * 15; i++) {
    std::vector<double> buckets_temp(NUM_PIXELS);
    double z_r = 0, z_i = 0;
    double z_temp;
    double c_r = 0, c_i = 0;
    int iter = 0;

    do {
      c_r = rng(-2.5, 1.); // Real (-2.5, 1)
      c_i = rng(-1., 1.); // Imaginary (-1, 1)
    } while (IsMandelbrot(c_r, c_i));

    for (int j = 0; j < NUM_PIXELS; j++) {
      buckets_temp[j] = buckets[j];
    }

    for (iter; iter < MAXITER; iter++) {
      z_temp = pow(z_r, 2.0) - pow(z_i, 2.0) + c_r;
      z_i = 2 * z_r * z_i + c_i;
      z_r = z_temp;

      if (pow(z_r, 2.0) + pow(z_i, 2.0) > 4) {
        for (int i = 0; i < NUM_PIXELS; i++) {
          buckets[i] = buckets_temp[i];
        }
        break;
      }

      int x = (int)(((z_r + 1.6) / 2.5) * (width));
      int y = (int)(((z_i + 1.03) / 2.05) * (height));

      if (0 <= x && x < width && 0 <= y && y < height) {
        buckets_temp[height * static_cast<std::vector<double, std::allocator<double>>::size_type>(x) + y]++;
      }
    }

  }

  for (int i = 0; i < NUM_PIXELS; i++) {
    if (buckets[i] > max_value) {
      max_value = int(buckets[i]);
    }
  }

  for (int i = 0; i < NUM_PIXELS; i++) {
    buckets[i] /= max_value;
    image[i * 4ull] = int(buckets[i] * 255);
    image[i * 4ull + 1] = int(buckets[i] * 51);
    image[i * 4ull + 2] = int(buckets[i] * 51);
    image[i * 4ull + 3] = 255;
  }
}


namespace {

  using Real = long double;
  using ImageReal = float;

  struct TargetProperties {
    Real x;
    Real y;
    Real zoom;
  };
  struct BitmapProperties {
    int height;
    int width;
  };

  constexpr TargetProperties target1 = { -0.4, 0, 0.32 };  // the full set
  constexpr TargetProperties target2 = { -1.25275, -0.343, 250 };
  constexpr TargetProperties target3 = { -0.1592, -1.0317, 80.5 };
  constexpr TargetProperties target4 = { -0.529854097, -0.667968575, 80.5 };
  constexpr TargetProperties target5 = { -0.657560793, 0.467732884, 70.5 };
  constexpr TargetProperties target6 = { -1.185768799, 0.302592593, 90.5 };
  constexpr TargetProperties target7 = { 0.443108035, 0.345012263, 4000 };
  constexpr TargetProperties target8 = { -0.647663050, 0.380700836, 1275 };
  constexpr TargetProperties target9 = { -0.0443594, -0.986749, 88.2 };  // steckles.com sample image

  constexpr TargetProperties target = target1;
  constexpr BitmapProperties bitmap = { 4096, 4096 };

  //convert color.ppm -normalize -level 0,90% -gamma 0.3 -scale 800x800 nebulabrot.jpg
  constexpr int max_orbit_length = 50000;

  mxws<uint32_t> prng;

  struct Complex {
    Real a = 0;
    Real b = 0;
    Real c = 0;
    Real d = 0;

    explicit Complex() = default;
    constexpr explicit Complex(Real a, Real b, Real c, Real d) : a(a), b(b), c(c), d(d) {}

    constexpr Real mag2() const { return c * c + d * d; }
    Real dist2(Real cp, Real dp) const { return (c - cp) * (c - cp) + (d - dp) * (d - dp); }
  };

  struct ImageCoord {
    int y;
    int x;
    constexpr explicit ImageCoord(int x, int y) : y(y), x(x) {}
    bool isOnScreen() const { return (0 <= y && y < bitmap.height && 0 <= x && x < bitmap.width); }
  };

  ImageCoord Project(const Complex& s)
  {
    // Notice that the image is rotated 90 degrees compared to the original complex plane.
    // Negative real components (target.x) turn into upward displacements (small coord.y).
    return ImageCoord(
      int((s.b - target.y) * target.zoom * bitmap.height + (bitmap.width / 2)),
      int((s.a - target.x) * target.zoom * bitmap.height + (bitmap.height / 2))
    );
  }

  constexpr Complex Unproject(ImageCoord p)
  {
    return Complex(
      0, 0,
      (p.y - (bitmap.height / 2)) / (target.zoom * bitmap.height) + target.x,
      (p.x - (bitmap.width / 2)) / (target.zoom * bitmap.height) + target.y
    );
  }

  Complex orbit[max_orbit_length];
  int olen = 0;

  ImageReal(&image_buffer)[3][bitmap.height][bitmap.width] = *new ImageReal[1][3][bitmap.height][bitmap.width]{};

  Complex RandomInRadiusAround(Real x, Real y, Real r)
  {
    const auto r2 = r * r;
    while (true) {
      Complex c = Complex(0, 0, prng(-r, r), prng(-r, r));
      if (c.mag2() < r2) {
        c.c += x;
        c.d += y;
        return c;
      }
    }
  }

  // This function is essentially gamma-correction at both ends of the range [0,1].
  // It applies the inverse of a sigmoid function to artificially increase the
  // brightness of small-but-not-zero values and decrease the brightness of
  // large-but-not-1 values, thus compressing the output into the middle of the range.
  // However, this still tends to make things too dark, so after that,
  // we multiply everything by two (saturating at 1.0).
  // Then, we multiply by 256 to get a final pixel value.
  //
  unsigned char GainCorrect(ImageReal x)
  {
    auto bend = [](ImageReal x) {
      assert(0 <= x && x <= 0.5);
      ImageReal y = ImageReal(pow(2 * x, 0.3219280948873623) / 2);
      assert(0 <= y && y <= 0.5);
      return y;
    };
    ImageReal gain_adjusted = (x <= 0.5) ? bend(x) : (1 - bend(1 - x));
    ImageReal brightened = gain_adjusted * 2;  // This factor seems unexplained.
    auto m = brightened * 256;
    return unsigned char((m < 0) ? 0 : (m < 255) ? m : 255);
  }

  void OutputPPM(const char* fname, ImageReal max0, ImageReal max1, ImageReal max2)
  {
    FILE* fp = fopen(fname, "wb");
    assert(fp != nullptr);
    fprintf(fp, "P6\n%d %d\n255\n", bitmap.width, bitmap.height);
    for (int y = 0; y < bitmap.height; ++y) {
      for (int x = 0; x < bitmap.width; ++x) {
        unsigned char rgb[3] = {
            GainCorrect(image_buffer[0][y][x] / max0),
            GainCorrect(image_buffer[1][y][x] / max1),
            GainCorrect(image_buffer[2][y][x] / max2),
        };
        fwrite(rgb, 1, 3, fp);
      }
    }
    fclose(fp);
  }

  void OutputPGM(const char* fname, int c, ImageReal cmax)
  {
    FILE* fp = fopen(fname, "wb");
    assert(fp != nullptr);
    fprintf(fp, "P5\n%d %d\n255\n", bitmap.width, bitmap.height);
    for (int y = 0; y < bitmap.height; ++y) {
      for (int x = 0; x < bitmap.width; ++x) {
        unsigned char rgb[1] = {
            GainCorrect(image_buffer[c][y][x] / cmax),
        };
        fwrite(rgb, 1, 1, fp);
      }
    }
    fclose(fp);
  }

  void DrawBuffer()
  {
    printf("DrawBuffer!\n");
    ImageReal max0 = 0.01f;  // avoid div-by-zero on a completely black image
    ImageReal max1 = 0.01f;
    ImageReal max2 = 0.01f;
    for (int y = 0; y < bitmap.height; ++y) {
      for (int x = 0; x < bitmap.width; ++x) {
        max0 = std::max(max0, image_buffer[0][y][x]);
        max1 = std::max(max1, image_buffer[1][y][x]);
        max2 = std::max(max2, image_buffer[2][y][x]);
      }
    }

    OutputPPM("color.ppm", max0, max1, max2);
    OutputPGM("red.pgm", 0, max0);
    OutputPGM("green.pgm", 1, max1);
    OutputPGM("blue.pgm", 2, max2);
  }

  constexpr Real computeEscapeMagnitude2()
  {
    Real m1 = Unproject(ImageCoord{ 0, 0 }).mag2();
    Real m2 = Unproject(ImageCoord{ bitmap.width - 1, 0 }).mag2();
    Real m3 = Unproject(ImageCoord{ 0, bitmap.height - 1 }).mag2();
    Real m4 = Unproject(ImageCoord{ bitmap.width - 1, bitmap.height - 1 }).mag2();
    Real result = 4.0;
    if (result < m1) result = m1;
    if (result < m2) result = m2;
    if (result < m3) result = m3;
    if (result < m4) result = m4;
    return result;
  }

  Real computeEscapeMagnitude()
  {
    return sqrt(computeEscapeMagnitude2());
  }

  bool Evaluate(const Complex& s, int it)
  {
    constexpr Real escape = computeEscapeMagnitude2();

    assert(it <= max_orbit_length);
    olen = 0;
    Real a = s.a;
    Real b = s.b;
    for (int i = 0; i < it; ++i) {
      Real t = a * a - b * b + s.c;
      b = 2 * a * b + s.d;
      a = t;
      if (a * a + b * b > escape) {
        return true;
      }
      orbit[olen++] = Complex(a, b, s.c, s.d);
    }
    return false;
  }

  void DrawOrbit(int c, const Complex* d_orbit, int len)
  {
    for (int i = 0; i < len; ++i) {
      ImageCoord p = Project(d_orbit[i]);
      if (p.isOnScreen()) {
        image_buffer[c][p.y][p.x] += ImageReal(1);
        if (target.y == 0) {
          // If the region being graphed is symmetrical, then make the image symmetrical too.
          image_buffer[c][p.y][bitmap.width - p.x - 1] += ImageReal(1);
        }
      }
    }
  }

  Real Contrib()
  {
    Real contrib = 0;
    for (int i = 0; i < olen; ++i) {
      ImageCoord p = Project(orbit[i]);
      if (p.isOnScreen()) {
        contrib += 1;
      }
    }
    return contrib / Real(olen);
  }

  Real TransitionProbability(Real q1, Real q2, Real olen1, Real olen2)
  {
    return (olen1 * q2) / (olen2 * q1);
  }

  // Most mutations should be small,
  // but occasionally throw in a completely new point
  // to maintain ergodicity.
  //
  Complex Mutate(const Complex& c)
  {
    static int counter = 5;
    if (--counter != 0) {
      Complex n = c;
      Real r2 = (1.L / target.zoom) * 0.1;
      Real phi = prng(2 * std::numbers::pi);
      Real r = r2 * exp(-prng(0, 6.9077));
      n.c += r * cos(phi);
      n.d += r * sin(phi);
      return n;
    }
    else {
      counter = 5;
      return RandomInRadiusAround(0, 0, computeEscapeMagnitude());
    }
  }

  /***********************************************************************/
  /*                                                                     */
  /* Recursively find a point whose orbit passes through the screen.     */
  /* Not required per-se, but it's a lot faster than just trying random  */
  /* samples over and over again until you find a good one, especially   */
  /* at higher zooms.                                                    */
  /*                                                                     */
  /***********************************************************************/

  bool FindInitialSample(Complex& c, Real x, Real y, Real rad, int f)
  {
    if (f > 500) {
      return false;
    }
    Complex seed;
    int m = -1;
    Real closest = 1e20;
    for (int i = 0; i < 200; ++i) {
      Complex tmp = RandomInRadiusAround(x, y, rad);
      if (Evaluate(tmp, 50000)) {
        if (Contrib() > 0) {
          c = tmp;
          return true;
        }
        for (int q = 0; q < olen; ++q) {
          Real d = orbit[q].dist2(target.x, target.y);
          if (d < closest) {
            m = q;
            closest = d;
            seed = tmp;
          }
        }
      }
    }
    return FindInitialSample(c, seed.c, seed.d, rad / 2, f + 1);
  }

  void RenderBuddhabrotClassic()
  {
    static constexpr int Iterations[] = { 5000, 500, 50 };  // TODO
    unsigned short m = 0;
    while (true) {
      Complex z = RandomInRadiusAround(0, 0, computeEscapeMagnitude());
      if (Evaluate(z, 50000)) {
        for (int c = 0; c < 3; ++c) {
          if (olen <= Iterations[c]) {
            DrawOrbit(c, orbit, olen);
            if (++m == 0) {
              DrawBuffer();
            }
          }
        }
      }
    }
  }

  /***********************************************************************/
  /*                                                                     */
  /* As many disparate parts of the Mandelbrot set can produce orbits    */
  /* that pass through the screen, an optimization is to run multiple    */
  /* copies of the algorithm from different starting values; that way    */
  /* you're less likely to miss important orbits early on.               */
  /*                                                                     */
  /***********************************************************************/

  constexpr int metro_threads = 30;

  Complex z_samples[metro_threads];
  Real c_samples[metro_threads];
  Real l[metro_threads][3];
  Real o[metro_threads][3];

  void BuildInitialSamplePoints()
  {
    for (int i = 0; i < metro_threads; ++i) {
      l[i][0] = o[i][0] = l[i][1] = o[i][1] = l[i][2] = o[i][2] = 1;

      Complex m;
      while (!FindInitialSample(m, 0, 0, 2, 0)) {
        printf("Couldn't find seed %d!\n", i);
      }
      Evaluate(m, 50000);
      z_samples[i] = m;
      c_samples[i] = Contrib();
    }
  }

  void RenderBuddhabrot()
  {
    static constexpr int Iterations[] = { 50000, 5000, 500 };
    unsigned short m = 0;
    while (true) {
      for (int ss = 0; ss < metro_threads; ++ss) {
        Complex n_sample = Mutate(z_samples[ss]);
        if (Evaluate(n_sample, 50000)) {  // it does eventually escape
          Real n_contrib = Contrib();
          if (n_contrib != 0) {
            for (int c = 0; c < 3; ++c) {
              int it_length = Iterations[c];
              if (olen <= it_length) {  // it escapes within it_length iterations
                Real T1 = TransitionProbability(it_length, l[ss][c], olen, o[ss][c]);
                Real T2 = TransitionProbability(l[ss][c], it_length, o[ss][c], olen);
                Real alpha = (n_contrib * T1) / (c_samples[ss] * T2);
                if (alpha > prng(1.)) {
                  c_samples[ss] = n_contrib;
                  z_samples[ss] = n_sample;
                  l[ss][c] = it_length;
                  o[ss][c] = olen;
                  DrawOrbit(c, orbit, olen);
                  if (++m == 0) {
                    DrawBuffer();
                  }
                }
              }
            }
          }
        }
      }
    }
  }

} // anonymous namespace

void Test_Buddhabrot()
{
  BuildInitialSamplePoints();
  RenderBuddhabrotClassic();
}