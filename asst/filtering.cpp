/* -----------------------------------------------------------------
 * File:    filtering.cpp
 * Created: 2015-09-22
 * -----------------------------------------------------------------
 *
 * Image convolution and filtering
 *
 * ---------------------------------------------------------------*/

#include "filtering.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include "basicImageManipulation.h"
using namespace std;

Image boxBlur(const Image &im, int k, bool clamp) { // Safe to assume k is odd
  // --------- HANDOUT  PS02 ------------------------------
  Image output(im.width(), im.height(), im.channels());
  for (int h = 0; h < im.height(); h++) { // Iterate all possible pixels given height and width
    for (int w = 0; w < im.width(); w++) {
      for (int c = 0; c < im.channels(); c++) {
  
      float sum = 0.0f; // Initialize sum for averaging
      for (int h0 = h - (k - 1)/2; h0 < h - (k - 1)/2 + k; h0++) {
        for (int w0 = w - (k - 1)/2; w0 < w - (k - 1)/2 + k; w0++) { // Iterate over kernel-contained coordinates
          sum += im.smartAccessor(w0, h0, c, clamp);
        }
      }
      output(w, h, c) = sum / (k * k); // Put averages in output pixels by dividing by k^2
      }
    }
  }
  return output; // Return output image
}

Image Filter::convolve(const Image &im, bool clamp) const {
  // --------- HANDOUT  PS02 ------------------------------
  // Write a convolution function for the filter class

  Image output(im.width(), im.height(), im.channels()); // Initialize output image of same size
  for (int h = 0; h < im.height(); h++) { // Iterate pixels in row-major iteration order
    for (int w = 0; w < im.width(); w++) {
      for (int c = 0; c < im.channels(); c++) {
  
      float sum = 0.0f; // Initialize sum for averaging
      int f_h = 0;      // Initialize index for kernel y value
      for (int h0 = h - (height() - 1)/2; h0 < h - (height() - 1)/2 + height(); h0++) { // Row-major local iteration
        int f_w = 0;    // Index for kernel x value
        for (int w0 = w - (width() - 1)/2; w0 < w - (width() - 1)/2 + width(); w0++) {
          sum += im.smartAccessor(w0, h0, c, clamp) * (*this)(f_w, f_h);
          f_w += 1; // Increment kernel x value
        }
        f_h += 1; // Increment kernel y value
      }

      output(w, h, c) = sum; // Put convolution sum into output image
      }
    }
  }
  return output; // change this
}

Image boxBlur_filterClass(const Image &im, int k, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Reimplement the box filter using the filter class.
  // check that your results match those in the previous function "boxBlur"
  
  vector<float> kernel;
	for (int n = 0; n < k * k; n++) {
		kernel.push_back(1.0f / (float(k * k))); // Generate averaging kernel of k * k size
	}

	Filter blur(kernel, k, k);       // Initialize filter then return convolution
	return blur.convolve(im, clamp); 
}

Image gradientMagnitude(const Image &im, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Uses a Sobel kernel to compute the horizontal and vertical components
  // of the gradient of an image and returns the gradient magnitude.
  return im; // change this
}

vector<float> gauss1DFilterValues(float sigma, float truncate) {
  // --------- HANDOUT  PS02 ------------------------------
  // Create a vector containing the normalized values in a 1D Gaussian filter
  // Truncate the gaussian at truncate*sigma.
  return vector<float>();
}

Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate,
                              bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Gaussian blur across the rows of an image
  return im;
}

vector<float> gauss2DFilterValues(float sigma, float truncate) {
  // --------- HANDOUT  PS02 ------------------------------
  // Create a vector containing the normalized values in a 2D Gaussian
  // filter. Truncate the gaussian at truncate*sigma.
  return vector<float>();
}

Image gaussianBlur_2D(const Image &im, float sigma, float truncate,
                      bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Blur an image with a full 2D rotationally symmetric Gaussian kernel
  return im;
}

Image gaussianBlur_separable(const Image &im, float sigma, float truncate,
                             bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Use principles of seperabiltity to blur an image using 2 1D Gaussian
  // Filters
  return im;
}

Image unsharpMask(const Image &im, float sigma, float truncate, float strength,
                  bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Sharpen an image
  return im;
}

Image bilateral(const Image &im, float sigmaRange, float sigmaDomain,
                float truncateDomain, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Denoise an image using the bilateral filter
  return im;
}

Image bilaYUV(const Image &im, float sigmaRange, float sigmaY, float sigmaUV,
              float truncateDomain, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // 6.865 only
  // Bilateral Filter an image seperatly for
  // the Y and UV components of an image
  return im;
}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
Image impulseImg(int k) {
  // initlize a kxkx1 image of all 0's
  Image impulse(k, k, 1);

  // set the center pixel to have intensity 1
  int center = floor(k / 2);
  impulse(center, center, 0) = 1.0f;

  return impulse;
}

// ------------- FILTER CLASS -----------------------
Filter::Filter() : width_(0), height_(0) {}
Filter::Filter(const vector<float> &fData, int fWidth, int fHeight)
    : kernel_(fData), width_(fWidth), height_(fHeight) {
  assert(fWidth * fHeight == (int)fData.size());
}

Filter::Filter(int fWidth, int fHeight)
    : kernel_(std::vector<float>(fWidth * fHeight, 0)), width_(fWidth),
      height_(fHeight) {}

const float &Filter::operator()(int x, int y) const {
  if (x < 0 || x >= width_)
    throw OutOfBoundsException();
  if (y < 0 || y >= height_)
    throw OutOfBoundsException();

  return kernel_[x + y * width_];
}

float &Filter::operator()(int x, int y) {
  if (x < 0 || x >= width_)
    throw OutOfBoundsException();
  if (y < 0 || y >= height_)
    throw OutOfBoundsException();

  return kernel_[x + y * width_];
}
// --------- END FILTER CLASS -----------------------
