/*
This program extracts data from level-3 HDF bin files and
writes them out as an ASCII table.

Regions of interest can be specified by latitude and longitude
boundaries or by a central coordinate and a radius in kilometers.

Norman Kuring		14-Feb-2003	Original development
Norman Kuring		11-Dec-2007	Fix memory-overrun bug and add a
					couple of calls to free().
Norman Kuring		21-Dec-2011	Give a precision when printing out
					bit strings to avoid unwanted printing
					of uninitialized memory.
Norman Kuring		21-Mar-2013	Change the latbin array from 32-bit
					floats to 64-bit to reduce rounding-
					error effects on bin numbers at smaller
					bin sizes.  Thanks to George White for
					pointing this one out.
*/

/*
 The following functions are based on the pseudocode found in Appendix A of:

 Campbell, J.W., J.M. Blaisdell, and M. Darzi, 1995:
 Level-3 SeaWiFS Data Products: Spatial and Temporal Binning Algorithms.
 NASA Tech. Memo. 104566, Vol. 32,
 S.B. Hooker, E.R. Firestone, and J.G. Acker, Eds.,
 NASA Goddard Space Flight Center, Greenbelt, Maryland
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

#include <Rcpp.h>
using namespace Rcpp;

#define NUMROWS		4320 //2160
#define EARTH_RADIUS	6371.229

static int	basebin[NUMROWS];
static int	numbin[NUMROWS];
static float	latbin[NUMROWS];
static int	totbins;

double	constrain_lat(double lat);
double	constrain_lon(double lon);

void initbin(void){
  int	row;

  basebin[0] = 1;

  for(row = 0; row < NUMROWS; row++) {

    latbin[row] = ((row + 0.5) * 180.0 / NUMROWS) - 90.0;

    numbin[row] = (int)(2.0 * NUMROWS * cos(latbin[row] * M_PI / 180.0) + 0.5);

    if(row > 0) {
      basebin[row] = basebin[row - 1] + numbin[row - 1];
    }
  }

  totbins = basebin[NUMROWS - 1] + numbin[NUMROWS - 1] - 1;
}

//’ The length of a string (in characters).
//’
//’ @param str input character vector
//’ @return characters in each element of the vector
// [[Rcpp::export]]
int lat2row(double lat){
  int	row;

  row = (int)((90 + lat)*NUMROWS/180.0);
  if(row >= NUMROWS) row = NUMROWS - 1;
  return(row);
}

// [[Rcpp::export]]
int rowlon2bin(int row, double lon){
  int	col;
  int	bin;

  lon = constrain_lon(lon);
  col = (int) ((lon + 180.0) * numbin[row] / 360.0);
  if(col >= numbin[row]) col = numbin[row] - 1;
  bin = basebin[row] + col;
  return(bin);
}

// [[Rcpp::export]]
Rcpp::DataFrame latlon2bin(NumericVector lat, NumericVector lon){

  int row, col;
  int n = lat.size();
  NumericVector bin(n);

  for(int i = 0; i < n; i++) {

    /* Constrain latitudes to [-90,90] and longitudes to [-180,180]. */
    lat[i] = constrain_lat(lat[i]);
    lon[i] = constrain_lon(lon[i]);

    row = lat2row(lat[i]);
    col = (int)((lon[i] + 180.0)*numbin[row]/360.0);
    if(col >= numbin[row]) col = numbin[row] - 1;
    bin[i] = basebin[row] + col;

  }

  Rcpp::DataFrame res = Rcpp::DataFrame::create(
    Rcpp::Named("bin_index") = bin);

  return(res);
}

int getRowIndex(double lat) {

  int rowIndex = ((lat + 90) * NUMROWS / 180);
  rowIndex = min(rowIndex, NUMROWS - 1);

  return rowIndex;

}

// [[Rcpp::export]]
Rcpp::DataFrame bin2latlon(NumericVector binIndex){

  initbin();

  int	row;

  int n = binIndex.size();

  NumericVector lat(n);
  NumericVector lon(n);
  // NumericVector ri(n);
  // NumericVector bi(n);

  // Rcpp::Rcout << binIndex;

  for(int i = 0; i < n; i++) {

    int	row;

    row = NUMROWS - 1;

    if(binIndex[i] < 1) binIndex[i] = 1;

    while(binIndex[i] < basebin[row]) row--;

    lat[i] = latbin[row];
    lon[i] = 360.0*(binIndex[i] - basebin[row] + 0.5)/numbin[row] - 180.0;

  }

  Rcpp::DataFrame df = Rcpp::DataFrame::create(
    Rcpp::Named("longitude") = lon,
    Rcpp::Named("latitude") = lat);

  return(df);

}

double constrain_lat(double lat){
  if(lat >  90) lat =  90;
  if(lat < -90) lat = -90;
  return(lat);
}

double constrain_lon(double lon){
  while(lon < -180) lon += 360;
  while(lon >  180) lon -= 360;
  return(lon);
}

/*** R
lonlat <- tibble::tribble(
  ~bin_num, ~nobs, ~nscenes, ~weights, ~time_rec,
  245535, 3, 1, 1.732051,  947813952,
  245536, 1, 1, 1.000000,  315937984,
  247290, 1, 1, 1.000000,  315943840,
  249046, 9, 2, 3.828427, 2843453696,
  249047, 11, 3, 5.700170, 3475382528,
  250809, 5, 1, 2.236068, 1579689984)

df <- bin2latlon(lonlat$bin_num)
df

latlon2bin(lat = df$latitude, lon = df$longitude)
*/
