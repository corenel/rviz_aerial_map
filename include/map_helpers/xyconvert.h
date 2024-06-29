#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

class XYConvert {
 public:
  std::vector<std::pair<double, double>> gcj2wgs(
      const std::vector<std::pair<double, double>>& gcj_coords) {
    std::vector<std::pair<double, double>> wgs_coords;
    for (const auto& coord : gcj_coords) {
      double lng = coord.first;
      double lat = coord.second;
      double dlat = transformlat(lng - 105.0, lat - 35.0);
      double dlng = transformlng(lng - 105.0, lat - 35.0);
      double radlat = lat / 180.0 * M_PI;
      double magic = sin(radlat);
      magic = 1 - ee * magic * magic;
      double sqrtmagic = sqrt(magic);
      dlat = (dlat * 180.0) / ((a * (1 - ee)) / (magic * sqrtmagic) * M_PI);
      dlng = (dlng * 180.0) / (a / sqrtmagic * cos(radlat) * M_PI);
      double mglat = lat + dlat;
      double mglng = lng + dlng;
      wgs_coords.emplace_back(lng * 2 - mglng, lat * 2 - mglat);
    }
    return wgs_coords;
  }

  std::vector<std::pair<double, double>> wgs2gcj(
      const std::vector<std::pair<double, double>>& wgs_coords) {
    std::vector<std::pair<double, double>> gcj_coords;
    for (const auto& coord : wgs_coords) {
      double lng = coord.first;
      double lat = coord.second;
      double dlat = transformlat(lng - 105.0, lat - 35.0);
      double dlng = transformlng(lng - 105.0, lat - 35.0);
      double radlat = lat / 180.0 * M_PI;
      double magic = sin(radlat);
      magic = 1 - ee * magic * magic;
      double sqrtmagic = sqrt(magic);
      dlat = (dlat * 180.0) / ((a * (1 - ee)) / (magic * sqrtmagic) * M_PI);
      dlng = (dlng * 180.0) / (a / sqrtmagic * cos(radlat) * M_PI);
      double gcjlat = lat + dlat;
      double gcjlng = lng + dlng;
      gcj_coords.emplace_back(gcjlng, gcjlat);
    }
    return gcj_coords;
  }

  std::vector<std::pair<double, double>> gcj2bd(
      const std::vector<std::pair<double, double>>& gcj_coords) {
    std::vector<std::pair<double, double>> bd_coords;
    for (const auto& coord : gcj_coords) {
      double gcjLon = coord.first;
      double gcjLat = coord.second;
      double z = sqrt(gcjLon * gcjLon + gcjLat * gcjLat) +
                 0.00002 * sin(gcjLat * M_PI * 3000.0 / 180.0);
      double theta = atan2(gcjLat, gcjLon) +
                     0.000003 * cos(gcjLon * M_PI * 3000.0 / 180.0);
      double bdLon = z * cos(theta) + 0.0065;
      double bdLat = z * sin(theta) + 0.006;
      bd_coords.emplace_back(bdLon, bdLat);
    }
    return bd_coords;
  }

  std::vector<std::pair<double, double>> bd2gcj(
      const std::vector<std::pair<double, double>>& bd_coords) {
    std::vector<std::pair<double, double>> gcj_coords;
    for (const auto& coord : bd_coords) {
      double bdLon = coord.first;
      double bdLat = coord.second;
      double x = bdLon - 0.0065;
      double y = bdLat - 0.006;
      double z = sqrt(x * x + y * y) - 0.00002 * sin(y * M_PI * 3000.0 / 180.0);
      double theta = atan2(y, x) - 0.000003 * cos(x * M_PI * 3000.0 / 180.0);
      double gcjLon = z * cos(theta);
      double gcjLat = z * sin(theta);
      gcj_coords.emplace_back(gcjLon, gcjLat);
    }
    return gcj_coords;
  }

  std::vector<std::pair<double, double>> wgs2bd(
      const std::vector<std::pair<double, double>>& wgs_coords) {
    return gcj2bd(wgs2gcj(wgs_coords));
  }

  std::vector<std::pair<double, double>> bd2wgs(
      const std::vector<std::pair<double, double>>& bd_coords) {
    return gcj2wgs(bd2gcj(bd_coords));
  }

  bool out_of_china(double lon, double lat) {
    if (lon > 73.66 && lon < 135.05 && lat > 3.86 && lat < 53.55) {
      return false;
    }

    return true;
  }

 private:
  const double a = 6378245.0;
  const double ee = 0.006693421883570923;

  double transformlat(double lng, double lat) {
    double ret = -100.0 + 2.0 * lng + 3.0 * lat + 0.2 * lat * lat +
                 0.1 * lng * lat + 0.2 * sqrt(fabs(lng));
    ret += (20.0 * sin(6.0 * lng * M_PI) + 20.0 * sin(2.0 * lng * M_PI)) * 2.0 /
           3.0;
    ret += (20.0 * sin(lat * M_PI) + 40.0 * sin(lat / 3.0 * M_PI)) * 2.0 / 3.0;
    ret += (160.0 * sin(lat / 12.0 * M_PI) + 320 * sin(lat * M_PI / 30.0)) *
           2.0 / 3.0;
    return ret;
  }

  double transformlng(double lng, double lat) {
    double ret = 300.0 + lng + 2.0 * lat + 0.1 * lng * lng + 0.1 * lng * lat +
                 0.1 * sqrt(fabs(lng));
    ret += (20.0 * sin(6.0 * lng * M_PI) + 20.0 * sin(2.0 * lng * M_PI)) * 2.0 /
           3.0;
    ret += (20.0 * sin(lng * M_PI) + 40.0 * sin(lng / 3.0 * M_PI)) * 2.0 / 3.0;
    ret += (150.0 * sin(lng / 12.0 * M_PI) + 300.0 * sin(lng * M_PI / 30.0)) *
           2.0 / 3.0;
    return ret;
  }
};
