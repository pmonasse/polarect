/**
 * @file match.hpp
 * @brief Match between two image points.
 * @author Pascal Monasse, Pierre Moulon
 *
 * Copyright (c) 2010-2011 Pascal Monasse, Pierre Moulon
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD License. You
 * should have received a copy of this license along this program. If
 * not, see <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef MATCH_H
#define MATCH_H

#include "libNumerics/matrix.h"
#include <vector>
#include <fstream>
#include <sstream>

/// Save matching position between two points.
struct Match {
    Match() {}
    Match(float ix1, float iy1, float ix2, float iy2)
    : x1(ix1), y1(iy1), x2(ix2), y2(iy2) {}
    float x1, y1, x2, y2;

    /**
    * Transform into matrix where each column is (x1 y1 x2 y2)^T.
    * \param m The matches that we transform.
    * \return Matrix 4xn where n is the size of vector \a m.
    */
    static libNumerics::matrix<double> toMat(const std::vector<Match>& m) {
      libNumerics::matrix<double> M(4, static_cast<int>(m.size()));
      std::vector<Match>::const_iterator it = m.begin();
      for(int j=0; it!=m.end(); ++it, ++j) {
        M(0,j) = it->x1;
        M(1,j) = it->y1;
        M(2,j) = it->x2;
        M(3,j) = it->y2;
      }
      return M;
    }

    /**
    * Load the corresponding matches from file.
    * \param nameFile   The file where matches were saved.
    * \param vec_match  The loaded corresponding points.
    * \return bool      True if everything was ok, otherwise false.
    */
    static bool loadMatch(const char* nameFile, std::vector<Match>& vec_match)
    {
      vec_match.clear();
      std::ifstream f(nameFile);
      while( f.good() ) {
        std::string str;
        std::getline(f, str);
        if( f.good() ) {
          std::istringstream s(str);
          Match m;
          s >> m;
          if(!s.fail() )
            vec_match.push_back(m);
        }
      }
      return !vec_match.empty();
    }
    
    /**
    * Save the corresponding matches to file.
    * \param nameFile   The file where matches will be saved.
    * \param vec_match  The matches that we want to save.
    * \return bool True if everything was ok, otherwise false.
    */
    static bool saveMatch(const char* nameFile, const std::vector<Match>& vec_match)
    {
      std::ofstream f(nameFile);
      if( f.is_open() ) {
        std::vector<Match>::const_iterator it = vec_match.begin();
        for(; it != vec_match.end(); ++it)
          f << *it;
      }
      return f.is_open();
    }

    /// Lexicographical ordering of matches. Used to remove duplicates.
    friend bool operator<(const Match& m1, const Match& m2)
    {
      if(m1.x1 < m2.x1) return true;
      if(m1.x1 > m2.x1) return false;

      if(m1.y1 < m2.y1) return true;
      if(m1.y1 > m2.y1) return false;

      if(m1.x2 < m2.x2) return true;
      if(m1.x2 > m2.x2) return false;

      return (m1.y2 < m2.y2);
    }

    /// Comparison Operator
    friend bool operator==(const Match& m1, const Match& m2)
    {
      return (m1.x1==m2.x1 && m1.y1==m2.y1 &&
              m1.x2==m2.x2 && m1.y2==m2.y2);
    }

    friend std::ostream& operator<<(std::ostream& os, const Match & m)
    {
      return os << m.x1 << " " << m.y1 << " "
        << m.x2 << " " << m.y2 << std::endl;
    }

    friend std::istream& operator>>(std::istream & in, Match & m)
    {
      return in >> m.x1 >> m.y1 >> m.x2 >> m.y2;
    }
};

/// Inverse match
inline Match inv(const Match& m) {
    return Match(m.x2,m.y2,m.x1,m.y1);
}

#endif // MATCH_H
