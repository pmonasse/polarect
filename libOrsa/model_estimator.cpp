/**
 * @file model_estimator.cpp
 * @brief Model estimation by ORSA (aka AC-RANSAC) algorithm
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 *
 * Copyright (c) 2010-2011,2020 Pascal Monasse
 * Copyright (c) 2010-2011 Pierre Moulon
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "model_estimator.hpp"

namespace orsa {

/// Matrix \a data is mxn, representing n datapoints of dimension m.
ModelEstimator::ModelEstimator(const Mat &data, bool symmetric)
: symError(symmetric), data_(data) {}

/// If multiple solutions are possible, return false.
bool ModelEstimator::ComputeModel(const std::vector<int> &indices,
                                  Model *model) const {
  std::vector<Model> models;
  Fit(indices, &models);
  if(models.size() != 1)
    return false;
  *model = models.front();
  return true;
}

} // namespace orsa
