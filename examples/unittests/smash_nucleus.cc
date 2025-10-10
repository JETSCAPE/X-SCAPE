/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "SMASHNucleusWrapper.h"

#include "gtest/gtest.h"

using namespace Jetscape;

// Test the constructor of the SmashNucleusWrapper
TEST(SMASHNucleusWrapperTest, ConstructorTest) {
  SMASHNucleusWrapper smash_nucleus_wrapper;

  EXPECT_EQ(smash_nucleus_wrapper.GetId(), "SMASHNucleusWrapper");
}