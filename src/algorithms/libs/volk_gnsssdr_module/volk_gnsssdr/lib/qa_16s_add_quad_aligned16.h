/* Copyright (C) 2010-2014 (see AUTHORS file for a list of contributors)
 *
 * This file is part of GNSS-SDR.
 *
 * GNSS-SDR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GNSS-SDR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNSS-SDR. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDED_QA_16S_ADD_QUAD_ALIGNED16_H
#define INCLUDED_QA_16S_ADD_QUAD_ALIGNED16_H

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

class qa_16s_add_quad_aligned16 : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE (qa_16s_add_quad_aligned16);
    CPPUNIT_TEST (t1);
    CPPUNIT_TEST_SUITE_END ();

private:
    void t1 ();
};


#endif /* INCLUDED_QA_16S_ADD_QUAD_ALIGNED16_H */
