#ifndef UTIL_IO_UTIL_H
#define UTIL_IO_UTIL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>
#include <sstream>

namespace Util
{

   /**
   * Return string representation of an integer.
   *
   * \param n integer to be converted.
   */
   std::string toString(int n);

   /**
   * Strip trailing whitespace from a string.
   *
   * \param  string string (stripped upon return).
   * \return length of stripped string.
   */
   int rStrip(std::string& string);

   /**
   * Extract string from stream, and compare to expected value. 
   *
   * \throw Exception if input value differs from expected value.
   *
   * \param in input stream
   * \param expected expected value of string read from stream
   */
   void checkString(std::istream& in, const std::string& expected);

   /**
   * Read the next line into a stringstream.
   *
   * \param in   input stream from which to read.
   * \param line stringstream containing line, on output.
   */
   bool getLine(std::istream& in, std::stringstream& line);

   /**
   * Get the next non-empty line of input, skipping any empty lines.
   *
   * Variant of std::getline() that skips empty lines.
   *
   * \param in input stream from which to read.
   * \param line string with next non-empty line, on output.
   * \return true if not end-of-file, false if end-of-file.
   */
   bool getNextLine(std::istream& in, std::string& line);

   /**
   * Read the next non-empty line into a stringstream.
   *
   * Variant of std::getline() that skips empty lines and uses stringstream.
   *
   * \param in   input stream from which to read.
   * \param line stringstream containing next non-empty line, on output.
   * \return true if not end-of-file, false if end-of-file.
   */
   bool getNextLine(std::istream& in, std::stringstream& line);

}
#endif
