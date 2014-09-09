#ifndef UTIL_PARSER_STRING_CPP
#define UTIL_PARSER_STRING_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ParserString.h"

namespace Util
{

   /*
   * Constructor.
   */
   ParserString::ParserString() 
    : stringPtr_(0),
      end_(0),
      cursor_(0),
      c_('\0')
   {}

   /*
   * Initialize string and cursor.
   */
   void ParserString::setString(const std::string& string, int cursor)
   {
      stringPtr_ = &string;
      end_ = string.length();
      setCursor(cursor);
   }

   /*
   * Set cursor. String must already be set.
   */
   void ParserString::setCursor(int cursor)
   {
      if (!stringPtr_) {
         UTIL_THROW("No string");
      }
      if (cursor >= end_) {
         UTIL_THROW("Error: cursor >= end_");
      }
      if (cursor < 0) {
         UTIL_THROW("Error: cursor < 0");
      }
      cursor_ = cursor;
      c_ = (*stringPtr_)[cursor];
   }

}
#endif
