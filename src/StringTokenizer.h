/*
 ***********************************************************************
 * Class: StringTokenizer                                              *
 * By Arash Partow - 2000                                              *
 * URL: http://www.partow.net/programming/stringtokenizer/index.html   *
 *                                                                     *
 * Copyright Notice:                                                   *
 * Free use of this library is permitted under the guidelines and      *
 * in accordance with the most current version of the Common Public    *
 * License.                                                            *
 * http://www.opensource.org/licenses/cpl.php                          *
 *                                                                     *
 ***********************************************************************
*/



#ifndef INCLUDE_STRINGTOKENIZER_H
#define INCLUDE_STRINGTOKENIZER_H


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>


class StringTokenizer
{

   public:

    StringTokenizer(const std::string& _str, const std::string& _delim);
   ~StringTokenizer(){};
    int         countTokens();
    std::string nextToken();

    inline bool hasMoreTokens() {
     return (token_str.length() > 0);
  };

    inline int nextIntToken(){
       return atoi(nextToken().c_str());
  };

    inline double nextFloatToken(){
       return atof(nextToken().c_str());
  };

    inline std::string remainingString(){
       return token_str;
  };

    std::string nextToken(const std::string& delim);
    std::string filterNextToken(const std::string& filterStr);

   private:

    std::string  token_str;
    std::string  delim;

};

#endif

