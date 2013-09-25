/*

Copyright (c) 2003, Cornell University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
    - Neither the name of Cornell University nor the names of its
        contributors may be used to endorse or promote products derived from
        this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

*/

/////////////////////////////////////////////////////////////////////
///
/// ItemsetOutput.cpp
/// File stream output for itemsets
///
/////////////////////////////////////////////////////////////////////

#include "ItemsetOutput.h"


using namespace std;
/////////////////////////////////////////////////////////////////////
/// @addtogroup InputOutput
/** @{ */

/////////////////////////////////////////////////////////////////////
/// Open file for outputting itemsets
///
/// @param filename           filename for output
/////////////////////////////////////////////////////////////////////
ItemsetOutput::ItemsetOutput(char *filename) {
    out = fopen(filename,"wt");
	 //out = fopen(filename,"a");
}

/////////////////////////////////////////////////////////////////////
/// Close file output stream
/////////////////////////////////////////////////////////////////////
ItemsetOutput::~ItemsetOutput() {
    if (out)
        fclose(out);
}

/////////////////////////////////////////////////////////////////////
/// Check if file is open
///
/// @return true if file is open
/////////////////////////////////////////////////////////////////////
int ItemsetOutput::isOpen() {
    if (out)
        return 1;
    else
        return 0;
}

/////////////////////////////////////////////////////////////////////
/// Print out an itemset
///
/// @param length          number of items in itemset
/// @param itemset         array of items
/// @param support         support of itemset (-1 to suppress)
/////////////////////////////////////////////////////////////////////
void ItemsetOutput::printSet(double score, double score2,int length, std::string *iset, std::string *tids,int support) {

    fprintf(out, "%d ", length);
    fprintf(out, "%d ", support);
	fprintf(out, "%f ", score);
	fprintf(out, "%f\n", score2);

	for (int index = 0; index < length; index++)
        fprintf(out, "%s ", iset[index].c_str());

	fprintf(out, "\n");
    for (int j = 0; j < support; j++){
         fprintf(out, "%s ", tids[j].c_str());
    } 

	fprintf(out, "\n");


}


/** @} */