/*
 
 debi.cpp: C++ implementation for the DeBi algorithm. The algorithm  is based on a
 well known data mining approach called frequent itemset. MAFIA algorithm is used
 for mining maximal frequent itemsets. The details of the DeBi algorithm can be found
 in the following paper.
 
 "Discovering Differentially Expressed Biclusters using a Frequent Itemset Approach" 
 by Serin et al.
 
 ==============
 
 USAGE: 
 
 ----Input----
 
 ./debi inFile outFilePath
 
 inFile: has to be tab deliminated file with the first row 'sample ids' and first column 'gene ids'. 
 outFilePath: is the path of the biclustering results folder
 
 -b binarization level (for binary data not needed otherwise it must be defined) 
 -o maximum overlap level of the biclusters default=0.5 (overlap=1: no overlap is allowed, overlap=0: all the biclusters are taken, no filtering applied)
 -p pattern of regulation up or down (u, d) default=u
 
 
 look at README.txt for sample commands. 
 
 ----Output----
 
 output file format
 <bicluster gene size> <bicluster sample size> <bicluster score> <normalized bicluster score>
 genes in bicluster 
 samples in bicluster
 
 ==============
 
 (c) 2010, Akdes Serin, Max Planck Institute for Molecular Genetics, Berlin
 Germany
 
 
 --- MAFIA ALGORITHM 
 
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

//#define DEBUG
#include <cstring>
#include <stdlib.h>
#include <climits>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include <list>
#include <vector>
#include <algorithm>
#include <functional>
#include <string>
#include <assert.h>
#include <cmath>
#include <map>
#include "Bitmap.h"
#include "TreeNode.h"
#include "ItemsetOutput.h"
#include "nmath.hpp"
#include "dpq.h"
#include "Rmath.h"

using namespace std;
using Rmath::phyper;

#define DESCRIPTION "find biclusters with DeBi algorithm"
#define MAX_NUM_ITEMS 100000


/// @defgroup GlobalVariables Global Variables
/// Global variables
/** @{ */
typedef vector<TreeNode *> NodeList;
typedef vector<TreeNode *> BranchList;
typedef vector<Bitmap *> BitmapList;
typedef vector<BaseBitmap *> BaseBitmapList;
typedef vector<int> ItemSet;
typedef map<long, ItemSet *> HashTable;
vector<int> Transactions;
vector <int> ItemList;


/// Simple class for storing subtree size estimates
class SubtreeEstimate {
public:
	int Count; ///< Count of actual subtrees counted
	int Sum; ///< Sum of subtree sizes
	SubtreeEstimate() {
		Count = 0;
		Sum = 0;
	}
	~SubtreeEstimate() {
	}
};

/// Simple class for storing tail elements of each node of the tree
class TailElement {
public:
	int Count; ///< Support of the 1-extension
	int Item; ///< Item-id for this1 -extension
	
	TailElement() {
		Count = 0;
		Item = 0;
	}
	
	~TailElement() {
	}
	bool operator<(const TailElement& rhs) const {
		return this->Count < rhs.Count;
	}
	;
};




/// @defgroup CommandLineParameters Command Line Parameters/Variables
/// Commmand line parameters from user or inferred
/** @{ */
double duration = -1;
char* outFilename; ///< filename for output
ItemsetOutput *outFile; ///< file for ouput
int ItemCount=0; ///<  of items in the file
int TransCount=0; ///< # of transactions in the file
int MS=0; ///< maximum sup as transaction count
/** @} */

/// @defgroup added Command Line Parameters/Variables
/** @} */
double overlap; ///< user-defined  min overlap size
double alpha; ///< user-defined min pvalue for gene sample association
int MSF; ///< user-defined min sample size of the biclusters
double level; ///< user-defined gene expression discretization threshold
char pattern; ///<  user-defined regulation pattern either up or down
double noise; ///< user-defined gene noise level in the biclusters (not used)
/** @} */


/// @defgroup added GlobalVariables Global Variables
/** @{ */
BitmapList TransTemp; 
vector<int> rowIds;
vector<string> rowToName; 
vector<string> columnToName; 
BaseBitmapList MFT; ///< List of transaction ids of maximally frequent itemsets
HashTable HT; ///< Hash table of transaction supports
vector<int> SupportCountList; ///< List that stores support count
vector<double> SignificanceScore;
vector<double> NSignificanceScore; 
vector<double> MaxNSignificanceScore; 

BaseBitmap* TempName; ///< Temporary buffer for one name bitmap
SubtreeEstimate* EstimateBuffy; ///< Buffer of subtree estimates
int *ItemsetBuffy; ///< Buffer for writing itemsets to file output
int *sim; ///< Vector for bicluster similarities, used in filtering step
int minLen = 10; ///< Final minumum bicluster gene size 
//
BitmapList MFI1; ///< Gene sets in the found biclusters after step 1
BitmapList MFT1; ///< Sample sets in the found biclusters after step 1
BitmapList MFI2; ///< Gene sets in the found biclusters after step 2
BitmapList MFT2; ///< Sample sets in the found biclusters after step 2
double seedLength; ///< Bicluster seed length
char type; ///< Type of the data either binary or continiuous
int setNum = 0; ///< # of biclusters
int *Counters;
double p_score;
double sum_pscore;
int setSize = 0;
double max_pscore=0;
int max_alpha=1;
int max_setSize=0;
/** @} */


/** @} */
int VerScale = 1; ///< Scaling factor for transactions
int HorScale = 1; ///< Scaling factor for items
bool GoFHUT = true; ///< FHUT flag
bool HUTMFI = true; ///< HUTMFI flag
bool PEPrune = true; ///< PEPrune flag -- parent equivalent pruning
bool Reorder = true; ///< Reorder flag
/** @} */

/// @defgroup CounterVariables Counter Variables
/// Counter variables for information gathering
/** @{ */
int CountFHUT = 0; ///< # of times FHUT was successful
int CountNodes = 0; ///< # of frequent nodes in the tree
int CountCounts = 0; ///< # of Counts or all nodes in the tree
int CountAnds = 0; ///< # of ANDs of normal bitmaps
int CountSmallAnds = 0; ///< # of compressed bitmap ANDs
int CountPEPrunes = 0; ///< # of PEPruning
int CountCheckPosition = 0; ///< # of CheckPosition calls
int CountHUTMFI = 0; ///< # of HUTMFI attempts
int CountHUTMFISuccess = 0; ///< # of HUTMFI successes
int CountRebuilds; ///< # of Rebuilds
/** @} */

/// @defgroup ProgramVariables Program Parameters/Variables
/// Useful program parameters/counters
/** @{ */
int maxtail = 0;
int FinalMFISize = 0;
int MFISize = 0; ///< MFI size before pruning
int MFIDepth = 0; ///< The aggregated depth of the all MFI elements
int F1size = 0; ///< # of frequent 1-itemsets after merging repeats
int FullF1size = 0; ///< # of frequent 1-itemsets
int k = 50; ///< # of items checked for a MFI lookup
int MAX_compID = 1; ///< max compression ID
int projectDepth = -1; ///< depth of the bitmap you're projecting from
int EstimateSize; ///< size of subtree estimation buffer
int EstimateDiv = 5; ///< bucket size by frequent tail length
int maxItemsetSize = 0; ///< max size of a frequent itemset
/** @} */

/// @defgroup DataVariables Data variables
/// Complex data structure variables
/** @{ */
NodeList F1; ///< List of frequent 1-itemsets
BitmapList TransBuffy; ///< Buffer of transaction bitmaps
BaseBitmapList NameBuffy; ///< Buffer of name bitmaps
NodeList NodeBuffy; ///< Buffer of tree nodess
TreeNode *Root; ///< The root (the nullset)
TailElement *gTail; ///< global tail pointer
TailElement *TailBuffy; ///< Common Buffer for tail elements
int *ItemMap; ///< For renaming items after sorting by support
BaseBitmapList MFI; ///< List of Maximally Frequent Itemsets
/** @} */

/// @defgroup TimingVariables Timing Variables
/// Variables for timing (and instrumenting the code)
/** @{ */
time_t total_start, total_finish;
double total_time;
time_t read_start, read_finish;
double read_time;
time_t algorithm_start, algorithm_finish;
double algorithm_time;
time_t print_start, print_finish;
double print_time;
/** @} */


/////////////////////////////////////////////////////////////////////
///  
///  start of added/modified functions for DeBi algorithm
///
/////////////////////////////////////////////////////////////////////

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  
///  deletes the variables, sets counter variables to zero
///
/////////////////////////////////////////////////////////////////////
void destroy() {
	
	
	
	if (FullF1size != 0) {
		
		delete []ItemsetBuffy;
		delete Root;
		delete []ItemMap;
		delete []gTail;
		delete []TailBuffy;
		delete []EstimateBuffy;	
		
		NodeList::iterator bli1 = NodeBuffy.begin();
		while (bli1 != NodeBuffy.end()) {
			delete (*bli1);
			bli1=NodeBuffy.erase(bli1);	
		}
		
		NodeBuffy.clear();
		
		
		NodeList::iterator bli2 = F1.begin();
		while (bli2 != F1.end()) {
			
			delete *bli2;
			bli2=F1.erase(bli2);
		}
		
		F1.clear();
		
		
		
		BaseBitmapList::iterator bli3 = MFI.begin();
		while (bli3 != MFI.end()) {
			delete *bli3;
			bli3=MFI.erase(bli3);
		}
		
		MFI.clear();
		
		
		BaseBitmapList::iterator bli4 = MFT.begin();
		while (bli4 != MFT.end()) {
			delete *bli4;
			bli4=MFT.erase(bli4);
			
		}
		
		MFT.clear();
		
		
		BaseBitmapList::iterator bli5 = NameBuffy.begin();
		while (bli5 != NameBuffy.end()) {
			delete *bli5;
			bli5=NameBuffy.erase(bli5);
			
		}
		
		
		BitmapList::iterator bli6 = TransBuffy.begin();
		while (bli6 != TransBuffy.end()) {
			delete *bli6;
			bli6=TransBuffy.erase(bli6);
			
		}
		
	}
	
	SupportCountList.clear(); ///< List that stores support count
	VerScale = 1; ///< Scaling factor for transactions
	HorScale = 1; ///< Scaling factor for items
	GoFHUT = true; ///< FHUT flag
	HUTMFI = true; ///< HUTMFI flag
	PEPrune = true; ///< PEPrune flag -- parent equivalent pruning
	Reorder = true; ///< Reorder flag
	/** @} */
	
	/// @defgroup CounterVariables Counter Variables
	/// Counter variables for information gathering
	/** @{ */
	CountFHUT = 0; ///< # of times FHUT was successful
	CountNodes = 0; ///< # of frequent nodes in the tree
	CountCounts = 0; ///< # of Counts or all nodes in the tree
	CountAnds = 0; ///< # of ANDs of normal bitmaps
	CountSmallAnds = 0; ///< # of compressed bitmap ANDs
	CountPEPrunes = 0; ///< # of PEPruning
	CountCheckPosition = 0; ///< # of CheckPosition calls
	CountHUTMFI = 0; ///< # of HUTMFI attempts
	CountHUTMFISuccess = 0; ///< # of HUTMFI successes
	
	/** @} */	
	/// @defgroup ProgramVariables Program Parameters/Variables
	/// Useful program parameters/counters
	/** @{ */
	maxtail = 0;
	MFISize = 0; ///< MFI size before pruning
	MFIDepth = 0; ///< The aggregated depth of the all MFI elements
	F1size = 0; ///< # of frequent 1-itemsets after merging repeats
	FullF1size = 0; ///< # of frequent 1-itemsets
	k = 50; ///< # of items checked for a MFI lookup
	MAX_compID = 1; ///< max compression ID
	projectDepth = -1; ///< depth of the bitmap you're projecting from
	EstimateSize = 0; ///< size of subtree estimation buffer
	EstimateDiv = 5; ///< bucket size by frequent tail length
	maxItemsetSize = 0; ///< max size of a frequent itemset
	/** @} */
	
	
	
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
/// Convert string to double
/////////////////////////////////////////////////////////////////////
double string2double(string in)
{
	stringstream ss;
	ss<<in;
	double out;
	ss>>out;
	return out;
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
/// Convert double to string
/////////////////////////////////////////////////////////////////////
string double2string(double in)
{
	std::ostringstream ss;
	ss << in;
	std::string str = ss.str();
	return str;
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  Return the name corresponding to a row index.
/////////////////////////////////////////////////////////////////////
string getRowName(unsigned int index)
{
	if (index >= rowToName.size())
    {
		cerr << "ERROR: row/probe index " << index << " is at least as large as "
		<< " number of rows/probes" << rowToName.size() << endl;
		return("");
    }
	
	return(rowToName[index]);
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  Return the column name corresponding to a column index.
/////////////////////////////////////////////////////////////////////
string getColumnName(unsigned int index)
{
	if (index >= columnToName.size())
    {
		cerr << "ERROR: sample index " << index << " is at least as large as "
		<< " number of columns" << columnToName.size() << endl;
		return("");
    }
	
	return(columnToName[index]);
}


/////////////////////////////////////////////////////////////////////
/// Insert pointer to node in order of increasing support
///
/// @param newNode item node to add to F1
/////////////////////////////////////////////////////////////////////

void AddToF1(TreeNode *newNode) {
	if (F1.empty())
		F1.push_back(newNode);
	else {
		// use insertion sort for increasing support ordering
		NodeList::iterator noli = F1.begin();
		while (noli != F1.end()) {
			if ((*noli)->Trans->_count >= newNode->Trans->_count) {
				F1.insert(noli, newNode);
				return;
			}
			noli++;
		}
		
		// Add to end of F1 list
		F1.push_back(newNode);
	}
}


/*
 modified function
 */
/////////////////////////////////////////////////////////////////////
/// Read transaction data from file and build item bitmaps.
///
/// @param filename      name of file to be read in
/// @param pattern       up regulated or down regulated
/////////////////////////////////////////////////////////////////////
void F1FromFile(string filename, char pattern) {
	
	
	
	int itemIndex = 0;
	
	
	time(&read_start);
	
	int cnt,total=0;
	int z=-2;
	
	ifstream in(filename.c_str());
	if (!in)
	{
		cerr << "Cannot open file " << in << " containing gene expression data." << endl;
		
	}
	
	stringstream topline;
	string line,token;
	
	//get the first line
	getline(in,line);
	topline<<line<<"\t";
	
	//get all the sample names
	while(getline(topline,token,'\t'))
		if(!(token=="NA" || token==""))
		{
			if(++z>=0)
			{
				columnToName.push_back(token);
			}
		}
	
	itemIndex=0;
	int ind=0;
	
	//get all the rows
	while(getline(in,line))
	{
		
		stringstream ss;
		ss<<line;
		string gene;
		getline(ss,gene,'\t');	
		rowToName.push_back(gene);
		
		if(gene=="")
		{
			break;
		}
		
		cnt=0;
		ItemList.push_back(ind);
		
		Bitmap *trans = new Bitmap(z+1);
		
		for(int i=0;i<=z;i++)
		{
			string token;
			
			getline(ss,token,'\t');
			
			
			if(!(token=="NA" || token==""))
			{
				
				if(type=='b'){
					if(!(string2double(token)==1 || string2double(token)==0)){
						cerr << "File containing gene expression data is not binary. You have to set the binarization threshold parameter." << endl;
						exit(1);
					}
				}
				
				//get the up regulated pattern
				if(pattern=='u'){
					
					if(string2double(token)>=level){
						
						trans->FillEmptyPosition(i);
						cnt++;
					}
				}
				
				//get the down regulated pattern
				if(pattern=='d'){
					if((string2double(token))<=level){
						
						trans->FillEmptyPosition(i);
						cnt++;
					}
				}
			}
		}
		
		//fill the transaction list
		if(cnt>0)
		{
			trans->_count=cnt;
			TransTemp.push_back(trans);
			Counters[itemIndex]=cnt;
			total=total+cnt;
			rowIds.push_back(ind);
			itemIndex++;
			
		}
		if(cnt<=0){
			delete trans;
		}
		
		if (cnt >= MS) {
			MS=cnt; //support value of the gene with the highest support
			
		}
		
		ind++;
	}
	
	
	ItemCount =itemIndex;
	TransCount =z+1;
	double size=(rowToName.size()*TransCount);
	double perc=total/size; //percentage of ones
	
	MSF=TransCount*0.01; //minumum support value
	printf("Number of genes:	%d",rowToName.size());
	printf("\t Number of samples:    %d\n",TransCount);
	printf("Percentage of ones:    %.3f\n",perc);
	
	// build the normal bitmaps -- Preallocated memory for the bitmaps
	// Create F1
	itemIndex = 0;
	for  (int itemIndex=0; itemIndex<ItemCount;itemIndex++) {
		
		if (Counters[itemIndex] >= MS) {
			
			
			TreeNode *node = new TreeNode(NULL, TransTemp[itemIndex], 1, -1, itemIndex, -1, 0);
			AddToF1(node);
			
		}
		
	}
	
	time(&read_finish);
	read_time = difftime(read_finish, read_start);
	printf("Reading data input time:    %.3f seconds.\n", read_time);
	
	
}


/*
 added function
 */
/////////////////////////////////////////////////////////////////////
/// Deleting the items from the search space
/// @param deleteItems items that we want to remove from the search space
/////////////////////////////////////////////////////////////////////
void UpdateF1(vector<int> deleteItems, char *filename) {
	
	
	MS--;
	int *F1IndexList = new int[ItemCount]; // the id's of F1 items
	int *InvF1IndexList = new int[ItemCount]; // the id's of F1 items
	int itemIndex = 0;
	
	//Initialize counters
	for (int ct = 0; ct < ItemCount; ++ct) {
		F1IndexList[ct] = -1;
		InvF1IndexList[ct] = -1;
	}
	
	
	//Deleting the items from the search space
	vector<int>::iterator itNum = ItemList.begin();
	if( deleteItems.size()>0){
		for (int i = 0; i < deleteItems.size(); i++) {
			for(; itNum < ItemList.end(); itNum++)
			{
				
				if(deleteItems[i]==(*itNum))
				{
					
					itNum=ItemList.erase(itNum);
				}
				
				
			}
			Counters[deleteItems[i]] = 0;
		}
	}
	
	
	//build the normal bitmaps -- Preallocated memory for the bitmaps
	//Create F1
	itemIndex = 0;
	int F1items = 0;
	for  (int itemIndex=0; itemIndex<ItemList.size();itemIndex++) {
		
		if (Counters[ItemList[itemIndex]] >= MS) {
			F1IndexList[F1items++] = ItemList[itemIndex];
			InvF1IndexList[ItemList[itemIndex]] = F1items - 1;
			TreeNode *node = new TreeNode(NULL, TransTemp[ItemList[itemIndex]], 1, -1,
										  F1IndexList[F1items-1], -1, 0);
			
			AddToF1(node);
			
		}
		
	}
	
	delete[] F1IndexList;
	delete[] InvF1IndexList;
	
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
/// apply Fisher exact test
/////////////////////////////////////////////////////////////////////


double fisher (double a, double b, double c, double d)
{
	if(a<0 ||b<0||c<0||d<0)
	{
		return 1;
	}
	else{
		return(phyper(a - 1, a+c, d+b, a+b,0,0));
	}
}
/*
 added function
 */
///////////////////////////////////////////////////////////////////////////////////////
/// permute the labels for the emprical distribution of cluster scores
///////////////////////////////////////////////////////////////////////////////////////

double Permute( int i,int len,int sample) {
	
	
	int a, b, c, d ;
	double p,score,sum_score=0;
	srand ( time(NULL) );
	Bitmap *trs = new Bitmap(TransCount);
	Bitmap *clOnes;
	
	for (int m=0;m<sample;m++)
	{
		
		int CountSmallAnds=0;
		trs->FillZeros();
		trs->Count(CountSmallAnds);
		
		while (trs->_count<=len){ //keep the num of ones in phenotype vector     
  		 	
			trs->FillEmptyPosition(rand() % (TransCount-1)); //permute the phenotype vector
			trs->Count(CountSmallAnds);	 	 
		}
	   	
		score=0;
		int itemIndex = 0;
		BitmapList::iterator bli = TransTemp.begin();
		
		//calculate the bicluster score for the permuted labels
		while (bli != TransTemp.end()) {
			
			if (MFI2[i]->CheckPosition(itemIndex, CountCheckPosition) > 0) {
				
				Bitmap *item =  * bli;	
				clOnes = new Bitmap(TransCount);
				clOnes->AndOnly(*item, *trs, CountSmallAnds);
				clOnes->Count(CountSmallAnds);			
				
				a = clOnes->_count;
				b = ((*bli)->_count) - a;
				c = trs->_count - a;
				d = (TransCount - (trs->_count)) - b;
				p=fisher(a,b,c,d);
				score = score -log(p);
				
				delete clOnes;
				
			}
			
			itemIndex++;
			bli++;
		}

		sum_score=sum_score+score;
	}
	
	delete trs;

	return(sum_score/sample); //return the estimated mean from the emprical distrubution
	
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
/// calculate the bicluster score
/////////////////////////////////////////////////////////////////////
void calcBicScore() {
	
	
	sum_pscore=0; //sum of biclustering score 
	double permuted_score=0; //estimated mean for the bicluster
	double n_score=0; //normalized biclustering score
	//calcSampleFreq(); //for weighted permutation
	
	setSize = 0;
	
	if (setNum > 0) {
		// translate bitmap to list of ids
		
		for (int i = 0; i < setNum; i++) {
			
			NSignificanceScore.push_back(0);
			if (sim[i] == 0) {
			
				if ((MFI2[i]->_count) >= minLen) {
				
					setSize++;				
					
					string* TID = new string[TransCount];
					
					int j2 = 0;
					for (int cc = 0; cc < TransCount; cc++) {
						if (MFT2[i]->CheckPosition(cc, CountCheckPosition) > 0) {
							TID[j2] = getColumnName(cc);
							j2++;
						}
					}
					
					//added small integer for zero permuted scores
					permuted_score=Permute(i,j2,500)+1;
					//permuted_score=PermuteWS(i,j2,500)+1;  //weighted sampling
					n_score=SignificanceScore[i]/permuted_score;
					NSignificanceScore[i]=n_score;
				    sum_pscore=sum_pscore+n_score; 
				    delete []TID;
				}
			}
		}
		
	} 
	delete []sim;
}

/*
 modified function
 */
/////////////////////////////////////////////////////////////////////
/// print the biclusters out to a file
/////////////////////////////////////////////////////////////////////

void PrintMFI() {
	
	
	sum_pscore=0; //sum of biclustering score 
	double n_score=0; //normalized biclustering score
	//calcSampleFreq(); //for weighted permutation
	
	if (!outFile->isOpen()) {
		cerr << "Output file not open!" << endl;
		exit(1);
	}
	setSize = 0;
	
	if (setNum > 0) {
		// translate bitmap to list of ids
		
		for (int i = 0; i < setNum; i++) {
			if (sim[i] == 0) {
				
				if ((MFI2[i]->_count) >= minLen) {
					setSize++;				
					int j1 = 0;
					string* ITEMS = new string[ItemCount];
					string* TID = new string[TransCount];
					for (int cc = 0; cc < ItemCount; cc++) {
						if (MFI2[i]->CheckPosition(cc, CountCheckPosition) > 0) {
							ITEMS[j1]= getRowName(rowIds[cc]);
							j1++;
						}
					}
					int j2 = 0;
					for (int cc = 0; cc < TransCount; cc++) {
						if (MFT2[i]->CheckPosition(cc, CountCheckPosition) > 0) {
							TID[j2] = getColumnName(cc);
							j2++;
						}
					}

					outFile->printSet(SignificanceScore[i],MaxNSignificanceScore[i],j1,ITEMS, TID, j2);
				    sum_pscore=sum_pscore+n_score; 
				    delete []ITEMS;
				    delete []TID;
				}
			}
		}
		
	} else {
		outFile->printSet(0,0, 0,NULL,NULL, 0);
	}
	
	delete []sim;
	delete outFile;
}


/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  filter out the similar sets
///  similarity is assesed by the user-defined overlap parameter
/////////////////////////////////////////////////////////////////////

void Filter() {
	
	int i, j;
	
	double common = 1;
	
	sim = new int[setNum];
	
	for (int k = 0; k < setNum; k++) {
		sim[k] = 0;
	}
	
	// sort the biclusters based on their size
	if (overlap>0){
		
		for (i = 0; i < setNum; i++) {
			
			for (j = 0; j < (setNum - i) - 1; j++) {
				
				if ((MFT2[j]->_count * MFI2[j]->_count) < (MFT2[j + 1]->_count
														   * MFI2[j + 1]->_count)) {
					
					Bitmap *temp1 = MFT2[j];
					Bitmap *temp2  = MFI2[j];
					
					temp1->_count = MFT2[j]->_count;
					temp2->_count = MFI2[j]->_count;
					
					MFT2[j] = MFT2[j + 1];
					MFI2[j] = MFI2[j + 1];
					
					MFT2[j]->_count = MFT2[j + 1]->_count;
					MFI2[j]->_count = MFI2[j + 1]->_count;
					
					MFT2[j + 1] = temp1;
					MFI2[j + 1] = temp2;
					
					double temp_score=SignificanceScore[j];
					SignificanceScore[j]=SignificanceScore[j+1];
					SignificanceScore[j+1]=temp_score;
				}
			}
			
		}
		
		//filter similar sets
		
		for (i = 1; i < setNum; i++) {
			
			int it = 0;
			
			while (it < i) {
				
				if (sim[it] == 0) {
					
					Bitmap *ItemsRes = new Bitmap(ItemCount);
					Bitmap *TransRes = new Bitmap(TransCount);
					
					ItemsRes->AndOnly(*MFI2[it], *MFI2[i], CountSmallAnds);
					TransRes->AndOnly(*MFT2[it], *MFT2[i], CountSmallAnds);
					
					ItemsRes->_count = ItemsRes->Count(CountSmallAnds);
					TransRes->_count = TransRes->Count(CountSmallAnds);
					
					int l2 = MFI2[i]->_count * MFT2[i]->_count;
					
					common = (ItemsRes->_count) * (TransRes->_count);
					
					common = common / l2;
					if(overlap==1 && common >0)
					{
						sim[i] = 1; //similar sets
					}
					else if (overlap<1  && common >= (1-overlap)) {
						
						sim[i] = 1; //similar sets
					}
				
					delete ItemsRes;
					delete TransRes;
					
					
				}
				it++;
				
			}
			
		}
	}
	
}



/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  in step 2 precompute cutoff values on the contingency table 
/////////////////////////////////////////////////////////////////////
int  PrecomputeCutoffIn(int oneoutcutoff,int trs_cnt,int ICount) {
	
	long double p;
	double a, b, c, d, n,onecutoff,found ;
	onecutoff = 0;
	found = 1;
	
	a = trs_cnt;
	b = (oneoutcutoff)/3;
	c = 0;
	d = (ICount - (trs_cnt)-((oneoutcutoff)/3));
	n = (a + b + c + d);
	p=fisher(a,b,c,d);
	
	// precompute cut-off values on the contingency table entries that yield a p-value just higher than alpha
	for (int cut = 0; cut < (trs_cnt ) && found == 1; cut++) {
		if (p > (alpha)) {
			onecutoff = a;
			found = 0;
		}
		
		c++;
		p=p*(double(c+d)/double(a+b))*double(a/c);
		a--;
		
		
	}
	return(onecutoff);
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  in step 2 precompute cutoff values on the contingency table 
/////////////////////////////////////////////////////////////////////
int PrecomputeCutoffOut(int onecutoff,int t,int trs_cnt,int ICount) {
	
	long double p;
	double a, b, c, d, n,oneoutcutoff,found ;
	
	a = onecutoff+((trs_cnt-onecutoff)*t/3);
	b = 0;
	c = trs_cnt-(onecutoff+((trs_cnt-onecutoff)*t/3));
	d = (ICount - (trs_cnt));
	n = (a + b + c + d);
	found = 1;
	p=fisher(a,b,c,d);
	
	// precompute cut-off values on the contingency table entries that yield a p-value just higher than alpha
	for (int cut = 0; cut < (ICount - (trs_cnt) ) && found == 1; cut++) {
		
		if (p > (alpha)) {
			oneoutcutoff = b;
			found = 0;
		}
		b++;
		p=p*((double(a+b)/double(c+d))*double(d/b));
		d--;
		
	}
	
	return (oneoutcutoff);
	
}



/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  extend the bicluster 
/////////////////////////////////////////////////////////////////////
void  ExtendGenes(Bitmap *conds ,Bitmap **ext ) {
	
	int a, b, check=1,c, d,oneoutcutoff=0,onecutoff=0,oneoutcutoff2=0,onecutoff2=0,onecutoff3=0,oneoutcutoff3=0 ;
	int itemIndex = 0;
	int cnt = 0;
	double thr=1;
	
	onecutoff=PrecomputeCutoffIn(0,conds->_count,TransCount);
	oneoutcutoff=PrecomputeCutoffOut(0,3,conds->_count,TransCount);
	onecutoff2=PrecomputeCutoffIn(oneoutcutoff*2,conds->_count,TransCount);
	oneoutcutoff2=PrecomputeCutoffOut(onecutoff,1,conds->_count,TransCount);
	onecutoff3=PrecomputeCutoffIn(oneoutcutoff,conds->_count,TransCount);
	oneoutcutoff3=PrecomputeCutoffOut(onecutoff,2,conds->_count,TransCount);
	
	BitmapList::iterator bli = TransTemp.begin();
	Bitmap *clOnes = new Bitmap(TransCount);
	while (bli != TransTemp.end()) {
		check = 1;
		int CountSmallAnds=0;
		Bitmap *item = *bli;
		
		clOnes->AndOnly(*item, *conds, CountSmallAnds);
		clOnes->Count(CountSmallAnds);
		
		
		//apply fisher's exact test	if smaller than precomputed values	
		if (clOnes->_count >= thr && (clOnes->_count) > onecutoff
			&& ((item->_count) - (clOnes->_count)) < oneoutcutoff) {
			
			
			if( (clOnes->_count < (onecutoff+((conds->_count-onecutoff)/3))) && ((item->_count) - (clOnes->_count)) > oneoutcutoff2){
			 	check=0;
			}
			
			if( (clOnes->_count) < (onecutoff2) && ((item->_count) - (clOnes->_count)) > (oneoutcutoff*2/3)){
			 	check=0;
			}
			
			
			if( (clOnes->_count < (onecutoff+((conds->_count-onecutoff)*2/3))) && ((item->_count) - (clOnes->_count)) > oneoutcutoff3){
			 	check=0;
			}
			
			if( (clOnes->_count) < (onecutoff3) && ((item->_count) - (clOnes->_count)) > (oneoutcutoff/3)){
			 	check=0;
			}
			
			if(check==1){
				a = clOnes->_count;
				b = (item->_count) - a;
				c = conds->_count - a;
				d = (TransCount - (conds->_count)) - b;
				
				double p = 0;
				p=fisher(a,b,c,d);
				
				if (p > 0) {
					
					if (p <= (alpha)) {
						
						p=fisher(a,b,c,d);
						(*ext)->FillEmptyPosition(itemIndex);
						p_score = p_score -log(p);
						cnt++;
						
					}
				}
			}
		}
		
		bli++;
		itemIndex++;
		clOnes->FillZeros();
		
	}
	delete clOnes;
	(*ext)->_count = cnt;
	
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  step 2 extending the biclusters using Fisher's exact test
/////////////////////////////////////////////////////////////////////

void StatTestF1Simple() {

	for (int w = 0; w < FinalMFISize; w++) {
		Bitmap *genes=new Bitmap(ItemCount);
		p_score=0;
		Bitmap *conds=new Bitmap(TransCount);
		
		int d=0;
		for (int cc = 0; cc < TransCount; cc++) {
			if (MFT1[w]->CheckPosition(cc, CountCheckPosition) > 0) {
				conds->FillEmptyPosition(cc);
				d++;
			}
		}
		
		conds->_count=d;
		ExtendGenes(conds,&genes);
		
		if ((genes->_count) >= seedLength){
			
			MFI2.push_back(genes);
			MFT2.push_back(conds);
			SignificanceScore.push_back(p_score);
			
			setNum++;
			
		}
		else{
			delete conds;
			delete genes;
		}
		
		
	}
	
	
}

/*
 added function
 */
/////////////////////////////////////////////////////////////////////
///  optimize the alpha value
/////////////////////////////////////////////////////////////////////

void optimizeAlpha(int min_range,int max_range,int it){

	for (int f=min_range;f<max_range;f=f++)
	{
		
		alpha=pow(double(10),double(it*f*-1));
		StatTestF1Simple();
		Filter();
		
		calcBicScore();	 
		
		if (setNum > 0){

			if(sum_pscore>max_pscore)
			{
				max_pscore=sum_pscore;
				max_alpha=f;
				max_setSize=setSize;
				
				MaxNSignificanceScore.clear();
				
				for (int i = 0; i < setNum; i++) {
				
					MaxNSignificanceScore.push_back(NSignificanceScore[i]);
				}
				
			}
			
			
			BitmapList::iterator bli = MFI2.begin();
			
			while (bli != MFI2.end()) {
				delete *bli; 	
				bli=MFI2.erase(bli); 
			}
			
			
			BitmapList::iterator bli2 = MFT2.begin();
			
			while (bli2 != MFT2.end()) {
				delete *bli2; 
				bli2=MFT2.erase(bli2); 
				
			}
			
			
			MFI2.clear();
			MFT2.clear();
			
			NSignificanceScore.clear();
			SignificanceScore.clear();
		}
		setNum=0;
		if(setSize==0)
			break;
		
		
	}

}

/////////////////////////////////////////////////////////////////////
///  
///  start of functions of MAFIA algorithm 
///
/////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////
/// Output itemset (don't need to save bitmaps for FI)
///
/// @param C                 the current node
/////////////////////////////////////////////////////////////////////

void AddToFI(TreeNode *C) {
	int itemsetIndex = 0;
	for (int cc = 0; cc < FullF1size; cc++) {
		if (C->Name->CheckPosition(cc, CountCheckPosition) > 0) {
			ItemsetBuffy[itemsetIndex] = ItemMap[cc];
			itemsetIndex++;
		}
	}
	
	
	if (C->Name->_count > maxItemsetSize)
		maxItemsetSize = C->Name->_count;
	
	// Update stat variables
	MFIDepth += C->Name->_count;
	MFISize++;
	
}

/** @} */

/*********************************************************************
 Algorithmic components
 *********************************************************************/
/// @defgroup AlgorithmicComponents Algorithmic Component Functions
/// Algorithm components (HUTMFI, PEP, etc.)
/** @{ */

/////////////////////////////////////////////////////////////////////
/// Check for an existing superset of name in the MFI
///
/// @param location          The node we're at.  Use this to examine only
///                          the relevant portion of the MFI
/// @return True - if superset found.
///         False - if no superset.
/////////////////////////////////////////////////////////////////////

bool LMFISuperSet(TreeNode* location) {
	return (location->rEnd > location->rBegin);
}

/*
 modified function
 */
/////////////////////////////////////////////////////////////////////
/// Add this node's name bitmap to the MFI list
///
/// @param C                 the current node
/////////////////////////////////////////////////////////////////////

void AddToMFI(TreeNode *C) {
	
	BaseBitmap *name = new BaseBitmap(FullF1size);
	int cnt = 0;
	for (int cc = 0; cc < FullF1size; cc++) {
		
		if (C->Name->CheckPosition(cc, CountCheckPosition) > 0) {
			name->FillEmptyPosition(cc);
			cnt++;
		}
		
	}
	name->_count = cnt;
	// add to MFI
	MFI.push_back(name);
	C->Name->_count = cnt;
	
	//	Bitmap *tr = new Bitmap(*C->Trans);
	BaseBitmap *tr = new BaseBitmap(TransCount);
	cnt = 0;
	for (int cc = 0; cc < TransCount; cc++) {
		
		if (C->Trans->CheckPosition(cc, CountCheckPosition) > 0) {
			tr->FillEmptyPosition(cc);
			cnt++;
		}
		
	}
	
	C->Trans->_count = cnt;
	tr->_count=cnt;
	MFT.push_back(tr);
	
	// update the end of the relevant
	C->rEnd++;
	
	// Update stat variables
	MFIDepth += C->Name->_count;
	MFISize++;
	
	SupportCountList.push_back(C->Trans->_count);
	if (name->_count > maxItemsetSize)
		maxItemsetSize = name->_count;
	
	
}

/*
 modified function
 */
/////////////////////////////////////////////////////////////////////
/// Add this node's name bitmap to the FCI list
///
/// @param C                 the current node
/////////////////////////////////////////////////////////////////////

void AddToFCI(TreeNode *C) {
	// Search HT
	HashTable::iterator h = HT.find(C->Trans->_count);
	
	// If the support of node C is NOT in HashSup
	if (h == HT.end()) {
		// Add a new list to the HashSup
		ItemSet* newList = new ItemSet();
		newList->reserve(500);
		newList->push_back(MFI.size());
		HT.insert(HashTable::value_type(C->Trans->_count, newList));
		
		// Else add pointer to last item in iName to HT entry
	} else {
		for (ItemSet::reverse_iterator goli = (*h).second->rbegin(); goli
			 != (*h).second->rend(); goli++)
			if (MFI[*goli]->Superset(C->Name))
				return;
		
		// search the table
		(*h).second->push_back(MFI.size());
	}
	
	// copy name bitmap
	BaseBitmap *name = new BaseBitmap(*C->Name);
	
	BaseBitmap *tr = new Bitmap(*C->Trans);
	// add to MFI
	MFI.push_back(name);
	MFT.push_back(tr);
	
	// Update stat variables
	MFIDepth += C->Name->_count;
	MFISize++;
	
	SupportCountList.push_back(C->Trans->_count);
	
	if (name->_count > maxItemsetSize)
		maxItemsetSize = name->_count;
}

/*
 modified function
 */
int SortLMFI(int rBegin, int rEnd, int sortBy) {
	int left = rBegin;
	int right = rEnd - 1;
	while (left <= right) {
		while (left < rEnd && !MFI[left]->CheckPosition(sortBy,
														CountCheckPosition))
			left++;
		while (right >= rBegin && MFI[right]->CheckPosition(sortBy,
															CountCheckPosition))
			right--;
		if (left < right) {
			// we are now at a point where MFI[left] is relevant
			// and MFI[right] is not since left < right, we swap the two
			BaseBitmap* tempBitmap = MFI[left];
			MFI[left] = MFI[right];
			MFI[right] = tempBitmap;
			tempBitmap = NULL;
			
			BaseBitmap* tempBit = MFT[left];
			MFT[left] = MFT[right];
			MFT[right] = tempBit;
			tempBit = NULL;
			
			int tempSupport = SupportCountList[left];
			SupportCountList[left] = SupportCountList[right];
			SupportCountList[right] = tempSupport;
			
			left++;
			right--;
		}
	}
	
	// the first relevant one for the next node is left
	return left;
}

/////////////////////////////////////////////////////////////////////
/// Determine whether a HUTMFI is true.
///     - if HUT is in MFI, then HUT is frequent
///     and the subtree rooted at this node can be pruned
///
/// @return True if HUT is in the MFI
/////////////////////////////////////////////////////////////////////

bool CheckHUTMFI(TreeNode *C, int iTAIL) {
	
	// for each element i in the tail form {head} U {i} and check for a
	//     superset in the MFI
	int rBegin = C->rBegin;
	int rEnd = C->rEnd;
	for (; iTAIL < C->tEnd; iTAIL++) {
		rBegin = SortLMFI(rBegin, rEnd, F1[gTail[iTAIL].Item]->Prefix);
		
		if (rEnd <= rBegin)
			return false;
	}
	
	return true;
}

/////////////////////////////////////////////////////////////////////
/// Dynamically reorder the elements in the tail by increasing support
///    - Expand all children and sort by increasing support
///    - Remove infrequent children
///
/// @param C                 current node
/// @param iTAIL             index into tail of current node
/// @param useComp           whether compressed bitmaps should be used
/// @param NoChild           whether C has any frequent children
/// @param AllFreq           whether all children are frequent
/////////////////////////////////////////////////////////////////////

void ReorderTail(TreeNode *C, int &iTAIL, bool useComp, bool &NoChild,
				 bool &AllFreq) {
	
	int tailIndex = 0;
	int lol = 0;
	// for each tail element
	for (lol = iTAIL; lol < C->tEnd; lol++) {
		Bitmap *trans = TransBuffy[C->Depth];
		int theCount = 0;
		
		// Compress the bitmaps
		if (useComp && (F1[gTail[lol].Item]->compID != C->compID)) {
			F1[gTail[lol].Item]->Trans->BuildRelComp(*TransBuffy[projectDepth]);
			F1[gTail[lol].Item]->compID = C->compID;
			CountRebuilds++;
		}
		
		// Use the compressed bitmaps
		if (useComp) {
			// AND the compressed bitmaps and count the result
			trans->AndCompOnly(*C->Trans, *F1[gTail[lol].Item]->Trans,
							   CountSmallAnds);
			theCount = trans->SmallCount(CountCounts);
			
			// use the full bitmaps
		} else {
			// AND & count the bitmaps
			trans->AndOnly(*C->Trans, *F1[gTail[lol].Item]->Trans, CountAnds);
			theCount = trans->Count(CountSmallAnds);
		}
		
		// If the results is frequent
		if (theCount >= MS) {
			// if PEP pruning holds
			if (PEPrune && (trans->_count == C->Trans->_count)) {
				// Move tail element from tail to head
				C->Name->Or(*C->Name, *F1[gTail[lol].Item]->Name);
				CountPEPrunes++;
				
				// add tail element to reordered tail
			} else {
				NoChild = false;
				
				// create new tail element
				TailBuffy[tailIndex].Count = theCount;
				TailBuffy[tailIndex].Item = gTail[lol].Item;
				tailIndex++;
			}
		} else
			AllFreq = false;
	}
	
	sort(TailBuffy, TailBuffy + tailIndex);
	
	// Set the begin and end values of the new tail
	iTAIL = C->tEnd;
	C->tEnd = iTAIL + tailIndex;
	C->tBegin = iTAIL;
	int r = 0;
	
	// Copy new tail into next slots
	for (lol = iTAIL; lol < C->tEnd; lol++) {
		gTail[lol] = TailBuffy[r];
		r++;
	}
}

/////////////////////////////////////////////////////////////////////
/// Simply copy over the tail without expanding any of the children
///    for pure DFS (no expansion of all children)
///
/// @param C                 current node
/// @param iTAIL             index into tail of current node
/////////////////////////////////////////////////////////////////////

void NoorderTail(TreeNode *C, int &iTAIL) {
	
	// set begin and end tail pointers
	iTAIL = C->tEnd;
	C->tEnd = C->tEnd - C->tBegin + C->tEnd;
	C->tBegin = iTAIL;
	int r = 0;
	
	// copy over old tail to new tail
	for (int lol = iTAIL; lol < C->tEnd; lol++) {
		gTail[lol] = gTail[lol - C->tEnd + C->tBegin];
		r++;
	}
}

/////////////////////////////////////////////////////////////////////
/// The main MAFIA algorithm function
///
/// @param C                 the current node
/// @param HUT               whether this is a HUT check (left most branch)
/// @param FHUT              [output] whether the HUT is frequent
/// @param useComp           if compression has been switched on
/////////////////////////////////////////////////////////////////////

void MAFIA(TreeNode *C, bool HUT, bool &FHUT, bool useComp) {
	
	CountNodes++;
	int iTAIL = C->tBegin; // index into the tail
	bool NoChild = true; // whether the node has any frequent children
	bool AllFreq = true; // whether all the children are frequent
	FHUT = false; // whether this is a FHUT
	int beforeCountNodes = CountNodes;
	int frequentTailSize = C->tEnd - iTAIL;
	
	if (iTAIL < C->tEnd) {
		
		if (C != Root) {
			if (Reorder) {
				ReorderTail(C, iTAIL, useComp, NoChild, AllFreq);
			} else {
				NoorderTail(C, iTAIL);
			}
		}
		
		frequentTailSize = C->tEnd - iTAIL;
		
		int estimateTail = (int) (frequentTailSize / (double) EstimateDiv);
		if (estimateTail > 0 && C->Trans->_count != TransCount) {
			double estimateSubTree = EstimateBuffy[estimateTail].Sum
			/ (double) EstimateBuffy[estimateTail].Count;
			double support = C->Trans->_count / (double) TransCount;
			double factor = 11.597 - 29.914 * (support - .52392) * (support
																	- .52392);
			double cost = abs(factor * frequentTailSize / (1 - 1.2 * support));
			//double cost = 5 * frequentTailSize / (1 - support);
			
			// check if relative comp should be performed
			if ((!useComp) && (estimateSubTree > cost)) {
				
				// build the relative bitmap for source [node bitmap] (all 1's)
				C->Trans->BuildSource();
				
				// remember the depth of the FULL bitmap your projecting from
				projectDepth = C->Depth - 1;
				
				// increment the ID
				C->compID = MAX_compID;
				MAX_compID++;
				useComp = true;
			}
		}
		
		// Candidate generation - extend the Head with the tail elements
		// We start from the end of the tail and move backwards
		// Therefore the tail is iterated through in increasing support,
		// but is stored in decreasing support.
		while (iTAIL < C->tEnd) {
			// form a one-extension
			Bitmap *trans = TransBuffy[C->Depth];
			BaseBitmap *name = NameBuffy[C->Depth];
			TreeNode *newNode = NodeBuffy[C->Depth];
			
			// create name for the new node
			name->Or(*C->Name, *F1[gTail[iTAIL].Item]->Name);
			
			// compress the bitmaps if warranted
			if (useComp && (F1[gTail[iTAIL].Item]->compID != C->compID)) {
				// build the relative for this node
				F1[gTail[iTAIL].Item]->Trans->BuildRelComp(
														   *TransBuffy[projectDepth]);
				F1[gTail[iTAIL].Item]->compID = C->compID;
				CountRebuilds++;
			}
			
			int theCount = 0;
			
			// use the compressed bitmaps for ANDing and counting
			if (useComp) {
				// AND and count small bitmaps
				trans->AndCompOnly(*C->Trans, *F1[gTail[iTAIL].Item]->Trans,
								   CountSmallAnds);
				if (Reorder)
					trans->_count = gTail[iTAIL].Count;
				else
					theCount = trans->SmallCount(CountCounts);
			} else {
				// AND and count the full bitmaps
				trans->AndOnly(*C->Trans, *F1[gTail[iTAIL].Item]->Trans,
							   CountAnds);
				
				if (Reorder)
					trans->_count = gTail[iTAIL].Count;
				else
					theCount = trans->Count(CountSmallAnds);
			}
			if (!Reorder && PEPrune && (theCount == C->Trans->_count)) {
				CountPEPrunes++;
				C->Name->Or(*C->Name, *F1[gTail[iTAIL].Item]->Name);
				iTAIL++;
				continue;
			}
			
			// Determine whether this candidate will be a HUT
			// Conceptually the leftmost branch of the tree is a HUT check
			if ((iTAIL != C->tBegin) && (C != Root))
				HUT = 0;
			else
				HUT = 1;
			
			if (!AllFreq)
				HUT = 0;
			
			if (Reorder || (theCount >= MS)) {
				// form the 1-extension node
				newNode->setTreeNode(name, trans, C->Depth + 1, C->compID,
									 F1[gTail[iTAIL].Item]->Prefix, iTAIL + 1, C->tEnd);
				
				// setup the LMFI for the next level; it contains all
				// itemsets in LMFI for this level that also include the
				// one we're extending the node with.  We do sort of a
				// quicksort thing to move the relevant itemsets for the
				// next node to the end of the portion of the MFI relevant
				// to this node
				
				newNode->rEnd = C->rEnd;
				newNode->rBegin = SortLMFI(C->rBegin, C->rEnd, newNode->Prefix);
				
				// Check for HUT in MFI for remaining tail
				if (HUTMFI && newNode->tBegin != newNode->tEnd && !HUT) {
					CountHUTMFI++;
					
					if (CheckHUTMFI(newNode, newNode->tBegin)) {
						// stop generation of extensions
						CountHUTMFISuccess++;
						AllFreq = false;
						break;
					}
				}
				
				NoChild = false;
				
				// recurse down the tree
				MAFIA(newNode, HUT, FHUT, useComp);
				
				// Add those discovered from lower levels to the current LMFI
				// LMFI_l = LMFI_l \union LMFI_{l+1}
				// all we need to do is to update the end pointer
				C->rEnd = newNode->rEnd;
			} else
				AllFreq = false;
			
			// if this was a successful HUT check
			if (FHUT) {
				// keep going up the tree
				if (HUT) {
					return;
					
					// reached start of HUT, so stop generation of subtree
					// rooted at this node
				} else {
					FHUT = false;
					break;
				}
			}
			
			// Move on the next tail element
			iTAIL++;
		}
	}
	
	// if this is a FHUT
	if (GoFHUT && HUT && AllFreq) {
		FHUT = true;
		CountFHUT++;
	}
	
	// if this node is childless and not in MFI
	if (NoChild && !LMFISuperSet(C)) {
		AddToMFI(C);
	}
	
	int subtreeSize = CountNodes - beforeCountNodes + 1;
	int estimateTail = (int) (frequentTailSize / (double) EstimateDiv);
	if (estimateTail > 0 && C->Trans->_count != TransCount) {
		EstimateBuffy[estimateTail].Count++;
		EstimateBuffy[estimateTail].Sum += subtreeSize;
	}
}

/////////////////////////////////////////////////////////////////////
/// Merge repeated itemsets into one combined itemset
///    - e.g. if (transaction set of item 4) AND
///    (transaction set item 5) = (transaction set item 5),
///    then item 5 is a duplicate of item 4
///    due to increasing support
/////////////////////////////////////////////////////////////////////

void MergeRepeatedItemsets() {
	NodeList::iterator bali = F1.begin();
	Bitmap out(*(*bali)->Trans);
	
	// for each frequent 1-itemset
	while (bali != F1.end()) {
		out.FillOnes();
		NodeList::iterator noli = bali;
		noli++;
		
		// search for a copy of the itemset's transaction set
		while (noli != F1.end()) {
			// stop when count is no longer the same
			if ((*bali)->Trans->_count != (*noli)->Trans->_count)
				break;
			else {
				// AND itemsets with the same count
				out.AndOnly(*(*bali)->Trans, *(*noli)->Trans, CountAnds);
				out.Count(CountSmallAnds);
				
				// check for a duplicate
				if (out._count == (*noli)->Trans->_count) {
					(*bali)->Name->Or(*(*bali)->Name, *(*noli)->Name);
					noli=F1.erase(noli);
				} else
					noli++;
			}
		}
		bali++;
	}
	
}

string getFileName(const string& s) {
	
	char sep = '/';
	
#ifdef _WIN32
	sep = '\\';
#endif
	
	size_t i = s.rfind(sep, s.length( ));
	if (i != string::npos) {
		return(s.substr(i+1, s.length( ) - i));
	}
	
	return("");
}






/*
 main function of DeBi
 */
/** @} */
/// @defgroup MainFunction Main Function
/// Main function for running the program
/** @{ */
/////////////////////////////////////////////////////////////////////
// main function for the program
/////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	
	// default paramaters
	int u = 0;
	string outFileNameS;
	overlap = 0.5;
	char **optarg = NULL;
	char *s; /* to traverse the options */
	pattern = 'u';
	type='b';
	seedLength = 3;
	level = 2;
	noise= -1;
	//outFilename ="./";
	// Check parameters
	if (argc > 3) { /* if arguments are given */
		
		fprintf(stderr, "\n%s - %s\n", argv[0], DESCRIPTION);
		
	}/* print a startup message */
	
	else { /* if no arguments given */
		printf("usage: %s [options] infile outfile \n", argv[0]);
		
		printf("-s#	minimum support of an item set"
			   "(default: %d) \n", MSF);
		
		printf("-o#	percentage of maximum overlap"
			   "(default: %g)\n", overlap);
		
		
		printf("-b#	binarization level"
			   "(default: %g)\n", level);
		
		printf("-p#	regulation pattern"
			   "(default: %c) \n", pattern);
		
		exit(1);
		
		
	}
	
	for (int i = 1; i < argc; i++) { /* traverse arguments */
		s = argv[i]; /* get option argument */
		if (optarg) {
			*optarg = s;
			optarg = NULL;
			continue;
		}
		if ((*s == '-') && *++s) { /* -- if argument is an option */
			while (*s) { /* traverse options */
				switch (*s++) { /* evaluate switches */
						break;
						
					case 'o':
						overlap =strtod(s, &s);
						break;
					case 'a':
						alpha = strtod(s, &s);
						break;
					case 'p':
						pattern = *s;
						if(pattern=='u'){
							outFileNameS="up";
						}
						else{
							outFileNameS="down";
						}
						break;
					case 'b':
						level = strtod(s, &s);
						type='r';
						break;
						
						
				}
				if (optarg && *s) {
					*optarg = s;
					optarg = NULL;
					break;
				}
			}
		} 
		
	}
	time(&total_start);

	// if the data is binarized, no need to define the pattern and discretization threshold parameters
	if(type=='b'){
		pattern='u';
		level=1;
	}

	outFilename = argv[2];
	Counters = new int[MAX_NUM_ITEMS];  //to count the support of items
	
	F1FromFile(argv[1], pattern);
	string path=argv[1];
	string inFileName=getFileName(path);
	
	time(&algorithm_start);
	clock_t start;
	
	start = clock();
	
	//minumum sample size is 5
	if(MSF<5)
	{
		MSF=5;
	}
	
	printf( "\n----Running DeBi algorithm with the following parameters----" );
	
	printf("\nminimum sample numbers in biclusters: "
		   "%d \n", MSF);
	
	printf("percentage of maximum overlap: "
		   "%g \n", overlap);
	
	printf("binarization level: "
		   "%g \n", level);
	
	printf("regulation pattern: "
		   "%s \n", outFileNameS.c_str());
	
	
	printf( "\n\nFinding seed biclusters...");
	
	vector<int> deleteItems; 
	bool check = 1;
	while (MS >= MSF & check) {
		
		
		// Set size of F1
		FullF1size = F1size = F1.size();
		// if F1 is not empty
		if (FullF1size != 0) {
			
			// Create a null node to begin the DFS tree
			Bitmap *NullTrans = new Bitmap(TransCount);
			NullTrans->FillOnes();
			NullTrans->_count = TransCount;
			ItemSet NullList;
			
			BaseBitmap *NullName = new BaseBitmap(FullF1size);
			//	MFI.reserve(100000);
			int p = 0;
			
			ItemMap = new int[FullF1size];
			ItemsetBuffy = new int[FullF1size];
			// Rename items in F1
			for (NodeList::iterator nli = F1.begin(); nli != F1.end(); nli++) {
				
				// store old itemid
				ItemMap[p] = (*nli)->Prefix;
				
				// assign new itemid
				(*nli)->Prefix = p;
				
				// assign name bitmaps
				
				(*nli)->Name = new BaseBitmap(FullF1size);
				(*nli)->Name->FillEmptyPosition(p);
				(*nli)->Name->_count = 1;
				p++;
			}
			
			MergeRepeatedItemsets();
			F1size = F1.size();			// Create global tail
			maxtail = F1size * (F1size + 1) / 2;
			gTail = new TailElement[maxtail];
			
			// Create buffer for sorting
			TailBuffy = new TailElement[F1size];
			
			
			// Create buffer for estimating size of each subtree
			EstimateSize = (int) ceil(F1size / (double) EstimateDiv);
			EstimateBuffy = new SubtreeEstimate[EstimateSize];
			for (int estimateIndex = 0; estimateIndex < EstimateSize; estimateIndex++) {
				EstimateBuffy[estimateIndex].Count = 1;
				EstimateBuffy[estimateIndex].Sum = estimateIndex * EstimateDiv
				* estimateIndex * EstimateDiv / 2;
			}
			
			// Initialize global tail
			int uu;
			for (uu = 0; uu < maxtail; uu++) {
				gTail[uu].Item = -1;
				gTail[uu].Count = 0;
			}
			
			// Fill global tail
			for (uu = 0; uu < F1size; uu++) {
				gTail[uu].Item = uu;
				
				gTail[uu].Count = F1[uu]->Trans->_count;
				
				// assign tail index
				F1[uu]->tBegin = uu + 1;
				if (uu == F1size - 1)
					F1[uu]->tBegin = -1;
				
				// add a buffer element for each item in F1
				BaseBitmap *name = new BaseBitmap(FullF1size);
				NameBuffy.push_back(name);			
				
				Bitmap *buff = new Bitmap(TransCount);
				TransBuffy.push_back(buff);			
				
				TreeNode *newNode = new TreeNode();
				NodeBuffy.push_back(newNode);
				
				
			}
			
			srand(666);
			bool FHUT;
			
			// create root node and its associated tail
			Root = new TreeNode(NullName, NullTrans, 0, 0, -1, 0, F1size);
			
			//Nothing is in MFI, so nothing is relevant
			Root->rBegin = 0;
			Root->rEnd = 0;
			GoFHUT = true; // FHUT flag
			HUTMFI = true; // HUTMFI flag
			PEPrune = true; // PEPrune flag
			Reorder = true; // Reorder flag
			MAFIA(Root, true, FHUT, false);
			
			
			if (MFISize != 0) {
				
				for (int i = 0; i < MFISize; i++) {
					
					if (MFI[i]->_count >= seedLength && FullF1size != 0 && MFT[i]->_count >= MS) {
						
						
						Bitmap* TID = new Bitmap(TransCount);	
						
						for (int cc = 0; cc < FullF1size; cc++) {
							
							if (MFI[i]->CheckPosition(cc, CountCheckPosition)> 0) {
								deleteItems.push_back(ItemMap[cc]);
								
							}
							
						}
						int d=0;
						for (int cc = 0; cc < TransCount; cc++) {
							if (MFT[i]->CheckPosition(cc, CountCheckPosition) > 0) {
								TID->FillEmptyPosition(cc);
								d++;
							}
						}
						
						// insert the the found sets
						TID->_count = MFT[i]->_count;
						MFT1.push_back(TID);
						FinalMFISize++;
						
						
					}
					
					
				}
				
			}
			
			
			if ((u + seedLength) > ItemCount | MS < MSF) {
				
				check = 0;
				continue;
			}
			
			delete NullTrans;
			delete NullName;
			destroy();
			UpdateF1(deleteItems, argv[1]);
	
		} 
		else{
			FullF1size = 0;
			UpdateF1(deleteItems, argv[1]);
		}
		
		
		deleteItems.clear();
		
	}
	
	printf( "\n\n");
	printf("Seed sets are found...\n");
	printf("Number of seed sets are: %.3d\n\n", FinalMFISize);
	
	delete[] Counters;
	
	ItemList.clear();
	alpha=0.01;
	
	printf("Searching for the optimum alpha value...");
	
	//find the interval for searching the optimum alpha 
	optimizeAlpha(1,10,10);
	max_pscore=0;
	
	//minumum and maximum intervals for searching optimum alpha 
	int min_range=(max_alpha*10)-9;
	int max_range=(max_alpha*10)+9;
	if(max_alpha==1)
	{
		if(TransCount>1000) //if the column number is more than 1000 then the optimum alpha is smaller than 1e-10
			min_range=10;
		else
			min_range=2;
	}
	
	//find the optimum alpha in the intervals
	optimizeAlpha(min_range,max_range+1,1);
	printf("\nAlpha value is optimized...");
	
	alpha=pow(double(10),double(max_alpha*-1));
	string overlap_val = double2string(overlap);
	string level_val = double2string(level);
	string val = double2string(alpha);
	string outFilename2=new char[1000];
	outFilename2=string(outFilename)+"/"+inFileName+"."+outFileNameS+".overlap"+overlap_val+".threshold"+level_val+".alpha"+val+".txt";
	outFile = new ItemsetOutput((char *)outFilename2.c_str());
	
	printf("\nPrinting the biclusters...");
	
	if(max_pscore>0)
	{
		StatTestF1Simple();
		Filter();
		PrintMFI();
	}
	
	
  	printf("\nSignficance Score of the Bicluster: " "%g \n", max_pscore);
  	printf("\nTotal number of biclusters: " "%i \n", max_setSize);
	
	
	time(&total_finish);
	total_time = difftime(total_finish, total_start);
	printf("Total time:            %.2f seconds.\n\n", total_time);
	return 0;
}
