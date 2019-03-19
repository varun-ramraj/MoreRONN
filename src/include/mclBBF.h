#ifndef __MCLBBF_H__
#define __MCLBBF_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <omp.h>
#include <sys/time.h>
#include <stdint.h>
#include <sys/timeb.h>
#include <deque>
#include <algorithm>

#define _USE_MATH_DEFINES //for M_PI
#include <math.h>

#define USE_MCL_CLUSTERING 1

using namespace std;

//how many cross-validation models do we have?
//currently statically set to 10
//TODO: make this dynamic??
static int MODELS = 10;


//alpha value for BBF
static double CONST_ALPHA = 4.0;

//experimentally determined
//this cutoff is for when alpha == 1
static double CLUSTER_WEIGHT_CUTOFF = 0.434;


//MoreRONN file loading
int load_moreronn_data(string datafile, string mString);

vector<double> new_detect(string curquery);

//takes as many strings as required
int callBBF_driver(vector<string> qry, vector<string> hdrs, string mod_fn, string pdf_fn1, double d_weight, string newdatafile, int printVerbosePredictionHeader = 1);

double weight(int win, int res);

int main(int argc,char **argv);

unsigned int split(string &txt, vector<string> &strs, string ch);


int printProgBar(int percent);

static string master_prob_output = "disorder.prb";

int write_output(string curquery, string header, vector<double> scores, int printVerbosePredictionHeader = 1);

#endif
