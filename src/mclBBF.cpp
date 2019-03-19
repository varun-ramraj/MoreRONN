/*

   The script generates 10 "estimate.rec" files (one per cross-validation model); in each we have the prediction for a sequence under testing using one trained cross-validation model

*/


#include "mclBBF.h"

using namespace std;


//for the weight function
//Then, if the optimizable parameter is x, 0 < x < 1;
//weight = 1 - x * cos (P)
//x = 0, all weights are 1
//x = 1, weights are 0 at each end and 2 in the middle

//RESULT: This function does not improve scores, it only
//makes the peaks and troughs sharper. Set to 0.0 so weights
//are even.
double OPT_PARAM = 0.5;

vector<string> query;
vector<string> headers;

double		disorder_weight;

//holds query's windows
//note that index in vector == start idx of window
//in original query sequence
//vector<string> query_seq_windowed;

//moreronn data holders
//map <int, vector<double> > clusters_weights;
vector<string> all_seqs; //all sequences, indexed in vector
vector<int> all_xvals; //all cross-validation indices
vector< vector<int> > iprotres; //BLOSUM indices of prototypes

//stores references to seq indices in all_seqs vector
map <int, map< int, vector<int> > > clusters_arrays_db_seqs;
vector<vector<int> > clust_to_seqs;
vector <vector<double> > clusters_weights;

vector<double> error_terms;

map<int, vector<double> > xval_pdfs;

int		nD;//number of database sequences
int		nW; //window length


// A b C D E F G H I j K L M  N  o P  Q  R  S  T  u V  W  x Y
int		INDEX[25]={0,0,1,2,3,4,5,6,7,0,8,9,10,11,0,12,13,14,15,16,0,17,18,0,19};

int		Dayhoff[20][20]={
    { 40, 24, 32, 32, 16, 36, 28, 28, 28, 24, 28, 32, 36, 32, 24, 36, 36, 32,  8, 20},
    { 24, 80, 12, 12, 16, 20, 20, 24, 12,  8, 12, 16, 20, 12, 16, 32, 24, 24,  0, 32},
    { 32, 12, 48, 44,  8, 36, 36, 24, 32, 16, 20, 40, 28, 40, 28, 32, 32, 24,  4, 16},
    { 32, 12, 44, 48, 12, 32, 36, 24, 32, 20, 24, 36, 28, 40, 28, 32, 32, 24,  4, 16},
    { 16, 16,  8, 12, 68, 12, 24, 36, 12, 40, 32, 16, 12, 12, 16, 20, 20, 28, 32, 60},
    { 36, 20, 36, 32, 12, 52, 24, 20, 24, 16, 20, 32, 28, 28, 20, 36, 32, 28,  4, 12},
    { 28, 20, 36, 36, 24, 24, 56, 24, 32, 24, 24, 40, 32, 44, 40, 28, 28, 24, 20, 32},
    { 28, 24, 24, 24, 36, 20, 24, 52, 24, 40, 40, 24, 24, 24, 24, 28, 32, 48, 12, 28},
    { 28, 12, 32, 32, 12, 24, 32, 24, 52, 20, 32, 36, 28, 36, 44, 32, 32, 24, 20, 16},
    { 24,  8, 16, 20, 40, 16, 24, 40, 20, 56, 48, 20, 20, 24, 20, 20, 24, 40, 24, 28},
    { 28, 12, 20, 24, 32, 20, 24, 40, 32, 48, 56, 24, 24, 28, 32, 24, 28, 40, 16, 24},
    { 32, 16, 40, 36, 16, 32, 40, 24, 36, 20, 24, 40, 28, 36, 32, 36, 32, 24, 16, 24},
    { 36, 20, 28, 28, 12, 28, 32, 24, 28, 20, 24, 28, 56, 32, 32, 36, 32, 28,  8, 12},
    { 32, 12, 40, 40, 12, 28, 44, 24, 36, 24, 28, 36, 32, 48, 36, 28, 28, 24, 12, 16},
    { 24, 16, 28, 28, 16, 20, 40, 24, 44, 20, 32, 32, 32, 36, 56, 32, 28, 24, 40, 16},
    { 36, 32, 32, 32, 20, 36, 28, 28, 32, 20, 24, 36, 36, 28, 32, 40, 36, 28, 24, 20},
    { 36, 24, 32, 32, 20, 32, 28, 32, 32, 24, 28, 32, 32, 28, 28, 36, 44, 32, 12, 20},
    { 32, 24, 24, 24, 28, 28, 24, 48, 24, 40, 40, 24, 28, 24, 24, 28, 32, 48,  8, 24},
    {  8,  0,  4,  4, 32,  4, 20, 12, 20, 24, 16, 16,  8, 12, 40, 24, 12,  8,100, 32},
    { 20, 32, 16, 16, 60, 12, 32, 28, 16, 28, 24, 24, 12, 16, 16, 20, 20, 24, 32, 72}};
int Blosum62[20][20]={
    { 4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-2,-3,-2},
    { 0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2},
    {-2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0, 1,-3,-4,-3},
    {-1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0, 0,-3,-3,-2},
    {-2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3},
    { 0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3,-2,-2,-2,-2, 0, 1, 0,-2,-3},
    {-2,-3, 1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1, 0,-2,-2, 2},
    {-1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-2, 1,-3,-1},
    {-1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0, 0,-3,-3,-2},
    {-1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-2, 3,-2,-1},
    {-1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1,-2,-1,-1},
    {-2,-3, 1, 0,-3, 0,-1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2},
    {-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-1, 7,-1,-2,-1, 1,-2,-4,-3},
    {-1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0, 0,-2,-2,-1},
    {-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2},
    { 1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2},
    {-1,-1, 1, 0,-2, 1, 0,-2, 0,-2,-1, 0, 1, 0,-1, 1, 4,-2,-3,-2},
    { 0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2,-2, 4,-3,-1},
    {-3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-3,-3,11, 2},
    {-2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7}};

//load MoreRONN clusters, weights , error terms and PDFs
//from the single database file
int load_moreronn_data(string datafile)
{
    stringstream *convert;

    string line; //line buffer
    int lines = 0; //line number counter

    int pdf_start_line = 0;

    string current_cluster = "";
    string current_pdf = "";

    int curclust = -1;

    int total_clust_count = 0;
    //find the first pdf header
    ifstream in1(datafile.c_str());

    if (in1.is_open())
    {
	while(in1.good())
	{
	    getline(in1, line);

	    if (line.length() > 0)
	    {
		if (line.compare("===PDF_0===") == 0)
		{
		    pdf_start_line = lines;
		}
		else if (line[0] == '>' && line[1] != 'E')
		{
		    total_clust_count++;
		}
	    }

	    lines++;
	}
    }
    else
    {
	cerr << "Check data file!" << endl;
	exit(-1);
    }

    in1.close();


    //reset line counter
    lines = 0;

    //initialise cluster_weights with known cluster counts
    clusters_weights.resize(total_clust_count);

    //now read in data again
    ifstream in(datafile.c_str());

    if (in.is_open())
    {
	while (in.good())
	{
	    getline(in, line);

	    if (line.length() > 0)
	    {
		//get the pdfs loaded in
		if (lines >= pdf_start_line)
		{
		    if (line[0] == '=')
		    {
			current_pdf = line;

		    }
		    else
		    {	

			int intid = atoi(current_pdf.substr(7,1).c_str());

			double pdfval = strtod(line.c_str(), NULL);

			xval_pdfs[intid].push_back(pdfval);
		    }

		}

		else //cluster and error data
		{
		    if (line[0] == '>' && line[1] != 'E') //cluster header, not error term
		    {

			//get weights
			vector<string> toptkn;
			split(line, toptkn, "|");

			current_cluster = toptkn[0];
			curclust++;

			string weights = toptkn[1];

			vector<string> weighttkns;

			split(weights, weighttkns, " ");

			for (int k = 0; k < weighttkns.size(); k++)
			{
			    double scr = strtod(weighttkns.at(k).c_str(), NULL);


			    clusters_weights[curclust].push_back(scr);
			}
		    }

		    else if (line[0] == '>' && line[1] == 'E')
		    {
			//load error term

			//get error terms
			vector<string> toptkn;
			split(line, toptkn, "|");

			string error_header = toptkn[0];

			string errors = toptkn[1];

			vector<string> errtkns;

			split(errors, errtkns, " ");

			for (int k = 0; k < errtkns.size(); k++)
			{
			    double scr = strtod(errtkns.at(k).c_str(), NULL);
			    error_terms.push_back(scr);
			}


		    }

		    else //populate cluster
		    {

			//first, tokenise each line
			vector<string> tkns;

			split(line, tkns, "\t");

			int tmp_xidx = atoi(tkns[1].c_str());

			//set window size
			//TODO: maybe DON'T hammer this poor variable
			//for every line in the file??!
			nW = tkns[0].length();

			//add to data structures
			all_seqs.push_back(tkns[0]);
			all_xvals.push_back(tmp_xidx);
			clusters_arrays_db_seqs[curclust][tmp_xidx].push_back(all_seqs.size() - 1);

		    }

		}

		lines++;	
	    }
	}
    }


    in.close();

    delete convert;

    //stores all sequence indices in a cluster
    clust_to_seqs.resize(clusters_arrays_db_seqs.size());

    //populate clust_to_seqs setup
    for (int yy = 0; yy < clusters_arrays_db_seqs.size(); yy++)
    {
	map<int, vector<int> >::iterator it;

	for (it = clusters_arrays_db_seqs[yy].begin(); it != clusters_arrays_db_seqs[yy].end(); it++)
	{
	    map<int, vector<string> > tmpmap;

	    //key...check that the cross-validation index
	    //is NOT the same as the current run (loop counter)
	    for (int k = 0; k < (it->second).size(); k++)
	    {
		int addidx = (it->second).at(k);

		clust_to_seqs[yy].push_back(addidx);
	    }

	}
    }

    //we can clear this data structure
    //because it's not needed in the big loop
    clusters_arrays_db_seqs.clear();

    //load prototype indices
    iprotres.resize(all_seqs.size());


    for (int i = 0; i < all_seqs.size(); i++)
    {
	for (int j = 0; j < all_seqs[i].length(); j++)
	{
	    iprotres[i].push_back(INDEX[(int)(all_seqs[i][j] - 'A')]);
	}
    }


    fprintf(stderr, "MoreRONN database loaded\n\n");

    fflush(stderr);

    return 0;
}

//weighting function
//currently unused, returns 1.0 when OPT_PARAM is 0.0
//which it is.
double weight(int win, int res)
{

    double P = (2.0 * M_PI) * ((double)res / (double)(win-1));

    //if global var OPT_PARAM is 0, all weights will be 1
    //reducing to original unweighted case

    return (1.0 - (OPT_PARAM * cos(P)));

}

//takes in query sequence
//returns vector of per-residue predictions
vector<double> new_detect(string curquery)
{

    cerr << "Calculating per-residue disorder probability..." << endl;
    //length of query
    int qlen = curquery.length();


    //indices of the sequence and prototype residues, precalculated
    int iseqres[qlen];


    for (int i = 0; i < qlen; i++)
    {
	iseqres[i] = INDEX[(int)(curquery[i] - 'A')];
    }

    int numclusts = clusters_weights.size();

    //end data setup

    //initialise data structures to store weights and probabilities
    vector< vector<double> > sumWeight(10, vector<double>(qlen, 0.0));
    vector< vector<double> > sumWeightProb(10, vector<double>(qlen, 0.0));

    //perform alignments
    double y[10], fOrder, fDisor, pOrder, pDisor;

    //legacy, when OpenMP was required due to slow file loading and cumbersome
    //data structures. Left in to remember how to get it working if a future
    //developer so desires. It will NOT work with this version the way it's written

    //Also, if MoreRONN is to be run on a server, it's wise to leave it single-threaded
    //to allow concurrent access by multiple users. The footprint is small enough
    //and it runs quickly (1s per 250aa, approx.) so OpenMP is probably not required
    //#pragma omp parallel for private(y, fOrder, fDisor, pOrder, pDisor, pWeight)	

    //window string holder...we could speed this up by using a deque
    //as suggested here:
    //http://stackoverflow.com/questions/12148296/c-speeding-up-multiple-substr-or-equivalent-function-calls-for-parsing-of-a
    deque<char> tmpseq(curquery.begin(), curquery.begin() + nW);


    timeval t1, t2;
    double elapsedTime;

    // start timer
    gettimeofday(&t1, NULL);

    //loop over query windows: B	
    for(int i = 0; i <= qlen - nW; i++)
    {
	for (int q = 0; q < MODELS; q++)
	{
	    y[q] = 0.0;

	}

	double rho1[MODELS], rho0[MODELS];

	//loop over clusters: C (each iter 0.002ms)
	for(int j = 0; j < numclusts; j++)
	{

	    int best_prot_idx[MODELS];

	    for (int q = 0; q < MODELS; q++)
	    {
		rho1[q] = -10000.0;
		best_prot_idx[q] = -1;
	    }

	    //tmprep is now the cluster name
	    //string tmprep = vector_cluster_number_to_rep_seq.at(j);

	    int num_elements_in_cluster = clust_to_seqs[j].size();

	    //loop over elements in a cluster: D
	    for (int v = 0; v < num_elements_in_cluster; v++)
	    {

		//tmpdat is the INDEX OF the
		//prototype sequence being examined
		int tmpdat = clust_to_seqs[j][v];

		double score = 0.0;
		double scoreDD = 0.0;

		//go through each residue
		for(int w = 0; w < nW; w++)
		{
		    //score += Blosum62[iseqres[i+w]][INDEX[(int)(all_seqs[tmpdat][w] - 'A')]];			  		
		    score += Blosum62[iseqres[i+w]][iprotres[tmpdat][w]];
		    scoreDD += Blosum62[iprotres[tmpdat][w]][iprotres[tmpdat][w]]; 		
		}

		//score /= scoreDD;

		//loop over each cross-validation fold
		for (int q = 0; q < MODELS; q++)
		{
		    if (all_xvals[clust_to_seqs[j][v]] != q && score > rho1[q])
		    {
			rho1[q] = score;
			best_prot_idx[q] = tmpdat;
		    }
		}

	    }

	    //loop over each cross-validation fold
	    for (int q = 0; q < MODELS; q++)
	    {
		rho0[q] = 0.0;

		//best_prot_idx will remain at -1 (its initial value) 
		//if the cluster was empty for this fold
		if (best_prot_idx[q] >= 0)
		{
		    for (int win = 0; win < nW; win++)
		    {
			//rho0[q] += Blosum62[INDEX[(int)(all_seqs[best_prot_idx[q]][win] - 'A')]][INDEX[(int)(all_seqs[best_prot_idx[q]][win] - 'A')]];
			//lookup for indices is pre-calculated
			rho0[q] += Blosum62[iprotres[best_prot_idx[q]][win]][iprotres[best_prot_idx[q]][win]];

		    }

		    //set y score per fold

		    y[q] += clusters_weights[j][q] * exp ( CONST_ALPHA * ((rho1[q] - rho0[q]) / rho0[q]));


		}

		//add error term for this fold.. just once
		//do it at the last iteration
		if (j == numclusts - 1)
		{
		    y[q] += error_terms[q];
		}

	    }


	} //end C loop (over clusters)

	//loop over each cross-validation fold
	for (int q = 0; q < MODELS; q++)
	{

	    //retrieving PDFs from static arrays
	    double mu[2];
	    double sigma[2];

	    mu[0] = xval_pdfs[q][0];
	    mu[1] = xval_pdfs[q][1];
	    sigma[0] = xval_pdfs[q][2];
	    sigma[1] = xval_pdfs[q][3];


	    //$VR$: bug fixed by Ron in Feb07
	    fOrder=exp(-0.5*pow(y[q]-mu[0],2.0)/sigma[0])/(sqrt((2.0 * M_PI) * sigma[0])); 


	    fDisor=exp(-0.5*pow(y[q]-mu[1],2.0)/sigma[1])/(sqrt((2.0 * M_PI) * sigma[1]));


	    //$VR$: Bernoulli Model, I think?
	    pDisor=disorder_weight * fDisor / ( (1.0 - disorder_weight) * fOrder + disorder_weight * fDisor);


	    for (int r = 0; r < nW; r++)
	    {
		//new weight structures
		//get weight
		//double temp_weight = weight(nW, r);
		//sumWeightProb[i+r] += temp_weight * pDisor;
		//sumWeight[i+r] += temp_weight;

		//every score is weighted equally
		sumWeightProb[q][i+r] += pDisor;
		sumWeight[q][i+r] += 1.0;

	    }


	}

	//adjust the window loop if not last element
	if (i != qlen - 1)
	{
	    tmpseq.pop_front();
	    tmpseq.push_back(curquery[i+nW]);
	}

	//output a progress bar
	double fraction = (double)(i+1) / (double)(qlen-nW + 1);
	printProgBar(fraction * 100.0);

    } //end outermost loop (over query window)

    // stop timer
    gettimeofday(&t2, NULL);

    // compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms


    vector<vector<double> > X_per_fold_pred(qlen, vector<double>(10, 0.0));

    vector<double> per_fold_mean(qlen, 0.0);

    //calculate probabilities per fold
    for (int q = 0; q < MODELS; q++)
    {		
	for(int i = 0; i < qlen; i++)
	{			
	    X_per_fold_pred[i][q] = sumWeightProb[q][i] / sumWeight[q][i];			
	}		
    }

    //calculate average probabilities (to be written to output file)
    for(int r = 0; r < qlen; r++)
    {	
	for(int m = 0; m < MODELS; m++)
	{
	    per_fold_mean[r] += X_per_fold_pred[r][m];
	}

	per_fold_mean[r] /= (double)MODELS;

    }

    fprintf(stderr, "\nStats:\nQuery length: %d\nWindow length: %d\nNumber of windows: %d\nPrediction time (s): %fs\n\n", qlen, nW, qlen-nW+1, elapsedTime/1000);

    return per_fold_mean;

}

//write to stdout by default
int write_output(string curquery, string header, vector<double> scores, int printVerbosePredictionHeader)
{

    int qlen = curquery.length();

    //write out master disorder probability file
    FILE *fp;

    //open it with APPEND access, so it keeps writing
    //allows multiple sequences and scores to be written
    //fp=fopen(master_prob_output.c_str(),"a");

    //write out the header 
    fprintf(stdout, ">%s\n", header.c_str());
    fflush(stdout);
    

    //write out metrics if desired

    if (printVerbosePredictionHeader == 1) {

	int npl = 70; //number of residues per line
	int nr = 0;
	string pred = "";
	for (int i = 0; i < qlen; i++)
	{
	    //bug fix: if it's not a number, it's likely due to the sequence being shorter
	    //than window length, and it should be annotated differently!
	    if (isnan(scores[i])) { pred += "?"; }
	    else if (scores[i] > 0.6) { pred += "#"; }
	    else if (scores[i] > 0.5) { pred += "="; }
	    else if (scores[i] > 0.4) { pred += "-"; }
	    else { pred += " "; }
	}

	string lineout = ">";
	string scoreout = ">";
	for (int j = 0; j < qlen; j++)
	{
	    lineout += curquery[j];
	    scoreout += pred[j];
	    nr++;

	    if (nr == npl)
	    {
		fprintf(stdout, "%s\n", lineout.c_str());
		fprintf(stdout, "%s\n", scoreout.c_str());
		fprintf(stdout, ">$\n"); //spacer
		fflush(stdout);

		lineout = ">";
		scoreout = ">";
		nr = 0;
	    }
	}
	//last bits, should definitely be less than npl
	if (nr != 0)
	{
	    fprintf(stdout, "%s\n", lineout.c_str());
	    fprintf(stdout, "%s\n", scoreout.c_str());
	    fprintf(stdout, ">$\n"); //spacer
	    fflush(stdout);
	}

    }

    //now write per-residue scores
    for(int r = 0; r < qlen; r++)
    {

	fprintf(stdout,"%c\t%lf\n", curquery[r], scores[r]);

    }
    fflush(stdout);

    //don't close stdout...
    //fclose(fp);

    fprintf(stderr,"\nFinished. Wrote output file(s)...\n----------------\n");
    fflush(stdout);

    return 0;
}

//split on tokens
unsigned int split(string &txt, vector<string> &strs, string ch)
{
    strs.clear();

    char *pch;

    pch = strtok((char *)txt.c_str(), (const char *)ch.c_str());

    while (pch != NULL)
    {
	//printf ("%s -->",pch);
	strs.push_back((string)pch);
	pch = strtok (NULL, (const char *)ch.c_str());
    }


    return strs.size();
}


int callBBF_driver(vector<string> qry, vector<string> hdrs, string mod_fn, string pdf_fn1, double d_weight, string newdatafile, int printVerbosePredictionHeader)
{

    //store variables and load data
    disorder_weight = d_weight;

    query = qry;

    headers = hdrs;

    load_moreronn_data(newdatafile);

    int numqueries = query.size();

    cerr << "Analysing " << numqueries << " sequences..." << endl;

    for (int i = 0; i < numqueries; i++)
    {
	fprintf(stdout,"(Sequence %d of %d): %s\n", i+1, numqueries, headers.at(i).c_str());

	//perform prediction
	vector<double> curscores = new_detect(query.at(i));

	//write prediction to file
	write_output(query.at(i), headers.at(i), curscores, printVerbosePredictionHeader);
    }

    return 0;

}

//adapted from:
//http://nakkaya.com/2009/11/08/command-line-progress-bar/

//prints a pretty progress bar in the terminal
int printProgBar( int percent )
{

    string bar;

    for(int i = 0; i < 50; i++)
    {

	if( i < (percent/2))
	{
	    bar.replace(i, 1, "=");
	}

	else if ( i == (percent / 2))
	{
	    bar.replace(i, 1, ">");
	}

	else
	{
	    bar.replace(i, 1, " ");
	}

    }

    cerr << "\r" "[" << bar << "] ";
    cerr.width( 3 );
    cerr << percent << "%     " << flush;

    return 0;
}
