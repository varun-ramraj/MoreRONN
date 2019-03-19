/*=========================================================
 **
 ** This program is designed for running MoreRONNv4.9 for
 ** detecting disordered residues in proteins.
 ** Date:	March 2012
 ** Designer: Varun Ramraj (University of Oxford)
 ** Programmer:	Varun Ramraj (University of Oxford)
 ** Contact:	varun@strubi.ox.ac.uk

 ** Based on RONN, a disorder predictor
 ** designed by Zheng Rong Yang (Exeter University)
 ==========================================================*/

#include <unistd.h>
#include "mclBBF.h"
#define MORERONN_VERSION 4.9

using namespace std;
    
//DATA_PATH variable set in Makefile.am at compile-time
string DEFAULT_DB_PATH = (string)DATA_PATH + "/moreronn_master_data.dat";

//delete unnecessary spaces, line breaks et c.
//taken from:
//http://stackoverflow.com/questions/83439/remove-spaces-from-stdstring-in-c

//final solution taken from
//http://stackoverflow.com/a/14233269
struct RemoveDelimiter
{
    bool operator()(char c)
    {
	return (c =='\r' || c =='\t' || c == ' ' || c == '\n');
    }
};

void delUnnecessary(string &str)
{
    /*
       int size = str.length();
       for(int j = 0; j<=size; j++)
       {
       for(int i = 0; i <= j; i++)
       {
       if(str[i] == ' ' && str[i+1] == ' ')
       {
       str.erase(str.begin() + i);
       }
       else if(str[0]== ' ')
       {
       str.erase(str.begin());
       }
       else if(str[i] == '\0' && str[i-1]== ' ')
       {
       str.erase(str.end() - 1);
       }
       }
       }

       return str;
       */

    //str.erase(remove(str.begin(), str.end(), ' '), str.end());
    str.erase( std::remove_if( str.begin(), str.end(), RemoveDelimiter()), str.end());
}

int main(int argc, char *argv[])
{
    //number of models, characters and the format
    int nM=10;
    int nChar;

    //writing out disorder.prb
    vector< vector<double> > X_per_fold_pred; //stores each residue's prediction per fold
    vector<double> per_fold_mean; //stores the mean for each fold

    stringstream *convert;
    FILE *fp;


    //$VR$: VERY IMPORTANT. CONTROLS COST FUNCTION
    //keep it at 0.5 for now, given that it gives best
    //results for MCL-clustered data
    double disorder_weight = 0.5;

    string fModel, fPDF, str;

    //MoreRONN master data file
    //DATA_PATH variable set in Makefile.am at compile-time
    string fNewData = DEFAULT_DB_PATH;
    bool dbPathChanged = false;

    //read in query sequences and FASTA headers (if they exist)
    string query;
    vector<string> queries;
    vector<string> headers;

    char* input_filename;
    int is_stdin = -1;
    ifstream ifile;
    string testoutputstring;

    int printVerbosePredictionHeader = 1;
    //argument check
    int opt;

    while ((opt = getopt(argc, argv, "sp:f:w:d:")) != -1)
    {
	switch (opt)
	{
	    case 's': //standard input
		is_stdin = 1;
		break;

	    case 'f': //input file name
		input_filename = optarg;
		is_stdin = 0;
		break;

	    case 'w': //disorder weight
		disorder_weight = atof(optarg);
		break;

	    case 'd': //database file
		fNewData = (string)optarg;
		dbPathChanged = true;
		break;

	    case 'p': //include prediction header if value is 1 (default)
		printVerbosePredictionHeader = atoi(optarg);
		break;

	    default: //display help and usage
		is_stdin = -1;
		break;

	}
    }

    if (is_stdin == -1)
    {
	fprintf(stderr, "MoreRONN Disorder Predictor version: %.1f\n", MORERONN_VERSION); 
	fprintf(stderr, "Usage: %s [-f = filename -- FASTA-format OR -s = Standard Input] [-w = disorder_weight (optional) = default 0.5] [-d = database file] [-p0/-p1 = omit/include symbolic prediction header (default is to include)]\n", argv[0]); 
	exit(EXIT_FAILURE);
    }

    if (disorder_weight != 0.5)
    {
	fprintf(stderr, "New disorder weight: %f\n", disorder_weight);
    }

    if (dbPathChanged)
    {
	fprintf(stderr, "Loading non-standard database: %s\n", fNewData.c_str());
    }

    string curseq = "";
    string curheader = "";
    
    if (is_stdin == 1)
    {
	//cerr << "I can take stdin now." << endl;

	string line;
	//bug reported by Stephen Graham, using getline() 
	//will avoid breaking on spaces in the FASTA header
	//also removed all the old stringstream code since
	//it's redundant. Might as well just getline() directly
	//into the sequence and header vectors

	while (getline(cin, line, '\n'))
	{
	    if (line.length() > 0)
	    {
		if (line[0] == '>') //header
		{

		    if (curseq.length() > 0)
		    {
			headers.push_back(curheader);

			queries.push_back(curseq);
		    }

		    curheader = line.substr(1);


		    curseq = "";
		}
		else
		{
		    curseq += line;
		}
	    }
	}
	//cout << headers.size() << endl;
	//cout << queries.size() << endl;

	//cout << headers[headers.size() - 1] << endl;

	//push last one
	headers.push_back(curheader);
	queries.push_back(curseq);


	//cout << headers[headers.size() - 1] << endl;
	//cout << queries[headers.size() - 1] << endl;
    
    } 

    else if (is_stdin == 0) //take a file name
    {

	ifstream qfp(input_filename);

	if (qfp.is_open())
	{
	    while (qfp.good())
	    {
		getline(qfp, query);

		if (query.length() > 0)
		{

		    if (query[0] == '>') //header
		    {

			if (curseq.length() > 0)
			{
			    headers.push_back(curheader);

			    queries.push_back(curseq);
			}

			curheader = query.substr(1);


			curseq = "";
		    }

		    else
		    {
			curseq += query;

		    }

		}
	    }

	    qfp.close();
	}
	else
	{
	    cerr << "Can't open " << argv[1] << endl;
	    return -1;
	}


	//push last one
	headers.push_back(curheader);
	queries.push_back(curseq);

    }

    //cerr << curheader << endl << curseq << endl;

    //remove unnecessary spaces in the query sequences	
    for (int x = 0; x < queries.size(); x++)
    {
	delUnnecessary(queries.at(x));

    }


    //run 10-fold prediction and combine results
    callBBF_driver(queries, headers, fModel, fPDF, disorder_weight, fNewData, printVerbosePredictionHeader);


    return 0;
}


