from cgi import parse_qs, escape
import json
import sys
import subprocess
import time
import os
import logging
from collections import OrderedDict

invalid_aa = ['B', 'Z', 'J', 'O', 'U', 'X']
res_limit = 3500 #don't exceed this!

STATIC_URL_PREFIX = '/html/'
STATIC_FILE_DIR = 'html/'

PREDICT_PREFIX = '/predict'

SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))
print SCRIPTDIR

MORERONN_BIN_DIR = os.path.join(SCRIPTDIR, 'bin')

MIME_TABLE = {'.txt': 'text/plain',
      '.html': 'text/html',
      '.css': 'text/css',
      '.js': 'application/javascript',
      '.map': 'application/json',
      '.png': 'image/png',
      '.gif': 'image/gif'
     }


def validate_input(inputseqs):

    curheaders = []
    num_aa = 0

    for line in inputseqs.split('\n'):
        if len(line) > 0:
            if line.startswith('>'):
                if line not in curheaders:
                    curheaders.append(line)

                else:
                    return 1 #duplicate header
            else:
                if len(line.strip()) < 15:
                    return 4 #sequence too short

                num_aa += len(line.strip())

                if (num_aa > res_limit):
                    return 3 #too many amino acids

                for s in line.strip().upper():
                    if s in invalid_aa: #invalid amino acid
                        return 2

    if len(curheaders) < 1:
        return 5 #improper FASTA format TODO: make this a better test

    return 0

#serve static files
def content_type(path):
    name, ext = os.path.splitext(path)

    if ext in MIME_TABLE:
        return MIME_TABLE[ext]

    else:
        return "application/octet-stream"

def static_middleware(environ, start_response):

    path = SCRIPTDIR + environ['PATH_INFO']

    if (environ['PATH_INFO'].startswith(STATIC_URL_PREFIX) or
        environ['PATH_INFO'].startswith(STATIC_FILE_DIR)) and os.path.exists(path):
        #return the file
        content = open(path, 'rb').read()

        headers = [('content-type', content_type(SCRIPTDIR + path))]
        start_response('200 OK', headers)
        return [content]

    else:
        pass #should be a 404


def moreronn_middleware(environ, start_response):

    fastaSequences = ''
    returnvalues = {}

    responseheaders = [('Content-Type', 'application/json')]

    #parse the POSTed data
    request_body_size = 0
    try:
        request_body_size = int(environ.get('CONTENT_LENGTH', 0))
    except:
        pass

    request_body = environ['wsgi.input'].read(request_body_size)

    parameters = parse_qs(request_body)

    if 'MoreRONNVersion' in parameters:
        #wait for a bit, for added UI effect
        #this sleep() does nothing else but make the page look good...
        time.sleep(1)

        mrn_proc = subprocess.Popen([os.path.join(MORERONN_BIN_DIR, 'moreRONN')], stderr = subprocess.PIPE)

        mrn_out = mrn_proc.communicate()[1]

        for line in mrn_out.split('\n'):
            if 'version' in line:
                returnvalues['version'] = line[line.find(':')+1:].strip()
                start_response('200 OK', responseheaders)
                return json.dumps(returnvalues)


    if 'fastaSequences' in parameters:
        fastaSequences = parameters['fastaSequences'][0]
    else:
        fastaSequences = ''

    if len(fastaSequences) > 0:
        #validate headers
        isInputValid = validate_input(inputseqs = fastaSequences)

        if isInputValid != 0:
            if isInputValid == 1:
                returnvalues['message'] = "FASTA headers must be unique!"
            elif isInputValid == 2:
                returnvalues['message'] = "Please check your sequences for invalid amino acids (B, J, O, U, X, Z)"
            elif isInputValid == 3:
                returnvalues['message'] = "This server is limited to %s residues. Please try fewer sequences pr download the standalone application." % (res_limit,)
            elif isInputValid == 4:
                returnvalues['message'] = "Please make sure input sequences are at least 15 residues long."
            elif isInputValid == 5:
                returnvalues['message'] = "Please provide sequences in FASTA format."


            #send the error response back right away
            start_response('200 OK', responseheaders)

            return json.dumps(returnvalues)

        #if we made it this far, we are probably okay to run MoreRONN

        #call moreRONN, redirect stderr to stdout to capture everything
        mrn_proc = subprocess.Popen([os.path.join(MORERONN_BIN_DIR, 'moreRONN'), '-s' , '-d', os.path.join(MORERONN_BIN_DIR, 'moreronn_master_data.dat')], stdin = subprocess.PIPE, stdout = subprocess.PIPE)

        mrn_out = mrn_proc.communicate(input = fastaSequences)[0]

        #parse the output nicely
        #MoreRONN gives us a nice output for its headers, but we will reconstruct
        #it ourselves using the per-residue scores, since we can then return a more
        #useful data structure for graphing, user interaction et c.
        predictions = OrderedDict()
        rawscores = OrderedDict()
        rawsequences = OrderedDict()

        raw_full_output = {}

        curseq = ''
        for line in mrn_out.split('\n'):
            if len(line) > 0:
                if line.startswith('(Sequence '): #get the sequence name out

                    curseq = line[line.find(':')+1:].strip()
                    predictions[curseq] = ''
                    rawsequences[curseq] = []
                    rawscores[curseq] = []
                    raw_full_output[curseq] = ''

                elif line.startswith('>') == False:

                    tkns = line.strip().split('\t')
                    if len(tkns) == 2:
                        rawsequences[curseq].append(tkns[0])

                        score = float(tkns[1])
                        rawscores[curseq].append(score)

                        if score < 0.4:
                            predictions[curseq] += '!' #convert this to a blank later
                        elif score >= 0.4 and score < 0.5:
                            predictions[curseq] += '-'
                        elif score >= 0.5 and score < 0.6:
                            predictions[curseq] += '='
                        elif score > 0.6:
                            predictions[curseq] += '#'
                        else:
                            predictions[curseq] += '?' #invalid score

                #and finally, we want to capture the full output anyway
                raw_full_output[curseq] += line.strip() + "\n"


        returnvalues['fastaSequences'] = fastaSequences
        returnvalues['message'] = "OK"
        returnvalues['predictions'] = predictions
        returnvalues['rawscores'] = rawscores
        returnvalues['rawsequences'] = rawsequences
        returnvalues['raw_full_output'] = raw_full_output

        start_response('200 OK', responseheaders)

        return json.dumps(returnvalues)

    else:
        start_response('200 OK', [('Content-Type', 'text/html')])
        return "<html></html>"

#global application
def application(environ, start_response):
    if environ['PATH_INFO'].startswith(STATIC_URL_PREFIX):
        return static_middleware(environ, start_response)
    elif environ['PATH_INFO'].startswith(PREDICT_PREFIX):
        return moreronn_middleware(environ, start_response)

    start_response('200 OK', [('Content-Type', 'text/html')])
    cont = open(os.path.join(SCRIPTDIR, os.path.join('html','index.html')), 'r').read()
    return cont

if __name__ == "__main__":
    from wsgiref.simple_server import make_server
    srv = make_server('localhost', 9090, application)
    srv.serve_forever()
