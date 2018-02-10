var app = angular.module('moreronnApp', ['chart.js', 'ui-notification', 'ui.bootstrap']);

//filter to collapse sequence and prediction together
//into the familiar MoreRONN header, with added metadata
//about the raw scores for user interaction
app.filter('format_headers', function($sce) {

  return function(seqinput, scoreinput, rawscores) {

      var outseq = '';
      var outscore = '';
      var finalval = '';
      var npl = 125; //number of residues per line; make this fluid later?

      var curcount = 0;
      for (var i = 0; i < seqinput.length; i++)
      {
          outseq += '<a title="' + seqinput[i] + ': ' + rawscores[i] + '" data-toggle="tooltip">' + seqinput[i] + '</a>';

          if (scoreinput[i] == '!') { //replace exclamation mark with a space
            //outscore += '<a title="' + seqinput[i] + ': ' + rawscores[i] + '" data-toggle="tooltip">' + '&nbsp;' + '</a>';
            outscore += '&nbsp;'
          }
          else if (scoreinput[i] == '=' || scoreinput[i] == '#') { //make disorder red
            outscore += '<font color="red">' + scoreinput[i] + '</font>';
            //outscore += scoreinput[i];
          }

          else {
            outscore += scoreinput[i];
          }

          curcount++;

          //format the string once we get past the npl
          if (curcount == npl)
          {
            finalval += '<br />' + outseq + '<br />' + outscore + '<br />';

            outseq = '';
            outscore = '';
            curcount = 0;

          }

      }

      //whatever's left
      finalval += '<br />' + outseq + '<br />' + outscore + '<br />';

      return $sce.trustAsHtml(finalval);

  };

});

app.controller('moreronnController', function($scope, $http, $sce, $location, $anchorScroll, $uibModal, Notification) {

    //current moreronn version on server
    $scope.moreronnVersion = 0.0; //seed value; will be populated by binary
    $scope.serverOnline = true;

    //submitted form data
    $scope.formData = {};

    $scope.message = ""; //to get notification messages from the server

    $scope.predictions = '';
    $scope.rawscores = '';
    $scope.rawsequences = '';

    $scope.isSubmittingFasta = false;
    $scope.isGettingVersion = false;
    $scope.isRenderingChart = false;

    //charting stuff
    $scope.chartlabels = [];
    $scope.chartseries = [];
    $scope.chartdata = [];
    $scope.rawResults = '';
    $scope.chartcolours = ['#FD1F5E', '#0000CF']; //for the O/D bdy. and pred.


    //submit fasta sequences
    var getVersion = function() {

        $scope.isGettingVersion = true;
      	$http({

      	    //call the WSGI driver here
      	    method: 'POST',
      	    url: ajaxURL,
      	    data: $.param({'MoreRONNVersion': '0'}),
      	    headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
            timeout: 10000 //'cos we really don't want to waste time with this

      	})
      	.success(function(data) {

            $scope.moreronnVersion = data.version;
            $scope.isGettingVersion = false;
            $scope.serverOnline = true;

      	})
        .error(function(data) {

            $scope.moreronnVersion = 0.0;
            $scope.isGettingVersion = false;
            $scope.serverOnline = false;

        });

    };

    //get moreRONN version
    getVersion();

    //download raw data as a blob
    $scope.downloadRawData = function(thekey) {
        var outBlob = new Blob([$scope.rawResults[thekey]], { type: 'text/plain; charset=utf-8' });
        var fname = thekey + ".txt";
        saveAs(outBlob, fname);

    };

    //submit fasta sequences
    $scope.submitSequences = function() {

        $scope.isSubmittingFasta = true;

        //make strings uppercase
        var myinput = '';

        var breaks = $scope.formData['fastaSequences'].split('\n');

        for (var j = 0; j < breaks.length; j++)
        {
          if (breaks[j][0] == '>')
          {
            myinput = myinput + breaks[j] + '\n';
          }
          else {
            myinput = myinput + breaks[j].toUpperCase() + '\n';
          }
        }

        $scope.formData['fastaSequences'] = myinput;


      	$http({

      	    //call the WSGI driver here
      	    method: 'POST',
      	    url: ajaxURL,
      	    data: $.param($scope.formData),
      	    headers: { 'Content-Type': 'application/x-www-form-urlencoded' }

      	})

      	.success(function(data) {

      	    $scope.message = data.message;

            if ($scope.message == "OK")
            {
              Notification.success({message: "<strong>Sequences processed!</strong><br /><em>Please scroll down for results.</em>", delay: 3000});
            }
            else
            {
              Notification.error({message: "<strong>Error!</strong><br /><em>" + $scope.message + "</em>", delay: 3000});
            }

      	    //retrieve the prediction headers
      	    //var mypreds = data.predictions.join(',').replace(/,/g, '<br />').replace(/!/g, '&nbsp;').replace(/SPACER/g, '<hr />').replace(/\>\$/g, '<br />');
      	    //$scope.predictionHeaders = $sce.trustAsHtml(mypreds);
            $scope.predictions = data.predictions;
            $scope.rawscores = data.rawscores;
            $scope.rawsequences = data.rawsequences;
            $scope.rawResults = data.raw_full_output;

            $scope.isSubmittingFasta = false;

      	})

        .error(function(data) {

            $scope.serverOnline = false;
            $scope.isSubmittingFasta = false;

            Notification.error({message: "<strong>The server is currently offline.</strong><br /> \
            <em>We apologise for the inconvenience. Please try again later.</em>", delay: 5000});

        });
    };

    //render the Chart.js line chart
    $scope.renderChart = function (predKey) {

      $scope.isRenderingChart = true;

      //clear the chart data
      $scope.chartdata = [];
      $scope.chartlabels = [];
      $scope.chartseries = [];

      //default charting options
      $scope.chartoptions = {

        animation: false,

        responsive: true,
        scaleShowVerticalLines: false,
        bezierCurveTension: 0.1,

        datasetFill: false,

        pointDot : false,

        scaleOverride: true,
        scaleSteps: 10,
        scaleStepWidth: 1/10,
        scaleStartValue: 0.0,

        scaleIntegersOnly: false,

      };

      $scope.chartseries = ['O/D Score Boundary', predKey];

      //push the data
      var boundaryArr = [];
      var xlabels = []; //populate a blank array of labels
      for (var j = 0; j < $scope.rawsequences[predKey].length; j++)
      {
        boundaryArr.push(0.5);
        xlabels.push("");
      }

      $scope.chartdata.push(boundaryArr);
      $scope.chartdata.push($scope.rawscores[predKey]);

      $scope.chartlabels = xlabels;

      //scroll to the chart
      $location.hash('disorderChart');
      $anchorScroll();

      $scope.isRenderingChart = false;

    };


    //
    // $scope.chartClick = function (points, evt) {
    //   console.log(points, evt);
    // };

});
