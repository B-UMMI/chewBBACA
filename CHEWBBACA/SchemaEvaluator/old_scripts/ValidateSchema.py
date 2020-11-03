#!/usr/bin/env python3

import os
import sys
import json
import time
import argparse
from operator import itemgetter

from Bio import SeqIO

try:
	from SchemaEvaluator import CheckCDS,alleleSizeStats
except ImportError:
	from CHEWBBACA.SchemaEvaluator import CheckCDS,alleleSizeStats

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True
    print(program+" not found")
    sys.exit()
    return "Not found"

def main(genes, cpuToUse, htmlFile, transTable, threshold, splited, title, logScale, OneBadGeneNotConserved, skipHeavy):

    outputpath = os.path.dirname(htmlFile)

    starttime = "\nStarting Script at : " + time.strftime("%H:%M:%S-%d/%m/%Y")
    print (starttime)

    print ("Checking all programs are installed")
    print ("Checking mafft installed... " + str(which('mafft')))
    print ("Checking clustalw2 installed... " + str(which('clustalw2')))

    try:
        f = open(genes, 'r')
        f.close()
    except IOError:
        listbasename = os.path.basename(os.path.normpath(genes))

        with open("listGenes" + listbasename + ".txt", "w") as f:
            for gene in os.listdir(genes):
                try:
                    genepath = os.path.join(genes, gene)
                    #gene_fp2 = HTSeq.FastaReader(genepath)
                    for allele in SeqIO.parse(genepath, "fasta"):
                        break
                    f.write(genepath + "\n")
                except Exception as e:
                    print(e)
                    pass

        genes = "listGenes" + listbasename + ".txt"

    genebasename = str(os.path.basename(genes))
    genebasename = genebasename.split(".")
    genebasename.pop()
    genebasename = ".".join(genebasename)

    notConservedgenes, totalgenes, genesWOneAllele, boxplot, histplot, allelenumberplot, listgenesBoxOrdered, totalnumberofgenes, boxListLink, allAllelesStats = alleleSizeStats.main(
        genes, threshold, OneBadGeneNotConserved, True, logScale, outputpath, splited)

    histplot = str(json.dumps(histplot))
    allelenumberplot = str(json.dumps(allelenumberplot))
    allAllelesStats = str(json.dumps(allAllelesStats))

    statsPerGene = CheckCDS.main(genes, transTable, True, outputpath, cpuToUse,skipHeavy)

    # stats values are ordered in a list allelesNotMultiple3,listStopcodonsInside,listnotStartCodon,numberOfAlleles
    htmlgenespath = os.path.join(outputpath, "genes_html/")
    relpath = os.path.relpath(htmlgenespath, outputpath)

    if not os.path.exists(htmlgenespath):
        os.makedirs(htmlgenespath)

    with open(htmlFile, "w") as f:
        f.write(
            "<!DOCTYPE html>\n<html>\n<head><script type='text/javascript' src='https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.5/d3.min.js'></script>\n")
        f.write(
            "<script type='text/javascript' src='https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js'></script>")
        f.write("""<!-- Latest compiled and minified JavaScript -->
		<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js" integrity="sha384-0mSbJDEHialfmuBBQP6A4Qrprq5OVfW37PRR3j5ELqxss1yVqOtnepnHVP9aJ7xS" crossorigin="anonymous"></script>""")
        f.write("""<style type="text/css">
      body {
        padding-bottom: 60px;
        padding-left: 20px;
        padding-right: 20px;
      }

      /* Custom container */
      .container {
        margin: 0 auto;
      }
      .container > hr {
        margin: 60px 0;
      }

      /* Main marketing message and sign up button */
      .jumbotron {
        text-align: center;
      }
      .jumbotron h1 {
        font-size: 100px;
        line-height: 1;
      }
      .jumbotron .lead {
        font-size: 24px;
        line-height: 1.25;
      }
      .jumbotron .btn {
        font-size: 21px;
        padding: 14px 24px;
      }

    </style>
		<!-- Latest compiled and minified CSS -->
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous">
<script>
function openList1(documentid) {
    var list = document.getElementById(documentid);

    if (list.style.display == "none"){
        list.style.display = "block";
    }else{
        list.style.display = "none";
    }
}
</script>""")

        f.write("""<style type="text/css">
		ul {
    /*min-height: 300px;*/
    -webkit-column-count: 4;
       -moz-column-count: 4;
            column-count: 4; /*4 is just placeholder -- can be anything*/
}
li {
    display: table;
    padding-bottom: 20px; 
    margin-right: 30px;
}
li a {
    color: rgb(0, 162, 232);
}

.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-vn4c{background-color:#D2E4FC;text-align:center}
.tg .tg-qpvr{background-color:#009901;vertical-align:top}
.tg .tg-14d4{background-color:#91df0a;text-align:center;vertical-align:top}

</style>\n</head>\n<body>""")

        f.write("""<div class="jumbotron">
  <h2>""" + title + """</h2>
  <p>Explore the analysis by clicking the analysis buttons</p></div>""")

        f.write("<h2>Allele size analysis using a mode +/- " + str(threshold) + "</h2>")
        f.write(
            "<p> Genes are considered not conserved if >1 allele are outside the mode +/-0.05 size. Genes with only 1 allele outside the threshold are considered conserved</p>\n")
        f.write("<h3>" + str(totalgenes) + " total loci</h3>")
        f.write("\n<h3>" + str(len(notConservedgenes)) + " loci with high length variability</h3>")
        f.write(
            """<button onclick = "openList1('ollist1')">Show/hide list</button><ol id='ollist1' style='display: none;'>""")

        for elem in notConservedgenes:
            f.write("<li><a href = '" + str(elem) + "' target='_blank'>" + (os.path.basename(elem)).replace(".html",
                                                                                                            "") + "</a></li>")

        f.write("</ol>\n<h3>" + str(len(genesWOneAllele)) + " loci with only one allele</h3>\n")

        f.write(
            """<button onclick = "openList1('ollist2')">Show/hide list</button><ol id='ollist2' style='display: none;'>""")

        for elem in genesWOneAllele:
            f.write("<li><a href = '" + str(elem) + "' target='_blank'>" + (os.path.basename(elem)).replace(".html",
                                                                                                            "") + "</a></li>\n")

        f.write("</ol>\n")

        f.write("""<div class='container'>
					<div class='row'>
						<div class="col-sm-3">
						<button id='button2' class="btn btn-success btn-block active">Allele length analysis</button>
						</div>
						<div class="col-sm-3">
						<button id='button3' class="btn btn-success btn-block active">Allele number analysis</button>
						</div>
						<div class="col-sm-3">
						<button id='button1' class="btn btn-success btn-block active">Loci size variation analysis</button>
						</div>
						<div class="col-sm-3">
						<button id='button4' class="btn btn-success btn-block active">CDS analysis</button>
						</div>
					</div>
				</div>""")

        f.write("""<div id="fig03" ><div id="fig03b" style="width: 1000px; height: 600px;"></div><div id="fig03a" style="width: 1300px; height: 700px;"></div></div>
		<div id="fig01" style="display:none; width: 100%;">
		<h2>Size boxplot for all loci</h2><p>Box plot for the sizes of each locus. Loci ordered by increasing median allele size</p>
		<p>Use the zoom button and hover the mouse over a box/median to see the locus name and points data</p>
		<p>-->Use the following buttons to navigate through all loci</p>
		<button id='buttonbackward' > Previous 500 loci < </button>
		<button id='buttonforward' > Next 500 loci > </button>
		""")

        f.write("""<div id="figbox0" classe="container" style="display:none;width: 1500px; height: 1000px;"></div>""")

        f.write(
            """</div><div id="fig02" classe="container" style="display:none;width: 1000px; height: 600px;"></div>""")

        f.write("""<script type="text/javascript">
					
					$("#buttonforward").click(function(){
					  boxplotPage=boxplotPage+1;
					  $.ajax({
									async: false,
									url: "jsonbox"+boxplotPage+".js",
									dataType: "script",
									
									error: function(){
										  boxplotPage=boxplotPage-1;
											},
								  success: function(){
						  										
										numberofgenes=numberofgenes+jsonboxList.length
										var listboxes=[];
										var i = 0;
										var len = jsonboxList.length;
										for (; i < len; i++) { 
											var trace = {
														  y: jsonboxList[i],
														  name: listgenes[i],
														  text: listlinks[i],
														  type: 'box',
														  marker: {size: 2},
														};
											listboxes.push(trace)
										}
										var layout = {
													  title: numberofgenes+'/'+totalGenes,
													  showlegend: false,
													  xaxis: {
																autotick: false,
																showgrid: true,
																showticklabels: false
															  },
													};
										Plotly.newPlot('figbox0', listboxes,layout);
						
						var myPlot = document.getElementById('figbox0');
							myPlot.on('plotly_click', function(data){
								link = '';
								for(var i=0; i < data.points.length; i++){
									link = data.points[i].data.text;
								}
								window.open(link, '_blank');
							});
						
						}
					    
					});
					return false;
					}); 
					</script>""")

        f.write("""<script type="text/javascript">
					$("#buttonbackward").click(function(){
					  boxplotPage=boxplotPage-1;
					  numberofgenes=numberofgenes-jsonboxList.length;
					  $.ajax({
									async: false,
									url: "jsonbox"+boxplotPage+".js",
									dataType: "script",
									
									error: function(){
										numberofgenes=numberofgenes+jsonboxList.length;
										  boxplotPage=boxplotPage+1;
											},
								  success: function(){
						  										
										var listboxes=[];
										var i = 0;
										var len = jsonboxList.length;
										for (; i < len; i++) { 
											var trace = {
														  y: jsonboxList[i],
														  name: listgenes[i],
														  text: listlinks[i],
														  type: 'box',
														  marker: {size: 2},
														};
											listboxes.push(trace)
										}
										var layout = {
													  title: numberofgenes+'/'+totalGenes,
													  showlegend: false,
													  xaxis: {
																autotick: false,
																showgrid: true,
																showticklabels: false
															  },
													};
										Plotly.newPlot('figbox0', listboxes,layout);
						
						var myPlot = document.getElementById('figbox0');
							myPlot.on('plotly_click', function(data){
								link = '';
								for(var i=0; i < data.points.length; i++){
									link = data.points[i].data.text;
								}
								window.open(link, '_blank');
							});
						
						}
					    
					});
					return false;
					}); 
					</script>""")

        f.write("""<script type="text/javascript">
					$("#button3").click(function(){
					  if ($("#fig03a").firstChild==undefined){
					  $.ajax({
								async: false,
								url: "json3.js",
								dataType: "script"
							});
						var listTraces=[];
						var i = 0;
						var len = jsonsScatterPlot.length;
						var listName=['Mode','Mean','Median'];
						for (; i < len; i++) { 
							var aux=jsonsScatterPlot[i];
							var listgenes=[];
							var j = 0;
							for (; j < aux[2].length; j++) { 
								genename=aux[2][j];
								listgenes.push(genename.replace(/^.*[\\\/]/, ''));
							}
							var trace = {
									x: aux[0],
									  y: aux[1],
									  name: listName[i],
									  customdata: aux[2],
									  text: listgenes,
									  mode: 'markers',
									  type: 'scattergl'
									};
							listTraces.push(trace)
								}
						var layout = {
									  yaxis: {title: "Number of alleles"},
									  xaxis: {
												title: "Allele size in bp",
											  },
									};
						Plotly.newPlot('fig03a', listTraces,layout);
							
						var trace = [{
									  x: listNumberDifAlleles,
									  type: 'histogram'
									}];
						var layout = {
									  barmode: "stack",
									  showlegend: false,
									  yaxis: {title: "Number of Loci"},
									  xaxis: {
												autotick: true,
												showgrid: true,
												title: "Number of Different Alleles",
												showticklabels: true
											  },
									};
						Plotly.newPlot('fig03b', trace,layout);
						
						var myPlot = document.getElementById('fig03a');
							myPlot.on('plotly_click', function(data){
								link = '';
								link = data.points[0].data.customdata[data.points[0].pointNumber];
								window.open(link, '_blank');
							});
						
						
					  }
					  $("#fig03").css({"display":"block"});
					  $("#fig01").css({"display":"none"});
					  $("#fig02").css({"display":"none"});
					  $("#fig04").css({"display":"none"});
					}); 
					</script>""")

        f.write("""<script type="text/javascript">
					$("#button2").click(function(){
					 if ($("#fig02").firstChild==undefined){
						$.ajax({
								async: false,
								url: "json2.js",
								dataType: "script"
							});
						
	
						var trace = [{
									  x: jsonHistPlot,
									  type: 'histogram'
									}];
						var layout = {
									  title: 'Distribution of allele mode sizes per locus',
									  barmode: "stack",
									  showlegend: false,
									  yaxis: {title: "Number of locus"},
									  xaxis: {
												autotick: true,
												showgrid: true,
												title: "Allele Size",
												showticklabels: true
											  },
									};
						Plotly.newPlot('fig02', trace,layout);
					}
					  $("#fig02").css({"display":"block"});
					  $("#fig03").css({"display":"none"});
					  $("#fig01").css({"display":"none"});
					  $("#fig04").css({"display":"none"});
					}); 
					</script>""")

        f.write("""<script type="text/javascript">
					var boxplotPage=0;
					var numberofgenes=0
					var totalGenes=""" + str(totalnumberofgenes) + """
					$("#button1").click(function(){
					 if ($("#figbox0").firstChild==undefined){
						$.ajax({
								async: false,
								url: "jsonbox0.js",
								dataType: "script"
							});
						var listboxes=[];
						var i = 0;
						numberofgenes=jsonboxList.length;
						var len = jsonboxList.length;
						for (; i < len; i++) { 
							var trace = {
										  y: jsonboxList[i],
										  name: listgenes[i],
										  text: listlinks[i],
										  type: 'box',
										  marker: {size: 2},
										};
							listboxes.push(trace);
						}
						var layout = {
									  title: numberofgenes+'/'+totalGenes,
									  showlegend: false,
									  xaxis: {
												autotick: false,
												showgrid: true,
												showticklabels: false
											  },
									};
						Plotly.newPlot('figbox0', listboxes,layout);
						
						var myPlot = document.getElementById('figbox0');
							myPlot.on('plotly_click', function(data){
								link = '';
								for(var i=0; i < data.points.length; i++){
									link = data.points[i].data.text;
								}
								window.open(link, '_blank');
							});
						}
						
					  $("#fig01").css({"display":"block"});
					  $("#figbox0").css({"display":"block"});
					  $("#fig02").css({"display":"none"});
					  $("#fig03").css({"display":"none"});
					  $("#fig04").css({"display":"none"});
					}); 
					</script>""")
        f.write("""<script type="text/javascript">
					$("#button4").click(function(){
					  $("#fig04").css({"display":"block"});
					  $("#fig02").css({"display":"none"});
					  $("#fig03").css({"display":"none"});
					  $("#fig01").css({"display":"none"});
					}); 
					</script>""")

        f.write("""<div id="fig04" classe="container" style="display:none">
		<title>Schema Validation Results</title>
		<h3>Summary of problematic alleles per locus using the <a href='http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11'>NCBI translation table 11</a></h3>
		<div id="fig04B" classe="container" style="width: 1200px; height: 700px;"></div>
		<p>Click the locus name to open the locus info page</p>
		<p>click the boxes with % to get the index of the alleles with problems</p>
		<table class="tg">
		  <tr>
			<th class="tg-031e">Gene</th>
			<th class="tg-031e">Number alleles not multiple of 3</th>
			<th class="tg-031e">Number alleles w/ >1 stop codons</th>
			<th class="tg-031e">Number alleles wo/ Start/Stop Codon</th>
			<th class="tg-qpvr" .background-color='#59b300'>Number of alleles (% alleles w/ issues) </th>
		  </tr>""")
        ordered = []
        orderedBySize = []
        for key, value in statsPerGene.items():
            numberMultip = float(len(value[0]))
            numberStop = float(len(value[1]))
            numberStart = float(len(value[2]))
            total = float(value[3])
            totalpercent = ((numberMultip + numberStart + numberStop) / total) * 100
            ordered.append([key, totalpercent])
            orderedBySize.append([key, total])
        ordered = sorted(ordered, key=itemgetter(-1))
        ordered.reverse()
        orderedBySize = sorted(orderedBySize, key=itemgetter(-1))

        newlist = []
        auxHistList = [[], [], [], [], []]
        i = 0
        while i < len(ordered):
            item = ordered[i]

            aux = []
            aux.append(item[0])

            value = statsPerGene[item[0]]
            numberMultip = float(len(value[0]))
            numberStop = float(len(value[1]))
            numberStart = float(len(value[2]))
            total = float(value[3])

            aux.append(value)
            newlist.append(aux)
            name = os.path.basename(str(item[0]))
            name = name.split(".")
            name = name[0]
            if (numberMultip > 0 or numberStop > 0 or numberStart > 0):
                locusHTML = os.path.join(relpath, (os.path.basename(str(item[0]))).replace(".fasta", ".html"))
                f.write("<tr id=" + str(item[
                                            0]) + """>\n<td class='tg-vn4c' onclick="window.open('""" + locusHTML + """')" style='cursor:pointer'>""" + name + "</td>\n<td class='tg-vn4c' onclick='a(this);'>" + str(
                    int(numberMultip)) + " (" + str('{0:.2f}'.format(
                    (numberMultip / total) * 100)) + "%)" + "</td>\n<td class='tg-vn4c' onclick='a(this);'>" + str(
                    int(numberStop)) + " (" + str('{0:.2f}'.format(
                    (numberStop / total) * 100)) + "%)" + "</td>\n<td class='tg-vn4c' onclick='a(this);'>" + str(
                    int(numberStart)) + " (" + str('{0:.2f}'.format(
                    (numberStart / total) * 100)) + "%)" + "</td>\n<td class='tg-14d4' onclick='a(this);'>" + str(
                    int(total)) + " (" + str('{0:.2f}'.format(
                    ((numberMultip + numberStart + numberStop) / total) * 100)) + "%)" + "</td>\n</tr>")
            i += 1
        i = 0
        while i < len(orderedBySize):
            item2 = orderedBySize[i]
            value2 = statsPerGene[item2[0]]

            if (len(value2[0]) > 0 or len(value2[1]) > 0 or len(value2[2]) > 0):
                locusHTML = os.path.join(relpath, (os.path.basename(str(item2[0]))).replace(".fasta", ".html"))
                auxHistList[0].append(locusHTML)
                auxHistList[1].append(float(len(value2[0])))
                auxHistList[2].append(float(len(value2[1])))
                auxHistList[3].append(float(len(value2[2])))
                auxHistList[4].append(
                    float(value2[3]) - float(len(value2[2])) - float(len(value2[1])) - float(len(value2[0])))
            i += 1

        auxHistList = str(json.dumps(auxHistList))
        f.write("</table>")
        f.write(
            """<div id='AllelesWissues'></div><button onclick="$('#AllelesWissues').empty();">clean</button></div>""")
        f.write("""\n<script type='text/javascript'>
		var histCDSanalysis=""" + auxHistList + """;
		var listTraces=[];
		var trace = {
					  x: histCDSanalysis[1],
					  name: 'Non multiple 3',
					  link: histCDSanalysis[0],
					  orientation: 'h',
					  type: 'bar'
					};
		var trace1 = {
					  x: histCDSanalysis[2],
					  name: '> 1 stop codon',
					  link: histCDSanalysis[0],
					  orientation: 'h',
					  type: 'bar'
					};
		var trace2 = {
					  x: histCDSanalysis[3],
					  name: 'No Start/Stop codon',
					  link: histCDSanalysis[0],
					  orientation: 'h',
					  type: 'bar'
					};
		var trace3 = {
					  x: histCDSanalysis[4],
					  name: 'CDS Alleles',
					  link: histCDSanalysis[0],
					  orientation: 'h',
					  type: 'bar'
					};			
		
					
		var layout = {
					  barmode: "stack",
					  yaxis: {
								autotick: false,
								showgrid: false,
								showticklabels: false
					  
							},
					  xaxis: {title: "Number of occurrences"},
								
					};
					
		listTraces.push(trace3,trace2,trace,trace1);			
		Plotly.newPlot('fig04B', listTraces,layout);
		
		var myPlot = document.getElementById('fig04B');
					myPlot.on('plotly_click', function(data){
						i= data.points[0].y
						link = data.points[0].data.link[i];
						window.open(link, '_blank');
					});
		
		</script>
		""")
        f.write("""\n<script type='text/javascript'>function a(element) {
	var id = $(element).closest("tr").attr("id");
	var badalleles=[];
	for (i = 0; i < alleles.length; i++) { 
		if((alleles[i])[0]==id){
			badalleles=(alleles[i])[1];
			break;
			}
		}
	var notmulti=(badalleles[0]).join('; ');
	var stopcodon=(badalleles[1]).join('; ');
	var startcodon=(badalleles[2]).join('; ');
	var name=(id.split("/")).slice(-1)[0]
	name=(name.split("."))[0]
	$('#AllelesWissues').append('<h2> Locus: '+name+'</h2>');
	$('#AllelesWissues').append('<p> Alleles not multiple of 3: '+notmulti+'</p><p> Alleles with >1 stop codon: '+stopcodon+'</p><p>Alleles without start or stop codon: '+startcodon+'</p>');

	$('html,body').animate({
        scrollTop: $('#AllelesWissues').offset().top},'slow');
	
	}
	</script>""")

        f.write("\n<script type='text/javascript'>var alleles=" + json.dumps(newlist) + "</script>")
        f.write("</body>\n</html>")

    i = 0

    filename = "jsonbox0.js"
    for elem in boxplot:
        boxplotElem = str(json.dumps(elem))
        listGenesJson = str(json.dumps(listgenesBoxOrdered[i]))
        listLinkGenesJson = str(json.dumps(boxListLink[i]))
        filename = "jsonbox" + str(i) + (".js")
        with open((os.path.join(outputpath, filename)), "w") as f:
            f.write("var jsonboxList =" + str(boxplotElem) + ";var listgenes=" + str(
                listGenesJson) + ";var listlinks=" + str(listLinkGenesJson))
        i += 1

    with open((os.path.join(outputpath, "json2.js")), "w") as f:

        f.write("var jsonHistPlot =" + str(histplot))

    with open((os.path.join(outputpath, "json3.js")), "w") as f:

        f.write("var jsonsScatterPlot =" + str(allelenumberplot) + ";var listNumberDifAlleles=" + allAllelesStats)

    try:
        os.remove("listGenes" + listbasename + ".txt")
    except:
        pass

    print (starttime)
    print ("Finished Script at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))


if __name__ == "__main__":
    main()
