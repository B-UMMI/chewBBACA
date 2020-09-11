#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
import os
import argparse
import json
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import ClustalwCommandline
import multiprocessing
import copy

def reverseComplement(strDNA):

    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:

        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]

def translateSeq(DNASeq,transTable):
    seq=DNASeq
    tableid=transTable
    reversedSeq=False
    try:
        myseq= Seq(seq)
        protseq=Seq.translate(myseq, table=tableid,cds=True)
    except:
        reversedSeq=True
        try:
            seq=reverseComplement(seq)
            myseq= Seq(seq)
            protseq=Seq.translate(myseq, table=tableid,cds=True)

        except:
            try:
                seq=seq[::-1]
                myseq= Seq(seq)
                protseq=Seq.translate(myseq, table=tableid,cds=True)
            except:
                reversedSeq=False
                try:
                    seq=seq[::-1]
                    seq=reverseComplement(seq)
                    myseq= Seq(seq)
                    protseq=Seq.translate(myseq, table=tableid,cds=True)
                except Exception as e:

                    raise ValueError(e)
    return protseq,seq,reversedSeq


def call_mafft(path_to_save,genefile):

    try:
        print("maffting "+os.path.basename(genefile))
        mafft_cline = MafftCommandline(input=genefile)
        stdout, stderr = mafft_cline()
        with open(path_to_save, "w") as handle:
            handle.write(stdout)
        return True

    except Exception as e:
        print(e)
        return False
def call_clustalw(genefile,htmlpath):

    try:
        print("clustaling "+os.path.basename(genefile))
        inputfile=os.path.join(htmlpath,(os.path.basename(genefile)).replace(".fasta","_aligned.fasta"))
        clw_cline = ClustalwCommandline("clustalw2", infile=inputfile, tree=True)
        clw_cline()

        return True

    except Exception as e:
        print(e)
        return False
#~ def main():
#~ 
    #~ parser = argparse.ArgumentParser(description="This program analyses cds")
    #~ parser.add_argument('-i', nargs='?', type=str, help='list genes', required=True)
    #~ parser.add_argument('-r', nargs='?', type=bool, help='Return values', required=False)
#~ 
    #~ args=parser.parse_args()
#~ 
#~ 
    #~ genes = args.i
#~ 
    #~ try:
        #~ ReturnValues=bool(args.r)
    #~ except:
        #~ ReturnValues=False
        #~ pass
#~ 
    #~ analyzeCDS(genes,ReturnValues)

def main(genes,transTable,ReturnValues,outputpath,cpu,skipClustalMafft):

    gene_fp = open( genes, 'r')

    stopc=0
    notStart=0
    notMultiple=0
    totalalleles=0

    htmlgenespath=os.path.join(outputpath,"genes_html/")

    if not os.path.exists(htmlgenespath):
        os.makedirs(htmlgenespath)

    listgenes=[]
    for gene in gene_fp:

        gene = gene.rstrip('\n')
        gene = gene.rstrip('\r')
        listgenes.append(gene)
    statsPerGene={}


    if not skipClustalMafft:
        pool = multiprocessing.Pool(cpu)

        for gene in listgenes:

            gene = gene.rstrip('\n')
            gene = gene.rstrip('\r')

            alignFileName=os.path.join(htmlgenespath,(os.path.basename(gene)).replace(".fasta","_aligned.fasta"))

            pool.apply_async(call_mafft,args=[alignFileName,gene])

        pool.close()
        pool.join()

        pool = multiprocessing.Pool(cpu)

        for gene in listgenes:

            gene = gene.rstrip('\n')
            gene = gene.rstrip('\r')

            alignFileName=os.path.join(htmlgenespath,(os.path.basename(gene)).replace(".fasta","_aligned.fasta"))

            pool.apply_async(call_clustalw,args=[gene,htmlgenespath])

        pool.close()
        pool.join()

    toPrintCDSStats="Locus\tFrameshift\tNo Start or Stop\tMore than 1 Stop\t Other\n"
    for gene in listgenes:

        gene = gene.rstrip('\n')
        gene = gene.rstrip('\r')

        alignFileName=os.path.join(htmlgenespath,(os.path.basename(gene)).replace(".fasta","_aligned.fasta"))


        listStopc=[]
        listnotStart=[]
        listnotMultiple=[]
        listOther=[]

        print("processing : "+str(os.path.basename(gene)))

        k=0

        multiple=True
        #gene_fp2 = HTSeq.FastaReader(gene)

        alleleSizes=[]
        alleleNames=[]
        alleleSizesNotMultiple=[]
        alleleSizesNotMultipleNames=[]
        alleleSizesTransError=[]
        alleleSizesTransErrorNames=[]
        # translate each allele and report the error if unable to translate
        for allele in SeqIO.parse(gene, "fasta"):

            k+=1
            realAlleleID=allele.id.split("_")[-1]
            # if allele is not multiple of 3 it's useless to try to translate
            if (len(str(allele.seq)) % 3 != 0):
                multiple=False
                listnotMultiple.append(realAlleleID)
                alleleSizesNotMultipleNames.append(realAlleleID)
                alleleSizesNotMultiple.append(len(str(allele.seq)))
                pass
            else:
                try:
                    protseq,seq,reversedSeq=translateSeq(allele.seq, transTable)
                    alleleNames.append(realAlleleID)
                    alleleSizes.append(len(str(allele.seq)))

                except Exception as err:
                    if "Extra in frame stop codon found" in str(err):
                        stopc+=1
                        listStopc.append(str(k))
                    elif "is not a start codon" in str(err):
                        notStart+=1
                        listnotStart.append(str(k))
                    elif "is not a stop codon" in str(err):
                        notStart+=1
                        listnotStart.append(realAlleleID)
                    else:
                        listOther.append(realAlleleID)
                        print(err)

                    alleleSizesTransErrorNames.append(realAlleleID)
                    alleleSizesTransError.append(len(str(allele.seq)))
                    #print "allele "+str(k)+" is not translating"
                    pass


        relpath=os.path.relpath(gene,outputpath)
        statsPerGene[relpath]=copy.deepcopy([listnotMultiple,listStopc,listnotStart,k])
        totalalleles+=k

        #create html per gene
        genename=(os.path.basename(gene)).split(".")
        genename.pop()
        genename=".".join(genename)
        with open(htmlgenespath+genename+".html", "w") as f:
            f.write("<meta charset='UTF-8'><!DOCTYPE html>\n<html>\n<head><script type='text/javascript' src='https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.5/d3.min.js'></script>\n")
            f.write("<script type='text/javascript' src='https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js'></script>")
            f.write("""<!-- Latest compiled and minified JavaScript -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js" integrity="sha384-0mSbJDEHialfmuBBQP6A4Qrprq5OVfW37PRR3j5ELqxss1yVqOtnepnHVP9aJ7xS" crossorigin="anonymous"></script>""")
            f.write("""<style type="text/css">
      body {
        padding-top: 10px;
        padding-bottom: 60px;
        padding-left: 20px;
      }

      /* Custom container */
      .container {
        margin: 0 auto;
        max-width: 1000px;
      }
      .container > hr {
        margin: 60px 0;
      }

      /* Main marketing message and sign up button */
      .jumbotron {
        margin: 40px 0;
        text-align: center;
        padding-top: 30px;
        padding-bottom: 30px;
      }
      .jumbotron h1 {
        font-size: 100px;
        line-height: 1;
      }
      
      phylocanvas {
          width: 100%;
          height: 30em;
        }
        </style>
      <!-- Latest compiled and minified CSS -->
      <script type="application/javascript" src="https://cdn.rawgit.com/phylocanvas/phylocanvas-quickstart/v2.3.0/phylocanvas-quickstart.js"></script>
            <script src="https://cdn.bio.sh/msa/latest/msa.min.gz.js"></script>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous">""")

            f.write("""<div class="jumbotron">
                      <h2>Analysis of the locus """+(os.path.basename(gene)).replace(".fasta","")+"""</h2>
                      <p>Explore the analysis by clicking the analysis buttons</p></div>""")
            relpath=os.path.relpath(gene,htmlgenespath)
            #print relpath
            f.write("<input type=button onClick=window.open('"+relpath+"') value='Get fasta file'>\n<p></p>""")
            f.write("<input type=button onClick=window.open('"+os.path.join((os.path.relpath(alignFileName,alignFileName)),(os.path.basename(gene)).replace(".fasta","_aligned.fasta"))+"') value='Get msa file'>\n""")

            f.write("""<p>MSA alignement done using mafft and visualized using msaviewer</p>
                <p>Tree built with clustalw2 and visualized using phylocanvas</p>
                <div class='container'>
                    <div class='row'>
                        
                        <div class="col-sm-3">
                        <button id='button3' class="btn btn-info btn-block active">Allele Size plot</button>
                        </div>
                        <div class="col-sm-3">
                        <button id='button1' class="btn btn-info btn-block active">NJ tree</button>
                        </div>
                        <div class="col-sm-3">
                        <button id='button2' class="btn btn-info btn-block active">MSA alignment</button>
                        </div>
                    </div>
                </div>""")


            f.write("""</head>\n<body><p></p><div id='snippetDiv' style='display:none'></div> <p></p><div id='phylocanvas' style='display:none;border:solid'>
            <div id="pc-buttons">
                <input type="text" id="searchbox" onkeyup="search(this);">
              <button id="rectangular" class="btn btn-default btn-sm">Rectangular</button>
              <button id="circular" class="btn btn-info btn-sm">Circular</button>
              <button id="radial" class="btn btn-default btn-sm">Radial</button>
              <button id="diagonal" class="btn btn-default btn-sm">Diagonal</button>
              <button id="hierarchical" class="btn btn-default btn-sm">Hierarchical</button>
            </div>
            
            </div><div id='histdiv'></div>\n""")


            f.write("""
            
            <script type="application/javascript">
                        var tree;
                        function doPhylocanvas() {
                          tree = Phylocanvas.createTree('phylocanvas', {
                          // other options
                          contextMenu : [{
                            text: 'Normal Menu',
                            handler: 'triggerNormal',
                            internal: false,
                            leaf: false
                          }, {
                            text: 'Internal Menu',
                            handler: 'triggerInternal',
                            internal: true,
                            leaf: false
                          }, {
                            text: 'Save as PNG',
                            handler: 'exportCurrentTreeView',
                            internal: false,
                            leaf: false
                          }]
                          });
                          var treestr=$.ajax({
                            url: '"""+os.path.join((os.path.relpath(alignFileName,alignFileName)),(os.path.basename(gene)).replace(".fasta","_aligned.ph"))+"""',
                            async: false
                         }).responseText;
                         tree.setTreeType('circular');
                          
                          tree.load(treestr);
                        }  
                          function search (ele) {
                            if(ele.value !== ""){
                              someleaves=tree.findLeaves(ele.value);
                              for (index = 0, len = someleaves.length; index < len; ++index) {
                                    if (someleaves[index]== null){
                                    someleaves.selected = true;
                                    
                                    }
                                    else{
                                    someleaves[index].selected = true;
                                    }
                                    
                                }
                            }
                            else {
                              tree.branches.E.cascadeFlag('selected', false);
                              tree.branches.E.cascadeFlag('highlighted', false);
                              tree.draw();
                                }
                            }
                        
                        $(document).on('click','#pc-buttons .btn', {} ,function(e){
                            $('#pc-buttons .btn').removeClass('btn-info');
                            $('#pc-buttons .btn').addClass('btn-default');
                            $(this).addClass('btn-info');
                            tree.setTreeType(this.id);
                        });
                        
                        
                        
                      </script>
            <script>
            

                    var yourDiv = document.getElementById('snippetDiv');
                    var menuDiv = document.createElement('div');
                    var msaDiv = document.createElement('div');
                    yourDiv.appendChild(menuDiv);
                    yourDiv.appendChild(msaDiv);

                    /* global yourDiv */
                    var opts = {
                      el: msaDiv,
                      importURL: '"""+os.path.join((os.path.relpath(alignFileName,alignFileName)),(os.path.basename(gene)).replace(".fasta","_aligned.fasta"))+"""',
                      colorscheme: {"scheme": "nucleotide"},
                    };
                    
                    opts.vis = {
                      scaleslider:true,
                    };
                    
                    var m = msa(opts);
                    m.render()

                    var defMenu = new msa.menu.defaultmenu({
                      el: menuDiv,
                      msa: m
                    });
                    defMenu.render();
                    
                    function addColumnFilter(menu){    
                        var msa = menu.msa;
                        var hidden = [];    threshold = 100 / 100;
                       var maxLen = msa.seqs.getMaxLength();
                       var hidden = [];
                       // TODO: cache this value
                       var conserv = msa.g.stats.scale(msa.g.stats.conservation());

                       var end = maxLen - 1;
                       for (var i = 0; 0 < end ? i <= end : i >= end; 0 < end ? i++ : i--) {
                         if (conserv[i] == threshold) {
                             hidden.push(i);
                         }
                       }
                       return msa.g.columns.set("hidden", hidden);}
                      
                     function addColumnFilter2(menu){    
                        var msa = menu.msa;
                        var seqToHideStr = prompt("Enter Seq to show", 'eg allele1,allele2');
                        var seqToHide = seqToHideStr.split(',');

                       // TODO: cache this value
                       var i = 0;
                       var len = (m.seqs.models.length)-1

                       for (; len >= i; ) {
                            var auxSeq=msa.seqs.models[len]
                            if ( $.inArray(auxSeq.attributes.name,seqToHide) < 0) {
                                msa.seqs.remove(msa.seqs.at(len));

                                }
                            len--;
                            }
                    
                       }
                      
                     $($('div').find('ul')[2]).prepend('<li id="removePoly">Hide Non Polymorphic Sites</li>');    $('#removePoly').click(function(){
                            addColumnFilter(defMenu);
                        });
                    
                    $($('div').find('ul')[2]).prepend('<li id="removeSpecSeq">Show specific set of seqs</li>');    $('#removeSpecSeq').click(function(){
                            addColumnFilter2(defMenu);
                        });	
                        
                    delete defMenu.views['10_import'];
                    </script>""")
            f.write("""<script type="text/javascript">
                    $("#button2").click(function(){
                      $("#snippetDiv").css({"display":"block"});
                      $("#histdiv").css({"display":"none"});
                      $("#phylocanvas").css({"display":"none"});
                    }); 
                    </script>
        
                    <script type="text/javascript">
                        $("#button3").click(function(){
                          $("#histdiv").css({"display":"block"});
                          $("#phylocanvas").css({"display":"none"});
                          $("#snippetDiv").css({"display":"none"});
                        }); 
                        </script>
            
                    <script type="text/javascript">
                        $("#button1").click(function(){
                          $("#phylocanvas").css({"display":"block"});
                          if ($("#phylocanvas__canvas").length < 1){
                            doPhylocanvas();}
                          $("#snippetDiv").css({"display":"none"});
                          $("#histdiv").css({"display":"none"});
                        }); 
                        </script>""")



            alleleSizes=str(json.dumps([alleleNames,alleleSizes]))
            alleleSizesNotMultiple=str(json.dumps([alleleSizesNotMultipleNames,alleleSizesNotMultiple]))
            alleleSizesTransError=str(json.dumps([alleleSizesTransErrorNames,alleleSizesTransError]))

            f.write("<script type='text/javascript'>var jsonScat1 ="+alleleSizes+""";
            
            
                        
                        var listTrace=[];
                        var trace1 = {
                                    x: jsonScat1[0],
                                      y: jsonScat1[1],
                                      text: jsonScat1[0],
                                      name: 'CDS alleles',
                                      mode: 'markers',
                                      type: 'scattergl'
                                    };
                        listTrace.push(trace1);""")
            if 	len(alleleSizesNotMultipleNames)>0:


                f.write("var jsonScat2 ="+alleleSizesNotMultiple+""";
                
                        var trace2 = {
                                    x: jsonScat2[0],
                                      y: jsonScat2[1],
                                      text: jsonScat2[0],
                                      name: 'Alleles not multiple of 3',
                                      mode: 'markers',
                                      type: 'scattergl'
                                    };
                        listTrace.push(trace2);
                        """)
            if 	len(alleleSizesTransErrorNames)>0:

                f.write("var jsonScat3 ="+alleleSizesTransError+""";
                            var trace3 = {
                                    x: jsonScat3[0],
                                      y: jsonScat3[1],
                                      text: jsonScat3[0],
                                      name: 'Non translatable alleles',
                                      mode: 'markers',
                                      type: 'scattergl'
                                    };
                        listTrace.push(trace3);
                        """)

            f.write("""var layout = {
                                      title: 'Allele size scatter plot',
                                      yaxis: {title: "DNA bp allele length"},
                                      xaxis: {
                                                title: "Allele number",
                                              },
                                    };
                        Plotly.newPlot('histdiv', listTrace,layout);
            
            
                    </script>"""+"</body></html>")


        if not(not listnotMultiple and not listnotStart and not listStopc and not listOther):
            if not listnotMultiple:
                listnotMultiple.append("-")
            if not listnotStart:
                listnotStart.append("-")
            if not listStopc:
                listStopc.append("-")
            if not listOther:
                listOther.append("-")
            toPrintCDSStats+= os.path.basename(gene)+"\t"+','.join(listnotMultiple)+"\t"+','.join(listnotStart)+"\t"+','.join(listStopc)+"\t"+','.join(listOther)+"\n"
    #print str(stopc) + " alleles have stop codons inside"
    #print str(notStart) + " alleles don't have start codons"
    print("total of alleles : " + str(totalalleles))

    with open('non_cds_alleles.tsv', "w") as f:
        f.write(toPrintCDSStats)
    print(toPrintCDSStats)

    return statsPerGene

if __name__ == "__main__":
    main()

