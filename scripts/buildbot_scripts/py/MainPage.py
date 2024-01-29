# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#================================== NEW ================================

from buildbot.status.web.base import HtmlResource
#import mainwork
import os

class HomePage(HtmlResource):
    title = "ICON Buildbot"
#    global test_Info=""
#    global change_Info=""
#    global exp_Info=""

    def head(self, request):
        h='''
<link href="archive.css" rel="stylesheet" type="text/css" />
        '''
        return h

    def main_page(self, request):
      global change_Info
      
# Create a emty Experiment List type

      ExpList = []
      ExpListInfo = []

# Read all files in the "icon_run" directory and write the filename which starts with "exp.test_" into the
# Experiment List type

      for filename in os.listdir("icon_run"):
	if filename.find("exp.test_") == 0:
          ExpDict = {"Description" : " ","Model" : " "}
          ExpDict["ExpName"] = filename.lstrip("exp.")
            
          file = open("icon_run/" + filename)
          while 1:
            line = file.readline()
            if not line:
              break
             
            if line.find("_Description_") >= 0:
	      tmp,Info = line.split("_bb_table_Description_",1)
              Info.replace("\n","")
              ExpDict["Description"] = Info.strip(" ")
                
            if line.find("_Model_") >= 0:
	      tmp,Info = line.split("_bb_table_Model_",1)
              Info.replace("\n","")
              ExpDict["Model"] = Info.strip(" ")
                
          file.close
          ExpList.append(filename.lstrip("exp."))
          ExpListInfo.append(ExpDict)
          
          file = open("icon_run/svn_info")
          while 1:
            line = file.readline()
            if not line:
              break
          
            if line.find("Revision:") >= 0:
	      tmp,tmp_Rev = line.split("Revision:",1)
              Rev = tmp_Rev.strip(" ").replace("\n","")
              
          file.close
           
            
# Sort the Experiment List type

      ExpList.sort()
	
# Write the upper part of the main page

      data  = '''
<h1>ICON Buildbot</h1>
Buildbot starts every night a test suite on a set of builders differing in machines, compilers 
and parallelization setup. A standard Buildbot test unpacks the model code and scripts from the  
repositors, configures and compiles the codes for the selected builder, and runs a series of 
standard test experiments and related post-processings. Figures resulting from the 
post-processings are archived and are uploaded to web pages for comparison with reference 
figures.<p>
      '''
	
# Write the Status result part to the main page

      data += '''
<h2>Status results</h2>
<table border="1" cellspacing="5px" cellpadding="2" >
   <thead>
     <tr valign="baseline">
       <td><i>Display</i></td>
       <td><i>What it shows</i></td>
     </tr>
   </thead>
   <tfoot></tfoot>
   <tbody>
     <tr valign="baseline">
       <td class=\"status_result_a\"><a href="waterfall?show_events=false&show_time=604800">Waterfall</a></td>
       <td class=\"status_result_b\">Status of each step of the Buildbot tests for the past 7 days</td>
     </tr>

     <tr valign="baseline">
       <td class=\"status_result_a\"><a href="grid?width=10">Grid</a></td>
       <td class=\"status_result_b\">Status of the complete Builbot test for the 10 latest tested revisions</td>
     </tr>

     <tr valign="baseline">
       <td class=\"status_result_a\"><a href="one_box_per_builder">Latest build</a></td>
       <td class=\"status_result_b\">Status of the complete Builbot test for the latest test on each builder</td>
     </tr>

     <tr valign="baseline">
       <td class=\"status_result_a\"><a href="one_line_per_build">List</a></td>
       <td class=\"status_result_b\">Status of the complete Builbot test for the latest 15 tests</td>
     </tr>

     <tr valign="baseline">
       <td class=\"status_result_a\"><a href="buildslaves">Build slaves</a></td>
       <td class=\"status_result_b\">Information on machines ("slaves"), builders, etc.</td>
     </tr>
   </tbody>
   </table>
      '''

# Create the experiment List as a table

      
      data += "<h2>Latest nightly test results</h2>\n"
      data += "<h3>Source : trunk/icon-dev &#64; &lt;r" + Rev + "&gt; </h3>"
      data += "<table border=\"1\" cellspacing=\"5px\" cellpadding=\"2\">\n"
      data += "  <thead>\n"
      data += "    <tr valign=\"baseline\">\n"
      data += "      <td width=\"250\"><i>Experiment (exp)</i></td>\n"
      data += "      <td width=\"380\"><i>Description</i></td>\n"
      data += "      <td width=\"200\"><i>Model</i></td>\n"
      data += "    </tr>\n  </thead>"
      data += "  <tfoot></tfoot>\n"

# Build HTML Table info

      data += "  <tbody>\n"
      for exp in ExpList:
        for e in ExpListInfo:
          if exp == e.get('ExpName'):
	    break
	      
 	data += "    <tr valign=\"baseline\">\n"
        data += "      <td class=\"last_night\"><a href=\"plot?exp=" 
	data += e.get('ExpName') + "&modus=nightly\">" + e.get('ExpName') + "</a></td>\n"
        data += "      <td class=\"last_night\">" + e.get('Description') + "</td>\n"
        data += "      <td>" + e.get('Model') + "</td>\n"
        data += "   </tr>"
      
      data += "  </tbody>\n</table>\n"
#============================================

# Archive

#============================================

      data += "<h2>Archived test results</h2>"
#      data += "<a href=\"plot?exp="+ExpList[1]+"&modus=archive\"><h3>Archive viewer</h3></a>"
      data += "<a href=\"plot?exp=test_hat_jww&modus=archive\"><h3>Archive viewer</h3></a>"
      
#============================================

# footer

#============================================

      data += "<h3>About BuildBot</h3>\n"
      data += "<ul>\n"
      data += "  <li><a href=\"about\">Buildbot version used here</a></li>\n"
      data += "  <li>Buildbot home page: <a href=\"http://trac.buildbot.net\" target=\"_blank\">http://trac.buildbot.net</a></li>\n"
      data += "</ul>\n"

      return data
  
   
    def get_info(self, request):
        global change_Info
	
        change_Info="NotSet"
        
	if "status" in request.args:
          try:
            change_Info = request.args["status"][0]
          except ValueError:
            pass
	
        return None
    
    def body(self, request):
        global change_Info
	self.get_info(request)
        data = self.main_page(request)	
        return data
#================================== NEW ================================
