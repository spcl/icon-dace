# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#================================== NEW ================================

from buildbot.status.web.base import HtmlResource
from twisted.web.util import Redirect
import os, time
from buildbot import version, util
from datetime import datetime, timedelta

class EXP_plot(HtmlResource):
    title = "Exp Plots"
    RevisionNr = "0000"

    def head(self, request):
        h='''
<link href="archive.css" rel="stylesheet" type="text/css" />

<script type='text/javascript'>
 function anzeigen(das,was){
 if (document.getElementById(das).style.display=='none') {
  document.getElementById(das).style.display='block'; 
  was.src="open.gif";
 }
 else {
  document.getElementById(das).style.display='none';
   was.src='closed.gif';
  }
}
 </script>
        '''
        return h
						
   
    def body(self, request):
        global exp_plot_Info
	global exp_file_Info
        global l_nightly
        global modus
        global RevisionNr
	global selected_builder
	global selected_build
	global selected_branch
	global selected_Date_from
	global selected_Date_to
	global ExpPlotError
        global DEBUG
	
	self.get_info(request)
	
        status = self.getStatus(request)	
        Archive_Button_Dict = {}
        ExpPlotError = 0
        DEBUG = False
# -----------------------------------------------------------
        def add_file(p,e,Dict):
	  ret = False
	  file_L = []
	  p += e + "/plots/"
	  for f in os.listdir(p):
            if DEBUG:
              print "== add_file =="
              print f

	    file_L.append(f)
       	    Archive_Button_Dict['file'].append(f)
	    ret = True
	    
# If ret == True we got a file in the list
	  if ret:
            Dict[e] = file_L
	    
          return ret

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
	def add_exp(p,b,Dict):
	  ret = False
	  exp_D = {}
	  p += b + "/"
	  for e in os.listdir(p):
            if DEBUG:
              print "== add_exp =="
              print e

	    if e.find("test_") == 0:
#	      if e == exp_plot_Info or exp_plot_Info == "all":
              if add_file(p,e,exp_D):
		Archive_Button_Dict['exp'].append(e)
	        ret = True

# If ret == True we found an existing experiment mit a file list
	  if ret:
            Dict[b] = exp_D
          
	  return ret

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
	def add_build(p,c,Dict):
	  ret = False
	  build_D = {}
	  p += c + "/"
	  for b in os.listdir(p):
            if DEBUG:
              print "== add_build =="
              print b
	    if add_exp(p,b,build_D):
	      Archive_Button_Dict['build'].append(b)
	      ret = True

# If ret == True we found a build with correct sub-information
	  if ret:
            Dict[c] = build_D
          
	  return ret

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
	def add_comp(p,r,r_Dict):
	  ret = False
	  comp_D = {}
	  if r != "":
	    p += r + "/"
	  else:
	    r = "trunk+icon-dev"
#          print "p " + p
#          print "r " + r
	    
	  for c in os.listdir(p):
            if DEBUG:
              print "== add_comp =="
              print c
	    if add_build(p,c,comp_D):
  	      Archive_Button_Dict['comp'].append(c)
	      ret = True

	  if ret:
            r_Dict[r] = comp_D
          
          return ret

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	def add_branch(p,D,d_Dict):
	  ret = False
	  rev_D = {}
	  p += D + "/"
	  
	  for b in os.listdir(p):
            if DEBUG:
              print "== add_branch =="
              print b

            B = b
	    if b.find("trunk") == 0 or b.find("tags") == 0 or b.find("branch") == 0:
              if DEBUG:
                print "== if case =="
	      if add_comp(p,b,rev_D):
	        Archive_Button_Dict['branch'].append(b)
	        ret = True
	    else:
              if DEBUG:
                print "== else  case =="

              B = "trunk+icon-dev"
	      if add_comp(p,"",rev_D):
	        Archive_Button_Dict['branch'].append(B)
	        ret = True
	      break
	
	  if ret:
            d_Dict[D] = rev_D

          return ret
	      
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	def add_rev(p,D,d_Dict):
	  ret = False
	  rev_D = {}
	  p += D + "/buildbot/"
	  for r in os.listdir(p):
            if DEBUG:
              print "== add_rev =="
              print r
            if (r >= Rev_from) and (r <= Rev_to):
              if DEBUG:
                print "== r OK =="
	      if add_branch(p,r,rev_D):
	        Archive_Button_Dict['rev'].append(r)
	        ret = True
	        
	  if ret:
            d_Dict[D] = rev_D
          
          return ret
	      
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        def create_Exp_Dict(d,r,br,c,b,eD):
	  build_D  = {}
	  comp_D   = {}
	  rev_D    = {}
	  branch_D = {}
	  date_D   = {}
	  
	  build_D[b]   = eD
	  comp_D[c]    = build_D
	  branch_D[br] = comp_D
	  rev_D[r]     = branch_D
	  date_D[d]    = rev_D
	  return date_D 
	  
        def create_date_Dict(r,br,c,b,eD):
	  rev_D    = {}
	  branch_D = {}
	  comp_D   = {}
	  build_D  = {}
	  
	  build_D[b]   = eD
	  comp_D[c]    = build_D
	  branch_D[br] = comp_D
	  rev_D[r]     = branch_D
	  return rev_D 
        
	def create_rev_Dict(br,c,b,eD):
	  branch_D = {}
	  build_D = {}
	  comp_D  = {}

          build_D[b] = eD
	  comp_D[c] = build_D
	  branch_D[br] = comp_D

	  return branch_D 

	def create_branch_Dict(c,b,eD):
	  build_D = {}
	  comp_D  = {}

          build_D[b] = eD
	  comp_D[c] = build_D

	  return comp_D 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	def create_comp_Dict(b,eD):
	  build_D = {}
	  build_D[b] = eD
	  return build_D 
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	def Exp_add(Exp_Dict,d,r,Br,c,b,eD):
	   
	   for e in eD:
             if e in Exp_Dict:
	       if d in Exp_Dict[e]:
	         if r in Exp_Dict[e][d]:
	           if Br in Exp_Dict[e][d][r]:
	             if c in Exp_Dict[e][d][r][Br]:
	               if b in Exp_Dict[e][d][r][Br][c]:
                         if DEBUG:
	                   print "build existiert " + c
	                 Exp_Dict[e][d][r][Br][c][b] = eD[e]
		       else:
	                 Exp_Dict[e][d][r][Br][c][b] = eD[e]
	             else:
	               Exp_Dict[e][d][r][Br][c] = create_comp_Dict(b,eD[e])
	           else:
	             Exp_Dict[e][d][r][Br] = create_branch_Dict(c,b,eD[e])
	         else:
	           Exp_Dict[e][d][r] = create_rev_Dict(Br,c,b,eD[e])
	       else:
	         Exp_Dict[e][d] = create_date_Dict(r,Br,c,b,eD[e])
             else:
	       Exp_Dict[e] = create_Exp_Dict(d,r,Br,c,b,eD[e])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
        def Comp_create_date(d,r,Br,eD):
	  branch_D  = {}
	  rev_D     = {}
	  date_D    = {}

	  branch_D[Br] = eD
          rev_D[r] = branch_D
	  date_D[d] = rev_D
	  return date_D 
	
        def Comp_create_rev(r,Br,eD):
	  branch_D  = {}
	  rev_D   = {}

	  branch_D[Br] = eD
          rev_D[r] = branch_D
	  return rev_D 

        def Comp_create_Branch(Br,eD):
	  branch_D  = {}
	  
	  branch_D[Br] = eD
	  return branch_D 
	
	def Comp_add(Comp_Dict,d,r,Br,cD):
	  for c in cD:
           if c in Comp_Dict:
             if d in Comp_Dict[c]:
               if r in Comp_Dict[c][d]:
                 if Br in Comp_Dict[c][d][r]:
                   Comp_Dict[c][d][r][Br] = cD[c]
	         else:
                   Comp_Dict[c][d][r][Br] = cD[c]
	       else:
	         Comp_Dict[c][d][r] = Comp_create_Branch(Br,cD[c])
	     else:
	       Comp_Dict[c][d] = Comp_create_rev(r,Br,cD[c])
	   else:
	     Comp_Dict[c] = Comp_create_date(d,r,Br,cD[c])
# -----------------------------------------------------------
	
	if exp_plot_Info == "NotSet":
	  data = self.footer(status,request)
          return data

	Archive_Button_Dict['date']   = []
	Archive_Button_Dict['rev']    = []
	Archive_Button_Dict['branch'] = []
	Archive_Button_Dict['comp']   = []
	Archive_Button_Dict['build']  = []
	Archive_Button_Dict['exp']    = []
	Archive_Button_Dict['file']   = []
	
#        ExpListInfo = []

	RevisionNr = "0000"
        file = open("icon_run/svn_info")
        while 1:
          line = file.readline()
          if not line:
            break
          
          if line.find("Revision:") >= 0:
            tmp,tmp_Rev = line.split("Revision:",1)
            RevisionNr = tmp_Rev.strip(" ").replace("\n","")
            
        file.close

        if l_nightly:
	  selected_branch = "trunk+icon-dev"
          Rev_from = RevisionNr
	  Rev_to = Rev_from
 
	  Date_to = time.strftime("%Y-%m-%d",time.localtime(util.now()))
#ws	  Date_to = "2012-06-21"
	  now_time = time.strftime("%H%M",time.localtime(util.now()))
	      
	  yesterday = datetime.now() - timedelta(days=1)
	      
	  if now_time > "2200":
#ws	    Date_from = "2012-06-20"
            Date_from = Date_to
	  else:
            Date_from = yesterday.strftime("%Y-%m-%d")
#ws	    Date_from = "2012-06-20"
	
	else:
          Rev_from  = str(int(RevisionNr)-100)
          Rev_to    = RevisionNr
          now_min_10d = datetime.now() - timedelta(days=30)
          Date_from = now_min_10d.strftime("%Y-%m-%d")
          Date_to = time.strftime("%Y-%m-%d",time.localtime(util.now()))
          
          if selected_Date_from != "NotSet":
	    Date_from = selected_Date_from
	    
          if selected_Date_to != "NotSet":
	    
	    if selected_Date_to < time.strftime("%Y-%m-%d",time.localtime(util.now())):
	      Date_to = selected_Date_to
            else:
	      Date_to = time.strftime("%Y-%m-%d",time.localtime(util.now()))
	      
	  if selected_Rev_from != "NotSet":
	    Rev_from = selected_Rev_from  
	    
          if selected_Rev_to != "NotSet":
	    Rev_to = selected_Rev_to
	    
# Check if from time greater to time. If yes change it.
          if Date_from > Date_to:
	    t         = Date_to
	    Date_to   = Date_from
	    Date_from = t

# Check if from revision greater to revision. If yes change it.
          if int(Rev_from) > int(Rev_to):
	    t         = Rev_to
	    Rev_to   = Rev_from
	    rev_from = t
#========================================================================================================
#
# Begin searching for plots
#
#========================================================================================================

# Build new emty date_Dict type

        

# Build new emty rev_Dict type and save it in the date_Dict under the correct Date
	date_Dict = {}
        
	p_date = "public_html/archive/"
	for DATE in os.listdir(p_date):
          if DEBUG:
             print DATE
	  if DATE.find("20") == 0:
            if DEBUG:
              print "DATE.find(\"20\") == 0:"
              print Date_from 
              print Date_to
              print DATE
            if (DATE >= Date_from) and (DATE <= Date_to):
              if DEBUG:
                print DATE
	      if add_rev(p_date,DATE,date_Dict):
	        Archive_Button_Dict['date'].append(DATE)
	
# Build a new List where the Experiment Name ist the hightest index
	Exp_Dict = {}
	for d in date_Dict:
	  for r in date_Dict[d]:
	    for B in date_Dict[d][r]:
	      for c in date_Dict[d][r][B]:
	        for b in date_Dict[d][r][B][c]:
	          Exp_add(Exp_Dict,d,r,B,c,b,date_Dict[d][r][B][c][b])
             
# Build a new List where the Experiment Name ist the hightest index
	Comp_Dict = {}
	for d in date_Dict:
	  for r in date_Dict[d]:
	    for b in date_Dict[d][r]:
	      Comp_add(Comp_Dict,d,r,b,date_Dict[d][r][b])

        Archive_Button_Dict['date']   = list(sorted(set(Archive_Button_Dict['date'])))
 	Archive_Button_Dict['rev']    = list(sorted(set(Archive_Button_Dict['rev'])))
 	Archive_Button_Dict['branch'] = list(sorted(set(Archive_Button_Dict['branch'])))
	Archive_Button_Dict['comp']   = list(sorted(set(Archive_Button_Dict['comp'])))
	Archive_Button_Dict['build']  = list(sorted(set(Archive_Button_Dict['build'])))
	Archive_Button_Dict['exp']    = list(sorted(set(Archive_Button_Dict['exp'])))
	Archive_Button_Dict['file']   = list(sorted(set(Archive_Button_Dict['file'])))
        if DEBUG:
          print "==========================================="
          print Archive_Button_Dict['date']
          print Archive_Button_Dict['rev']
          print Archive_Button_Dict['branch']
          print Archive_Button_Dict['comp']
          print Archive_Button_Dict['build']  
          print Archive_Button_Dict['exp']
          print Archive_Button_Dict['file']
          print "==========================================="

	
        if not Archive_Button_Dict['rev']:
	  ExpPlotError = 1
#err	  data = "Revision List for this timeperode ist emty"
#err	  return data
	  
        if not l_nightly and ExpPlotError == 0:
	  min_Revision = Archive_Button_Dict['rev'][0]
	  max_Revision = Archive_Button_Dict['rev'][-1]
	  
	  if Rev_from <= min_Revision:
	    Rev_from = min_Revision
	    
	  if Rev_to >= max_Revision:
	    Rev_to   = max_Revision
	    
	

# Build a plot list which contains only the plots of the acual experiment

	t_file = []
	for Fi in Archive_Button_Dict['file']:
	  if Fi.find(exp_plot_Info+"_") >= 0:
	    t_file.append(Fi) 
            if DEBUG:
	      print Fi + " (" + exp_plot_Info + "_) included in the list"
	  else:
            if DEBUG:
	      print Fi + " (" + exp_plot_Info + "_) not included in the list"
	    
        if len(t_file) != 0:
	  Archive_Button_Dict['file'] = t_file
	else:
	  if ExpPlotError == 0:
	    ExpPlotError = 2
#          date = "<h2> No Plots available for experiment: <h1>" + exp_plot_Info + "</h1></h2>" 
#	  return date
	
	
	
	Plot_name_length = 290.
	for s in Archive_Button_Dict['file']:
	  if len(s)*6.9 > Plot_name_length:
	    Plot_name_length = len(s)*6.9
          if DEBUG:
	    print "String: " + s + " " + str(Plot_name_length) + " " + str(len(s)*7.0)
	    

#========================================================================================================
#
# End searching for plots
#
#========================================================================================================

#==============================================
# Begin of whole div
        

        data = ""
        
	
        data +=  "<div id=\"page_menu\" >\n"

# Begin of left side  div
	
	data += "<div style=\"float:left; padding:3px; margin:5px width:380px;\">\n"

# Begin of ICON Buildbot Archive div
        
        data += "<div id=\"arch_menu\" >\n"
	data += "<h1>ICON Buildbot Archive</h1>\n"
	

# Begin of form 
	
	data += "<form name=\"date_form\" method=\"POST\" action=\"plot/select\" class=\"command selectfile\">\n"

# Begin of table 
        
	data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\""+str(Plot_name_length+90)+"\">\n"

# Define cell widths

        data += "<colgroup><col width=\"90\"><col width=\""+str(Plot_name_length)+"\"></colgroup>\n"
	
#     ----------   Date     ----------   

        data += "<tr>\n"
#	data += "  <td style=\"text-align:left;\"><b>Date:</b> <br> from/to</td>\n"
	data += "  <td style=\"text-align:left;\"><b>Date:</b></td>\n"
	data += "  <td style=\"text-align:left;\">"
	
# Date Year from
        if l_nightly:
	  data += Date_from 
	else:
	  data += "    <input name=\"date_from\" type=\"text\" value=\"" + Date_from.split(' ')[0] + "\" size=\"10\" maxlength=\"10\" >\n"
	  

# Date Year to
        if l_nightly:
	  if Date_to != Date_from:
	    data += "&nbsp; to &nbsp" + Date_to
	else:
	  data += "&nbsp; to &nbsp"
	  data += "    <input name=\"date_to\" type=\"text\" value=\"" + Date_to.split(' ')[0] + "\" size=\"10\" maxlength=\"10\">\n"
	  if ExpPlotError == 0:
	    data += "    <input type=\"submit\" name=\"Date_button\" value=\">\">"
	
	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Revision     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\"><b>Revision:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
# Revision from 

        if l_nightly:
	  data += Rev_from
	else:
	  data += "    <input name=\"rev_from\" type=\"text\" value=\"" + Rev_from.split(' ')[0] + "\" size=\"10\" maxlength=\"10\">\n"
        
# Revision to 

        if l_nightly:
	  if Rev_to != Rev_from:
	    data += "&nbsp; to &nbsp" + Rev_to
	else:
	  data += "&nbsp; to &nbsp"
	  data += "    <input name=\"rev_to\" type=\"text\" value=\"" + Rev_to.split(' ')[0] + "\" size=\"10\" maxlength=\"10\">\n"
	  if ExpPlotError == 0:
	    data += "    <input type=\"submit\" name=\"Rev_button\" value=\">\">"

	data += "  </td>\n"
	data += "</tr>\n"
	
# Branch 

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Branch:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"

        if l_nightly:
	  data += selected_branch.replace('+', '/')
	else:
	  if ExpPlotError == 0:
	    data += "    <select  name=\"branch\" onChange=\"document.date_form.submit()\">\n"
	  else:
	    data += "    <select  disabled=\"disabled\" name=\"branch\" onChange=\"document.date_form.submit()\">\n"
	    
	  if selected_branch == "all":
	    data += "    <option selected> all </option>\n"
	  else:
	    data += "    <option> all </option>\n"
	    
	  for Br in Archive_Button_Dict['branch']:
	    if selected_branch == Br:
	      data += "    <option selected>"+ Br.replace('+', '/') + "</option>\n"
	    else:
	      data += "    <option>"+ Br.replace('+', '/') + "</option>\n"
	      
	  data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"

#     ----------    Computer     ----------   

	if not Comp_Dict.has_key(selected_builder):
	  selected_builder = "all"
	  
        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Builder:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
        if l_nightly:
	  data += "all"
	else:
	  if ExpPlotError == 0:
	    data += "    <select  name=\"builder\" onChange=\"document.date_form.submit()\">\n"
	  else:
	    data += "    <select   disabled=\"disabled\" name=\"builder\" onChange=\"document.date_form.submit()\">\n"
	  if selected_builder == "all":
	    data += "    <option selected> all </option>\n"
	  else:
	    data += "    <option> all </option>\n"
	
	  for Co in Archive_Button_Dict['comp']:
	    if selected_builder == Co:
	      data += "    <option selected>"+ Co + "</option>\n"
	    else:
	      data += "    <option>"+ Co + "</option>\n"
	  data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Build Nr.     ----------   
        if selected_builder != "all":
	  t_build = []
          for Da in Comp_Dict[selected_builder]:
            for Re in Comp_Dict[selected_builder][Da]:
	      for Br in Comp_Dict[selected_builder][Da][Re]:
	        for Bu in Comp_Dict[selected_builder][Da][Re][Br]:
		  t_build.append(Bu)
		    
          Archive_Button_Dict['build'] = list(sorted(t_build))

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Build:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"

        if l_nightly:
	  data += "all"
	else:
	  if ExpPlotError == 0:
	    data += "    <select  name=\"build\" onChange=\"document.date_form.submit()\">\n"
	  else:
	    data += "    <select  disabled=\"disabled\" name=\"build\" onChange=\"document.date_form.submit()\">\n"
	  if selected_build == "all":
	    data += "    <option selected> all </option>\n"
	  else:
	    data += "    <option> all </option>\n"
	    
	  for Bu in Archive_Button_Dict['build']:
	    if selected_build == Bu:
	      data += "    <option selected>"+ Bu + "</option>\n"
	    else:
	      data += "    <option>"+ Bu + "</option>\n"
	      
	  data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Exp. name     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Experiment:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"

#old bug begin	
	if ExpPlotError == 0:
	  data += "    <select  name=\"exp\" onChange=\"document.date_form.submit()\">\n"
	else:
	  data += "    <select  disabled=\"disabled\" name=\"exp\" onChange=\"document.date_form.submit()\">\n"
# old ug end
# new  begin
#        data += "    <select  name=\"exp\" onChange=\"document.date_form.submit()\">\n"
# new end
#ToDo   Include the correct showing of '_' in experment name and plot-name

	for Ex in Archive_Button_Dict['exp']:
	  if Ex == exp_plot_Info:
	    data += "      <option selected>" + Ex + "</option>\n"
	  else:
	    data += "      <option>" + Ex + "</option>\n"
	
	
	data += "  </select>\n"
	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Exp. file     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Plot:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
	if ExpPlotError == 0:
          data += "    <select  name=\"file\" onChange=\"document.date_form.submit()\">\n"
	else:
          data += "    <select disabled=\"disabled\"  name=\"file\" onChange=\"document.date_form.submit()\">\n"
	  
        if exp_file_Info == "NotSet" and len(Archive_Button_Dict['file']) > 0:
	  exp_file_Info = Archive_Button_Dict['file'][0].strip(".png")

        for f in Archive_Button_Dict['file']:
	  if f.strip(".png") == exp_file_Info:
	    data += "      <option selected>" + f.strip(".png") + "</option>\n"
	  else:
	    data += "      <option>" + f.strip(".png") + "</option>\n"

	data += "    </select>\n"
	data += "  </td>\n"
	data += "</tr>\n"
        
#     ----------    Search OK Button     ----------   

	data += "    <input type=\"hidden\" name=\"modus\" value=\""+ modus + "\">\n"
	data += "    <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\">\n"
	data += "    <input type=\"hidden\" name=\"builder\" value=\"" + selected_builder + "\">\n"
	data += "    <input type=\"hidden\" name=\"build\" value=\"" + selected_build + "\">\n"
	data += "    <input type=\"hidden\" name=\"branch\" value=\"" + selected_branch.replace('+', '/') + "\">\n"

# End of ICON Buildbot Archive and Form area

        data += "</table>\n"
        data += "</form>\n"

# End of ICON Buildbot Archive div
        
	data += "</div>\n"

#==================================================
#
# Begin of Available Plot List div
#
#==================================================

        data += "<div id=\"arch_file_name\">\n"
        data += "<h1>Available Plot List</h1>\n"
        data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
	data += "<tr><td style=\"text-align:left; valign:top;\">\n"
	
	e = 0
	r = 0
	z = 0
	d = 0
	c = 0
	b = 0
	f = 0
	
        for Ex in Archive_Button_Dict['exp']:
#        for Ex in Exp_Dict:
	  e += 1
          data += "\n<!--- " + str(e) + ". Experiment " + Ex + " -->\n"
          data += "\n<a href=\"javascript:anzeigen('date_" + str(e) 
	  data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	  data += "  <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + Ex + "</span>\n"
	  data += "</a><br>\n"
	  
	  data += "<div id=\"date_" + str(e) + "\" class=\"Rand\" style=\"display: none;\">\n"
	  for Da in sorted(Exp_Dict[Ex], reverse=True):
	    d += 1 
            data += "  <!--- " + str(d) + ". Date " + Da + " -->\n"
          
	    data += "  <a href=\"javascript:anzeigen('rev_" + str(d) 
	    data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	    data += "    <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" 
	    data += Da + "</span>\n"
	    data += "  </a><br>\n"
#------------------------------
	    data += "  <div id=\"rev_" + str(d) + "\" class=\"Rand\" style=\"display: none;\">\n"
	    for Re in sorted(Exp_Dict[Ex][Da], reverse=True):
	      r += 1 
              data += "    <!--- " + str(r) + ". Revision " + Re + " -->\n"
          
	      data += "    <a href=\"javascript:anzeigen('branch_" + str(r) 
	      data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	      data += "      <img name=\"bild_1\" src=\"closed.gif\" Brclass=\"Bild\"> <span style=\"color: #000000;\">" 
	      data += Re + "</span>\n"
	      data += "    </a><br>\n"
#------------------------------
	      data += "    <div id=\"branch_" + str(r) + "\" class=\"Rand\" style=\"display: none;\">\n"
#new begin
	      for Br in sorted(Exp_Dict[Ex][Da][Re]):
	        z += 1 
                data += "      <!--- " + str(z) + ". Branch " + Br + " -->\n" 
	        data += "      <a href=\"javascript:anzeigen('comp_" + str(z) 
		data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	        data += "        <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" 
		data += Br.replace('+', '/') + "</span>\n"
	        data += "      </a><br>\n"
#------------------------------
	        data += "      <div id=\"comp_" + str(z) + "\" class=\"Rand\" style=\"display: none;\">\n"
	        for Co in sorted(Exp_Dict[Ex][Da][Re][Br]):
	          c += 1 
                  data += "      <!--- " + str(c) + ". Computer " + Co + " -->\n" 
	          data += "      <a href=\"javascript:anzeigen('build_" + str(c) 
		  data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	          data += "        <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" 
		  data += Co + "</span>\n"
	          data += "      </a><br>\n"
#------------------------------
	          data += "      <div id=\"build_" + str(c) + "\" class=\"Rand\" style=\"display: none;\">\n"
	          for Bu in sorted(Exp_Dict[Ex][Da][Re][Br][Co], reverse=True):
	            b += 1 
                    data += "        <!--- " + str(b) + ". Build " + Bu + " -->\n" 
	            data += "        <a href=\"javascript:anzeigen('file_" + str(b) 
		    data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	            data += "          <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" 
		    data +=  Bu + "</span>\n"
	            data += "        </a><br>\n"
#------------------------------
	            data += "      <div id=\"file_" + str(b) + "\" class=\"Rand\" style=\"display: none;\">\n"
	            for Fi in sorted(Exp_Dict[Ex][Da][Re][Br][Co][Bu]):
	              f += 1 
                      data += "          <!--- " + str(f) + ". File " + Co + " -->\n" 
                    
                      data += "<input type=\"checkbox\" name=\"plot_check\" value=\"" + Fi.strip(".png") + "\" id=\"check1\""
		      if Fi.strip(".png") == exp_file_Info:
		        data += "checked=\"checked\""
		      else:
			data += "disabled=\"disabled\""

		      
		      if l_nightly:
		        data += "disabled=\"disabled\""
		      data += ">"
		      data += "<label for=\"check1\">"+ Fi.strip(".png") + "</label><br />\n"
                    data += "        </div>\n"
#------------------------------

                  data += "      </div>\n"
#------------------------------
	        data += "    </div>\n"
	      data += "    </div>\n"
#new end
#------------------------------
	  
            data += "  </div>\n"
#------------------------------
	  
          data += "</div>\n"
          

        data += "</table>\n"
	  
# End of Available Plot List  div

        data += "</div>\n"

# End of left side  div
        
	data += "</div>\n"
	
#===================================================================================


#===================================================================================

	data += "<div style=\"float:left; padding:3px; margin:5px  width:200px;\">\n"

#===================================================================================

# div vor Reference Plot

#===================================================================================
	
        if l_nightly:
	  data += "  <div id=\"arch_ref\">\n"
#ToDo Ref plot mit branch
	  p = "public_html/reference_plot/" + exp_plot_Info
          if os.path.isdir(p):
	    save_date = os.listdir(p)
	    save_date = sorted(save_date)[-1]
            p += "/" + save_date
	    ref_date = os.listdir(p)[-1]
	  
            p2 = p + "/" + ref_date + "/buildbot/"
	    ref_rev = os.listdir(p2)[-1]

            p2 += "/" + ref_rev
	    ref_branch = os.listdir(p2)[-1]
            
	    p3 = p2 + "/" + ref_branch
	    ref_comp = os.listdir(p3)[-1]

            p4 = p3 + "/" + ref_comp
	    ref_build = os.listdir(p4)[-1]
	  
            data += "  <h1>Reference Plot </h1>\n"
            
	    data += "  <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"670\">\n"
	    data += "    <tr>\n    <td>\n"
            data += "      <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"670\">\n"
            data += "        <tr>\n"
            data += "            <td style=\"text-align:left;\"><b>Date:</b>&nbsp;" + ref_date
	    
            data += "            &nbsp;&nbsp;<b>Revision:</b>&nbsp;" 
            data += "<a href=\"https://code.zmaw.de/projects/icon/repository/revisions/" + ref_rev + "\">"
            data += ref_rev + "</a>"
	    
            data += "           &nbsp;&nbsp;<b>Branch:</b>&nbsp;"
            data += "<a href=\"https://code.zmaw.de/projects/icon/repository/show/" + ref_branch.replace('+', '/') + "\">"
            data += "trunk/icon-dev</a>"            
            data += "            &nbsp;&nbsp;<b>Builder:</b>&nbsp;" 
            data += "<a href=\"builders/" + ref_comp + "\">"
            data +=              ref_comp +"</a>"
            
	    data += "            &nbsp;&nbsp;<b>Build:</b>&nbsp;<a href=\"builders/" + ref_comp + "/builds/"+ ref_build + "\">" + ref_build + "</a>\n"
	    data += "</td>\n"
	    data += "<td style=\"text-align:left;\">\n"
	    data += "  <form name=\"replace_form\" method=\"POST\" action=\"plot/replace\" class=\"command replace\">\n"
            data += "    <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\" >\n"
            data += "    <input type=\"hidden\" name=\"rev\" value=\""+ ref_rev + "\" >\n"
            data += "    <input type=\"submit\" name=\"Replace\" value=\"Replace\"  maxlength=\"10\">\n"
	    data += "  </form>\n"            
	    data += "            </td>\n"

	    data += "        </tr>\n      </table>\n"
            data += "    </td></tr>\n    <tr><td style=\"text-align:left; valign:top;\">\n"
          
	    data += "      <img src=\"reference_plot/"
	    data +=   exp_plot_Info + "/" 
	    data += save_date + "/" 
	    data += ref_date + "/buildbot/" 
	    data += ref_rev + "/"
	    data += ref_branch +"/" 
	    data += ref_comp +"/" 
	    data += ref_build +"/" 
	    data += exp_plot_Info + "/plots/"
	    data += exp_file_Info + ".png\" />\n"
	  
	    data += "    </td></tr>\n"
	    data += "  </table>"
          else:
            data += "<h1>Reference Plot</h1>\n"
#	    data += "<form name=\"replace_form\" method=\"POST\" action=\"plot/replace\" class=\"command replace\">"
#            data += "  <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\">\n"
#            data += "  <input type=\"submit\" name=\"Replace_2\" value=\"Replace\" >\n"
#	    data += "</form>"
            data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"670\">\n"
	    data += "  <tr>\n    <td>"
            data += "      <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"670\">\n"
#            data += "        <colgroup>"
#	    data += "          <col width=\"670\">"
#	    data += "        </colgroup>\n<re>\n"
            data += "        <tr>\n"
            data += "            <td style=\"text-align:left;\"><b>Date:</b>&nbsp; xxxx-xx-xx \n"
            data += "            &nbsp;&nbsp;&nbsp;<b>Revision:</b>&nbsp;  ???? \n"
            data += "            &nbsp;&nbsp;&nbsp;<b>Branch:</b>&nbsp;  ???? \n"
            data += "            &nbsp;&nbsp;&nbsp;<b>Builder:</b>&nbsp;  ???? \n"
            data += "            &nbsp;&nbsp;&nbsp;<b>Build:</b>&nbsp;  ???? </td>\n"
            data += "        <tr>\n"
#   	    data += "<td style=\"text-align:left;\">\n"
#	    data += "  <form name=\"replace_form\" method=\"POST\" action=\"plot/replace\" class=\"command replace\">\n"
#            data += "    <input type=\"hidden\" name=\"exp\" value=\"all\" >\n"
#            data += "    <input type=\"hidden\" name=\"rev\" value=\"all\" >\n"
#            data += "    <input type=\"submit\" name=\"Replace\" value=\"Replace\"  maxlength=\"10\">\n"
#	    data += "  </form>\n"
            
	    data += "            </td>\n"

            data += "</tr>\n      </table>\n"
            data += "    </td></tr>\n    <tr><td style=\"text-align:left; valign:top;\">\n"
	    data += "    </td>"
	    data += "  </tr>"
	    data += "</table>"
	    
	  data += "</div>\n "
	  
	
#===================================================================================

# div vor Exp Plots

#===================================================================================
        if l_nightly:
	  data += "<div id=\"ap_ref\">\n"
	else:
	  data += "<div id=\"arch_plots\">\n"
	  
	if ExpPlotError == 0:
	  data += "<h1>Plot Area</h1>\n"
	else:
	  print "Error: " + str(ExpPlotError)
	  data += "<h1><font color=\"#ff0000\">Error</font></h1>\n"

	if ExpPlotError > 0:
	  if ExpPlotError == 1:
	    data += "<h2>Revision List for this timeperode ist emty</h2>"
	    
	  if ExpPlotError == 2:
	    data += "No Plots available for experiment: <br /><b><i>" + exp_plot_Info + "</i><b>" 
	else:
	  p = "public_html/archive/" + Date_from + "/buildbot/" + Rev_from + "/"
	  data += "  <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"660\">\n"
	  
	  file_plotted = False
	  for comp in Archive_Button_Dict['comp']:
	    if comp == selected_builder or selected_builder == "all":
	      for Da in sorted(Comp_Dict[comp], reverse=True):
		for Re in sorted(Comp_Dict[comp][Da], reverse=True):
		  for Br in sorted(Comp_Dict[comp][Da][Re]):
		    if Br == selected_branch or selected_branch == "all":
		      for Bu in sorted(Comp_Dict[comp][Da][Re][Br], reverse=True):
			if Bu == selected_build or selected_build == "all":
			  for Ex in Comp_Dict[comp][Da][Re][Br][Bu]:
			    for Fi in Comp_Dict[comp][Da][Re][Br][Bu][Ex]:
			      if Fi.strip(".png") == exp_file_Info:
				file_plotted = True
				data += "  <tr><td>\n"
				data += "    <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"660\">\n"
	      
				data += "        <colgroup>"
				data += "          <col width=\"670\">"
				data += "        </colgroup>\n"
		
				data += "      <tr>\n"
				data += "        <td style=\"text-align:left;\"><b>Date:</b>&nbsp;" + Da + "\n"
				data += "        &nbsp;&nbsp;&nbsp;<b>Revision:</b>&nbsp;" 
				data += "<a href=\"https://code.zmaw.de/projects/icon/repository/revisions/" + Re + "\">"
				data += Re + "</a> "
				
				data += "        &nbsp;&nbsp;&nbsp;<b>Branch:</b>&nbsp;" 
				data += "<a href=\"https://code.zmaw.de/projects/icon/repository/show/" + Br.replace('+', '/') + "\">"
				data +=  Br.replace('+', '/') + "</a> "
				
				data += "        &nbsp;&nbsp;&nbsp;<b>Builder:</b>&nbsp;" 
				data += "<a href=\"builders/" + comp + "\">"
				data +=  comp + "</a>\n"
				
				data += "        &nbsp;&nbsp;&nbsp;<b>Build:</b>&nbsp;<a href=\"builders/" + comp + "/builds/"+ Bu + "\">" 
				data += Bu  + "</a>"
				
				data += "</a></td>\n"
				data += "        <td><b>&nbsp</b></td>\n"
				data += "</tr>\n"
		  
				data += "    </table>\n"
				data += "  </td></tr>\n"
				data += "  <tr><td style=\"text-align:left;\">\n"        
				br_path = "public_html/archive/" + Da + "/buildbot/"+ Re + "/" + Br + "/" + comp 
				if os.path.isdir(br_path):
				  data += "    <img src=\"archive/" + Da + "/buildbot/"+ Re + "/" + Br + "/" + comp +"/" +  Bu + "/" + Ex + "/plots/" + Fi + "\" />\n"
				else:
				  data += "    <img src=\"archive/" + Da + "/buildbot/"+ Re + "/" + comp +"/" +  Bu + "/" + Ex + "/plots/" + Fi + "\" />\n"
				data += "  </td></tr>\n"
          if not file_plotted:
            data += "<br />No plost available for this archive setting"
	  data += "</table>\n"
        data += "</div>\n"
        data += "</div>\n"
        data += "</div>\n"

#        data += '<hr /><div class="footer">\n'
#        data += '</div>\n'
	data += self.footer(status,request)

        return data

#==========================================================================================
#
#
#
#==========================================================================================
    
    def select(self, req):
      global selected_builder
#------------------------------------------------------------------------------------------
# This routine is used to fined the start- and end-date of the revision range which is given
# by the to values rf and rt.
# If no date is found
#------------------------------------------------------------------------------------------
      def get_Dates (rf,rt):
        if DEBUG:
          print "get_Dates " + rf + " " + rt
	df    = "9999-12-31"
	dt    = "0000-01-01"
	left_rev  = "00000000000"
	right_rev = "99999999999"
	
	p = "public_html/archive/"
	for d in os.listdir(p):
	  pd = p + d + "/buildbot/"
	  for r in os.listdir(pd):
	    if (r >= rf) and (r <= rt):
	      if df > d:
	        df = d        
	      if dt < d:
	        dt = d

            if left_rev < r and r < rf :
	      left_rev = r        
	    if right_rev > r and r > rt:
	      right_rev = r


        if DEBUG:
          print "right_rev " + left_rev + " " + right_rev      
	if dt == "0000-01-01" or df == "9999-12-31":
#	  return get_Dates (left_rev,right_rev)
	  return (df,dt,left_rev,right_rev)
	        
        if DEBUG:
          print "df,dt " + df + " " + dt
        return (df,dt,rf,rt)        
#------------------------------------------------------------------------------------------
      def get_Rev (df,dt):
	rf = "99999999"
	rt = "00000000"
	
	p = "public_html/archive/"
	for d in os.listdir(p):
	  if (d >= df) and (d <= dt):
	    pd = p + d + "/buildbot/"
	    for r in os.listdir(pd):
	      if rf > r:
	        rf = r
	      if rt < r:
	        rt = r
        if DEBUG:
          print "rf,rt " + rf + " " + rt
        return (rf,rt)        
#------------------------------------------------------------------------------------------
      if DEBUG:
        print "======  select ====="
        print req
        print req.args
        print "======  select ====="
      m       = req.args.get("modus",[None])[0]
      e       = req.args.get("exp",[None])[0]
      f       = req.args.get("file",[None])[0]
      builder = req.args.get("builder",[None])[0]
      build   = req.args.get("build",[None])[0]
      fRev    = "0000"
      tRev    = "0000"
      lRev    = False
      lDate   = False
     
      Rev_button_pushed = False
      
      Redirect_Info = "../plot?exp=" + e  + "&modus=" + m
      
      if f.find(e + "_") >= 0:
        Redirect_Info += "&file=" + f

      Redirect_Info += "&build=" + build

#     Check if information all is given in the html-call
      if builder != "all":
	Redirect_Info += "&builder=" + builder
	if selected_builder != builder:
          build = "all"
          
#     Check if the start revision number is given in the html-call
      if "rev_from" in req.args:
        try:
          fRev = req.args["rev_from"][0]
          if int(fRev) < 4774:
              fRev = "4774"
          lRev = True
        except ValueError:
          pass
      
#     Check if the end revision number is given in the html-call
      if "rev_to" in req.args:
        try:
          tRev = req.args["rev_to"][0]
          if int(tRev) < 4774:
              tRev = "4774"
          lRev = True
        except ValueError:
          pass

#     Check if the start date is given in the html-call
      if "date_from" in req.args:
        try:
          fDate = req.args["date_from"][0]
          if fDate < "2011-06-19":
              fDate = "2011-06-19"
          if fDate > time.strftime("%Y-%m-%d",time.localtime(util.now())):
              fDate = time.strftime("%Y-%m-%d",time.localtime(util.now()))
          lDate = True
        except ValueError:
          pass
      
#     Check if the end date is given in the html-call
      if "date_to" in req.args:
        try:
          tDate = req.args["date_to"][0]
          if tDate < "2011-06-19":
              tDate = "2011-06-19"
          if tDate > time.strftime("%Y-%m-%d",time.localtime(util.now())):
              tDate = time.strftime("%Y-%m-%d",time.localtime(util.now()))

          lDate = True
        except ValueError:
          pass

#     Check if the revision- or date-button was pressed
      if "Rev_button" in req.args:
	Rev_button_pushed = True
        lDate = True
        fDate,tDate,fRev,tRev = get_Dates(fRev,tRev)
        
        if fDate == "9999-12-31" or tDate == "0000-01-01":
          fDate,tDate,fRev,tRev = get_Dates(fRev,tRev)
          
        if DEBUG:
          print "fDate,tDate " + fDate + " " + tDate
          print "fRev,tRev " + fRev + " " + tRev
        
      if "Date_button" in req.args:
	Date_button_pushed = True
        lRev = True
        fRev,tRev = get_Rev(fDate,tDate)
            
#     Include the revision start/end values into the new html-call
      if lRev:
        Redirect_Info += "&rev_from=" + fRev      
        Redirect_Info += "&rev_to="   + tRev      
	
#     Include the new start/end date into the new html-call
      if lDate:
        Redirect_Info += "&date_from=" + fDate      
        Redirect_Info += "&date_to=" + tDate
        
#     Include the selected branch info  into the new html-call
      if "branch" in req.args:
        try:
          tBranch = req.args["branch"][0]
          Redirect_Info += "&branch=" + tBranch.replace('+', '/')
        except ValueError:
          pass
	

#     Redirect to the new build html-call
      return Redirect(Redirect_Info)
#      return Redirect("../plot?exp=" + e + "&file=" + f + "&modus=")

#==========================================================================================
#
#
#
#==========================================================================================
    
    def replace(self, req):
      e = req.args.get("exp",[None])[0]
      r = req.args.get("rev",[None])[0]
#      f = req.args.get("file",[None])[0]
      return Redirect("../reference?exp="+e+"&start="+r)

#==========================================================================================
#
#
#
#==========================================================================================

    def default(self, req):
      return Redirect("../plot")
    
#==========================================================================================
#
#
#
#==========================================================================================
    
    def getChild(self, path, req):
      if path == "select":
        return self.select(req)
      elif path == "replace":
        return self.replace(req)
      else:
        return self.default(req)

#==========================================================================================
#
#
#
#==========================================================================================
    
    def get_info(self, request):
        global exp_plot_Info
	global exp_file_Info
        global l_nightly
        global modus
	global selected_builder
	global selected_build
	global selected_branch

	global selected_Date_from
	global selected_Date_to
	
	global selected_Rev_from
	global selected_Rev_to

	exp_plot_Info      = "NotSet"
	exp_file_Info      = "NotSet"
	modus              = "nightly"
	selected_builder   = "all"
	selected_build     = "all"
	selected_branch    = "all"
	selected_Date_from = "NotSet"
	selected_Date_to   = "NotSet"
	selected_Rev_from  = "NotSet"
	selected_Rev_to    = "NotSet"
	
        l_nightly = True
        
        if "exp" in request.args:
          try:
            exp_plot_Info = request.args["exp"][0]
          except ValueError:
            pass

        if "branch" in request.args:
          try:
            selected_branch = request.args["branch"][0]
	    selected_branch = selected_branch.replace('/', '+')
#            print "branch_if " + selected_branch
          except ValueError:
            pass
	  
#        print "branch " + selected_branch
	
        if "builder" in request.args:
          try:
            selected_builder = request.args["builder"][0]
          except ValueError:
            pass
        
        if "build" in request.args:
          try:
            selected_build = request.args["build"][0]
          except ValueError:
            pass

        if "modus" in request.args:
          try:
	    if request.args["modus"][0] == "archive":
	      modus = "archive"
              l_nightly = False
          except ValueError:
            pass
        
	if "file" in request.args:
	  try:
	    exp_file_Info = request.args["file"][0]
	  except ValueError:
	    pass
							     
	if "date_from" in request.args:
	  try:
	    selected_Date_from = request.args["date_from"][0]
	  except ValueError:
	    pass

	if "date_to" in request.args:
	  try:
	    selected_Date_to = request.args["date_to"][0]

	  except ValueError:
	    pass

	if "rev_from" in request.args:
	  try:
	    selected_Rev_from = request.args["rev_from"][0]
	  except ValueError:
	    pass

	if "rev_to" in request.args:
	  try:
	    selected_Rev_to = request.args["rev_to"][0]
	  except ValueError:
	    pass

	return None

#==========================================================================================
#
#
#
#==========================================================================================
    
    def footer(self, status, req):
        global RevisionNr
        # TODO: this stuff should be generated by a template of some sort
        projectURL = status.getProjectURL()
        projectName = status.getProjectName()
        data = '<div class="footer">\n'

        welcomeurl = self.path_to_root(req) + "index.html"
        data += '  [<a href="%s">mainpage</a>]\n' % welcomeurl
        data += "  <br />\n"

#        data += "  <br />\n"
        data += ("  Page built: " +
                 time.strftime("%a %d %b %Y %H:%M:%S",
                               time.localtime(util.now()))
                 + "\n")
        data += '</div>\n'

        return data
#================================== NEW ================================
