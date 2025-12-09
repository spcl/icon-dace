# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

require 'thread'
#= config =====================================================================#
dtimes                       = [600, 900, 1800]
k_veloc_hs                   = [1.0E+5,1.0E+4,1.0E+3]
k_veloc_vs                   = [1.0E-3,1.0E-2,1.0E-1]
k_pot_temp_hs                = [1.0E+2,1.0E+3,1.0E+4]
k_pot_temp_vs                = [1.0E-2,1.0E-3,1.0E-4]
expl_vertical_velocity_diffs = [1]
testcases                    = [33]
#= internals ==================================================================#
$configFile   = "variations.log"
$configString = ''
$inputScript  = 'exp.pk_lsmfl_4variations'
$inputContent = File.open($inputScript).read
$outputScript = 'exp.pk_lsmfl'
$outputs      = []
$threads      = []
THREADEDRUN   = false
#= lib ========================================================================#
def createRunScriptSetting( dtime,
                            k_veloc_h,
                            k_veloc_v,
                            k_pot_temp_h,
                            k_pot_temp_v,
                            expl_vertical_velocity_diff,
                            testcase,
                            index)

  fileContent = $inputContent.dup
  fileContent.sub!(/DTIME/                      , dtime.to_s)
  fileContent.sub!(/K_VELOC_H/                  , k_veloc_h.to_s)
  fileContent.sub!(/K_VELOC_V/                  , k_veloc_v.to_s)
  fileContent.sub!(/K_POT_TEMP_H/               , k_pot_temp_h.to_s)
  fileContent.sub!(/K_POT_TEMP_V/               , k_pot_temp_v.to_s)
  fileContent.sub!(/EXPL_VERTICAL_VELOCITY_DIFF/, expl_vertical_velocity_diff.to_s)
  fileContent.sub!(/TESTCASEOCE/                , testcase.to_s)

  $outputs <<  fileContent
end

def createVariationDatabase(dtime,
                            k_veloc_h,
                            k_veloc_v,
                            k_pot_temp_h,
                            k_pot_temp_v,
                            expl_vertical_velocity_diff,
                            testcase,
                            index)

    if index == 0 then
      $configString << '|runtime|' << %w[index dtime k_veloc_h k_veloc_v k_pot_temp_h k_pot_temp_v expl_vertical_velocity_diff testcase].join("|") << '|' << "\n"
    end
    $configString << '||'<< [index,dtime,k_veloc_h,k_veloc_v,k_pot_temp_h,k_pot_temp_v,expl_vertical_velocity_diff,testcase].join('|') << '|'<<"\n"
end

def writeLogFile(filename)
  File.open(filename,"w") {|f|
    f << $configString
  }
end
def writeRunScript(index,content)
  puts index
  File.open($outputScript+"_#{index}","w") {|f| f << content}
end
#================================================================================#
# MAIN SCRIPT
runscriptIndex=0
dtimes.each {|dtime|
  k_veloc_hs.each{|k_veloc_h|
    k_veloc_vs.each {|k_veloc_v|
      k_pot_temp_hs.each {|k_pot_temp_h|
        k_pot_temp_vs.each {|k_pot_temp_v|
          expl_vertical_velocity_diffs.each {|expl_vertical_velocity_diffs|
            testcases.each {|testcase|

              createRunScriptSetting(dtime,
                                     k_veloc_h,
                                     k_veloc_v,
                                     k_pot_temp_h,
                                     k_pot_temp_v,
                                     expl_vertical_velocity_diffs,
                                     testcase,
                                     runscriptIndex)
              createVariationDatabase(dtime,
                                     k_veloc_h,
                                     k_veloc_v,
                                     k_pot_temp_h,
                                     k_pot_temp_v,
                                     expl_vertical_velocity_diffs,
                                     testcase,
                                     runscriptIndex)
              runscriptIndex += 1
            }}}}}}}

writeLogFile($configFile)

if THREADEDRUN
  $outputs.each_with_index {|o,i|
    $threads << Thread.new(o,i) {|filecontent,index|
      writeRunScript(index, filecontent)
    }
  }
  $threads.each {|th| th.join}
else
  $outputs.each_with_index {|filecontent,index|
    writeRunScript(index, filecontent)
  }
end
