<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
 <xsl:template match="/">
   <![CDATA[
       === XSL Utility Script: Mark Files for Deletion in DWD release =================== 
                                                                                          
       Mark a file to be deleted from the DWD release:       				

          svn propset dwd_release remove "file"

       Get a list of all files/directories                                            
       with attribute "dwd_release", value=remove		    

         svn -R proplist --verbose --xml | xsltproc scripts/postprocessing/dwd_release.xslt - 
      										   
       2014-05-28 : F. Prill, DWD                                                         
                                                                                          
       ================================================================================== 

Files to be deleted in DWD release:   ]]>

   <xsl:apply-templates select="properties/target/property[@name = 'dwd_release']"/>
 </xsl:template>

 <xsl:template match="property[@name = 'dwd_release']">
   <xsl:text>&#xa;</xsl:text>
   <xsl:if test="text()='remove'">
     <xsl:value-of select="../@path"/>
   </xsl:if>
 </xsl:template>

</xsl:stylesheet>


