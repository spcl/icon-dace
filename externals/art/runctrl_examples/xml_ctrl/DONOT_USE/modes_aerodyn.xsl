<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template match="/">
    <html>
      <body>
        <h2>The single aerosol definition</h2>
        <table border="2">
          <tr bgcolor="9acd32">
            <th>aerosol</th>            
            <th>kind</th>            
            <th>d_gn</th>            
            <th>sigma_g</th>
            <th>shift2larger</th>
            <th>shift_diam</th>
            <th>condensation</th>
            <th>icoag</th>
          </tr>
          <xsl:for-each select="modes/aerosol">
            <tr>
              <td><xsl:value-of select="@id" /></td>
              <td><xsl:value-of select="kind" /></td>
              <td><xsl:value-of select="d_gn" /></td>
              <td><xsl:value-of select="sigma_g" /></td>
              <td><xsl:value-of select="shift2larger" /></td>
              <td><xsl:value-of select="shift_diam" /></td>
              <td><xsl:value-of select="condensation" /></td>
              <td><xsl:value-of select="icoag" /></td>              
            </tr>
          </xsl:for-each>
        </table>
      </body>
    </html>
  </xsl:template>
</xsl:stylesheet>
