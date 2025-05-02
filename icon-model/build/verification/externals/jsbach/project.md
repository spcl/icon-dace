---
project: ICON-Land
author: MPI-M, MPI-BGC
author_description: [AUTHORS](page/1-general/02-authors.html)
website: https://icon-model.org
email: reiner.schnur@mpimet.mpg.de
project_dir: .
src_dir: ./src
include: ./doc/include
exclude: cdipio.inc
         jsb4_driver_echam.f90
page_dir: ./doc/guide
media_dir: ./doc/media
output_dir: ./doc/html
docmark: <
docmark_alt: ^
predocmark: >
predocmark_alt: |
summary: <img src="media/MPIM.png" width="100%" alt="MPI/M Logo">
         <p style="text-align:center;font-weight:bold;font-size:300%;">ICON-Land Documentation</p>
display: public
         protected
         private
css: ./doc/resources/user.css
max_frontpage_items: 0
source: true
incl_src: false
proc_internals: true
search: true
warn: false
graph: false
graph_maxdepth: 2
graph_dir: ./doc/graphs
coloured_edges: true
md_extensions: markdown.extensions.toc
               markdown.extensions.extra
               markdown.extensions.footnotes
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
---

------------------------------

@warning This documentation is work in progress!

@note
Use the navigation bar at the top of the screen to browse the guide as well as inline source
code documentation (modules, procedures, interfaces and derived types). Source files are not
listed as a whole but the views of subroutines and functions in the modules show their source
at the bottom of the page.

------------------------------

To get startet, visit the guide from the top navigation, or jump directly to

- [General Guide](page/1-general/index.html)
- [Science Guide](page/2-science-guide/index.html)
- [User Guide](page/3-user-guide/index.html)
- [Developer Guide](page/4-developer-guide/index.html)

------------------------------

{!LICENSE.md!}

